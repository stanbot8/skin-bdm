#ifndef FIBROBLAST_BEHAVIOR_H_
#define FIBROBLAST_BEHAVIOR_H_

#include <cmath>
#include <iostream>

#include "fibroblast/fibroblast.h"
#include "immune/immune_helpers.h"
#include "core/continuum_handoff.h"
#include "core/field_names.h"
#include "infra/sim_param.h"
#include "infra/util.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// Behavior: fibroblast lifecycle, TGF-beta-driven state transitions,
// collagen deposition, and migration.
//
// States: quiescent -> activated -> myofibroblast -> removed
// Transitions driven by local TGF-beta concentration.
// Only myofibroblasts produce TGF-beta (positive feedback) and deposit
// collagen (= mechanistic scar output).
// ---------------------------------------------------------------------------
struct FibroblastBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(FibroblastBehavior, Behavior, 1);

  FibroblastBehavior() {}
  virtual ~FibroblastBehavior() {}

  void Run(Agent* agent) override {
    auto* cell = dynamic_cast<Fibroblast*>(agent);
    if (!cell) return;

    auto* sim = Simulation::GetActive();
    auto* sp = sim->GetParam()->Get<SimParam>();
    auto* rm = sim->GetResourceManager();

    cell->IncrementAge();
    cell->IncrementStateAge();
    int age = cell->GetAge();

    // --- Hard lifespan limit ---
    int eff_lifespan = sp->fibroblast.lifespan;
    if (sp->diabetic.mode) {
      eff_lifespan = static_cast<int>(
          sp->fibroblast.lifespan * sp->diabetic.fibroblast_lifespan_factor);
    }
    if (age > eff_lifespan) {
      ContinuumHandoff(cell);
      cell->RemoveFromSimulation();
      return;
    }

    // --- Read local TGF-beta ---
    Real3 pos = cell->GetPosition();
    Real3 qpos = ClampToBounds(pos, sim->GetParam());

    auto* tgfb_grid = rm->GetDiffusionGrid(fields::kTGFBetaId);
    real_t local_tgfb = tgfb_grid->GetValue(qpos);

    // --- Receptor-mediated TGF-beta endocytosis ---
    // TGF-beta binds TbRII, recruits TbRI, complex internalized via
    // clathrin-coated pits and degraded in lysosomes. Clearance scales
    // with cell density. (Vilar et al. 2006, doi:10.1016/j.jtbi.2006.03.024)
    // Myofibroblasts upregulate TbRI/TbRII expression ~5x as part of the
    // autocrine loop maintaining the contractile phenotype, internalizing
    // proportionally more ligand. (Tomasek et al. 2002, doi:10.1038/nrm809)
    {
      real_t eff_tgfb = local_tgfb;
      if (cell->GetFibroblastState() == kMyofibroblast)
        eff_tgfb *= 5.0;
      else if (cell->GetFibroblastState() == kFibroActivated)
        eff_tgfb *= 2.0;
      immune::ConsumeTGFBeta(qpos, sim, sp, eff_tgfb);
    }

    // --- State transitions ---
    auto state = cell->GetFibroblastState();
    int state_age = cell->GetStateAge();

    if (state == kFibroQuiescent) {
      int eff_activation_delay = sp->fibroblast.activation_delay;
      if (sp->diabetic.mode) {
        eff_activation_delay = static_cast<int>(
            sp->fibroblast.activation_delay *
            sp->diabetic.fibroblast_activation_factor);
      }
      if (local_tgfb > sp->fibroblast.activation_threshold ||
          state_age > eff_activation_delay) {
        cell->SetFibroblastState(kFibroActivated);
        if (sp->debug_fibroblast) {
          std::cout << "[fibroblast] quiescent->activated tgfb=" << local_tgfb
                    << " age=" << state_age << std::endl;
        }
      }
    } else if (state == kFibroActivated) {
      // YAP/TAZ mechanosensing: stiff substrate promotes myofibroblast
      // differentiation by lowering TGF-beta threshold. Nuclear YAP/TAZ
      // drives alpha-SMA expression on stiff matrices.
      // Dupont et al. 2011 (doi:10.1038/nature10137)
      real_t myofib_thresh = sp->fibroblast.myofibroblast_threshold;
      if (sp->mechanotransduction.enabled) {
        auto* stiff_grid = rm->GetDiffusionGrid(fields::kStiffnessId);
        real_t stiffness = stiff_grid->GetValue(qpos);
        if (stiffness > sp->stiffness_yap_threshold) {
          real_t yap_factor = (stiffness - sp->stiffness_yap_threshold) /
                              (1.0 - sp->stiffness_yap_threshold);
          myofib_thresh *= (1.0 - 0.5 * yap_factor);  // up to 50% lower
        }
      }
      if (local_tgfb > myofib_thresh &&
          state_age > sp->fibroblast.myofibroblast_delay) {
        cell->SetFibroblastState(kMyofibroblast);
        if (sp->debug_fibroblast) {
          std::cout << "[fibroblast] activated->myofibroblast tgfb=" << local_tgfb
                    << " age=" << state_age << std::endl;
        }
      }
      // Stochastic apoptosis/emigration for activated fibroblasts: not all
      // recruited fibroblasts differentiate; some undergo apoptosis or
      // emigrate back to surrounding tissue (Grinnell 1994).
      if (state_age > sp->fibroblast.activated_apoptosis_onset &&
          age > sp->fibroblast.min_lifespan) {
        auto* random = sim->GetRandom();
        if (random->Uniform(0, 1) < sp->fibroblast.activated_apoptosis_rate) {
          ContinuumHandoff(cell);
          cell->RemoveFromSimulation();
          return;
        }
      }
    } else if (state == kMyofibroblast) {
      // Immediate apoptosis if TGF-beta drops to near-zero
      if (local_tgfb < sp->fibroblast.apoptosis_threshold &&
          age > sp->fibroblast.min_lifespan) {
        ContinuumHandoff(cell);
        cell->RemoveFromSimulation();
        return;
      }
      // Stochastic apoptosis program (Desmouliere 1995): after onset delay,
      // myofibroblasts undergo apoptosis at a per-step probability
      if (state_age > sp->fibroblast.apoptosis_onset) {
        auto* random = sim->GetRandom();
        if (random->Uniform(0, 1) < sp->fibroblast.apoptosis_rate) {
          ContinuumHandoff(cell);
          cell->RemoveFromSimulation();
          return;
        }
      }
    }

    // --- TGF-beta production (myofibroblasts only) ---
    // Production tapers with state age (mechanical tension resolves as wound
    // contracts; Tomasek et al. 2002, doi:10.1038/nrm809).
    if (cell->GetFibroblastState() == kMyofibroblast) {
      real_t tgfb_rate = sp->fibroblast.tgfb_rate;
      // Exponential taper: halves every ~580 steps (~2.4 days)
      real_t taper = std::exp(-sp->fibroblast.tgfb_taper * cell->GetStateAge());
      tgfb_rate *= taper;
      ScaledGrid sg(tgfb_grid, sp);
      sg.AgentDeposit(sg.Index(qpos), tgfb_rate);
      if (sp->debug_fibroblast) {
        std::cout << "[fibroblast] tgfb pos=(" << pos[0] << "," << pos[1]
                  << "," << pos[2] << ") rate=" << tgfb_rate
                  << " scaled=" << tgfb_rate * sg.agent_factor << std::endl;
      }
    }

    // --- Collagen deposition (myofibroblasts only) ---
    // Collagen is non-diffusing (structural ECM) -- no ScaledGrid.
    if (cell->GetFibroblastState() == kMyofibroblast) {
      auto* col_grid = rm->GetDiffusionGrid(fields::kCollagenId);
      size_t col_idx = col_grid->GetBoxIndex(qpos);
      real_t deposit;
      if (sp->mech_collagen_deposition) {
        // Mechanistic: constitutive + TGF-b-responsive collagen synthesis.
        // Basal: epigenetically locked myofibroblast program (constitutive).
        // Responsive: Michaelis-Menten TGF-b receptor occupancy modulation.
        deposit = sp->mech_collagen_basal +
                  sp->mech_collagen_vmax *
                  local_tgfb / (sp->mech_collagen_tgfb_km + local_tgfb);
      } else {
        // Parametric: constant rate (Murphy et al. 2012)
        deposit = sp->fibroblast.collagen_deposition_rate;
      }
      if (sp->diabetic.mode) {
        deposit *= sp->diabetic.collagen_factor;
      }
      // Lactate boost (Hunt et al. 2007, doi:10.1089/ten.2007.0115)
      if (sp->lactate.enabled) {
        auto* lac_grid = rm->GetDiffusionGrid(fields::kLactateId);
        if (lac_grid) {
          real_t lac = lac_grid->GetValue(qpos);
          deposit *= (1.0 + sp->lactate.collagen_boost * lac);
        }
      }
      // NO suppression (Schaffer et al. 1996, doi:10.1016/S0022-4804(96)80068-5)
      if (sp->nitric_oxide.enabled) {
        auto* no_grid = rm->GetDiffusionGrid(fields::kNitricOxideId);
        if (no_grid) {
          real_t no_val = no_grid->GetValue(qpos);
          deposit *= std::max(0.0, 1.0 - sp->no_collagen_suppression * no_val);
        }
      }
      // O2-dependent hydroxylation (Myllyharju 2003)
      if (sp->fibroblast.collagen_o2_half_max > 0) {
        auto* o2_grid = rm->GetDiffusionGrid(fields::kOxygenId);
        if (o2_grid) {
          real_t o2 = std::max(static_cast<real_t>(0), o2_grid->GetValue(qpos));
          deposit *= o2 / (sp->fibroblast.collagen_o2_half_max + o2);
        }
      }
      col_grid->ChangeConcentrationBy(col_idx, deposit);
    }

    // --- MMP production (activated and myofibroblast states) ---
    // Fibroblasts produce pro-MMP-1 (collagenase) and pro-MMP-3 (stromelysin)
    // as zymogens; activation occurs extracellularly (Visse & Nagase 2003).
    if (sp->mmp.enabled &&
        cell->GetFibroblastState() != kFibroQuiescent) {
      ScaledGrid sg(rm->GetDiffusionGrid(fields::kProMMPId), sp);
      real_t mmp_rate = sp->mmp.fibroblast_rate;
      if (sp->diabetic.mode) {
        mmp_rate *= sp->diabetic.mmp_factor;
      }
      sg.AgentDeposit(sg.Index(qpos), mmp_rate);
    }

    // --- TIMP production (activated and myofibroblast states) ---
    // Fibroblasts are the primary TIMP-1 source during the proliferative
    // phase, providing delayed negative feedback on MMP activity.
    // Vaalamo et al. 1999 (doi:10.1046/j.1523-1747.1999.00508.x)
    if (sp->mmp.enabled &&
        cell->GetFibroblastState() != kFibroQuiescent) {
      immune::ProduceTIMP(qpos, sim, sp, sp->mmp.timp_fibroblast_rate);
    }

    // --- Fibronectin deposition (activated and myofibroblast states) ---
    // Fibroblasts deposit cellular fibronectin as provisional matrix
    // scaffold for keratinocyte migration (Clark 1990).
    // TGF-beta stimulates FN production via Smad3 signaling
    // (Ignotz & Massague 1986, doi:10.1016/0092-8674(86)90438-9).
    // Myofibroblasts produce EDA-FN splice variant at reduced rate
    // (Serini et al. 1998, doi:10.1083/jcb.142.3.873).
    // Fibronectin is non-diffusing ECM, no ScaledGrid.
    if (sp->fibronectin.enabled &&
        cell->GetFibroblastState() != kFibroQuiescent) {
      auto* fn_grid = rm->GetDiffusionGrid(fields::kFibronectinId);
      size_t fn_idx = fn_grid->GetBoxIndex(qpos);
      real_t fn_cur = fn_grid->GetConcentration(fn_idx);
      real_t fn_cap = sp->fibronectin.carrying_capacity;
      real_t deposit = sp->fibronectin.deposition_rate;
      // Myofibroblasts produce EDA-FN at reduced rate (collagen-dominant)
      if (cell->GetFibroblastState() == kMyofibroblast) {
        deposit *= sp->myofibroblast_fn_fraction;
      }
      // TGF-beta coupling: Michaelis-Menten upregulation
      real_t tgfb_boost = 1.0 + sp->fn_tgfb_max_boost *
          local_tgfb / (sp->fn_tgfb_half_maximal + local_tgfb);
      deposit *= tgfb_boost;
      // Carrying capacity saturation
      if (fn_cap > 0 && fn_cur > 0) {
        deposit *= std::max(static_cast<real_t>(0), 1.0 - fn_cur / fn_cap);
      }
      // Collagen-dependent suppression: as collagen cross-links mature,
      // fibroblasts shift from FN secretion to collagen synthesis
      real_t col_fn_switch = sp->fibroblast.collagen_fn_transition;
      if (col_fn_switch > 0) {
        auto* col_grid = rm->GetDiffusionGrid(fields::kCollagenId);
        real_t col_local = col_grid->GetConcentration(col_grid->GetBoxIndex(qpos));
        deposit *= std::max(static_cast<real_t>(0), 1.0 - col_local / col_fn_switch);
      }
      if (deposit > 1e-10) {
        fn_grid->ChangeConcentrationBy(fn_idx, deposit);
      }
    }

    // --- Tropoelastin production (activated and myofibroblast states) ---
    // Fibroblasts secrete tropoelastin that self-assembles into elastic
    // fibers via lysyl oxidase cross-linking (Sage 1982, Kielty 2002).
    // Elastin is non-diffusing structural ECM -- no ScaledGrid.
    if (sp->elastin.enabled &&
        cell->GetFibroblastState() != kFibroQuiescent) {
      auto* el_grid = rm->GetDiffusionGrid(fields::kElastinId);
      el_grid->ChangeConcentrationBy(el_grid->GetBoxIndex(qpos),
                                     sp->elastin.production_rate);
    }

    // --- Hyaluronan synthesis (activated and myofibroblast states) ---
    // Fibroblasts produce HA via HAS2 synthase. HA provides hydrated
    // ground substance for cell migration (Toole 2004, Stern 2006).
    if (sp->hyaluronan.enabled &&
        cell->GetFibroblastState() != kFibroQuiescent) {
      ScaledGrid sg(rm->GetDiffusionGrid(fields::kHyaluronanId), sp);
      sg.AgentDeposit(sg.Index(qpos), sp->hyaluronan.production_rate);
    }

    // --- Migration (activated and myofibroblast states) ---
    if (cell->GetFibroblastState() == kFibroQuiescent) return;
    if (!sp->wound.enabled) return;

    // AGE-mediated migration impairment: glycated ECM reduces integrin
    // binding and focal adhesion turnover, slowing fibroblast crawling.
    // Rana et al. 2007 (doi:10.1109/nebc.2007.4413357)
    real_t age_speed_factor = 1.0;
    if (sp->diabetic.mode && sp->glucose_mod.enabled &&
        sp->age_fibroblast_migration_impair > 0) {
      auto* age_grid = rm->GetDiffusionGrid(fields::kAGEId);
      if (age_grid) {
        real_t age_val = age_grid->GetValue(qpos);
        age_speed_factor = std::max(static_cast<real_t>(0),
            1.0 - sp->age_fibroblast_migration_impair * age_val);
      }
    }

    // Edema migration impairment: interstitial fluid accumulation creates
    // hydrostatic resistance that slows cell crawling through the ECM.
    // Rutkowski & Swartz 2007 (doi:10.1016/j.tcb.2006.11.007)
    real_t edema_speed_factor = 1.0;
    if (sp->lymphatic.enabled) {
      auto* edema_grid = rm->GetDiffusionGrid(fields::kEdemaId);
      if (edema_grid) {
        real_t edema_val = edema_grid->GetValue(qpos);
        edema_speed_factor = std::max(static_cast<real_t>(0),
            1.0 - sp->edema_migration_impairment * edema_val);
      }
    }

    real_t combined_speed_factor = age_speed_factor * edema_speed_factor;

    // Wound contraction: myofibroblast traction force pulls wound edges
    // inward, modulated by substrate stiffness (stiffer = more force).
    // Hinz 2007 (doi:10.1038/sj.jid.5700613)
    if (cell->GetFibroblastState() == kMyofibroblast &&
        sp->mechanotransduction.enabled) {
      auto* stiff_grid = rm->GetDiffusionGrid(fields::kStiffnessId);
      real_t stiffness = stiff_grid->GetValue(qpos);
      real_t contraction = sp->stiffness_contraction_rate * stiffness;
      if (contraction > 1e-10) {
        real_t cdx = pos[0] - sp->wound.center_x;
        real_t cdy = pos[1] - sp->wound.center_y;
        real_t cdist = std::sqrt(cdx * cdx + cdy * cdy);
        if (cdist > 1e-6) {
          real_t cinv = 1.0 / cdist;
          cell->SetTractorForce({-cdx * cinv * contraction,
                                 -cdy * cinv * contraction, 0});
          return;
        }
      }
    }

    // Galvanotaxis: voltage gradient biases migration toward wound center.
    // The lateral electric field from collapsed TEP at wound edge guides
    // cell migration via electrophoresis of membrane receptors.
    // Zhao et al. 2006 (doi:10.1038/nature04925)
    Real3 galvano_dir = {0, 0, 0};
    if (sp->bioelectric.enabled) {
      auto* volt_grid = rm->GetDiffusionGrid(fields::kVoltageId);
      Real3 vgrad = {0, 0, 0};
      volt_grid->GetGradient(qpos, &vgrad, false);
      real_t vmag = std::sqrt(vgrad[0]*vgrad[0] + vgrad[1]*vgrad[1] +
                              vgrad[2]*vgrad[2]);
      if (vmag > 1e-8) {
        real_t vinv = 1.0 / vmag;
        // Cells migrate toward cathode (low voltage = wound center),
        // so follow negative gradient
        galvano_dir = {-vgrad[0] * vinv, -vgrad[1] * vinv, 0};
      }
    }

    // Chemotaxis: follow TGF-beta gradient
    if (sp->chemotaxis_enabled) {
      Real3 gradient = {0, 0, 0};
      tgfb_grid->GetGradient(qpos, &gradient, false);
      real_t grad_mag = std::sqrt(
          gradient[0] * gradient[0] +
          gradient[1] * gradient[1] +
          gradient[2] * gradient[2]);

      if (grad_mag > 1e-8) {
        real_t inv_mag = 1.0 / grad_mag;
        real_t chem_x = gradient[0] * inv_mag;
        real_t chem_y = gradient[1] * inv_mag;
        // Blend chemotaxis with galvanotaxis
        if (sp->bioelectric.enabled && (galvano_dir[0] != 0 || galvano_dir[1] != 0)) {
          real_t gw = sp->galvanotaxis_strength;
          chem_x = (1.0 - gw) * chem_x + gw * galvano_dir[0];
          chem_y = (1.0 - gw) * chem_y + gw * galvano_dir[1];
          real_t blend_mag = std::sqrt(chem_x*chem_x + chem_y*chem_y);
          if (blend_mag > 1e-8) {
            chem_x /= blend_mag;
            chem_y /= blend_mag;
          }
        }
        real_t speed = sp->fibroblast.migration_speed *
                       sp->chemotaxis_speed_scale * grad_mag * combined_speed_factor;
        cell->SetTractorForce({chem_x * speed, chem_y * speed, 0});
        return;
      }
    }

    // Geometric fallback: migrate toward wound center
    real_t dx = pos[0] - sp->wound.center_x;
    real_t dy = pos[1] - sp->wound.center_y;
    real_t dist = std::sqrt(dx * dx + dy * dy);
    if (dist < 1e-6) return;

    real_t inv_dist = 1.0 / dist;
    real_t dir_x = -dx * inv_dist;
    real_t dir_y = -dy * inv_dist;

    // Blend with galvanotaxis
    if (sp->bioelectric.enabled && (galvano_dir[0] != 0 || galvano_dir[1] != 0)) {
      real_t gw = sp->galvanotaxis_strength;
      dir_x = (1.0 - gw) * dir_x + gw * galvano_dir[0];
      dir_y = (1.0 - gw) * dir_y + gw * galvano_dir[1];
      real_t blend_mag = std::sqrt(dir_x*dir_x + dir_y*dir_y);
      if (blend_mag > 1e-8) {
        dir_x /= blend_mag;
        dir_y /= blend_mag;
      }
    }

    real_t r = sp->wound.radius;
    real_t speed = sp->fibroblast.migration_speed * (dist / r) * combined_speed_factor;

    cell->SetTractorForce({dir_x * speed, dir_y * speed, 0});
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIBROBLAST_BEHAVIOR_H_
