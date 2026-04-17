#ifndef PHOTON_SOURCE_HOOK_H_
#define PHOTON_SOURCE_HOOK_H_

#include "core/grid_registry.h"
#include "core/voxel_snapshot.h"
#include "core/signal_board.h"

namespace bdm {
namespace skibidy {

// Photon transport source hook: diffusion approximation for light in tissue.
//
// Models SLM-shaped light delivery for optogenetic stimulation or phototherapy.
// The fluence rate field is driven by a surface irradiance boundary condition
// (beam profile) and attenuated by tissue absorption and scattering via the
// diffusion approximation: D = 1 / (3 * (mu_a + mu_s')).
//
// Produces:
//   - Fluence field: photon density (light intensity) throughout tissue
//   - Opsin activation field: channelrhodopsin open-state fraction
// Modifies:
//   - ROS: phototoxicity from excessive light dose
//   - Temperature: absorbed photon energy converted to heat
//
// Jacques 2013 (doi:10.1088/0031-9155/58/11/R37)
// Mardinly et al. 2018 (doi:10.1038/s41593-018-0139-8)
// Fenno et al. 2011 (doi:10.1146/annurev-neuro-061010-113817)
struct PhotonSourceHook {
  DiffusionGrid* fluence_grid = nullptr;
  DiffusionGrid* opsin_grid = nullptr;
  DiffusionGrid* ros_grid = nullptr;
  DiffusionGrid* temp_grid = nullptr;
  const SimParam* sp_ = nullptr;
  bool active = false;

  inline void Init(const GridRegistry& reg, SignalBoard& sig) {
    sp_ = reg.Params();
    active = sp_->photon.enabled;
    if (!active) return;
    fluence_grid = reg.Get(fields::kFluenceId);
    opsin_grid = reg.Get(fields::kOpsinId);
    if (sp_->ros.enabled)
      ros_grid = reg.Get(fields::kROSId);
    if (sp_->temperature.enabled)
      temp_grid = reg.Get(fields::kTemperatureId);
  }

  // Dermal: light attenuation and opsin activation in tissue volume.
  inline void ApplyDermal(const VoxelSnapshot& snap, SignalBoard& sig) {
    if (!snap.in_wound) return;

    real_t dt = snap.dt;

    // Beer-Lambert source term: inject photons at surface, attenuate with depth.
    // Surface voxels (z near wound surface) receive irradiance scaled by beam
    // profile. Deeper voxels rely on diffusion grid transport.
    if (fluence_grid) {
      real_t z_surface = sp_->volume_z_cornified;
      real_t depth = z_surface - snap.z;  // distance from surface into tissue
      if (depth >= 0 && depth < 1.0) {
        // Surface layer: apply beam profile as source term
        real_t dx = snap.x - sp_->photon.beam_center_x;
        real_t dy = snap.y - sp_->photon.beam_center_y;
        real_t r2 = dx * dx + dy * dy;
        real_t beam_r2 = sp_->photon.beam_radius * sp_->photon.beam_radius;
        if (beam_r2 > 0 && r2 < beam_r2 * 4.0) {
          // Gaussian beam profile
          real_t intensity = sp_->photon.irradiance *
                             std::exp(-2.0 * r2 / beam_r2);
          fluence_grid->ChangeConcentrationBy(snap.idx, intensity * dt);
        }
      }

      // Absorption sink: mu_a removes photons from fluence field.
      // Uses per-voxel material absorption (tissue-appropriate attenuation).
      real_t fluence = fluence_grid->GetConcentration(snap.idx);
      if (fluence > 1e-10) {
        real_t mu_a = snap.mat ? snap.mat->absorption
                               : sp_->photon.absorption_coeff;
        real_t absorbed = mu_a * fluence * dt;
        fluence_grid->ChangeConcentrationBy(snap.idx, -absorbed);

        // Thermal coupling: absorbed photons heat tissue.
        // Scaled by material conductivity for tissue-appropriate heating.
        if (temp_grid) {
          real_t coupling = sp_->photon.thermal_coupling;
          if (snap.mat) {
            // Lower conductivity = more local heating (inverse scaling)
            coupling *= 0.5 / std::max(snap.mat->conductivity, 0.01);
          }
          temp_grid->ChangeConcentrationBy(snap.idx, coupling * absorbed);
        }

        // Phototoxicity: high fluence generates ROS
        if (ros_grid && fluence > sp_->photon.phototoxicity_threshold) {
          real_t excess = fluence - sp_->photon.phototoxicity_threshold;
          ros_grid->ChangeConcentrationBy(snap.idx,
              sp_->photon.phototoxicity_rate * excess * dt);
        }
      }
    }

    UpdateOpsin(snap.idx, dt);
  }

  // Two-state opsin kinetics: dO/dt = k_on * L * (sat - O) - k_off * O.
  // Clamped to [0, saturation].
  inline void UpdateOpsin(size_t idx, real_t dt) {
    if (!opsin_grid || !fluence_grid) return;
    real_t fluence = fluence_grid->GetConcentration(idx);
    real_t opsin = opsin_grid->GetConcentration(idx);

    real_t k_on = sp_->photon.opsin_activation_rate;
    real_t k_off = sp_->photon.opsin_deactivation_rate;
    real_t sat = sp_->photon.opsin_saturation;

    real_t delta = k_on * fluence * (sat - opsin) * dt
                 - k_off * opsin * dt;
    if (std::abs(delta) > 1e-12) {
      opsin_grid->ChangeConcentrationBy(idx, delta);
      real_t new_val = opsin_grid->GetConcentration(idx);
      if (new_val < 0) {
        opsin_grid->ChangeConcentrationBy(idx, -new_val);
      } else if (new_val > sat) {
        opsin_grid->ChangeConcentrationBy(idx, sat - new_val);
      }
    }
  }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // PHOTON_SOURCE_HOOK_H_
