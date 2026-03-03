#ifndef MATERIAL_H_
#define MATERIAL_H_

#include <array>
#include <string>
#include <unordered_map>
#include "biodynamo.h"

namespace bdm {
namespace skibidy {

// ---------------------------------------------------------------------------
// Material system: tissue type definitions with physical properties.
//
// Three-layer architecture:
//   Material (what tissue IS)  -> defines optical/mechanical/thermal properties
//   Module   (what biology CAN happen) -> field registration, hook code
//   Profile  (patient modifiers) -> diabetic, aged, etc.
//
// Each voxel carries a material_id (uint8_t). Hooks query per-voxel material
// properties through the VoxelSnapshot, enabling tissue-appropriate behavior
// without per-module branching.
//
// Material IDs:
//   0 = void (outside tissue domain)
//   1-9 = skin layers (epidermis, dermis, hypodermis)
//   10-19 = neural tissue (grey, white, meninges)
//   20-29 = bone/cartilage
//   30-39 = muscle
//   40+ = user-defined
// ---------------------------------------------------------------------------

struct MaterialProperties {
  // Optical (for photon transport module)
  real_t absorption = 0.0;    // mu_a (mm^-1)
  real_t scattering = 0.0;    // mu_s' reduced scattering (mm^-1)
  real_t anisotropy = 0.9;    // g factor

  // Mechanical
  real_t stiffness = 1.0;     // kPa
  real_t density = 1.0;       // g/cm^3

  // Thermal
  real_t conductivity = 0.5;  // W/(m*K)
  real_t specific_heat = 3500; // J/(kg*K)

  // Biological flags
  bool vascularized = true;
  bool innervated = true;

  // Derived: effective photon diffusion coefficient
  // D = 1 / (3 * (mu_a + mu_s'))
  real_t PhotonDiffusion() const {
    real_t mu_t = absorption + scattering;
    return (mu_t > 1e-10) ? 1.0 / (3.0 * mu_t) : 10.0;
  }

  // Derived: penetration depth (1/e attenuation)
  real_t PenetrationDepth() const {
    real_t mu_eff = std::sqrt(3.0 * absorption * (absorption + scattering));
    return (mu_eff > 1e-10) ? 1.0 / mu_eff : 100.0;
  }
};

// Material ID constants for built-in tissue types.
namespace material {
  constexpr uint8_t kVoid = 0;
  constexpr uint8_t kSkinEpidermis = 1;
  constexpr uint8_t kSkinDermis = 2;
  constexpr uint8_t kSkinHypodermis = 3;
  constexpr uint8_t kBrainGrey = 10;
  constexpr uint8_t kBrainWhite = 11;
  constexpr uint8_t kBoneCortical = 20;
  constexpr uint8_t kCartilage = 21;
  constexpr uint8_t kMuscle = 30;
  constexpr uint8_t kMaxMaterials = 64;
}  // namespace material

// ---------------------------------------------------------------------------
// MaterialRegistry: maps material IDs to property structs.
// Initialized from materials/*.toml at simulation startup.
// Compact array storage for cache-friendly per-voxel lookup.
// ---------------------------------------------------------------------------
class MaterialRegistry {
 public:
  MaterialRegistry() {
    // Default: void material at index 0
    props_[0] = MaterialProperties{};
  }

  // Register a material by ID.
  void Register(uint8_t id, const std::string& name,
                const MaterialProperties& props) {
    props_[id] = props;
    names_[id] = name;
  }

  // Fast lookup by ID (array indexed, no hash).
  const MaterialProperties& Get(uint8_t id) const {
    return props_[id];
  }

  const std::string& Name(uint8_t id) const {
    static const std::string empty;
    auto it = names_.find(id);
    return (it != names_.end()) ? it->second : empty;
  }

  // Register default skin materials (used when no custom materials loaded).
  void RegisterSkinDefaults() {
    MaterialProperties epi;
    epi.absorption = 0.4;
    epi.scattering = 2.0;
    epi.anisotropy = 0.8;
    epi.stiffness = 1.0;
    epi.density = 1.1;
    epi.conductivity = 0.21;
    epi.specific_heat = 3590;
    epi.vascularized = false;
    epi.innervated = true;
    Register(material::kSkinEpidermis, "skin_epidermis", epi);

    MaterialProperties derm;
    derm.absorption = 0.2;
    derm.scattering = 1.5;
    derm.anisotropy = 0.85;
    derm.stiffness = 5.0;
    derm.density = 1.1;
    derm.conductivity = 0.37;
    derm.specific_heat = 3300;
    derm.vascularized = true;
    derm.innervated = true;
    Register(material::kSkinDermis, "skin_dermis", derm);

    MaterialProperties hypo;
    hypo.absorption = 0.15;
    hypo.scattering = 1.0;
    hypo.anisotropy = 0.9;
    hypo.stiffness = 0.5;
    hypo.density = 0.95;
    hypo.conductivity = 0.20;
    hypo.specific_heat = 2500;
    hypo.vascularized = true;
    hypo.innervated = true;
    Register(material::kSkinHypodermis, "skin_hypodermis", hypo);
  }

  // Register default brain materials.
  void RegisterBrainDefaults() {
    MaterialProperties grey;
    grey.absorption = 0.1;
    grey.scattering = 1.0;
    grey.anisotropy = 0.9;
    grey.stiffness = 0.5;
    grey.density = 1.04;
    grey.conductivity = 0.51;
    grey.specific_heat = 3630;
    grey.vascularized = true;
    grey.innervated = false;
    Register(material::kBrainGrey, "brain_grey", grey);

    MaterialProperties white;
    white.absorption = 0.08;
    white.scattering = 4.0;
    white.anisotropy = 0.85;
    white.stiffness = 1.0;
    white.density = 1.04;
    white.conductivity = 0.50;
    white.specific_heat = 3600;
    white.vascularized = true;
    white.innervated = false;
    Register(material::kBrainWhite, "brain_white", white);
  }

  // Register default joint materials (bone + cartilage for RA studies).
  void RegisterJointDefaults() {
    MaterialProperties bone;
    bone.absorption = 0.02;
    bone.scattering = 2.5;
    bone.anisotropy = 0.9;
    bone.stiffness = 15000.0;  // kPa, cortical bone
    bone.density = 1.9;
    bone.conductivity = 0.32;
    bone.specific_heat = 1260;
    bone.vascularized = true;
    bone.innervated = true;
    Register(material::kBoneCortical, "bone_cortical", bone);

    MaterialProperties cart;
    cart.absorption = 0.05;
    cart.scattering = 1.2;
    cart.anisotropy = 0.85;
    cart.stiffness = 800.0;    // kPa, articular cartilage
    cart.density = 1.1;
    cart.conductivity = 0.21;
    cart.specific_heat = 3400;
    cart.vascularized = false;  // avascular tissue
    cart.innervated = false;    // aneural
    Register(material::kCartilage, "cartilage", cart);
  }

 private:
  std::array<MaterialProperties, material::kMaxMaterials> props_{};
  std::unordered_map<uint8_t, std::string> names_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // MATERIAL_H_
