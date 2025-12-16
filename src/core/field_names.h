#ifndef FIELD_NAMES_H_
#define FIELD_NAMES_H_

namespace bdm {
namespace skibidy {
namespace fields {

constexpr const char* kCalcium = "Calcium";
constexpr const char* kKGF = "KGF";
constexpr const char* kOxygen = "Oxygen";
constexpr const char* kStratum = "Stratum";
constexpr const char* kScar = "Scar";
constexpr const char* kWater = "Water";
constexpr const char* kInflammation = "Inflammation";
constexpr const char* kProInflammatory = "ProInflammatory";
constexpr const char* kAntiInflammatory = "AntiInflammatory";
constexpr const char* kTGFBeta = "TGFBeta";
constexpr const char* kCollagen = "Collagen";
constexpr const char* kTumor = "Tumor";
constexpr const char* kVascular = "Vascular";
constexpr const char* kBiofilm = "Biofilm";
constexpr const char* kVEGF = "VEGF";
constexpr const char* kMMP = "MMP";
constexpr const char* kFibronectin = "Fibronectin";
constexpr const char* kElastin = "Elastin";
constexpr const char* kHyaluronan = "Hyaluronan";
constexpr const char* kDermis = "Dermis";
constexpr const char* kImmunePressure = "ImmunePressure";
constexpr const char* kPH = "pH";
constexpr const char* kFibrin = "Fibrin";

// Derived composite field names (not DiffusionGrid IDs).
constexpr const char* kECMQuality = "ECMQuality";
constexpr const char* kTissueViability = "TissueViability";
constexpr const char* kWoundMicroenv = "WoundMicroenv";

// Substance IDs for BioDynaMo DiffusionGrid registration.
enum : int {
  kCalciumId = 0,
  kKGFId = 1,
  kOxygenId = 2,
  kStratumId = 3,
  kScarId = 4,
  kWaterId = 5,
  kInflammationId = 6,
  kProInflammatoryId = 7,
  kAntiInflammatoryId = 8,
  kTGFBetaId = 9,
  kCollagenId = 10,
  kTumorId = 11,
  kVascularId = 12,
  kBiofilmId = 13,
  kVEGFId = 14,
  kMMPId = 15,
  kFibronectinId = 16,
  kElastinId = 17,
  kHyaluronanId = 18,
  kDermisId = 19,
  kImmunePressureId = 20,
  kPHId = 21,
  kFibrinId = 22,
};

}  // namespace fields
}  // namespace skibidy
}  // namespace bdm

#endif  // FIELD_NAMES_H_
