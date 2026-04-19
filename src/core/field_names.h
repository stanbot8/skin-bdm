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
constexpr const char* kTemperature = "Temperature";
constexpr const char* kGlucose = "Glucose";
constexpr const char* kLactate = "Lactate";
constexpr const char* kNitricOxide = "NitricOxide";
constexpr const char* kBasalDensity = "BasalDensity";
constexpr const char* kTIMP = "TIMP";
constexpr const char* kProMMP = "ProMMP";
constexpr const char* kAGE = "AGE";
constexpr const char* kSenescence = "Senescence";
constexpr const char* kNerve = "Nerve";
constexpr const char* kROS = "ROS";
constexpr const char* kStiffness = "Stiffness";
constexpr const char* kLymphatic = "Lymphatic";
constexpr const char* kEdema = "Edema";
constexpr const char* kVoltage = "Voltage";
constexpr const char* kTNFAlpha = "TNFAlpha";
constexpr const char* kCartilage = "Cartilage";
constexpr const char* kIL6 = "IL6";
constexpr const char* kSynovialFluid = "SynovialFluid";
constexpr const char* kTCellDensity = "TCellDensity";
constexpr const char* kBone = "Bone";
constexpr const char* kScab = "Scab";
constexpr const char* kFluence = "Fluence";
constexpr const char* kOpsin = "Opsin";
constexpr const char* kScarMaturity = "ScarMaturity";

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
  kTemperatureId = 23,
  kGlucoseId = 24,
  kLactateId = 25,
  kNitricOxideId = 26,
  kBasalDensityId = 27,
  kTIMPId = 28,
  kProMMPId = 29,
  kAGEId = 30,
  kSenescenceId = 31,
  kNerveId = 32,
  kROSId = 33,
  kStiffnessId = 34,
  kLymphaticId = 35,
  kEdemaId = 36,
  kVoltageId = 37,
  kTNFAlphaId = 38,
  kCartilageId = 39,
  kIL6Id = 40,
  kSynovialFluidId = 41,
  kTCellDensityId = 42,
  kBoneId = 43,
  kScabId = 44,
  kFluenceId = 45,
  kOpsinId = 46,
  kScarMaturityId = 47,
};

}  // namespace fields
}  // namespace skibidy
}  // namespace bdm

#endif  // FIELD_NAMES_H_
