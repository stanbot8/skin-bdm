#ifndef TIME_MODEL_H_
#define TIME_MODEL_H_

#include "biodynamo.h"

namespace bdm {
namespace skibidy {

/// Time conversion helpers.  All internal bookkeeping uses integer steps;
/// TOML configs express durations in hours (or days for simulation length).
struct TimeModel {
  static constexpr real_t kDt = 0.1;  // hours per step

  static int HoursToSteps(double h) { return static_cast<int>(h / kDt); }
  static int DaysToSteps(double d) { return static_cast<int>(d * 24.0 / kDt); }
  static real_t StepsToHours(int s) { return s * kDt; }
  static real_t StepsToDays(int s) { return s * kDt / 24.0; }
};

}  // namespace skibidy
}  // namespace bdm

#endif  // TIME_MODEL_H_
