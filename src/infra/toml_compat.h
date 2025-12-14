// Bundled toml++ compatibility layer for BDM version independence.
// Provides TomlConfig typedef and config macros using toml++ API,
// bypassing BDM's config chain (which differs between v1.04 and v1.05).

#ifndef TOML_COMPAT_H_
#define TOML_COMPAT_H_

#include "tomlplusplus/toml.hpp"

namespace bdm {
namespace skibidy {

using TomlConfig = toml::table;

}  // namespace skibidy
}  // namespace bdm

// Config macros using toml++ API (independent of BDM version).
// These intentionally shadow any BDM-provided macros of the same name.

#ifdef BDM_ASSIGN_CONFIG_VALUE
#undef BDM_ASSIGN_CONFIG_VALUE
#endif
#define BDM_ASSIGN_CONFIG_VALUE(variable, config_key)           \
  {                                                             \
    auto bdm_toml_val_ =                                        \
        config.at_path(config_key).value<decltype(variable)>(); \
    if (bdm_toml_val_) {                                        \
      variable = *bdm_toml_val_;                                \
    }                                                           \
  }

#define BDM_ASSIGN_HOURS(var, key)                               \
  {                                                              \
    auto v_ = config.at_path(key).value<double>();               \
    if (v_) var = static_cast<int>(*v_ * 10.0);                  \
  }

#define BDM_ASSIGN_DAYS(var, key)                                \
  {                                                              \
    auto v_ = config.at_path(key).value<double>();               \
    if (v_) var = static_cast<int>(*v_ * 240.0);                 \
  }

#endif  // TOML_COMPAT_H_
