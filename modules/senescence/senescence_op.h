#ifndef SENESCENCE_OP_H_
#define SENESCENCE_OP_H_

// Senescence accumulation and SASP output are handled in fused_post.h
// (single-pass epidermal wound voxel loop). This header exists for
// organizational completeness; no standalone operation needed.
//
// Accumulation sources (fused_post.h):
//   1. Wound-induced DNA damage (basal rate in wound voxels)
//   2. Inflammation-driven (ROS, NF-kB activation)
//   3. AGE-driven (diabetic mode: glycative stress)
//
// SASP outputs (fused_post.h):
//   - Pro-inflammatory cytokines -> kInflammation field
//   - MMP-3/MMP-9 -> kProMMP field (zymogen)
//   - TGF-beta1 -> kTGFBeta field (fibrotic)
//
// Demaria et al. 2014 (doi:10.1016/j.devcel.2014.11.012)
// Coppe et al. 2008 (doi:10.1371/journal.pbio.0060301)

#endif  // SENESCENCE_OP_H_
