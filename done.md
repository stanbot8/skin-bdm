# Wiggum progress log

## Session 2026-04-16 (manual, pre-loop)
- Infrastructure cleanup: wound state, coarse map, dt bug fix, fused-op dedup (commit ddd882e)
- Round 2 infra: registration helpers, Q10, WoundGate, test bug fix (a2ab55e)
- Fix dangling GridContext crash (fdddc3c)
- 8 smoke tests for untested modules (e2080ee)
- PASS/FAIL validation gate with exit code (3ab78e3)
- Misc cleanup: deprecated API, --quick mode, ProMMP registration (6b70715)
- Decompose angiogenesis/pressure/photon source hooks (09a70f1)
- 4 mechanistic improvements for RMSE (closure apoptosis, TGF-beta sequestration, FN replacement, neutrophil degranulation) (83b4123)

## Wiggum loop started 2026-04-17
- Iteration 1: Kicked off 3-replicate validation (seeds 42/43/44) in WSL background. Cleaned 2 verbose static_cast<real_t>(1) patterns in diabetic/baseline_inflammation.h and fibroblast/post_hook.h. Attempted to add mechanism_test.cc with 3 integration tests for the 4 new mechanisms, but they hit the OperationRegistry singleton limitation (ops persist across test runs) and were removed. Next iteration: analyze replicate results.
- Iteration 2: Consolidated 3-include block (grid_registry+voxel_snapshot+signal_board) into src/core/hook_api.h across 35 hook files. Fixed original replicate run failure (bash $seed not expanding) by rewriting as scripts/wiggum_replicate.py. Re-kicked off 3-replicate validation in background (job bzg17m1lh). 9 tests pass. Commit dc728a4.
- Iteration 3: Added scripts/wiggum_compare.py for multi-replicate RMSE analysis. 3-replicate validation (seeds 42/43/44) WITH new mechanisms shows 4 consistent failures: macrophages 17.57% (stddev 1.52 = real), fibronectin 15.99% (stddev 1.32), myofibroblasts 18.11% (stddev 6.73 = seed-dependent), tgfb 15.05% (stddev 2.34). Macrophage trajectory analysis: day 1 sim 0% vs ref 40% (early arrival missing), day 10 sim 68% vs ref 28% (late rebound, needs M2 emigration). Kicked off baseline (mechanisms OFF) 3-replicate in background for delta analysis. Commit 75acc47.
- Iteration 4: Computed WITH vs WITHOUT delta across 3 replicates. 6 of 12 observables regressed (myofib +5.00%, mmp +2.48%, tgfb +1.89%, ph +1.09%, mac +0.94%, fibro +0.47%). Only closure (-0.33%) and neutrophils (-1.23%) improved. Disabled all 4 mechanisms by default (zeroed param values) while keeping code paths intact. Users can opt in via config. scripts/wiggum_delta.py preserved for future calibration sweeps. Commit afe26a5.
- Iteration 5: Corrected stale README claims: 130 -> 540 tests (verified via --gtest_list_tests), 41 -> 47 diffusion grids (verified via field_names.h). Kicked off diabetic-wound sim in background for next iteration's validation. Commit f9f4115.
- Iteration 6: Found + fixed bogus pH validation for non-normal conditions. Diabetic sim previously reported pH RMSE = 49.07% (seed 42) because validate_ph always used the normal-wound Schneider 2007 reference; diabetic/burn/pressure/surgical have different pH kinetics with no reference data. Gated pH validation to condition == "normal". Diabetic study now reports 6 real failures instead of 7 bogus. Commit 0928670.
- Iteration 7: Extended iteration 6's bug fix. validate_fibroblast had no condition awareness; validate_microenvironment only handled condition for TGF-b and MMP but not VEGF/FN. Both silently compared non-normal runs against normal-wound references. Gated fibroblast validation to condition=="normal" (no diabetic/burn/pressure/surgical fibroblast refs exist) and VEGF/FN to condition=="normal". print_summary and threshold gate now skip None RMSE values. Diabetic seed-42 fails down from 7 bogus to 5 real. Commit e324a28.
- Iteration 8 (final): Full test suite verification (29/29 pass across 9 suites). Audit found no further bugs: all modules have READMEs, no unused params, no column-name mismatches between metrics.csv and validators, no stale magic-number duplicates. Sim wall time 37s matches paper's 35s target. Remaining work is user-gated (parameter tuning they forbid), externally-gated (BDM OperationRegistry upstream fix), or speculative (new biology modules). Mission complete. Loop 18d1cd6d cancelled.
