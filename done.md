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
