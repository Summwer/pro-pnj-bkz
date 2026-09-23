# Supplementary materials

This package accompanies **Optimizing Two-Step Primal LWE Attacks with PnJBKZ Simulation**. Figure and table numbers below follow the current `main.tex`, including its appendices. Figures and tables have separate numbering sequences. Entries are ordered by their first appearance in the paper; a shared folder can support several later items.

## Material index

| Paper item | Material in this package | Description |
| --- | --- | --- |
| Table 2 | [Table_02_NIST_Estimates](Table_02_NIST_Estimates/) | Concrete primal-attack estimates for the selected NIST parameter sets. |
| Table 3; Tables 6–8 | [Tables_03_06-08_Slope_and_Cost](Tables_03_06-08_Slope_and_Cost/) | Simulated and observed reduction slopes and costs. |
| Table 4 | [40-035.log](40-035.log) | One supplied execution log for the LWE Challenge records; this file is not the complete set of record logs. |
| Table 5; Figures 9 and 29 | [Figures_09_29_Table_05_Strategy_Comparison](Figures_09_29_Table_05_Strategy_Comparison/) | Solver comparisons, strategy-generation information, and execution logs. |
| Figure 10 | [Figure_10_Reduction_Speedup](Figure_10_Reduction_Speedup/) | Reduction-step speedup on the 222-dimensional lattice for `(n, alpha) = (60, 0.010)`. |
| Figure 12 | [Figure_12_Strategy_Acceleration](Figure_12_Strategy_Acceleration/) | Estimated acceleration from the optimized strategy. |
| Figure 14 | External experiment entry point described below | Failure probability of the estimated terminal SVP dimension; no separate raw-data folder for this experiment is included here. |
| Figures 15–19 | [Figures_15-19_Practical_Cost_Model](Figures_15-19_Practical_Cost_Model/) | Pump and PnJBKZ practical cost measurements, fitting code, and plots. |
| Figures 21–28 | [Figures_21-28_Simulator_Verification](Figures_21-28_Simulator_Verification/) | PnJBKZ simulator verification code and GS-profile data. |
| Table 9 | External experiment entry point described below | Generation of strategies for the LWE Challenge records; no dedicated Table 9 folder is included here. |
| Additional material, without a current table number | [Additional_Sample_Size_Optimization](Additional_Sample_Size_Optimization/) | Sample-size optimization workbook retained from earlier experiments. |
| Additional simulator data | [simulator-test/gs-lengths-simulator_backup](simulator-test/gs-lengths-simulator_backup/) | Additional stored profiles, ratios, and plots. |

Items not listed in this index do not have a separately numbered data directory in this package. In particular, the folder for Figure 12 contains acceleration estimates, not the failure-probability experiment for Figure 14.

## Software and experimental environment

The research implementations are [pro-pnj-bkz](https://github.com/Summwer/pro-pnj-bkz) and [lwe-estimator-with-PnJBKZ](https://github.com/Summwer/lwe-estimator-with-PnJBKZ). Commands referring to scripts not included here should be run in the corresponding external repository after following its build instructions. The implementation option `enumbs` corresponds to the strategy-search implementation called PSSearch in the paper; older logs and plotting code also use the name EnumBS.

For example, in the configured `pro-pnj-bkz` repository:

```bash
source ./activate
python ProPnJBKZ_for_lwe.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --float_type "dd" --strategy_method "enumbs" --load_lwe "lwe_challenge" --max_RAM 43
```

This example selects `(n, alpha) = (40, 0.025)`, 32 CPU threads, two GPUs, `dd` floating-point arithmetic, and a memory bound of 43 in log2(bits). The recorded strategy choices include `enumbs`, `bssav1`, and `bssav2`. To generate a random LWE instance instead of loading an LWE Challenge instance, use `--load_lwe "lwe_instance"`.

The practical cost measurements supplied with the original experiment package use Intel Xeon 5128 CPUs (16 cores / 32 threads at 2.3 GHz), 1.48 TB RAM, two NVIDIA GeForce RTX 3090 GPUs, Ubuntu 20.04, GCC 8, CUDA 11.4, and NVIDIA driver 470.199.02. The paper specifies the machine for each benchmark; the record table also contains results from other machines.

For a different machine, measure its Pump and PnJBKZ costs and update `refined_practical_pump_cost_dd` and `refined_PnJBKZ_cost_model_dd` in the implementation's `cost.cpp`. The implementation also supports its theoretical model through `cost_model = 1`.

## Table 2 — NIST estimates

The logs are in [Table_02_NIST_Estimates/nist-est-log](Table_02_NIST_Estimates/nist-est-log/).

In the external `lwe-estimator-with-PnJBKZ` repository, build the C++ implementation using its `rebuild.sh`, then run the NIST experiment entry point:

```bash
./implement_all_NIST_schemes.sh
```

The experiment records the strategy and estimated cost for each parameter set, including the MATZOV22 list-decoding cost model. The original experiment workflow uses `enumbs(cumprob+prob)+list_decoding[MATZOV22].log` in `nist-round3-est-result`.

The comparison labeled Previous uses the cumulative-cost adaptation of [leaky-LWE-Estimator](https://github.com/lducas/leaky-LWE-Estimator), with the original experiment script:

```bash
sage NIST-pro-bkz-avgbeta-matzov.sage
```

For the Two-step baseline, use the external repository's `sage/NIST-round3` experiment in Sage:

```python
load("NIST-two-step-matzov.sage")
```

## Table 3 and Tables 6–8 — Reduction slopes and costs

The workbook is [simulate_slope_and_cost.xlsx](Tables_03_06-08_Slope_and_Cost/simulate_slope_and_cost.xlsx), sheet `simulate-slope-cost-lwechal`.

| Paper table | LWE Challenge `(n, alpha)` | Workbook rows |
| --- | --- | --- |
| 3 | `(40, 0.030)` | 6–11 |
| 6 | `(45, 0.020)` | 30–34 |
| 7 | `(50, 0.015)` | 42–47 |
| 8 | `(40, 0.025)` | 18–21 |

The Practical columns use only the five measured trials. The simulated columns are excluded from these averages. The time entries describe individual BKZ calls, rather than elapsed time accumulated over earlier calls. The reported log-time column retains the workbook's convention of averaging the five log2(time) measurements. With the four-decimal slope values printed in the paper, the relative slope error over these 21 calls is at most 7.26%.

The corresponding execution logs are under [test-result/lwechal-test](Figures_09_29_Table_05_Strategy_Comparison/test-result/lwechal-test/), particularly `lwechal-test-new-cost-modelv2/d4f-default-g6k/enumbs(32+2gpus)`. The three appendix comparisons use `dd` precision. The original external entry point for obtaining the simulation output is `strategy_simulation.py`.

## Table 4 — LWE Challenge records

The file [40-035.log](40-035.log) supplies one execution record. Refer to the paper for the full record table and its machine settings. The command for generating the record strategies is described under Table 9 below.

## Table 5 and Figures 9 and 29 — Strategy and solver comparisons

[Figures_09_29_Table_05_Strategy_Comparison](Figures_09_29_Table_05_Strategy_Comparison/) contains the cost-comparison workbooks, `draw_strategy_comparison.py`, and the `test-result` directory. Figures 9 and 29 use the same solver-comparison plots, including the jump-parameter ablations. Table 5 reports selected strategies and their generation costs from these experiments.

In the external solver repository, the original experiment entry points are `implement_lwechal_forall.sh` for LWE Challenges and `implement_lwe_instance_forall.sh` for randomly generated LWE instances. The recorded tests use `dd` precision. The local logs are organized under `test-result/lwechal-test` and `test-result/lwe-instance-test`.

To regenerate the comparison plots from the numerical arrays supplied in the plotting script, run it from its directory:

```bash
cd Figures_09_29_Table_05_Strategy_Comparison
python draw_strategy_comparison.py
```

The script requires NumPy, Matplotlib, SciPy, and adjustText.

## Figure 10 — Speedup in the reduction step

[Figure_10_Reduction_Speedup/lwechal-concrete-cost.xlsx](Figure_10_Reduction_Speedup/lwechal-concrete-cost.xlsx) contains the comparison for `(n, alpha) = (60, 0.010)` and lattice dimension 222. The original simulation entry point is `cpp/strategy_simulation.cpp` in the external estimator repository.

## Figure 12 — Acceleration from the optimized strategy

[Figure_12_Strategy_Acceleration/lwechal-practical-cost-est.xlsx](Figure_12_Strategy_Acceleration/lwechal-practical-cost-est.xlsx) contains the estimated costs and differences for the baseline and optimized strategies. These estimates use the PnJBKZ simulator and practical time-cost model. The original simulation entry point is also `cpp/strategy_simulation.cpp`.

## Figure 14 — Failure probability of SVP-dimension estimation

The experiment generates 100 random LWE instances for each tested `(n, alpha)` and compares the terminal-dimension estimates. Its original entry point, in the external `pro-pnj-bkz` repository, is:

```bash
python PumpDimEst_comparison.py
```

This script and its raw trials are not included in this supplementary package. The figure's two panels concern `alpha = 0.010` and `alpha = 0.015`; this experiment is separate from the Figure 12 workbook.

## Figures 15–19 — Practical cost models

The measurements and fitting scripts are in [Figures_15-19_Practical_Cost_Model](Figures_15-19_Practical_Cost_Model/).

| Paper figures | Local subdirectory | Relevant scripts |
| --- | --- | --- |
| 15, 18, 19 | `Practical-PnjBKZ-Cost-Model` | `draw_pump_cost_in_each_pnjBKZ.py`, `pnjbkz_cost_model.py` |
| 16, 17 | `Practical-Pump-Cost-Model` | `practical_pump_cost.py` |

Run each local script from its own subdirectory because it reads or writes relative paths. To collect new measurements, use the external solver's experiment entry points:

```bash
python practical_pump_cost_test.py | tee practical_pump_cost.log
python practical_PnJBKZ_cost_test.py
```

The latter records per-Pump costs in `each_pump_cost_in_pnjbkz`. Use those measurements with the fitting scripts to calibrate the target machine's cost model.

## Figures 21–28 — PnJBKZ simulator verification

The local script is [Figures_21-28_PnJBKZ_Simulator_Verification.py](Figures_21-28_Simulator_Verification/PnJBKZ_Simulator/Figures_21-28_PnJBKZ_Simulator_Verification.py). Run it with `PnJBKZ_Simulator` as the working directory:

```bash
cd Figures_21-28_Simulator_Verification/PnJBKZ_Simulator
python Figures_21-28_PnJBKZ_Simulator_Verification.py
```

It requires NumPy and Matplotlib. It reads the relative directory `gs-lengths-simulator` and writes profile plots, Real/Sim GS-norm ratios, and ratio plots into the corresponding experiment directory. The default configuration is `(n, alpha) = (70, 0.005)`, dimension 235, `beta = 90`, `J = 12`, 13 tours, and 20 trials. To reproduce the other parameter settings in Figures 21–28, select their settings in the script's `__main__` configuration block. One invocation processes the selected configuration, not every figure automatically.

| Paper figures | Experiment |
| --- | --- |
| 21 | Ideal-Pump verification, `(n, alpha) = (60, 0.010)`, 10 trials. |
| 22–23 | Practical settings, `(60, 0.010)`, 20 trials. |
| 24–25 | Practical settings, `(75, 0.005)`, 20 trials. |
| 26–28 | Practical settings, `(70, 0.005)`, 20 trials. |

The plotted ratio is `Real(||b_i*||) / Sim(||b_i*||)`. Observed profiles are averaged in logarithmic coordinates, so the observed norm in the ratio is a geometric mean. The practical GS-coordinate ratio statistics are independent of the slope/time workbook for Tables 3 and 6–8.

Additional archived profiles and plots are supplied in [simulator-test/gs-lengths-simulator_backup](simulator-test/gs-lengths-simulator_backup/). The local script does not automatically read that separate archive.

To collect fresh profiles, the original external solver workflow uses `lwe_challenge_gen_rr.py`, for example:

```bash
python lwe_challenge_gen_rr.py 75 --lwe/alpha 0.005 --bkz/jump 9 --pump/down_sieve True --bkz/blocksizes "[95,95,95,95,95,95,95,95,95,95,95,95]" --gpus 2 --threads 32
```

Match the block sizes, jump, tour count, and number of repeated trials to the intended experiment before using the verification script.

## Table 9 — Strategies for the LWE Challenge records

In the external implementation, the original generation command is:

```bash
./implement_unsolved_lwechal.sh
```

That workflow stores its results in `strategy-for-unsolved-lwechal`. This is an external-repository path, not a dedicated directory included in the present package.

## Additional sample-size optimization data

[Additional_Sample_Size_Optimization/optimize-m-estimation.xlsx](Additional_Sample_Size_Optimization/optimize-m-estimation.xlsx) retains the sample-size optimization experiments. The corresponding table from an earlier draft is not present in the current paper, so this directory has no figure or table number.

The original external entry point is `./optimize_m.sh`, with results stored under `optimize-m` in the implementation repository.
