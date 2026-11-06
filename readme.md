# Multi-EEG analysis code

MATLAB and R code for Gheza et al. (submitted to Nature Communications). Figure numbers below follow the submitted manuscript and SI.

The data aren't in this folder, and most scripts use absolute paths from my machines, so you'll need to adapt those before running anything. Plotting scripts need the `.mat` files from the earlier steps.

## Which script makes which figure

| Figure | Script | Notes |
|---|---|---|
| 1D-G, 2A-J, 3A | `Multi_EEG__behavioral_analyses.rmd` | Behavior |
| 3B-E | `dynamic_LCA_sims.m` and plotting helpers | Model simulations |
| 4B | `RSA_distances_classprob_SIDEnCue_plot_classificationaccuracy.m`, then `Multi_EEG__plots.Rmd` | Decoding accuracy, from the output of `RSA_distances_..._HPC.m` |
| 4C-D, F | `plot_RSA_betas.m`, then `Multi_EEG__plots.Rmd` | Model 1. The RT distribution in 4D comes from the behavioral notebook |
| 4E | `plot_RT_regress_RSAm1.m` | Fits from `RT_regress_RSAxCNG_m1.m` |
| 5B-C | `plot_RSA_betas_m3.m`, then `Multi_EEG__plots.Rmd` | Model 3 |
| 5D | `plot_RT_regress_RSAm3orEncoding.m` | Fits from `RT_regress_RSAm3xCNG_or_EncodingxCNG.m` with `RSA = 1` |
| 5E-G, S2 | `plot_RSA_betas_CSE_m3.m`, then `Multi_EEG__plots.Rmd` | CSE |
| 6A-E | `plot_RSA_betas_m3.m`, then `Multi_EEG__plots.Rmd` | Split by proportion congruency and distractor dimension |
| 6F | `Multi_EEG__behavioral_analyses.rmd` | Switch-trial TSPC |
| 7A-C, S5 (all panels) | `plot_conflict_geometry_panels.m` | |
| S1 | `aa6_tf_plot_cse.m`, `aa7_tf_RTs_regression_plot_2025.m` | Response-locked TF |
| S3-S4 | `plot_encoding_fullmodel_CSE.m`, then `Multi_EEG__plots.Rmd` | |

## General notes

- The dataset is experiment version 3: participants 22 and 24-48 (26 in total). Participant 23 is excluded in the analysis scripts.
- Model 3 in the scripts (`m3`) is model 2 in the paper. The `m2` scripts and exports are an unreported model.
- `onset = true` runs stimulus-locked analyses and `onset = false` response-locked ones. Check this setting before running a script.
- Files ending in `_HPC`, and functions with `local = 0`, expect cluster paths (`/scratch`, `/mnt/beegfs`).

## MATLAB exports read by R

`Multi_EEG__plots.Rmd` loads exports by file name with `read_export()`. Names start with the figure they feed (e.g. `Fig4_m1_RSA_betas.mat`). `_bst0` means no bootstrapping, and `notreported_` means the export isn't in the paper.

## Pipeline

All MATLAB scripts are in `EEG_analyses__MATLAB/`.

### 1. Preprocessing

Run in this order:

1. `pp1_preprocessing_filter_ICA.m`: import the BrainVision EEG, filter, estimate ICA weights.
2. `pp2_epoching_ICApruning.m`: epoch (stimulus- or response-locked), remove ICA components, baseline-correct.
3. `pp3_trial_rej.m`: reject trials, rereference, save the cleaned datasets and logs.
4. `decoding_prepare_sessions.m`: prepare the V3 session and recover task/stimulus markers. Saves the behavior for correct, unrejected trials and the task-sorted EEG (`sessions_data_[response_]sbj_N.mat`).

`FWHM_report.m` just reports the temporal FWHM of the wavelets.

### 2. Time-frequency power and decoding

1. `decoding_tf_on_sessions.m`: Morlet decomposition. Saves baseline-normalized single-trial band power (`sessions_data_[response_]tf_decoding_BN_sbj_N.mat`), which the RSA decoding and encoding scripts read. Run stimulus-locked first, because response-locked uses the same pre-stimulus baseline.
2. `decoding_tf_on_sessions_singletrial.m`: TF power with single-trial baselines, for section 6.
3. `RSA_distances_classprob_kfold_SIDEnCue_multiclass_repCV_HPC.m`: repeated cross-validated multiclass decoding of side and task. Its class probabilities are the input to RSA. It saves `distances_mltclass_sbj_N.mat`, but for legacy reasons the scripts that read it locally expect a renamed copy, `RSA_pc/multiclass/singletrial/distances_mltclass_repCV_effHPT_bal_sbj_N.mat` (under `rsa_data/`), so copy and rename the files after running it.
4. `RSA_distances_classprob_SIDEnCue_plot_classificationaccuracy.m`: group decoding accuracy. Exports `Fig4_classification_accuracy.mat`. `pt` picks the inspection plots (1 accuracy, 2 probability over time, 3 probability structure).

### 3. RSA fits

- `fit_RSAmodels_multiclassRDM_sngltrl_m1__dim_gen_side__std.m`: model 1 (dimension-specific stimulus-response, rule, and dimension-general side).
- `fit_RSAmodels_multiclassRDM_model3_sngltrl__std.m`: model 3 (conjunctions for cued and non-cued dimensions).

Helpers:

- `cov1para.m`: Ledoit-Wolf covariance shrinkage (third-party, BSD; license comments kept).
- `cong_control_ST_precur_function.m`: bootstrap/decimation indices matched on current and previous congruency, for the CSE analyses.
- `cong_control_ST_function.m`, `cong_control_ST_singlecued_function.m`, `cong_control_ST_singledstr_function.m`: bootstrap indices that balance congruency within PC conditions (overall, per cued task, per distractor). `plot_RSA_betas.m` and `plot_RSA_betas_m3.m` use them when `bst = 1`, so they affect the averages exported for Figures 4-6.
- `RSA_permutest.m`, `CSE_permutest.m`: cluster tests for the MATLAB inspection plots only. The time-series clusters in the paper are computed in `Multi_EEG__plots.Rmd` with `permuco`.

### 4. RSA results

- `plot_RSA_betas.m`: model 1 time courses, by dimension and by previous congruency. Exports `Fig4_m1_RSA_betas.mat` and `Fig4_m1_RSA_betas_bydim_bst0.mat`.
- `plot_RSA_betas_m3.m`: model 3 time courses, including the PC and distractor splits. Exports `Fig5_Fig6_m3_RSA_betas.mat` and `Fig6_m3_RSA_betas_bydim.mat` (both `bst = 1`) and `Fig5_m3_RSA_betas_bydim_bst0.mat`.
- `plot_RSA_betas_CSE_m3.m`: parametric, within/across-dimension, by-dimension and bidirectional CSE. Exports `Fig5_SFig2_CSE_m3_RSA_betas.mat`.
- `plot_RSA_betas_CSE_m1.m`: the same CSE checks for model 1.

MATLAB mostly prepares arrays here. The multi-panel figures are put together in `Multi_EEG__plots.Rmd`.

### 5. Encoding model and geometry

Run the matching model before its plotting script:

- `encoding_fullmodel_repCV__no_tempgen_HPC.m`: trial-level, time-resolved full model.
- `encoding_fullmodel_repCV__tempgeneralization_HPC.m`: train-time by test-time version.
- `encoding_fullmodel_repCV__similarity_HPC.m`: similarity/geometry outputs.

Plotting:

- `plot_encoding_fullmodel_CSE.m`: stimulus- or response-locked projections and CSE. With `basic_saveout = 1` (bst = 0) it writes `SFig3_basic_..._bst0` or `SFig4_resplocked_..._bst0`. With `basic_saveout = 0` (bst = 1) it writes `SFig3_CSE_encoding_dimGen` or `_dimSpec`.
- `plot_conflict_geometry_panels.m`: temporal-generalization maps, covariance-weighted similarity matrices, orthogonality time courses. It also has some inspection panels that didn't make it into the paper. 2-D clusters come from `cluster_perm_test_TG_2D.m`.
- `plot_encoding_tempgeneralization.m`: quick inspection only. It saves nothing and no panel in the paper comes from it.
- `RT_regress_RSAm3xCNG_or_EncodingxCNG.m` with `RSA = 0` fits the encoding-by-congruency RT models, and `plot_RT_regress_RSAm3orEncoding.m` plots them.

### 6. Time-frequency analyses

These use the single-trial power from `decoding_tf_on_sessions_singletrial.m` (`data_[response_]tf_singletrial_NB_sbj_N.mat`). The `aa6`/`aa7` prefixes are just old names.

- `tf_RT_regression.m`: single-trial TF/RT regression, stimulus- or response-locked. Saves `rt_power_regression_[stim|resp]locked2025_slope.mat`.
- `aa7_tf_RTs_regression_plot_2025.m`: group plot of the response-locked regression.
- `aa6_tf_plot_cse.m`: TF power CSE plots, topographies and cluster summaries, stimulus- or response-locked.

### 7. Model simulations

The code is in `Model_simulations__MATLAB/dynamic_LCA_model/`.

Run the dynamic leaky competing accumulator model, with dimension-specific time constants, like this:

```matlab
output = dynamic_LCA_sims(160000,-5,4)
```

- `generate_task.m`, `init_model_Hopfield.m` and `logistic.m` are required helpers.
- `make_cse_plot_rt.m`, `make_cse_plot_acc.m`, `make_dyn_plot.m` and `make_heatmap.m` plot the outputs.
- `ddm_sims_distractorspecific_Hopfield.m` is the older DDM-like implementation from the NHB 2025 model. It's here for reference and isn't the model in this paper.

## R code

In `Behavioral_and_figures__R/`:

- `Multi_EEG__behavioral_analyses.rmd`: imports the V3 behavior, builds the congruency and PC variables, fits the Bayesian models, makes the behavioral panels.
- `distributional_analysis_functions.R`: RT-quantile, delta-plot and delta-delta helpers used by the behavioral notebook.
- `Multi_EEG__plots.Rmd`: loads the MATLAB exports, runs the time-series cluster tests, makes the EEG panels.

## Requirements

- MATLAB with the Statistics and Machine Learning Toolbox and the Image Processing Toolbox.
- EEGLAB, ERPLAB and the BVA-IO plugin for preprocessing.
- FieldTrip (`ft_prepare_neighbours`, used by `aa6_tf_plot_cse.m`).
- `permutest.m` (Gerber, MATLAB File Exchange, not included). Only `RSA_permutest.m` and `CSE_permutest.m` call it.
- The R packages loaded in the notebooks, including `permuco`.
- Not included: settings files (`preprocessing_settings.mat` or `preprocessing_settings_mac.mat`), behavioral files, cleaned EEG, intermediate decoding/RSA/encoding results, and cluster job files.
