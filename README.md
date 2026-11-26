# How reliable is the link between effort and nucleus accumbens connectivity at rest? A psychometric evaluation of repeated assessments

## Authors

* Zsofia Karlocai
* Anne Kühnel
* Corinna Schulz
* Melina Grahlow
* Kristin Kaduk
* Andreas Frick
* Birgit Derntl
* Nils B. Kroemer

## Overview

This repository contains the data and the R and MATLAB code to reproduce the statistics and figures of the manuscript. Participants completed an effort allocation task (EAT) in five behavioral sessions (S1–S5) and two task-free fMRI sessions (S6, S7). The analyses estimate participant-specific effort parameters with linear mixed-effects models, assess the test–retest reliability of these parameters and of nucleus accumbens (NAcc) functional connectivity (FC), relate the effort parameters to NAcc FC and test how reproducible these brain–behavior associations are across sessions, and relate the effort parameters to wanting ratings and questionnaire measures (SHAPS, LARS, BIS/BAS).

* Unthresholded statistical maps (NeuroVault): https://neurovault.org/collections/23461/
* Study preregistration (parent study): https://osf.io/gpr52
* fmreli toolbox (voxel/ROI ICC): https://github.com/nkroemer/reliability

### Session naming

|Manuscript|Content|
|-|-|
|S1|intake session, EAT|
|S2–S5|weekly sessions, EAT with certain and uncertain trials|
|S6, S7|task-free fMRI sessions|

\---

## Repository structure

```
code/      analysis scripts
data/      anonymised input data
output/    results and figures (created by the scripts)
```

### Data files

Participant IDs (`sub-001`, …) are random and consistent across all files. Only the variables used in the analyses are included.

* `eat_trials.csv` – trial-level EAT data, S1–S5 (certain trials only)
* `mri_sample.csv` – MRI sample (n = 84): SPM covariates, z-scored S1 effort coefficients, sex, and NAcc FC extracted from the right IFG and subcallosal cortex ROIs
* `questionnaires.csv` – SHAPS, LARS, BIS and BAS sum scores
* `session_timing.csv` – start time of each behavioral session and days to the MRI sessions (no dates)
* `roi_brain_behavior_partial_r.csv` – per atlas ROI: partial r between NAcc FC and the S1–S5 effort coefficients, ROI NAcc FC reliability (ICC), network
* `roi_S6_S7_tmap_similarity.csv` – per atlas ROI: correlation of the S6 and S7 T maps of the S1 effort coefficients
* `fingerprinting_similarity_matrix.csv` – correlation of gray-matter NAcc FC maps between all participants, S6 (rows) × S7 (columns)
* `group_fc_overlap_S6_S7.csv` – Dice and Jaccard overlap of the S6 and S7 group NAcc FC maps (p < .01 uncorrected; Figure 2B)

Raw and preprocessed MRI data are not shared.

### Scripts

* `code/analyses.R` – behavioral and brain–behavior analyses, in manuscript order:

  * session intervals (Methods)
  * linear mixed-effects models of relative effort and wanting in S1
  * participant-specific coefficients (z-scored for SPM)
  * test–retest reliability of the effort coefficients across S1–S5 (ICC, Pearson correlations; Figure 2A, Table S1)
  * NAcc FC and effort (Figure 4)
  * reproducibility of brain–behavior associations across behavioral sessions S1–S5 (Figure 5) and MRI sessions S6–S7 (Figure 6)
  * effort–wanting and effort–questionnaire correlations (Figures S4–S6)
  * supplementary analyses (Figure S1, Table S2, Figures S2, S3)
* `code/fingerprinting.m` – functional connectome fingerprinting: identification accuracy, permutation and binomial tests, within- vs between-participant similarity (Figure 3B, 3C)

### Output

Created in `output/`:

* `session_intervals.csv` – intervals between the last behavioral session and the MRI sessions
* `coef_relEff.csv`, `coef_wanting.csv` – participant-specific S1 coefficients of the effort and wanting models
* `coef_relEff_ZScore.csv`, `coef_wanting_ZScore.csv` – z-scored coefficients (SPM predictors)
* `Figure2A_behavioral_reliability.png` – session-by-session correlations and ICCs of the effort coefficients
* `Figure3B_fingerprinting_matrix.png`, `Figure3C_ECDF.png`, `Figure3C_mean_similarity.png`, `fingerprinting_results.mat`
* `Figure4A_intercept_IFG.png`, `Figure4B_difficulty_SubCal.png` – covariate-adjusted scatterplots
* `Figure5A_profile_similarity.png`, `Figure5B_focal_ROI_trajectories.png`
* `Figure6A_spatial_reproducibility_distribution.png`, `Figure6B_spatial_reproducibility_extreme_ROIs.png`
* `effort_wanting_correlations.csv`, `effort_questionnaire_correlations.csv`
* `TableS1_behavioral_reliability.csv`, `TableS2_timing_sensitivity.csv`
* `FigureS1_sex_brain_behavior.png`, `FigureS1_sex_interactions.csv`
* `FigureS2.png` – `FigureS6.png`
* `sessionInfo.txt` – R and package versions

\---

## Variables

`eat\\\_trials.csv`

* `ID` – participant
* `Session` – behavioral session (1–5)
* `Trial` – trial number
* `RewType` – reward type (0 = food, 1 = money)
* `RewMag` – reward magnitude (low / high)
* `Difficulty` – difficulty (S1: 75 % / 85 % of maximum frequency; S2–S5: continuous)
* `RelEffort` – mean button-press rate relative to the individual maximum
* `Wanting` – wanting rating (VAS 0–100)

`mri_sample.csv`

* `Female` – 1 = women, 0 = men
* `cAge`, `cSex`, `cBMI` – centered covariates of the SPM models
* `FD` – log mean framewise displacement
* `Effort_intercept`, `Effortrewmag`, `Effort_diff` – z-scored S1 effort coefficients (SPM predictors)
* `FC_rIFGoper`, `FC_subcallosal` – NAcc FC (Fisher z, mean of S6 and S7) in the right IFG pars opercularis and subcallosal cortex atlas ROIs

`questionnaires.csv`

* `SHAPS` – sum of 14 items scored 1–4 (rescaled to 0–42 in the script)
* `LARS`, `BIS`, `BAS` – sum scores

`roi_brain_behavior_partial_r.csv`

* `Intercept_S1` … `Difficulty_S5` – partial correlations controlling for age, sex, BMI, FD and the other two coefficients
* `ICC_con` – ICC(3,1) (consistency) of ROI NAcc FC across S6 and S7 (fmreli)

\---

## Analyses

### Variable centering

All task predictors are grand-mean centered before model fitting (prefix `c`): `cTrial`, `cDifficulty`, `cRewMag`, `cRewType`.

### Linear mixed-effects models

**Model 1: S1 behavioral models.** Task conditions predicting relative effort and wanting:

```r
RelEffort ~ cRewMag * cDifficulty * cRewType + cTrial +
  (1 + cDifficulty + cRewMag + cRewType + cTrial | ID)

Wanting ~ cRewMag * cDifficulty * cRewType + cTrial +
  (1 + cDifficulty + cRewMag + cRewType + cTrial | ID)
```

Participant-specific coefficients (fixed effect + random deviation, `coef()`) for the intercept, reward magnitude and difficulty are extracted and z-scored across participants before entry into SPM.

**Model 2: cross-session reliability (Figure 2A, Table S1).** The effort model is re-fitted separately in each session S1–S5, using the participants who completed all five sessions. The z-scored coefficients are compared across sessions with pairwise Pearson correlations and the ICC(3,1) (two-way, consistency; `irr`).

### Brain–behavior associations

* **Figure 4:** partial correlations between the S1 effort coefficients and NAcc FC in the right IFG and subcallosal cortex, controlling for age, sex, BMI, FD and the other two coefficients.
* **Figure 5:** ROI-wise partial correlations recomputed with the coefficients of each behavioral session S1–S5, and the similarity of the ROI-wise association profiles across sessions.
* **Figure 6:** within-ROI correlation of the S6 and S7 T maps of the S1 coefficients.
* **Figure S1:** sex × behavior interactions for the two focal associations.
* **Table S2:** the focal associations of Figure 5B, additionally controlling for session start time and the interval to the MRI sessions.

### Fingerprinting (Figure 3)

Gray-matter NAcc FC maps are correlated between all participants across S6 and S7. A participant is identified correctly if their S6 map correlates most strongly with their own S7 map. Accuracy is tested against chance with a permutation test (1,000 permutations) and a binomial test.

### Correlational analyses

The effort coefficients are correlated (Pearson) with the wanting coefficients and with SHAPS, LARS, BIS and BAS scores (Figures S4–S6).

\---

## Visualizations

`analyses.R` generates the figures with `ggplot2`, `cowplot`, `ggcorrplot` and `patchwork`; `fingerprinting.m` generates Figure 3 in MATLAB. Panels are saved separately and were assembled in BioRender.

### Figure 2A: Behavioral test–retest reliability

* Type: correlation heatmaps, one per effort parameter (intercept, reward magnitude, difficulty)
* Cells: Pearson correlations of the z-scored coefficients between sessions S1–S5
* Annotation: ICC(3,1) across S1–S5
* Output: `Figure2A_behavioral_reliability.png`

### Figure 3: Functional connectome fingerprinting (`fingerprinting.m`)

* Panel B: participant-by-participant similarity matrix (S6 × S7) with identification accuracy
* Panel C (left): empirical cumulative distributions of within- and between-participant similarity (Fisher z) with 95% CIs
* Panel C (right): mean within- vs between-participant similarity (± SEM)
* Output: `Figure3B_fingerprinting_matrix.png`, `Figure3C_ECDF.png`, `Figure3C_mean_similarity.png`

### Figure 4: NAcc FC and effort

* Type: covariate-adjusted scatterplots with linear fit and 95% CI
* Panel A: intercept vs NAcc–right IFG FC
* Panel B: difficulty vs NAcc–subcallosal FC
* Both axes adjusted for age, sex, BMI, FD and the other two coefficients; partial r and p annotated
* Output: `Figure4A_intercept_IFG.png`, `Figure4B_difficulty_SubCal.png`

### Figure 5: Reproducibility across behavioral sessions S1–S5

* Panel A: heatmaps of the correlation between ROI-wise association profiles for each pair of sessions, per parameter
* Panel B: partial r across S1–S5 for the two focal associations; all other ROIs in grey
* Output: `Figure5A_profile_similarity.png`, `Figure5B_focal_ROI_trajectories.png`

### Figure 6: Reproducibility across MRI sessions S6–S7

* Panel A: violin, box and point plot of the within-ROI correlation between the S6 and S7 T maps (intercept, difficulty)
* Panel B: the five ROIs with the lowest and highest S6–S7 correlation per effect (NAcc seed excluded)
* Output: `Figure6A_spatial_reproducibility_distribution.png`, `Figure6B_spatial_reproducibility_extreme_ROIs.png`

### Supplementary figures

* Figure S1: focal associations by sex, with sex-specific fitted lines (`FigureS1_sex_brain_behavior.png`)
* Figure S2: partial r across S1–S5 in the two focal ROIs for all three parameters (`FigureS2_focal_ROI_behavioral_specificity.png`)
* Figure S3: ROI NAcc FC reliability (ICC) vs brain–behavior association, per ROI (A) and per network (B) (`FigureS3A.png`, `FigureS3B.png`)
* Figure S4: effort vs wanting coefficients (`FigureS4_effort_wanting.png`)
* Figure S5: correlation matrices of effort and wanting coefficients with questionnaire scores (`FigureS5_correlation_matrices.png`)
* Figure S6: effort coefficients vs BAS and BIS (`FigureS6_BIS_BAS_effort.png`)

\---

## Requirements

### R packages (R 4.4.1)

```r
# Data handling and visualization
library(tidyverse)     # ggplot2 >= 3.5.0
library(cowplot)       # plot composition, theme
library(patchwork)     # plot composition
library(ggcorrplot)    # correlation heatmaps (Figure 2A)
library(ggrepel)       # ROI labels (Figure S3)
library(reshape2)      # correlation matrices (Figure S5)

# Statistical modeling
library(lme4)          # linear mixed-effects models
library(lmerTest)      # p-values for mixed-effects models
library(irr)           # intraclass correlations
library(psych)         # correlation tests
```

### MATLAB

* MATLAB R2021b with the Statistics and Machine Learning Toolbox (`fingerprinting.m`)
* Optional: `plasma` colormap on the MATLAB path (otherwise `parula` is used)

### Upstream software (not included)

fMRIPrep 20.2.7, CONN 22.a, SPM12 (v7771), fmreli (beta 1.0)

### Theme settings

The R script sets a global plotting theme:

```r
theme_set(theme_cowplot(font_size = 12))
```

### Running the analyses

```r
setwd("path/to/repository")
source("code/analyses.R")
```

```matlab
run('path/to/repository/code/fingerprinting.m')
```

\---

## Citation

If you use this code or data, please cite:

> Karlocai, Zs., Kühnel, A., Schulz, C., Grahlow, M., Kaduk, K., Frick, A., Derntl, B., & Kroemer, N. B. How reliable is the link between effort and nucleus accumbens connectivity at rest? A psychometric evaluation of repeated assessments. (manuscript submitted)

\---

## Contact

For questions or issues, please contact the corresponding authors or open an issue in this repository.

\---

## Notes

* Wanting ratings are on a VAS from 0 to 100.
* FC values are Fisher z-transformed correlations (mean of left and right NAcc seed, averaged across S6 and S7 where stated).
* Uncertain trials (S2–S5) are excluded from all analyses and are not included in the data.
* Participants are kept in their original order in the mixed-effects models, because `lme4` estimates depend slightly on it.
* Missing values are excluded per analysis (complete cases or pairwise correlations).
* Not included: steps that require the individual MRI data, preprocessing (fMRIPrep, CONN), SPM second-level models (whole-brain results, Table S3), fmreli (Figure 2C), Figure 3A, and the extraction of the tables in `data/`.

