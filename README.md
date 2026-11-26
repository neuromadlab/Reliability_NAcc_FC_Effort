# How reliable is the link between effort and nucleus accumbens connectivity at rest?

Code and data for:

> Karlocai, Zs., Kühnel, A., Schulz, C., Grahlow, M., Kaduk, K., Frick, A., Derntl, B., & Kroemer, N. B. *How reliable is the link between effort and nucleus accumbens connectivity at rest? A psychometric evaluation of repeated assessments.* [Journal, year, DOI — add after acceptance]

- Study preregistration (parent study): https://osf.io/gpr52
- Unthresholded statistical maps: https://neurovault.org/collections/23461/
- fmreli toolbox (voxel/ROI ICC): https://github.com/nkroemer/reliability

## Repository structure

```
code/
  analyses.R          behavioral models, reliability, brain-behavior analyses,
                      figures and supplementary tables (manuscript order)
  fingerprinting.m    functional connectome fingerprinting (Figure 3B, 3C)
data/                 anonymised input data (see "Data")
output/               created by the scripts
```

## Session naming

| Manuscript | Content |
|---|---|
| S1 | intake session, effort allocation task (EAT) |
| S2–S5 | weekly sessions, EAT with certain and uncertain trials |
| S6, S7 | task-free fMRI sessions |

## Data

Participant IDs (`sub-001`, …) are random and consistent across all files. Only the variables used in the analyses are included.

| File | Content | Used for |
|---|---|---|
| `eat_trials.csv` | Trial-level EAT data, S1–S5 (certain trials only) | S1 mixed-effects models, behavioral reliability (Fig. 2A, Table S1), all participant coefficients |
| `mri_sample.csv` | MRI sample (n = 84): SPM covariates, S1 effort coefficients, sex, NAcc FC of the two focal ROIs | Fig. 4, Fig. S1, Table S2 |
| `questionnaires.csv` | SHAPS, LARS, BIS and BAS sum scores | Figs. S4–S6 |
| `session_timing.csv` | Start time of each behavioral session and days to the MRI sessions (no dates) | Session intervals (Methods), Table S2 |
| `roi_brain_behavior_partial_r.csv` | Per atlas ROI: partial r between NAcc FC and the S1–S5 coefficients, NAcc FC reliability (ICC), network | Fig. 5, Figs. S2, S3 |
| `roi_S6_S7_tmap_similarity.csv` | Per atlas ROI: correlation of the S6 and S7 T maps of the S1 coefficients | Fig. 6 |
| `fingerprinting_similarity_matrix.csv` | Correlation of whole-brain (gray-matter) NAcc FC maps between all participants, S6 (rows) × S7 (columns) | Fig. 3B, 3C |
| `group_fc_overlap_S6_S7.csv` | Dice and Jaccard overlap of the S6 and S7 group NAcc FC maps (p < .01 unc.) | Fig. 2B |

### Variables

`eat_trials.csv`

| Variable | Meaning |
|---|---|
| `ID` | participant |
| `Session` | behavioral session (1–5) |
| `Trial` | trial number |
| `RewType` | reward type (0 = food, 1 = money) |
| `RewMag` | reward magnitude (low / high) |
| `Difficulty` | difficulty (S1: 75 % / 85 % of maximum frequency; S2–S5: continuous) |
| `RelEffort` | mean button-press rate relative to the individual maximum |
| `Wanting` | wanting rating (VAS 0–100) |

`mri_sample.csv`

| Variable | Meaning |
|---|---|
| `Female` | 1 = women, 0 = men |
| `cAge`, `cSex`, `cBMI` | centered covariates of the SPM models |
| `FD` | log mean framewise displacement |
| `Effort_intercept`, `Effort_rewmag`, `Effort_diff` | z-scored S1 effort coefficients (SPM predictors) |
| `FC_rIFGoper`, `FC_subcallosal` | NAcc FC (Fisher z, mean of S6 and S7) in the right IFG pars opercularis and subcallosal cortex atlas ROIs |

`questionnaires.csv`: `SHAPS` (sum of 14 items scored 1–4; rescaled to 0–42 in the script), `LARS`, `BIS`, `BAS` (sum scores).

`roi_brain_behavior_partial_r.csv`: `Intercept_S1` … `Difficulty_S5` are partial correlations controlling for age, sex, BMI, FD and the other two coefficients; `ICC_con` is the ICC(3,1) (consistency) of ROI NAcc FC across S6 and S7 from fmreli.

Raw and preprocessed MRI data are not shared [state reason, e.g. participants did not consent to public sharing; access on reasonable request]. Group-level unthresholded maps are available on NeuroVault (link above).

## Software

- R 4.4.1 with tidyverse, lme4, lmerTest, irr, psych, cowplot, ggcorrplot, patchwork, ggrepel, reshape2 (ggplot2 ≥ 3.5.0)
- MATLAB R2021b with the Statistics and Machine Learning Toolbox (`fingerprinting.m`)
- Upstream (not included): fMRIPrep 20.2.7, CONN 22.a, SPM12 (v7771), fmreli (beta 1.0)

## Running the analyses

Run both scripts from the repository root; outputs are written to `output/`.

```r
setwd("path/to/repository")
source("code/analyses.R")
```

```matlab
cd path/to/repository
run code/fingerprinting.m
```

## Script-to-result map

| Result | Script (section) | Output |
|---|---|---|
| Session intervals (Methods) | analyses.R A2 | `session_intervals.csv` |
| S1 mixed-effects models (Methods) | analyses.R A3–A4 | console, `coef_relEff*.csv`, `coef_wanting*.csv` |
| Behavioral reliability, Fig. 2A | analyses.R B1 | `Figure2A_behavioral_reliability.png` |
| Dice/Jaccard, Fig. 2B | — | `data/group_fc_overlap_S6_S7.csv` |
| Fingerprinting, Fig. 3B, 3C | fingerprinting.m | `fingerprinting_results.mat`, `Figure3B/3C_*.png` |
| NAcc FC and effort, Fig. 4 | analyses.R B2 | `Figure4A/B_*.png` |
| Behavioral sessions S1–S5, Fig. 5 | analyses.R B3 | `Figure5A/B_*.png` |
| MRI sessions S6–S7, Fig. 6 | analyses.R B4 | `Figure6A/B_*.png` |
| Effort, wanting and questionnaires | analyses.R B5 | `effort_wanting_correlations.csv`, `effort_questionnaire_correlations.csv` |
| Table S1 | analyses.R C | `TableS1_behavioral_reliability.csv` |
| Fig. S1 | analyses.R C | `FigureS1_sex_brain_behavior.png`, `FigureS1_sex_interactions.csv` |
| Table S2 | analyses.R C | `TableS2_timing_sensitivity.csv` |
| Figs. S2, S3 | analyses.R C | `FigureS2_*.png`, `FigureS3A/B_*.png` |
| Figs. S4–S6 | analyses.R C | `FigureS4_*` – `FigureS6_*.png` |
| Whole-brain results, Table S3; Fig. 2C (voxel/ROI ICC); Fig. 3A | SPM, fmreli, MRI data | not included |

Figure panels are saved separately; the manuscript figures were assembled in BioRender.

## Not included

Preprocessing (fMRIPrep, CONN), SPM second-level models, fmreli, and the extraction of the tables in `data/` from the MRI data, as these steps require the individual MRI data. [Add: available on request / in a separate repository.]

## License

Code: [e.g. MIT]. Data: [e.g. CC BY 4.0].

## Contact

Nils B. Kroemer — nkroemer@uni-bonn.de
