#############################################################################
# How reliable is the link between effort and nucleus accumbens connectivity
# at rest? A psychometric evaluation of repeated assessments
#
# R analyses, in the order of the manuscript:
#
#   0  Setup
#   A  Methods
#      A1  Data
#      A2  Session intervals
#      A3  S1 mixed-effects models (effort, wanting)
#      A4  Participant-specific coefficients for SPM
#   B  Results
#      B1  Behavioral test-retest reliability (Figure 2A)
#      B2  NAcc FC and effort (Figure 4)
#      B3  Reproducibility across behavioral sessions S1-S5 (Figure 5)
#      B4  Reproducibility across MRI sessions S6-S7 (Figure 6)
#      B5  Effort, wanting and questionnaires
#   C  Supplementary material
#      Table S1   Behavioral reliability with and without S1
#      Figure S1  Sex x behavior interactions
#      Table S2   Timing sensitivity of the focal associations
#      Figure S2  Focal ROIs, all behavioral parameters
#      Figure S3  ROI FC reliability vs brain-behavior associations
#      Figure S4  Effort vs wanting estimates
#      Figure S5  Correlations with questionnaires
#      Figure S6  Effort estimates vs BIS/BAS
#
# All inputs are in data/ (see README). B2-B4, Table S2 and Figures S1-S3 use
# the MRI-derived tables (mri_sample.csv, roi_*.csv), which were extracted
# from the SPM second-level models. Figure 3 is made by fingerprinting.m.
#
# Run from the repository root (e.g. setwd("path/to/repository")).
##############################################################################


##############################################################################
# 0  SETUP
##############################################################################

library(tidyverse)
library(lme4)
library(lmerTest)
library(irr)
library(psych)
library(cowplot)
library(ggcorrplot)
library(patchwork)
library(ggrepel)
library(reshape2)

theme_set(theme_cowplot(font_size = 12))


# Paths ----------------------------------------------------------------------

data_dir   <- "data"
output_dir <- "output"

trial_file  <- file.path(data_dir, "eat_trials.csv")
mri_file    <- file.path(data_dir, "mri_sample.csv")
quest_file  <- file.path(data_dir, "questionnaires.csv")
timing_file <- file.path(data_dir, "session_timing.csv")
roi_behavior_file <- file.path(data_dir, "roi_brain_behavior_partial_r.csv")
roi_tmap_file     <- file.path(data_dir, "roi_S6_S7_tmap_similarity.csv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


read_data <- function(file) {
  read_csv(file, show_col_types = FALSE, col_types = cols(ID = col_character()))
}

save_png <- function(plot, filename, width, height, dpi = 600) {
  ggsave(file.path(output_dir, paste0(filename, ".png")), plot,
         width = width, height = height, units = "in", dpi = dpi, bg = "white")
}



drop_zero <- function(x) sub("^(-?)0\\.", "\\1.", x)
fmt_r <- function(r) {
  drop_zero(if (abs(r) >= 0.1) sprintf("%.2f", r) else format(signif(r, 2), scientific = FALSE))
}
fmt_p <- function(p) {
  if (p < 0.001) return("italic(p) < '.001'")
  paste0("italic(p) == '", drop_zero(sprintf(if (p < 0.1) "%.3f" else "%.2f", p)), "'")
}
cor_labels <- function(data, x, y, groups = character(0)) {
  data %>%
    group_by(across(all_of(groups))) %>%
    summarise(r = cor.test(.data[[x]], .data[[y]])$estimate,
              p = cor.test(.data[[x]], .data[[y]])$p.value,
              .groups = "drop") %>%
    mutate(r_label = map_chr(r, ~ paste0("italic(r) == '", fmt_r(.x), "'")),
           p_label = map_chr(p, fmt_p),
           label   = paste0(r_label, "*','~", p_label))
}


##############################################################################
# A  METHODS
##############################################################################

# A1  Data -------------------------------------------------------------------

d_EAT <- read_data(trial_file)

# A2  Session intervals ------------------------------------------------------

has_timing <- file.exists(timing_file)

if (has_timing) {

  session_timing <- read_data(timing_file)
  stopifnot(!anyDuplicated(session_timing[c("ID", "Session")]))

  intervals <- session_timing %>%
    group_by(ID) %>%
    summarise(S5_to_S6 = first(days_to_S6[Session == 5]),
              S5_to_S7 = first(days_to_S7[Session == 5]),
              S6_to_S7 = first(na.omit(days_to_S7 - days_to_S6)),
              .groups = "drop")

  # One participant whose last behavioral session took place after the MRI
  # session (negative S5-MRI interval) is excluded from the interval summary.
  s5_after_mri <- intervals$ID[intervals$S5_to_S6 < 0 | intervals$S5_to_S7 < 0]
  s5_after_mri <- s5_after_mri[!is.na(s5_after_mri)]
  cat("Excluded from the session-interval summary (S5 after MRI):",
      if (length(s5_after_mri)) paste(s5_after_mri, collapse = ", ") else "none", "\n")

  interval_summary <- intervals %>%
    filter(!ID %in% s5_after_mri) %>%
    pivot_longer(-ID, names_to = "Interval", values_to = "Days") %>%
    group_by(Interval) %>%
    summarise(N = sum(!is.na(Days)), N_missing = sum(is.na(Days)),
              Mean = mean(Days, na.rm = TRUE), SD = sd(Days, na.rm = TRUE),
              Min = min(Days, na.rm = TRUE), Max = max(Days, na.rm = TRUE),
              .groups = "drop")

  print(interval_summary)
  write.csv(interval_summary, file.path(output_dir, "session_intervals.csv"), row.names = FALSE)

} else {
  message("No ", timing_file, ": skipping session intervals (A2) and Table S2.")
}


# A3  S1 mixed-effects models: relative effort and wanting -------------------

d_session1 <- d_EAT %>%
  filter(Session == 1) %>%
  mutate(
    cTrial      = Trial - mean(Trial, na.rm = TRUE),
    cDifficulty = Difficulty - mean(Difficulty, na.rm = TRUE),
    cRewMag     = RewMag - mean(RewMag, na.rm = TRUE),
    cRewType    = RewType - mean(RewType, na.rm = TRUE),
    fID         = factor(ID, levels = unique(ID))
  )

cat("Session 1 participants:", n_distinct(d_session1$ID), "\n")

lm_wanting <- lmer(
  Wanting ~ cRewMag * cDifficulty * cRewType + cTrial +
    (1 + cDifficulty + cRewMag + cRewType + cTrial | fID),
  data = d_session1
)
print(summary(lm_wanting))

lm_Effort <- lmer(
  RelEffort ~ cRewMag * cDifficulty * cRewType + cTrial +
    (1 + cDifficulty + cRewMag + cRewType + cTrial | fID),
  data = d_session1
)
print(summary(lm_Effort))


# A4  Participant-specific coefficients for SPM ------------------------------

coef_wanting_raw <- coef(lm_wanting)$fID
coef_relEff_raw  <- coef(lm_Effort)$fID
coef_wanting_raw$ID <- rownames(coef_wanting_raw)
coef_relEff_raw$ID  <- rownames(coef_relEff_raw)

write.csv(coef_wanting_raw, file.path(output_dir, "coef_wanting.csv"), row.names = FALSE)
write.csv(coef_relEff_raw,  file.path(output_dir, "coef_relEff.csv"),  row.names = FALSE)

zscore_coefs <- function(coefs) {
  coefs %>%
    select(`(Intercept)`, cRewMag, cDifficulty, ID) %>%
    mutate(across(c(`(Intercept)`, cRewMag, cDifficulty), ~ as.numeric(scale(.x))))
}

write.csv(zscore_coefs(coef_wanting_raw),
          file.path(output_dir, "coef_wanting_ZScore.csv"), row.names = FALSE)
write.csv(zscore_coefs(coef_relEff_raw),
          file.path(output_dir, "coef_relEff_ZScore.csv"), row.names = FALSE)


##############################################################################
# B  RESULTS
##############################################################################


# B1  Behavioral test-retest reliability (Figure 2A) -------------------------

common_ids <- Reduce(intersect, lapply(1:5, function(s) {
  unique(d_EAT$ID[d_EAT$Session == s])
}))
cat("Participants with behavioral data in S1-S5:", length(common_ids), "\n")

coef_list <- lapply(1:5, function(session) {

  d_session <- d_EAT %>%
    filter(Session == session, ID %in% common_ids) %>%
    mutate(
      cTrial      = as.numeric(scale(Trial, scale = FALSE)),
      cDifficulty = as.numeric(scale(Difficulty, scale = FALSE)),
      cRewMag     = as.numeric(scale(RewMag, scale = FALSE)),
      cRewType    = as.numeric(scale(RewType, scale = FALSE)),
      fID         = factor(ID, levels = unique(ID))
    )

  model <- lmer(
    RelEffort ~ cRewMag * cDifficulty * cRewType + cTrial +
      (1 + cDifficulty + cRewMag + cRewType + cTrial | fID),
    data = d_session
  )

  coef(model)$fID %>%
    select(`(Intercept)`, cRewMag, cDifficulty) %>%
    mutate(across(everything(), ~ as.numeric(scale(.x))),
           ID = rownames(.),
           Session = session)
})

coef_long <- bind_rows(coef_list)

coef_wide <- coef_long %>%
  pivot_wider(
    names_from  = Session,
    values_from = c(`(Intercept)`, cRewMag, cDifficulty),
    names_sep   = "_S"
  )

stopifnot(!anyDuplicated(coef_wide$ID), nrow(coef_wide) == length(common_ids))

parameters <- c(
  "Intercept"        = "(Intercept)",
  "Reward magnitude" = "cRewMag",
  "Difficulty"       = "cDifficulty"
)

get_icc <- function(data, cols) {
  irr::icc(data[, cols, drop = FALSE],
           model = "twoway", type = "consistency", unit = "single")$value
}

for (label in names(parameters)) {
  cat("\nPairwise correlations:", label, "\n")
  print(psych::corr.test(coef_wide[, paste0(parameters[[label]], "_S", 1:5)]))
}

icc_S1_S5 <- sapply(parameters, function(p) get_icc(coef_wide, paste0(p, "_S", 1:5)))
cat("\nICC(3,1), S1-S5\n")
print(icc_S1_S5)

session_cor <- function(prefix) {
  m <- cor(coef_wide[, paste0(prefix, "_S", 1:5)], use = "pairwise.complete.obs")
  dimnames(m) <- list(paste0("S", 1:5), paste0("S", 1:5))
  m
}

mat_to_long <- function(mat) {
  expand.grid(Row = rownames(mat), Col = colnames(mat), stringsAsFactors = FALSE) %>%
    mutate(
      Row   = factor(Row, levels = rownames(mat)),
      Col   = factor(Col, levels = colnames(mat)),
      value = mat[cbind(as.numeric(Row), as.numeric(Col))]
    ) %>%
    filter(as.numeric(Row) > as.numeric(Col))
}

plot_reliability <- function(prefix, title, icc_value, show_legend = FALSE) {
  mat <- session_cor(prefix)
  p <- ggcorrplot::ggcorrplot(
    mat, type = "upper", lab = FALSE, title = title,
    ggtheme = theme_cowplot(), show.legend = show_legend,
    outline.color = "white", colors = c("white", "#4A90E2")
  ) +
    scale_fill_gradient(limits = c(0, 1), low = "white", high = "#4A90E2") +
    theme(
      plot.title  = element_text(hjust = 0.5, size = 18),
      axis.text.x = element_text(face = "bold", size = 18),
      axis.text.y = element_text(face = "bold", size = 18)
    ) +
    geom_text(data = mat_to_long(mat),
              aes(x = Col, y = Row, label = sprintf("%.2f", value)), size = 6) +
    annotate("text", y = -Inf, x = 2,
             label = paste0("ICC = ", sprintf("%.2f", icc_value)),
             hjust = -0.2, vjust = -1.5, fontface = "bold", size = 6)
  if (show_legend) p <- p + guides(fill = guide_colorbar(title = "Correlation"))
  p
}

Figure2A <- cowplot::plot_grid(
  plot_reliability("(Intercept)", "Intercept",        icc_S1_S5[["Intercept"]]),
  plot_reliability("cRewMag",     "Reward magnitude", icc_S1_S5[["Reward magnitude"]]),
  plot_reliability("cDifficulty", "Difficulty",       icc_S1_S5[["Difficulty"]],
                   show_legend = TRUE),
  ncol = 3, rel_widths = c(1, 1, 1.4)
)
save_png(Figure2A, "Figure2A_behavioral_reliability", width = 10, height = 5)


# B2  NAcc FC and effort (Figure 4) ------------------------------------------

# Whole-brain inference is performed in SPM. FC values were extracted from the
# atlas regions overlapping the significant clusters. The partial correlation
# shown equals the partial r of the focal predictor in the full model.

# mri_sample.csv: one row per participant of the MRI sample with the SPM
# design covariates (cAge, cSex, cBMI centered; FD = log mean framewise
# displacement), the z-scored S1 effort coefficients used in SPM
# (Effort_*), sex (Female: 1 = women, 0 = men) and NAcc FC extracted from the
# right IFG pars opercularis and subcallosal cortex atlas ROIs.

mri_sample <- read_data(mri_file) %>%
  left_join(
    coef_relEff_raw %>%
      transmute(ID,
                Intercept_raw  = `(Intercept)`,
                Difficulty_raw = cDifficulty),
    by = "ID"
  )


cat("Mean NAcc-IFG FC:", mean(mri_sample$FC_rIFGoper, na.rm = TRUE), "\n")
cat("Mean NAcc-subcallosal FC:", mean(mri_sample$FC_subcallosal, na.rm = TRUE), "\n")

plot_adjusted <- function(y, x_model, x_plot, controls, xlab, ylab) {

  d <- mri_sample %>%
    select(all_of(c(y, x_model, x_plot, controls))) %>%
    drop_na()

  fit <- lm(reformulate(c(x_model, controls), response = y), data = d)
  print(summary(fit))

  resid_on_controls <- function(v) {
    residuals(lm(reformulate(controls, response = v), data = d)) + mean(d[[v]])
  }
  d$x_adjusted <- resid_on_controls(x_plot)
  d$y_adjusted <- resid_on_controls(y)

  ggplot(d, aes(x = x_adjusted, y = y_adjusted)) +
    geom_point(color = "black", size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, color = "darkred") +
    geom_text(data = cor_labels(d, "x_adjusted", "y_adjusted"), aes(label = label),
              x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 4, parse = TRUE,
              inherit.aes = FALSE) +
    labs(x = xlab, y = ylab) +
    theme_cowplot(12) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12))
}

Figure4A <- plot_adjusted(
  y = "FC_rIFGoper", x_model = "Effort_intercept", x_plot = "Intercept_raw",
  controls = c("cAge", "cSex", "cBMI", "FD", "Effort_rewmag", "Effort_diff"),
  xlab = "Intercept (adjusted)", ylab = "NAcc-IFG FC"
)
save_png(Figure4A, "Figure4A_intercept_IFG", width = 4.2, height = 3.4, dpi = 300)

Figure4B <- plot_adjusted(
  y = "FC_subcallosal", x_model = "Effort_diff", x_plot = "Difficulty_raw",
  controls = c("cAge", "cSex", "cBMI", "FD", "Effort_rewmag", "Effort_intercept"),
  xlab = "Difficulty (adjusted)", ylab = "NAcc-SubCal FC"
)
save_png(Figure4B, "Figure4B_difficulty_SubCal", width = 4.2, height = 3.4, dpi = 300)


# B3  Reproducibility across behavioral sessions S1-S5 (Figure 5) ------------

sessions  <- paste0("S", 1:5)
behaviors <- c("Intercept", "Reward magnitude", "Difficulty")
behavior_colors <- c("Intercept" = "#D73027",
                     "Reward magnitude" = "#7B3294",
                     "Difficulty" = "#2166AC")

theme_paper <- theme_classic(base_size = 13) +
  theme(axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

partial_r_label <- expression("Brain-behavior association (partial " * italic(r) * ")")

# ROI-wise partial r between averaged S6/S7 NAcc FC and the S1-S5 behavioral
# coefficients (controlling for age, sex, BMI, FD and the other two
# coefficients), ROI NAcc FC reliability (ICC_con, fmreli ICC(3,1)) and
# network labels.
master <- read_csv(roi_behavior_file, show_col_types = FALSE)

required_master <- c("ROI", "ROI_short", "Network", "ICC_con",
                     paste0(rep(c("Intercept", "RewMag", "Difficulty"), each = 5),
                            "_S", 1:5))
missing_master <- setdiff(required_master, names(master))
if (length(missing_master)) {
  stop("Missing columns in ", roi_behavior_file, ": ", paste(missing_master, collapse = ", "))
}

roi_long <- master %>%
  select(ROI, ROI_short, Network, ICC_con,
         matches("^(Intercept|RewMag|Difficulty)_S[1-5]$")) %>%
  pivot_longer(
    cols          = matches("^(Intercept|RewMag|Difficulty)_S[1-5]$"),
    names_to      = c("Behavior", "Session"),
    names_pattern = "(Intercept|RewMag|Difficulty)_(S[1-5])",
    values_to     = "partial_r"
  ) %>%
  mutate(Behavior = factor(recode(Behavior, RewMag = "Reward magnitude"),
                           levels = behaviors),
         Session  = factor(Session, levels = sessions))

IFG_ROI    <- "IFG oper r"
SUBCAL_ROI <- "SubCalC"
stopifnot(IFG_ROI %in% roi_long$ROI_short, SUBCAL_ROI %in% roi_long$ROI_short)


# Figure 5A: similarity of ROI-wise association profiles across sessions

profile_one <- function(prefix, label) {
  R  <- cor(as.matrix(master[, paste0(prefix, "_S", 1:5)]),
            use = "pairwise.complete.obs")
  ij <- which(lower.tri(R), arr.ind = TRUE)
  tibble(Behavior = label,
         Earlier  = sessions[ij[, 2]],
         Later    = sessions[ij[, 1]],
         r        = R[ij])
}

profile_data <- bind_rows(
  profile_one("Intercept",  "Intercept"),
  profile_one("RewMag",     "Reward magnitude"),
  profile_one("Difficulty", "Difficulty")
) %>%
  mutate(Behavior = factor(Behavior, levels = behaviors),
         Earlier  = factor(Earlier, levels = c("S1", "S2", "S3", "S4")),
         Later    = factor(Later,   levels = c("S5", "S4", "S3", "S2")))

Figure5A <- ggplot(profile_data, aes(x = Earlier, y = Later, fill = r)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4.2) +
  facet_wrap(~Behavior, nrow = 1, axes = "all", axis.labels = "all") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 13) +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(size = 15, face = "bold", margin = margin(b = 4)),
    axis.text.x       = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
    axis.text.y       = element_text(size = 12, color = "black"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    axis.ticks        = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.08, "cm"),
    panel.spacing     = unit(0.45, "cm"),
    legend.title      = element_text(size = 13, face = "bold"),
    legend.text       = element_text(size = 11)
  )
save_png(Figure5A, "Figure5A_profile_similarity", width = 8.5, height = 3.5, dpi = 400)


# Figure 5B: focal associations across sessions

panel_levels <- c("Intercept – right IFG", "Difficulty – Subcallosal cortex")

fig5b_background <- bind_rows(
  roi_long %>% filter(Behavior == "Intercept")  %>% mutate(Panel = panel_levels[1]),
  roi_long %>% filter(Behavior == "Difficulty") %>% mutate(Panel = panel_levels[2])
) %>%
  mutate(Panel = factor(Panel, levels = panel_levels))

fig5b_data <- bind_rows(
  roi_long %>% filter(ROI_short == IFG_ROI, Behavior == "Intercept") %>%
    mutate(Panel = panel_levels[1]),
  roi_long %>% filter(ROI_short == SUBCAL_ROI, Behavior == "Difficulty") %>%
    mutate(Panel = panel_levels[2])
) %>%
  mutate(Panel = factor(Panel, levels = panel_levels))

Figure5B <- ggplot() +
  geom_line(data = fig5b_background, aes(Session, partial_r, group = ROI),
            color = "grey65", alpha = 0.20, linewidth = 0.35) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(data = fig5b_data, aes(Session, partial_r, group = 1, color = Panel),
            linewidth = 1.2) +
  geom_point(data = fig5b_data, aes(Session, partial_r, color = Panel), size = 2.3) +
  facet_wrap(~Panel, nrow = 1, axes = "all", axis.labels = "all") +
  scale_color_manual(values = setNames(c("#D73027", "#2166AC"), panel_levels)) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  labs(x = "Behavioral session", y = partial_r_label) +
  theme_paper +
  theme(legend.position = "none",
        strip.text  = element_text(size = 12.5, face = "bold"),
        axis.title  = element_text(size = 11.5),
        axis.text   = element_text(size = 10.5, color = "black"),
        plot.margin = margin(3, 5, 3, 5))
save_png(Figure5B, "Figure5B_focal_ROI_trajectories", width = 9.5, height = 4.5)


# B4  Reproducibility across MRI sessions S6-S7 (Figure 6) -------------------

# Within-ROI Pearson r between the unthresholded S6 and S7 T maps of the S1
# effort coefficients (Intercept, Difficulty).
mri <- read_csv(roi_tmap_file, show_col_types = FALSE) %>%
  mutate(Effect = factor(Effect, levels = c("Intercept", "Difficulty"))) %>%
  filter(!is.na(Effect), is.finite(Pearson_r))

stopifnot(all(c("Intercept", "Difficulty") %in% mri$Effect))

mri_colors <- c("Intercept" = "#2C7FB8", "Difficulty" = "#F16913")


# Figure 6A: distribution of within-ROI S6-S7 spatial correlations

set.seed(2026)

Figure6A <- ggplot(mri, aes(x = Effect, y = Pearson_r, fill = Effect, color = Effect)) +
  geom_violin(width = 0.72, trim = FALSE, alpha = 0.20, linewidth = 1.0) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.58, size = 2.15, stroke = 0) +
  geom_boxplot(width = 0.16, outlier.shape = NA, color = "black",
               alpha = 0.38, linewidth = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  scale_fill_manual(values = mri_colors) +
  scale_color_manual(values = mri_colors) +
  scale_y_continuous(breaks = seq(-1, 1, 0.2), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = NULL, y = "Within-ROI S6-S7 Pearson's r") +
  theme_classic(base_size = 13) +
  theme(
    legend.position    = "none",
    axis.text.x        = element_text(size = 14, face = "bold", color = "black",
                                      margin = margin(b = 4)),
    axis.text.y        = element_text(size = 11, color = "black"),
    axis.title.y       = element_text(size = 13),
    plot.margin        = margin(5, 5, 3, 8)
  )
save_png(Figure6A, "Figure6A_spatial_reproducibility_distribution",
         width = 8.5, height = 4.6, dpi = 400)


# Figure 6B: five lowest and five highest ROIs

short_roi <- function(x) {
  z <- tolower(x)
  side <- case_when(
    str_detect(z, "right|\\br\\b") ~ "r",
    str_detect(z, "left|\\bl\\b")  ~ "l",
    TRUE ~ ""
  )
  case_when(
    str_detect(z, "amygdala")                       ~ paste("Amy", side),
    str_detect(z, "inferior frontal.*operc|ifg oper") ~ paste("IFG", side),
    str_detect(z, "substantia nigra")               ~ paste("SN", side),
    str_detect(z, "ventral tegmental|\\bvta\\b")    ~ paste("VTA", side),
    str_detect(z, "nucleus ruber|red nucleus")      ~ paste("NR", side),
    str_detect(z, "middle temporal|\\bmtg\\b")      ~ paste("MTG", side),
    str_detect(z, "mammil")                         ~ paste("MB", side),
    str_detect(z, "supracalcarine")                 ~ paste("SCC", side),
    str_detect(z, "ventral pallid")                 ~ paste("VP", side),
    str_detect(z, "inferior temporal|\\bitg\\b")    ~ paste("ITG", side),
    str_detect(z, "subcallosal|subcalc")            ~ "SubCal",
    str_detect(z, "pallidum")                       ~ paste("Pallidum", side),
    str_detect(z, "vermis.*1.*2")                   ~ "Vermis 1 2",
    TRUE ~ str_remove(x, "\\s*\\(.*$")
  ) %>%
    str_squish()
}

fig6b_data <- mri %>%
  filter(!str_detect(str_to_lower(ROI_name), "accumbens")) %>%
  group_by(Effect) %>%
  group_modify(~ bind_rows(
    slice_min(.x, Pearson_r, n = 5, with_ties = FALSE),
    slice_max(.x, Pearson_r, n = 5, with_ties = FALSE)
  )) %>%
  ungroup() %>%
  distinct(Effect, Atlas_label, .keep_all = TRUE) %>%
  mutate(Label = short_roi(ROI_name))

make_ranked_mri_panel <- function(effect_name) {

  plot_data <- fig6b_data %>%
    filter(Effect == effect_name) %>%
    arrange(Pearson_r) %>%
    mutate(Plot_ID = factor(paste0("ROI_", Atlas_label),
                            levels = paste0("ROI_", Atlas_label)))

  ggplot(plot_data, aes(x = Pearson_r, y = Plot_ID)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.5) +
    geom_segment(aes(x = 0, xend = Pearson_r, yend = Plot_ID),
                 color = "grey82", linewidth = 0.45) +
    geom_point(aes(fill = Effect), shape = 21, size = 4.0, color = "black",
               stroke = 0.45, alpha = 0.95) +
    scale_fill_manual(values = mri_colors, guide = "none") +
    scale_y_discrete(labels = setNames(plot_data$Label, as.character(plot_data$Plot_ID))) +
    scale_x_continuous(limits = c(-1.05, 1.05), breaks = c(-1, -0.5, 0, 0.5, 1),
                       labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
                       expand = expansion(mult = c(0.01, 0.01))) +
    labs(x = "Within-ROI S6-S7 Pearson's r", y = NULL, title = effect_name) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title   = element_text(face = "bold", size = 14, hjust = 0.5,
                                  margin = margin(b = 4)),
      axis.title.x = element_text(size = 12),
      axis.text.x  = element_text(size = 10.5, color = "black"),
      axis.text.y  = element_text(size = 11, color = "black"),
      plot.margin  = margin(3, 5, 5, 8)
    )
}

Figure6B <- make_ranked_mri_panel("Intercept") | make_ranked_mri_panel("Difficulty")
save_png(Figure6B, "Figure6B_spatial_reproducibility_extreme_ROIs",
         width = 10, height = 5.2, dpi = 400)


# B5  Effort, wanting and questionnaires -------------------------------------

# Uses the non-z-scored S1 coefficients of the primary MRI sample.
# questionnaires.csv: SHAPS (sum of 14 items, 1-4), LARS, BIS and BAS sum scores.

questionnaires <- read_data(quest_file) %>%
  mutate(SHAPS = SHAPS - 14) %>%
  filter(is.finite(SHAPS), SHAPS >= 0)

add_questionnaires <- function(coefs) {
  coefs %>%
    transmute(ID,
              Intercept = `(Intercept)`,
              RewMag    = cRewMag,
              RewDiff   = cDifficulty) %>%
    inner_join(questionnaires, by = "ID")
}

QUEST_Estimates   <- add_questionnaires(coef_relEff_raw)
QUEST_Estimates_W <- add_questionnaires(coef_wanting_raw)

cat("Questionnaire sample, effort:",  n_distinct(QUEST_Estimates$ID), "\n")
cat("Questionnaire sample, wanting:", n_distinct(QUEST_Estimates_W$ID), "\n")

to_long <- function(d) {
  d %>%
    pivot_longer(c(Intercept, RewMag, RewDiff),
                 names_to = "Parameter", values_to = "Estimate") %>%
    mutate(Parameter = factor(
      recode(Parameter, RewMag = "Reward magnitude", RewDiff = "Difficulty"),
      levels = c("Intercept", "Reward magnitude", "Difficulty")
    ))
}

QUEST_Estimates_long   <- to_long(QUEST_Estimates)
QUEST_Estimates_long_W <- to_long(QUEST_Estimates_W)

theme_scatter <- theme(
  legend.position = "none",
  text        = element_text(face = "bold", size = 16),
  axis.text   = element_text(face = "plain", size = 12),
  axis.title  = element_text(size = 16),
  strip.text  = element_text(size = 16, face = "bold")
)

Estimates_corr <- inner_join(
  QUEST_Estimates_long   %>% select(ID, Parameter, Effort  = Estimate),
  QUEST_Estimates_long_W %>% select(ID, Parameter, Wanting = Estimate),
  by = c("ID", "Parameter")
)

# Correlations reported in the text (plotted in Figures S4-S6)
effort_wanting_r <- Estimates_corr %>%
  group_by(Parameter) %>%
  summarise(N = n(),
            r = cor.test(Effort, Wanting)$estimate,
            p = cor.test(Effort, Wanting)$p.value,
            .groups = "drop")
print(effort_wanting_r)
write.csv(effort_wanting_r, file.path(output_dir, "effort_wanting_correlations.csv"),
          row.names = FALSE)

effort_questionnaire_r <- QUEST_Estimates_long %>%
  pivot_longer(c(BAS, BIS, LARS, SHAPS), names_to = "Questionnaire", values_to = "Score") %>%
  group_by(Parameter, Questionnaire) %>%
  summarise(N = sum(complete.cases(Estimate, Score)),
            r = cor.test(Estimate, Score)$estimate,
            p = cor.test(Estimate, Score)$p.value,
            .groups = "drop")
print(effort_questionnaire_r, n = Inf)
write.csv(effort_questionnaire_r, file.path(output_dir, "effort_questionnaire_correlations.csv"),
          row.names = FALSE)



##############################################################################
# C  SUPPLEMENTARY MATERIAL
##############################################################################


# Table S1: behavioral reliability with and without S1 -----------------------

table_S1 <- tibble(
  Parameter     = names(parameters),
  `ICC (S1-S5)` = sprintf("%.3f", icc_S1_S5),
  `ICC (S2-S5)` = sprintf("%.3f", sapply(parameters, function(p) {
    get_icc(coef_wide, paste0(p, "_S", 2:5))
  }))
)

print(table_S1)
cat("Participants included:", nrow(coef_wide), "\n")
write.csv(table_S1, file.path(output_dir, "TableS1_behavioral_reliability.csv"),
          row.names = FALSE)


# Figure S1: sex x behavior interactions -------------------------------------


stopifnot(!anyNA(mri_sample$Female), all(mri_sample$Female %in% c(0, 1)))
fc_sex <- mri_sample %>%
  mutate(Sex = factor(if_else(Female == 1, "Women", "Men"), levels = c("Men", "Women")))

run_interaction <- function(y, x, x_raw, other_behaviors, label, xlab, ylab) {

  nuisance <- c("cAge", "cBMI", "FD", other_behaviors)
  d <- fc_sex %>% select(all_of(c(y, x, x_raw, nuisance, "Sex"))) %>% drop_na()

  fit  <- lm(reformulate(c(paste0(x, " * Sex"), nuisance), response = y), data = d)
  term <- paste0(x, ":SexWomen")
  ci   <- confint(fit)[term, ]

  interaction <- tibble(
    Finding  = label,
    Beta_Women_minus_Men = coef(fit)[[term]],
    CI95_lower = ci[[1]],
    CI95_upper = ci[[2]],
    p_value  = summary(fit)$coefficients[term, "Pr(>|t|)"]
  )
  print(summary(fit))

  scale_fit <- lm(reformulate(x_raw, response = x), data = d)
  raw_to_model <- function(v) unname(coef(scale_fit)[1] + coef(scale_fit)[2] * v)

  at_means <- function(newdata) {
    for (v in nuisance) newdata[[v]] <- mean(d[[v]])
    newdata
  }

  # x: remove covariate-related variance, keeping sex-specific group means.
  x_fit <- lm(reformulate(c("Sex", nuisance), response = x_raw), data = d)
  d$x_adjusted <- d[[x_raw]] - predict(x_fit, d) + predict(x_fit, at_means(d))

  # y: predicted at the adjusted x and mean covariates, plus the residual.
  d_adj <- at_means(d)
  d_adj[[x]] <- raw_to_model(d$x_adjusted)
  d$y_adjusted <- d[[y]] - predict(fit, d) + predict(fit, d_adj)

  line_data <- d %>%
    group_by(Sex) %>%
    reframe(x_adjusted = seq(min(x_adjusted), max(x_adjusted), length.out = 100)) %>%
    mutate(!!x := raw_to_model(x_adjusted)) %>%
    at_means()
  line_data <- bind_cols(line_data,
                         as.data.frame(predict(fit, line_data, interval = "confidence")))

  cols <- c("Men" = "#7F9779", "Women" = "#D67C1C")
  plot <- ggplot(d, aes(x = x_adjusted, y = y_adjusted, color = Sex)) +
    geom_point(size = 2.2, alpha = 0.8) +
    geom_ribbon(data = line_data,
                aes(x = x_adjusted, ymin = lwr, ymax = upr, fill = Sex),
                inherit.aes = FALSE, alpha = 0.12) +
    geom_line(data = line_data, aes(x = x_adjusted, y = fit, color = Sex),
              linewidth = 0.95) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(title = label, x = xlab, y = ylab, color = "Sex", fill = "Sex") +
    theme_classic(base_size = 17) +
    theme(axis.text    = element_text(color = "black", size = 15),
          axis.title   = element_text(size = 17),
          plot.title   = element_text(face = "bold", size = 18),
          legend.title = element_text(size = 16),
          legend.text  = element_text(size = 15),
          legend.position = "bottom")

  list(plot = plot, interaction = interaction)
}

ifg <- run_interaction(
  y = "FC_rIFGoper", x = "Effort_intercept", x_raw = "Intercept_raw",
  other_behaviors = c("Effort_rewmag", "Effort_diff"),
  label = "A  Intercept - right IFG",
  xlab = "Intercept (adjusted)", ylab = "NAcc-right IFG FC"
)

subcal <- run_interaction(
  y = "FC_subcallosal", x = "Effort_diff", x_raw = "Difficulty_raw",
  other_behaviors = c("Effort_rewmag", "Effort_intercept"),
  label = "B  Difficulty - subcallosal gyrus",
  xlab = "Difficulty (adjusted)", ylab = "NAcc-subcallosal FC"
)

interaction_table <- bind_rows(ifg$interaction, subcal$interaction)
print(interaction_table)
write.csv(interaction_table, file.path(output_dir, "FigureS1_sex_interactions.csv"),
          row.names = FALSE)

FigureS1 <- (ifg$plot | subcal$plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
save_png(FigureS1, "FigureS1_sex_brain_behavior", width = 11, height = 4.7, dpi = 400)


# Table S2: timing sensitivity of the focal associations ---------------------

# Session-specific partial correlations (Fig. 5B) repeated with behavioral visit
# start time and the mean absolute interval (days) to the two MRI visits.

if (has_timing) {

  timing <- session_timing %>%
    transmute(ID, Session, start_hour,
              gap_days = (abs(days_to_S6) + abs(days_to_S7)) / 2)

  timing_data <- coef_long %>%
    rename(Intercept = `(Intercept)`, RewMag = cRewMag, Difficulty = cDifficulty) %>%
    inner_join(mri_sample %>% select(ID, FC_rIFGoper, FC_subcallosal, cAge, cSex, cBMI, FD),
               by = "ID") %>%
    left_join(timing, by = c("ID", "Session"))

  partial_cor <- function(d, y, x, covars) {
    e <- qr.resid(qr(model.matrix(reformulate(covars), d)), as.matrix(d[c(y, x)]))
    cor(e[, 1], e[, 2])
  }

  timing_r <- function(d, roi, behavior) {
    controls <- c("cAge", "cSex", "cBMI", "FD", "RewMag",
                  setdiff(c("Intercept", "Difficulty"), behavior))
    timing_vars <- c("start_hour", "gap_days")
    full    <- drop_na(d, all_of(c(roi, behavior, controls)))
    matched <- drop_na(full, all_of(timing_vars))
    tibble(ROI = roi, Behavior = behavior,
           N_full = nrow(full), N = nrow(matched),
           r_full     = partial_cor(full, roi, behavior, controls),
           r_original = partial_cor(matched, roi, behavior, controls),
           r_adjusted = partial_cor(matched, roi, behavior, c(controls, timing_vars)),
           delta_r    = r_adjusted - r_original)
  }

  table_S2 <- timing_data %>%
    group_by(Session) %>%
    group_modify(~ bind_rows(timing_r(.x, "FC_rIFGoper",    "Intercept"),
                             timing_r(.x, "FC_subcallosal", "Difficulty"))) %>%
    ungroup() %>%
    arrange(ROI, Session)


  # The baseline should reproduce the partial r of Figure 5B (MATLAB).
  fig5_reference <- roi_long %>%
    filter((ROI_short == IFG_ROI & Behavior == "Intercept") |
             (ROI_short == SUBCAL_ROI & Behavior == "Difficulty")) %>%
    transmute(Session = as.integer(str_remove(Session, "S")),
              ROI = if_else(ROI_short == IFG_ROI, "FC_rIFGoper", "FC_subcallosal"),
              r_fig5 = partial_r)

  table_S2 <- left_join(table_S2, fig5_reference, by = c("Session", "ROI"))
 

  print(table_S2, n = Inf)
  write.csv(table_S2, file.path(output_dir, "TableS2_timing_sensitivity.csv"), row.names = FALSE)

}


# Figure S2: focal ROIs, all three behavioral parameters ---------------------

region_levels <- c("Right IFG", "Subcallosal cortex")

figS2_background <- bind_rows(
  roi_long %>% mutate(Region = region_levels[1]),
  roi_long %>% mutate(Region = region_levels[2])
) %>%
  mutate(Region = factor(Region, levels = region_levels))

figS2_data <- bind_rows(
  roi_long %>% filter(ROI_short == IFG_ROI)    %>% mutate(Region = region_levels[1]),
  roi_long %>% filter(ROI_short == SUBCAL_ROI) %>% mutate(Region = region_levels[2])
) %>%
  mutate(Region = factor(Region, levels = region_levels))

FigureS2 <- ggplot() +
  geom_line(data = figS2_background, aes(Session, partial_r, group = ROI),
            color = "grey65", alpha = 0.20, linewidth = 0.35) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(data = figS2_data, aes(Session, partial_r, group = 1, color = Behavior),
            linewidth = 1.5) +
  geom_point(data = figS2_data, aes(Session, partial_r, color = Behavior), size = 3) +
  facet_grid(Region ~ Behavior) +
  scale_color_manual(values = behavior_colors) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  labs(x = "Behavioral session", y = partial_r_label) +
  theme_paper +
  theme(legend.position = "none", strip.text = element_text(size = 13))
save_png(FigureS2, "FigureS2_focal_ROI_behavioral_specificity", width = 12, height = 8)


# Figure S3: ROI FC reliability vs brain-behavior associations ---------------

# A: ROI-wise ICC vs partial r

figS3a_data <- roi_long %>% filter(is.finite(ICC_con), is.finite(partial_r))

figS3a_labels <- figS3a_data %>%
  group_by(Behavior, Session) %>%
  summarise(Pearson_r = cor(ICC_con, partial_r), .groups = "drop") %>%
  mutate(label = sprintf("r = %.2f", Pearson_r),
         x = min(figS3a_data$ICC_con),
         y = 0.35)
print(figS3a_labels)

figS3a_extremes <- bind_rows(
  figS3a_data %>% group_by(Behavior, Session) %>%
    slice_min(partial_r, n = 5, with_ties = FALSE) %>% ungroup() %>%
    mutate(Group = "5 most negative associations"),
  figS3a_data %>% group_by(Behavior, Session) %>%
    slice_max(partial_r, n = 5, with_ties = FALSE) %>% ungroup() %>%
    mutate(Group = "5 most positive associations")
)

FigureS3A <- ggplot(figS3a_data, aes(ICC_con, partial_r)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey75") +
  geom_point(color = "grey75", size = 1.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              color = "grey45", linewidth = 0.7) +
  geom_point(data = figS3a_extremes, aes(color = Group), size = 2.3) +
  geom_text_repel(data = figS3a_extremes, aes(label = ROI_short, color = Group),
                  size = 2.5, max.overlaps = Inf, min.segment.length = 0,
                  show.legend = FALSE) +
  geom_text(data = figS3a_labels, aes(x, y, label = label), inherit.aes = FALSE,
            hjust = 0, fontface = "bold", size = 3.2) +
  facet_grid(Behavior ~ Session) +
  scale_color_manual(values = c("5 most negative associations" = "#2166AC",
                                "5 most positive associations" = "#D73027")) +
  coord_cartesian(ylim = c(-0.4, 0.4), clip = "off") +
  labs(x = "ROI NAcc FC reliability (ICC)", y = partial_r_label, color = NULL) +
  theme_paper +
  theme(legend.position = "bottom",
        axis.text     = element_text(size = 8),
        strip.text    = element_text(size = 10),
        panel.spacing = unit(0.6, "lines"))
save_png(FigureS3A, "FigureS3A_FC_reliability_vs_brain_behavior", width = 16, height = 9)


# B: network-wise relationship

network_colors <- c(
  "visual"            = "#781286",
  "somato-motor"      = "#4682B4",
  "dorsal attention"  = "#00760E",
  "ventral attention" = "#C43AFA",
  "limbic"            = "#A9C95C",
  "frontoparietal"    = "#E69422",
  "default mode"      = "#CD3E4E",
  "diencephalon"      = "#4D4D4D",
  "midbrain"          = "#00A6D6",
  "brainstem"         = "#8C564B"
)
network_order <- names(network_colors)

network_n <- roi_long %>%
  distinct(ROI, Network) %>%
  filter(!is.na(Network), Network != "") %>%
  count(Network, name = "N")

network_labels <- setNames(
  paste0(network_order, " (n=", network_n$N[match(network_order, network_n$Network)], ")"),
  network_order
)

figS3b_data <- roi_long %>%
  filter(!is.na(Network), Network != "", is.finite(ICC_con), is.finite(partial_r)) %>%
  group_by(Behavior, Network, Session) %>%
  summarise(N_ROIs = n(),
            Pearson_r = if (n() >= 4) cor(ICC_con, partial_r) else NA_real_,
            .groups = "drop") %>%
  mutate(Network = factor(Network, levels = network_order))

FigureS3B <- ggplot(figS3b_data,
                    aes(Session, Pearson_r, group = Network, color = Network)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  facet_wrap(~Behavior, nrow = 1, axes = "all", axis.labels = "all") +
  scale_color_manual(values = network_colors, breaks = network_order,
                     labels = network_labels, drop = TRUE) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  labs(x = "Behavioral session", y = "Correlation with ROI FC reliability",
       color = "Network / system") +
  theme_paper +
  theme(legend.position = "right",
        legend.text = element_text(size = 8),
        strip.text  = element_text(size = 11))
save_png(FigureS3B, "FigureS3B_network_trajectories", width = 12, height = 5.5)


# Figure S4: effort vs wanting estimates -------------------------------------

figS4_labels <- cor_labels(Estimates_corr, "Effort", "Wanting", "Parameter")

FigureS4 <- ggplot(Estimates_corr, aes(x = Effort, y = Wanting)) +
  geom_point(size = 4.5) +
  geom_text(data = figS4_labels, aes(label = label), x = -Inf, y = Inf,
            hjust = -0.1, vjust = 1.4, size = 5, parse = TRUE, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  geom_smooth(method = "lm", formula = y ~ x, linewidth = 0.3,
              alpha = 0.5, color = "black") +
  theme_scatter +
  labs(x = "Individual estimate (Effort)", y = "Individual estimate (Wanting)") +
  facet_wrap(~Parameter, scales = "free", ncol = 3)

save_png(FigureS4, "FigureS4_effort_wanting", width = 10, height = 3.75)


# Figure S5: correlations with questionnaires --------------------------------

var_order  <- c("Intercept", "RewMag", "RewDiff", "BAS", "BIS", "LARS", "SHAPS")
var_labels <- c(Intercept = "Intercept", RewMag = "Reward magnitude",
                RewDiff = "Difficulty", BAS = "BAS", BIS = "BIS",
                LARS = "LARS", SHAPS = "SHAPS")

make_corr_heatmap <- function(data, title, show_legend = TRUE) {
  corr_long <- round(cor(data[, var_order], use = "pairwise.complete.obs"), 2) %>%
    reshape2::melt() %>%
    mutate(Var1 = factor(Var1, levels = var_order),
           Var2 = factor(Var2, levels = var_order)) %>%
    filter(as.integer(Var1) <= as.integer(Var2))

  ggplot(corr_long, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = value), size = 4) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1)) +
    scale_x_discrete(labels = var_labels) +
    scale_y_discrete(limits = rev(var_order), labels = var_labels) +
    labs(title = title, x = "", y = "", fill = "r") +
    theme(
      legend.position = if (show_legend) "right" else "none",
      text        = element_text(face = "bold", size = 16),
      axis.text   = element_text(face = "plain", size = 12),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1)
    )
}

FigureS5 <- cowplot::plot_grid(
  make_corr_heatmap(QUEST_Estimates,   "Effort estimates",  show_legend = FALSE),
  make_corr_heatmap(QUEST_Estimates_W, "Wanting estimates", show_legend = TRUE),
  ncol = 2, rel_widths = c(1, 1.2)
)
save_png(FigureS5, "FigureS5_correlation_matrices", width = 12, height = 5)


# Figure S6: effort estimates vs BIS/BAS -------------------------------------

plot_bisbas <- function(scale_name, xlab) {
  ggplot(QUEST_Estimates_long, aes(x = Estimate, y = .data[[scale_name]])) +
    geom_point(size = 4.5) +
    geom_text(data = cor_labels(QUEST_Estimates_long, "Estimate", scale_name, "Parameter"),
              aes(label = label), x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.6,
              size = 5, parse = TRUE, inherit.aes = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.05))) +
    geom_smooth(method = "lm", formula = y ~ x, linewidth = 0.3,
                alpha = 0.5, color = "black") +
    theme_scatter +
    labs(x = xlab, y = scale_name) +
    facet_wrap(~Parameter, scales = "free_x", ncol = 3)
}

FigureS6 <- plot_bisbas("BAS", "") /
  plot_bisbas("BIS", "Individual estimate (Effort)")
save_png(FigureS6, "FigureS6_BIS_BAS_effort", width = 10, height = 8)


writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
