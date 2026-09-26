########################################################################
# T2 biomarkers and CRSwNP - Part C - Prognostic value of T2 biomarkers
########################################################################


########################################################################
# SECTION 1 — Packages
########################################################################

library(tidyverse)
library(readxl)
library(metafor)
library(gtsummary)
library(writexl)
library(broom)
library(patchwork)
library(cowplot)

user_name <- Sys.info()[["user"]]

# Load the MIAnalysis package (provides IPD_two_stage)
MIAnalysis_path <- file.path(
  "/Users", user_name,
  "Library/CloudStorage/OneDrive-USherbrooke/Recherche/R_packages/MIAnalysis/R"
)
devtools::load_all(MIAnalysis_path)

########################################################################
# SECTION 2 — Paths
########################################################################

Loading_path <- file.path(
  "/Users", user_name,
  "Library/CloudStorage/OneDrive-USherbrooke/Recherche/Projets/Pneumologie/ORACLE/ORACLE_NP_IgE/Data/After_part_A"
)

Saving_path <- file.path(
  "/Users", user_name,
  "Library/CloudStorage/OneDrive-USherbrooke/Recherche/Projets/Pneumologie/ORACLE/ORACLE_NP_IgE/Data/After_part_C"
)
dir.create(Saving_path, showWarnings = FALSE, recursive = TRUE)

Figure_path <- file.path(
  "/Users", user_name,
  "Library/CloudStorage/OneDrive-USherbrooke/Recherche/Projets/Pneumologie/ORACLE/ORACLE_NP_IgE/Figures_and_Tables/Part_C"
)
dir.create(Figure_path, showWarnings = FALSE, recursive = TRUE)

########################################################################
# SECTION 3 — Loading the dataset
########################################################################

load(file.path(Loading_path, "DATA_Raw.RData"))      # DATA_Raw
load(file.path(Loading_path, "DATA_imputed.RData"))  # DATA_imputed

########################################################################
# SECTION 4 — Trial selection
#
# Part A already excluded trials with no CRSwNP data collected and trials
# where all subjects have CRSwNP = No, so all remaining trials are included.
# Trials with few CRSwNP cases may produce extreme interaction estimates;
# these are handled by x_limits in IPD_two_stage (over-range CIs shown as arrows).
########################################################################

qualifying_trials <- unique(as.character(DATA_Raw$Enrolled_Trial_name))

cat("Trials included (n =", length(qualifying_trials), "):",
    paste(qualifying_trials, collapse = ", "), "\n")



########################################################################
# SECTION 5 — Two-stage IPD meta-analysis mirroring the one-stage model
#
# Outcome : Attack_number_during_followup (count, negative binomial)
# Offset  : Follow_up_duration_days
# Model   : count ~ biomarker_log10 * CRSwNP
#           + GINA_step_numeric + ACQ_score_0W + Attack_history_num
#           + FEV1_preBD_PCT_0W + FeNO_log10 + BEC_log10
#           + offset(log(Follow_up_duration_days))
# Strategy: IPD_two_stage() — fits within each trial × imputation,
#           Rubin-pools within trials, REML meta-analyzes across trials.
########################################################################

# ── 5.1  Prepare data (log10-scale biomarkers, all trials) ───────────
data_s5 <- DATA_imputed |>
  dplyr::filter(as.character(Enrolled_Trial_name) %in% qualifying_trials,
                .imp >= 1) |>
  dplyr::mutate(
    FeNO_log10 = log10(FeNO),
    BEC_log10  = log10(BEC),
    IgE_log10  = log10(IgE),
    CRSwNP     = factor(as.character(CRSwNP), levels = c("No", "Yes"))
  )


# IQR of each biomarker (imp == 1 as reference; used for per-IQR scaling)
imp1_iqr      <- dplyr::filter(data_s5, .imp == 1)
iqr_FeNO_orig <- IQR(imp1_iqr$FeNO,       na.rm = TRUE)
iqr_BEC_orig  <- IQR(imp1_iqr$BEC,        na.rm = TRUE)
iqr_IgE_orig  <- IQR(imp1_iqr$IgE,        na.rm = TRUE)
iqr_FeNO_log  <- IQR(imp1_iqr$FeNO_log10, na.rm = TRUE)
iqr_BEC_log   <- IQR(imp1_iqr$BEC_log10,  na.rm = TRUE)
iqr_IgE_log   <- IQR(imp1_iqr$IgE_log10,  na.rm = TRUE)
covars_s5 <- c(
  "GINA_step_numeric", "ACQ_score_0W", "Attack_history_num",
  "FEV1_preBD_PCT_0W", "FeNO_log10", "BEC_log10", "CRSwNP", "IgE_log10"
)

# ── 5.2  FeNO × CRSwNP (NB, count outcome) ───────────────────────────
covars_FeNO <- c(
  "GINA_step_numeric", "ACQ_score_0W", "Attack_history_num",
  "FEV1_preBD_PCT_0W", "BEC_log10", "IgE_log10", "CRSwNP", "FeNO_log10"
)

int2_FeNO <- IPD_two_stage(
  data            = data_s5,
  outcome_var     = "Attack_number_during_followup",
  predictor_vars  = c("FeNO_log10 * CRSwNP"),
  covariables     = covars_FeNO[covars_FeNO != "FeNO_log10"],
  imp_col         = ".imp",
  intercept_var   = "Enrolled_Trial_name",
  followup_offset = "Yes",
  followup_col    = "Follow_up_duration_days",
  model_type      = "nb",
  meta_method     = "DL",
  main_term       = "FeNO_log10:CRSwNPYes",
  xlab            = "Interaction term: log10(FeNO) \u00d7 CRSwNP",
  ref_line        = 1,
  x_limits        = c(0.1, 10),
  x_breaks        = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
  diamond_color   = "steelblue",
  pooled_color    = "steelblue",
  het_y           = 0.02,
  het_line_gap    = 0.025
)

cat("\n── FeNO × CRSwNP (NB, two-stage) ───────────────────────────────\n")
print(int2_FeNO$table |> dplyr::filter(grepl(":", term)))
int2_FeNO$forest_plot

ggsave(file.path(Figure_path, "Forest_int2_FeNO_CRSwNP.pdf"),
       plot = int2_FeNO$forest_plot, width = 12, height = 7,
       units = "in", device = "pdf")

# ── 5.3  BEC × CRSwNP (NB) ───────────────────────────────────────────
covars_BEC <- c(
  "GINA_step_numeric", "ACQ_score_0W", "Attack_history_num",
  "FEV1_preBD_PCT_0W", "FeNO_log10", "IgE_log10", "CRSwNP", "BEC_log10"
)

int2_BEC <- IPD_two_stage(
  data            = data_s5,
  outcome_var     = "Attack_number_during_followup",
  predictor_vars  = c("BEC_log10 * CRSwNP"),
  covariables     = covars_BEC[covars_BEC != "BEC_log10"],
  imp_col         = ".imp",
  intercept_var   = "Enrolled_Trial_name",
  followup_offset = "Yes",
  followup_col    = "Follow_up_duration_days",
  model_type      = "nb",
  meta_method     = "DL",
  main_term       = "BEC_log10:CRSwNPYes",
  xlab            = "Interaction term: log10(BEC) \u00d7 CRSwNP",
  ref_line        = 1,
  x_limits        = c(0.1, 10),
  x_breaks        = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
  diamond_color   = "darkviolet",
  pooled_color    = "darkviolet",
  het_y           = 0.02,
  het_line_gap    = 0.025
)

cat("\n── BEC × CRSwNP (NB, two-stage) ────────────────────────────────\n")
print(int2_BEC$table |> dplyr::filter(grepl(":", term)))
int2_BEC$forest_plot

ggsave(file.path(Figure_path, "Forest_int2_BEC_CRSwNP.pdf"),
       plot = int2_BEC$forest_plot, width = 12, height = 7,
       units = "in", device = "pdf")

# ── 5.4  IgE × CRSwNP (NB) ───────────────────────────────────────────
covars_IgE <- c(
  "GINA_step_numeric", "ACQ_score_0W", "Attack_history_num",
  "FEV1_preBD_PCT_0W", "FeNO_log10", "BEC_log10", "CRSwNP", "IgE_log10"
)

int2_IgE <- IPD_two_stage(
  data            = data_s5,
  outcome_var     = "Attack_number_during_followup",
  predictor_vars  = c("IgE_log10 * CRSwNP"),
  covariables     = covars_IgE[covars_IgE != "IgE_log10"],
  imp_col         = ".imp",
  intercept_var   = "Enrolled_Trial_name",
  followup_offset = "Yes",
  followup_col    = "Follow_up_duration_days",
  model_type      = "nb",
  meta_method     = "DL",
  main_term       = "IgE_log10:CRSwNPYes",
  xlab            = "Interaction term: log10(IgE) \u00d7 CRSwNP",
  ref_line        = 1,
  x_limits        = c(0.1, 10),
  x_breaks        = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
  diamond_color   = "darkorange",
  pooled_color    = "darkorange",
  het_y           = 0.02,
  het_line_gap    = 0.025
)

cat("\n── IgE × CRSwNP (NB, two-stage) ────────────────────────────────\n")
print(int2_IgE$table |> dplyr::filter(grepl(":", term)))
int2_IgE$forest_plot

ggsave(file.path(Figure_path, "Forest_int2_IgE_CRSwNP.pdf"),
       plot = int2_IgE$forest_plot, width = 12, height = 7,
       units = "in", device = "pdf")

# ── 5.4b  IgE × CRSwNP — per-IQR alternative ─────────────────────────────────
# Predictor rescaled: IgE_log10_iqr = IgE_log10 / IQR(IgE_log10)
# => 1 unit = 1 IQR of log10(IgE) = IQR of original IgE in kU/L

data_s5_iqr <- data_s5 |>
  dplyr::mutate(IgE_log10_iqr = IgE_log10 / iqr_IgE_log)

covars_IgE_base <- c(
  "GINA_step_numeric", "ACQ_score_0W", "Attack_history_num",
  "FEV1_preBD_PCT_0W", "FeNO_log10", "BEC_log10", "CRSwNP"
)

int2_IgE_iqr <- IPD_two_stage(
  data            = data_s5_iqr,
  outcome_var     = "Attack_number_during_followup",
  predictor_vars  = c("IgE_log10_iqr * CRSwNP"),
  covariables     = covars_IgE_base,
  imp_col         = ".imp",
  intercept_var   = "Enrolled_Trial_name",
  followup_offset = "Yes",
  followup_col    = "Follow_up_duration_days",
  model_type      = "nb",
  meta_method     = "DL",
  main_term       = "IgE_log10_iqr:CRSwNPYes",
  xlab            = sprintf("Interaction term: IgE (per %g kU/L) \u00d7 CRSwNP",
                            round(iqr_IgE_orig)),
  ref_line        = 1,
  x_limits        = c(0.1, 10),
  x_breaks        = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
  diamond_color   = "darkorange",
  pooled_color    = "darkorange",
  het_y           = 0.02,
  het_line_gap    = 0.025
)

int2_IgE_iqr$forest_plot

ggsave(file.path(Figure_path, "Forest_int2_IgE_perIQR_CRSwNP.pdf"),
       plot = int2_IgE_iqr$forest_plot, width = 12, height = 7,
       units = "in", device = "pdf")
# ── 5.5  Summary ─────────────────────────────────────────────────────
int2_summary <- dplyr::bind_rows(
  int2_FeNO$table |> dplyr::filter(grepl(":", term)) |> dplyr::mutate(biomarker = "FeNO (log10)"),
  int2_BEC$table  |> dplyr::filter(grepl(":", term)) |> dplyr::mutate(biomarker = "BEC (log10)"),
  int2_IgE$table  |> dplyr::filter(grepl(":", term)) |> dplyr::mutate(biomarker = "IgE (log10)")
) |>
  dplyr::mutate(
    IRR_CI  = sprintf("%.2f (%.2f - %.2f)", pooled_exp, ci_exp_low, ci_exp_high),
    p_text  = dplyr::case_when(p_value < 0.001 ~ "<0.001",
                               p_value < 0.01  ~ sprintf("%.3f", p_value),
                               TRUE            ~ sprintf("%.2f",  p_value)),
    I2_text = sprintf("%.1f%%", I2)
  ) |>
  dplyr::select(biomarker, IRR_CI, p_text, I2_text, k)

cat("\n── Section 5 pooled IRR (NB, fully adjusted, two-stage) ─────────\n")
print(int2_summary)


########################################################################
# SECTION 6 — Subgroup analysis: prognostic value of T2 biomarkers
#             separately within CRSwNP = No and CRSwNP = Yes
#
# For each biomarker (FeNO, BEC, IgE), fit IPD_two_stage within each
# CRSwNP stratum (NB count outcome, follow-up offset, fully adjusted).
# Forest plots are shown side by side: CRSwNP No | CRSwNP Yes.
########################################################################

# ── 6.1  Helper: run one stratified IPD_two_stage ────────────────────
run_subgroup <- function(biomarker_var, crswNP_level, colour,
                         covars_base, data = data_s5) {

  d_sub <- data |> dplyr::filter(CRSwNP == crswNP_level)

  # Covariates: all except CRSwNP and the biomarker of interest
  covars <- covars_base[!covars_base %in% c("CRSwNP", biomarker_var)]

  IPD_two_stage(
    data            = d_sub,
    outcome_var     = "Attack_number_during_followup",
    predictor_vars  = biomarker_var,
    covariables     = covars,
    imp_col         = ".imp",
    intercept_var   = "Enrolled_Trial_name",
    followup_offset = "Yes",
    followup_col    = "Follow_up_duration_days",
    model_type      = "nb",
    meta_method     = "DL",
    main_term       = biomarker_var,
    xlab            = "aRR",
    ref_line        = 1,
    x_limits        = c(0.1, 10),
    x_breaks        = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
    diamond_color   = colour,
    pooled_color    = colour,
    show_headers     = TRUE,
    bold_pooled      = TRUE,
    het_y           = 0.02,
    het_line_gap    = 0.025
  )
}

covars_base <- c(
  "GINA_step_numeric", "ACQ_score_0W", "Attack_history_num",
  "FEV1_preBD_PCT_0W", "FeNO_log10", "BEC_log10", "IgE_log10", "CRSwNP"
)

# ── 6.2  FeNO subgroup ────────────────────────────────────────────────
sub_FeNO_no  <- run_subgroup("FeNO_log10", "No",  "steelblue",  covars_base)
sub_FeNO_yes <- run_subgroup("FeNO_log10", "Yes", "steelblue",  covars_base)

cat("\nFeNO — CRSwNP No :\n");  print(sub_FeNO_no$table  |> dplyr::filter(term == "FeNO_log10"))
cat("\nFeNO — CRSwNP Yes:\n");  print(sub_FeNO_yes$table |> dplyr::filter(term == "FeNO_log10"))

panel_FeNO <- cowplot::plot_grid(
  cowplot::ggdraw() + cowplot::draw_label("CRSwNP = No",  fontface = "bold", size = 13),
  cowplot::ggdraw() + cowplot::draw_label("CRSwNP = Yes", fontface = "bold", size = 13),
  sub_FeNO_no$forest_plot,
  sub_FeNO_yes$forest_plot,
  ncol = 2, rel_heights = c(0.08, 1)
)

ggsave(file.path(Figure_path, "Subgroup_FeNO_CRSwNP.pdf"),
       plot = panel_FeNO, width = 18, height = 7, units = "in", device = "pdf")

# ── 6.3  BEC subgroup ────────────────────────────────────────────────
sub_BEC_no  <- run_subgroup("BEC_log10", "No",  "darkviolet", covars_base)
sub_BEC_yes <- run_subgroup("BEC_log10", "Yes", "darkviolet", covars_base)

cat("\nBEC — CRSwNP No :\n");  print(sub_BEC_no$table  |> dplyr::filter(term == "BEC_log10"))
cat("\nBEC — CRSwNP Yes:\n");  print(sub_BEC_yes$table |> dplyr::filter(term == "BEC_log10"))

panel_BEC <- cowplot::plot_grid(
  cowplot::ggdraw() + cowplot::draw_label("CRSwNP = No",  fontface = "bold", size = 13),
  cowplot::ggdraw() + cowplot::draw_label("CRSwNP = Yes", fontface = "bold", size = 13),
  sub_BEC_no$forest_plot,
  sub_BEC_yes$forest_plot,
  ncol = 2, rel_heights = c(0.08, 1)
)

ggsave(file.path(Figure_path, "Subgroup_BEC_CRSwNP.pdf"),
       plot = panel_BEC, width = 18, height = 7, units = "in", device = "pdf")

# ── 6.4  IgE subgroup ────────────────────────────────────────────────
sub_IgE_no  <- run_subgroup("IgE_log10", "No",  "darkorange", covars_base)
sub_IgE_yes <- run_subgroup("IgE_log10", "Yes", "darkorange", covars_base)

cat("\nIgE — CRSwNP No :\n");  print(sub_IgE_no$table  |> dplyr::filter(term == "IgE_log10"))
cat("\nIgE — CRSwNP Yes:\n");  print(sub_IgE_yes$table |> dplyr::filter(term == "IgE_log10"))

panel_IgE <- cowplot::plot_grid(
  cowplot::ggdraw() + cowplot::draw_label("CRSwNP = No",  fontface = "bold", size = 13),
  cowplot::ggdraw() + cowplot::draw_label("CRSwNP = Yes", fontface = "bold", size = 13),
  sub_IgE_no$forest_plot,
  sub_IgE_yes$forest_plot,
  ncol = 2, rel_heights = c(0.08, 1)
)

ggsave(file.path(Figure_path, "Subgroup_IgE_CRSwNP.pdf"),
       plot = panel_IgE, width = 18, height = 7, units = "in", device = "pdf")

# ── 6.5  Summary table ────────────────────────────────────────────────
extract_row <- function(obj, biomarker, var_name, crswNP_level) {
  obj$table |>
    dplyr::filter(term == var_name) |>
    dplyr::mutate(
      biomarker  = biomarker,
      CRSwNP     = crswNP_level,
      IRR_CI     = sprintf("%.2f (%.2f - %.2f)", pooled_exp, ci_exp_low, ci_exp_high),
      p_text     = dplyr::case_when(p_value < 0.001 ~ "<0.001",
                                    p_value < 0.01  ~ sprintf("%.3f", p_value),
                                    TRUE            ~ sprintf("%.2f", p_value)),
      I2_text    = sprintf("%.1f%%", I2)
    ) |>
    dplyr::select(biomarker, CRSwNP, IRR_CI, p_text, I2_text, k)
}

subgroup_summary <- dplyr::bind_rows(
  extract_row(sub_FeNO_no,  "FeNO (log10)", "FeNO_log10", "No"),
  extract_row(sub_FeNO_yes, "FeNO (log10)", "FeNO_log10", "Yes"),
  extract_row(sub_BEC_no,   "BEC (log10)",  "BEC_log10",  "No"),
  extract_row(sub_BEC_yes,  "BEC (log10)",  "BEC_log10",  "Yes"),
  extract_row(sub_IgE_no,   "IgE (log10)",  "IgE_log10",  "No"),
  extract_row(sub_IgE_yes,  "IgE (log10)",  "IgE_log10",  "Yes")
)

cat("\n── Section 6: prognostic IRR by CRSwNP subgroup ─────────────────\n")
print(subgroup_summary)

# ── 6.6  Summary forest plot: all 6 biomarker × subgroup combinations ─
df_plot_s <- dplyr::bind_rows(
  dplyr::mutate(sub_FeNO_no$table  |> dplyr::filter(term=="FeNO_log10"),
                biomarker="FeNO (log10)", subgroup="Without CRSwNP", colour="steelblue"),
  dplyr::mutate(sub_FeNO_yes$table |> dplyr::filter(term=="FeNO_log10"),
                biomarker="FeNO (log10)", subgroup="With CRSwNP",    colour="steelblue"),
  dplyr::mutate(sub_BEC_no$table   |> dplyr::filter(term=="BEC_log10"),
                biomarker="BEC (log10)",  subgroup="Without CRSwNP", colour="darkviolet"),
  dplyr::mutate(sub_BEC_yes$table  |> dplyr::filter(term=="BEC_log10"),
                biomarker="BEC (log10)",  subgroup="With CRSwNP",    colour="darkviolet"),
  dplyr::mutate(sub_IgE_no$table   |> dplyr::filter(term=="IgE_log10"),
                biomarker="IgE (log10)",  subgroup="Without CRSwNP", colour="darkorange"),
  dplyr::mutate(sub_IgE_yes$table  |> dplyr::filter(term=="IgE_log10"),
                biomarker="IgE (log10)",  subgroup="With CRSwNP",    colour="darkorange")
) |>
  dplyr::mutate(
    y         = c(8,7, 5,4, 2,1),
    row_label = paste0("    ", subgroup),
    est       = pooled_exp,
    lo        = ci_exp_low,
    hi        = ci_exp_high,
    sig       = p_value < 0.05,
    est_txt   = sprintf("%.2f (%.2f - %.2f)", pooled_exp, ci_exp_low, ci_exp_high),
    p_txt     = dplyr::case_when(p_value < 0.001 ~ "<0.001",
                                 p_value < 0.01  ~ sprintf("%.3f", p_value),
                                 TRUE            ~ sprintf("%.2f", p_value))
  )

new_ylim <- c(0.3, 9.5)

hdr_df <- tibble::tibble(y = c(9,6,3), label = c("FeNO (log10)","BEC (log10)","IgE (log10)"))

df_left <- data.frame(y   = c(hdr_df$y,  df_plot_s$y),
                      lab = c(hdr_df$label, df_plot_s$row_label),
                      face= c(rep("bold", 3), rep("plain", 6)))

p_left <- ggplot(df_left, aes(y=y, x=1, label=lab)) +
  geom_text(hjust=1, size=3.5, fontface=df_left$face, colour="black") +
  scale_x_continuous(limits=c(0,1), expand=c(0.01,0)) +
  scale_y_continuous(limits=new_ylim, breaks=NULL) +
  labs(x=NULL, y=NULL, title=NULL) +   # "Subgroup" header removed
  theme_void(base_size=12) +
  theme(plot.title  = element_text(hjust=1, face="bold", size=11,
                                   margin=margin(t=5.5, b=4)),
        plot.margin = margin(0, 0, 5.5, 5.5))

p_forest_s <- ggplot(df_plot_s, aes(y=y)) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey50") +
  geom_errorbarh(aes(xmin=lo, xmax=hi, colour=colour), height=0.28, linewidth=0.8) +
  geom_point(aes(x=est, colour=colour, shape=subgroup), size=3.8) +
  # Arrow from x=1 rightward + "Higher risk" label
  annotate("segment", x=1, xend=3.3, y=-2, yend=-2,
           arrow=arrow(type="closed", length=unit(0.18,"cm")),
           colour="grey30", linewidth=0.7) +
  annotate("text", x=2.1, y=-2.4, label="Higher risk",
           hjust=0.5, size=3.2, colour="grey30", fontface="italic") +
  scale_colour_identity() +
  scale_shape_manual(values=c("Without CRSwNP"=16,"With CRSwNP"=17)) +
  scale_x_log10(breaks=c(0.5,1,1.5,2,3), labels=c("0.5","1.0","1.5","2.0","3.0")) +
  scale_y_continuous(breaks=NULL) +   # limits handled by coord_cartesian
  coord_cartesian(xlim=c(0.4,3.5), ylim=new_ylim, clip="off") +
  labs(x="Adjusted Rate Ratio (aRR)", y=NULL) +
  theme_bw(base_size=11) +
  theme(axis.title.x=element_text(face="bold"),
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.grid=element_blank(), legend.position="none",
        plot.margin=margin(0, 2, 38, 2))

make_txt_col <- function(labs, title, cols, faces) {
  d <- data.frame(y=df_plot_s$y, lab=labs, col=cols, face=faces)
  ggplot(d, aes(y=y, x=0.02, label=lab, colour=col)) +
    geom_text(hjust=0, size=3.3, fontface=d$face) +
    scale_colour_identity() +
    scale_x_continuous(limits=c(0,1), expand=c(0.02,0)) +
    scale_y_continuous(limits=c(0.3,10.7), breaks=NULL) +
    labs(x=NULL, y=NULL, title=title) +
    theme_void(base_size=11) +
    theme(plot.title=element_text(hjust=0, face="bold", size=10, margin=margin(b=4)),
          plot.margin=margin(5.5,5,5.5,0))
}

# Centre-aligned, black text helper
make_txt_col_bw <- function(labs, title, faces) {
  d <- data.frame(y=df_plot_s$y, lab=labs, face=faces)
  ggplot(d, aes(y=y, x=0.5, label=lab)) +
    geom_text(hjust=0.5, size=3.3, fontface=d$face, colour="black") +
    scale_x_continuous(limits=c(0,1), expand=c(0.02,0)) +
    scale_y_continuous(limits=c(0.3,10.7), breaks=NULL) +
    labs(x=NULL, y=NULL, title=title) +
    theme_void(base_size=11) +
    theme(plot.title=element_text(hjust=0.5, face="bold", size=10,
                                margin=margin(t=5.5, b=4)),
          plot.margin=margin(0, 5, 5.5, 0))
}

make_txt_col_bw <- function(labs, title, faces) {
  d <- data.frame(y=df_plot_s$y, lab=labs, face=faces)
  ggplot(d, aes(y=y, x=0.5, label=lab)) +
    geom_text(hjust=0.5, size=3.3, fontface=d$face, colour="black") +
    scale_x_continuous(limits=c(0,1), expand=c(0.02,0)) +
    scale_y_continuous(limits=new_ylim, breaks=NULL) +
    labs(x=NULL, y=NULL, title=title) +
    theme_void(base_size=11) +
    theme(plot.title  = element_text(hjust=0.5, face="bold", size=10,
                                     margin=margin(t=5.5, b=4)),
          plot.margin = margin(0, 5, 5.5, 0))
}

p_ci <- make_txt_col_bw(df_plot_s$est_txt, "aRR (95% CI)",
                        ifelse(df_plot_s$sig, "bold", "plain"))

# Replace p-value column with biomarker × CRSwNP interaction p-value
# (from Section 5 int2 models; shown once per biomarker on the first subgroup row)
get_int_p <- function(obj) {
  obj$table |> dplyr::filter(grepl(":", term)) |> dplyr::pull(p_value) |> dplyr::first()
}
fmt_p <- function(p) dplyr::case_when(p < 0.001 ~ "<0.001",
                                       p < 0.01  ~ sprintf("%.3f", p),
                                       TRUE      ~ sprintf("%.2f",  p))
p_int_FeNO <- get_int_p(int2_FeNO)
p_int_BEC  <- get_int_p(int2_BEC)
p_int_IgE  <- get_int_p(int2_IgE)

p_int_vec   <- c(fmt_p(p_int_FeNO), "",
                 fmt_p(p_int_BEC),  "",
                 fmt_p(p_int_IgE),  "")
sig_int_vec <- c(p_int_FeNO < 0.05, FALSE,
                  p_int_BEC  < 0.05, FALSE,
                  p_int_IgE  < 0.05, FALSE)

p_int_col <- make_txt_col_bw(p_int_vec, "p interaction",
                             ifelse(sig_int_vec, "bold", "plain"))

p_summary <- cowplot::plot_grid(
  p_left, p_forest_s,
  p_ci      + theme(plot.margin = margin(0, 5, 5.5, -15)),
  p_int_col + theme(plot.margin = margin(0, 5, 5.5, -15)),
  nrow=1, rel_widths=c(0.8, 1.6, 0.7, 0.40),
  align="hv", axis="tb"
)
p_summary

ggsave(file.path(Figure_path, "Summary_Forest_Subgroup.pdf"),
       plot=p_summary, width=8, height=4, units="in", device="pdf")

# ── 6.7  Alternative: per-IQR version ──────────────────────────────────────────────────
# IRR are rescaled: est^IQR_log10 (equivalent to exp(beta * IQR_log10))
# Labels show the IQR in original units (ppb, x10^9/L, kU/L)
# iqr_* variables are defined right after data_s5 creation (Section 5.1)

lbl_FeNO <- sprintf("FeNO (per %g ppb)",         round(iqr_FeNO_orig))
lbl_BEC  <- sprintf("BEC (per %.2f \u00d710\u2079/L)", round(iqr_BEC_orig, 2))
lbl_IgE  <- sprintf("IgE (per %g kU/L)",          round(iqr_IgE_orig))

df_plot_iqr <- dplyr::bind_rows(
  dplyr::mutate(sub_FeNO_no$table  |> dplyr::filter(term=="FeNO_log10"),
                biomarker=lbl_FeNO, subgroup="Without CRSwNP", colour="steelblue",  iqr_log=iqr_FeNO_log),
  dplyr::mutate(sub_FeNO_yes$table |> dplyr::filter(term=="FeNO_log10"),
                biomarker=lbl_FeNO, subgroup="With CRSwNP",    colour="steelblue",  iqr_log=iqr_FeNO_log),
  dplyr::mutate(sub_BEC_no$table   |> dplyr::filter(term=="BEC_log10"),
                biomarker=lbl_BEC,  subgroup="Without CRSwNP", colour="darkviolet", iqr_log=iqr_BEC_log),
  dplyr::mutate(sub_BEC_yes$table  |> dplyr::filter(term=="BEC_log10"),
                biomarker=lbl_BEC,  subgroup="With CRSwNP",    colour="darkviolet", iqr_log=iqr_BEC_log),
  dplyr::mutate(sub_IgE_no$table   |> dplyr::filter(term=="IgE_log10"),
                biomarker=lbl_IgE,  subgroup="Without CRSwNP", colour="darkorange", iqr_log=iqr_IgE_log),
  dplyr::mutate(sub_IgE_yes$table  |> dplyr::filter(term=="IgE_log10"),
                biomarker=lbl_IgE,  subgroup="With CRSwNP",    colour="darkorange", iqr_log=iqr_IgE_log)
) |>
  dplyr::mutate(
    est_iqr   = pooled_exp^iqr_log,
    lo_iqr    = ci_exp_low^iqr_log,
    hi_iqr    = ci_exp_high^iqr_log,
    y         = c(8, 7, 5, 4, 2, 1),
    row_label = paste0("    ", subgroup),
    sig       = p_value < 0.05,
    est_txt   = sprintf("%.2f (%.2f\u2013%.2f)", est_iqr, lo_iqr, hi_iqr),
    p_txt     = dplyr::case_when(p_value < 0.001 ~ "<0.001",
                                 p_value < 0.01  ~ sprintf("%.3f", p_value),
                                 TRUE            ~ sprintf("%.2f",  p_value))
  )

hdr_df_iqr <- tibble::tibble(y = c(9, 6, 3), label = c(lbl_FeNO, lbl_BEC, lbl_IgE))

df_left_iqr <- data.frame(
  y    = c(hdr_df_iqr$y,     df_plot_iqr$y),
  lab  = c(hdr_df_iqr$label, df_plot_iqr$row_label),
  face = c(rep("bold", 3),   rep("plain", 6))
)

p_left_iqr <- ggplot(df_left_iqr, aes(y = y, x = 1, label = lab)) +
  geom_text(hjust = 1, size = 3.5, fontface = df_left_iqr$face, colour = "black") +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
  scale_y_continuous(limits = new_ylim, breaks = NULL) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_void(base_size = 12) +
  theme(plot.margin = margin(0, 0, 5.5, 5.5))

p_forest_iqr <- ggplot(df_plot_iqr, aes(y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lo_iqr, xmax = hi_iqr, colour = colour),
                 height = 0.28, linewidth = 0.8) +
  geom_point(aes(x = est_iqr, colour = colour, shape = subgroup), size = 3.8) +
  annotate("segment", x = 1, xend = 3.3, y = -2, yend = -2,
           arrow = arrow(type = "closed", length = unit(0.18, "cm")),
           colour = "grey30", linewidth = 0.7) +
  annotate("text", x = 2.1, y = -2.4, label = "Higher risk",
           hjust = 0.5, size = 3.2, colour = "grey30", fontface = "italic") +
  scale_colour_identity() +
  scale_shape_manual(values = c("Without CRSwNP" = 16, "With CRSwNP" = 17)) +
  scale_x_log10(breaks = c(0.5, 1, 1.5, 2, 3),
                labels  = c("0.5", "1.0", "1.5", "2.0", "3.0")) +
  scale_y_continuous(breaks = NULL) +
  coord_cartesian(xlim = c(0.4, 3.5), ylim = new_ylim, clip = "off") +
  labs(x = "Adjusted Rate Ratio per IQR (aRR)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.title.x    = element_text(face = "bold"),
        axis.text.y     = element_blank(),
        axis.ticks.y    = element_blank(),
        panel.grid      = element_blank(),
        legend.position = "none",
        plot.margin     = margin(0, 2, 38, 2))

make_txt_col_bw2 <- function(labs, title, faces) {
  d <- data.frame(y = df_plot_iqr$y, lab = labs, face = faces)
  ggplot(d, aes(y = y, x = 0.5, label = lab)) +
    geom_text(hjust = 0.5, size = 3.3, fontface = d$face, colour = "black") +
    scale_x_continuous(limits = c(0, 1), expand = c(0.02, 0)) +
    scale_y_continuous(limits = new_ylim, breaks = NULL) +
    labs(x = NULL, y = NULL, title = title) +
    theme_void(base_size = 11) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold", size = 10,
                                     margin = margin(t = 5.5, b = 4)),
          plot.margin = margin(0, 5, 5.5, 0))
}

p_ci_iqr      <- make_txt_col_bw2(df_plot_iqr$est_txt, "aRR (95% CI)",
                                   ifelse(df_plot_iqr$sig, "bold", "plain"))
p_int_col_iqr <- make_txt_col_bw2(p_int_vec, "p interaction",
                                   ifelse(sig_int_vec, "bold", "plain"))

p_summary_iqr <- cowplot::plot_grid(
  p_left_iqr, p_forest_iqr,
  p_ci_iqr      + theme(plot.margin = margin(0, 5, 5.5, -15)),
  p_int_col_iqr + theme(plot.margin = margin(0, 5, 5.5, -15)),
  nrow = 1, rel_widths = c(0.9, 1.6, 0.7, 0.40),
  align = "hv", axis = "tb"
)
p_summary_iqr

ggsave(file.path(Figure_path, "Summary_Forest_Subgroup_per_IQR.pdf"),
       plot = p_summary_iqr, width = 8, height = 4, units = "in", device = "pdf")
# SECTION 7 — One-stage RCS spline: log10(IgE) × CRSwNP
#
# All patients pooled in a single NB model with trial as a fixed
# factor (absorbs between-trial intercept differences).
# Model: rcs(IgE_log10, 4 knots) * CRSwNP + covars +
#        Enrolled_Trial_name + offset(log(days)).
# MI: Rubin pooling across imputations.
# Prediction: linear predictors averaged over all trial fixed effects
#   to marginalise the trial baseline rate.
########################################################################

library(rms)

# ── 7.1  Global knots (from first imputation, all trials combined) ─────
imp1      <- dplyr::filter(data_s5, .imp == 1)
knots_IgE <- as.numeric(quantile(imp1$IgE_log10,
                                  probs = c(0.05, 0.35, 0.65, 0.95),
                                  na.rm = TRUE))
cat("Knots (log10 IgE):", round(knots_IgE, 3), "\n")
cat("Knots (IgE IU/mL):", round(10^knots_IgE), "\n")

# ── 7.2  Prediction grid (covariates at overall imp-1 medians) ─────────
x_seq <- seq(
  quantile(imp1$IgE_log10, 0.05, na.rm = TRUE),
  quantile(imp1$IgE_log10, 0.95, na.rm = TRUE),
  length.out = 100
)

pred_base <- expand.grid(
  IgE_log10 = x_seq,
  CRSwNP    = factor(c("No", "Yes"), levels = c("No", "Yes")),
  KEEP.OUT.ATTRS = FALSE
) |>
  dplyr::mutate(
    GINA_step_numeric       = median(imp1$GINA_step_numeric,  na.rm = TRUE),
    ACQ_score_0W            = median(imp1$ACQ_score_0W,       na.rm = TRUE),
    Attack_history_num      = median(imp1$Attack_history_num, na.rm = TRUE),
    FEV1_preBD_PCT_0W       = median(imp1$FEV1_preBD_PCT_0W, na.rm = TRUE),
    FeNO_log10              = median(imp1$FeNO_log10,         na.rm = TRUE),
    BEC_log10               = median(imp1$BEC_log10,          na.rm = TRUE),
    Follow_up_duration_days = 365
  )

all_trials <- sort(unique(as.character(imp1$Enrolled_Trial_name)))
n_pred     <- nrow(pred_base)
imps       <- sort(unique(data_s5$.imp[data_s5$.imp >= 1]))
M          <- length(imps)

# ── 7.3  Fit one NB model per imputation; average predictions over trials
# For each imputation:
#   1. Fit glm.nb with rcs() * CRSwNP + covars + trial_factor + offset
#   2. Predict at each trial level; average over trials on link scale
#      → marginalises the trial-specific baseline attack rate
#   3. Store mean lp and mean SE across trials
# Then apply Rubin's rules across imputations.

all_lp <- matrix(NA_real_, nrow = n_pred, ncol = M)
all_se <- matrix(NA_real_, nrow = n_pred, ncol = M)

for (j in seq_along(imps)) {
  dat_j <- dplyr::filter(data_s5, .imp == imps[j]) |>
    dplyr::mutate(
      CRSwNP              = factor(CRSwNP, levels = c("No", "Yes")),
      Enrolled_Trial_name = factor(Enrolled_Trial_name)
    )
  trial_lvls <- levels(dat_j$Enrolled_Trial_name)

  fit_j <- tryCatch(
    MASS::glm.nb(
      Attack_number_during_followup ~
        rcs(IgE_log10, knots_IgE) * CRSwNP +
        GINA_step_numeric + ACQ_score_0W + Attack_history_num +
        FEV1_preBD_PCT_0W + FeNO_log10 + BEC_log10 +
        Enrolled_Trial_name +
        offset(log(Follow_up_duration_days)),
      data = dat_j
    ),
    error = function(e) { message("Imp ", imps[j], " failed: ", e$message); NULL }
  )
  if (is.null(fit_j)) next

  pr_list <- lapply(trial_lvls, function(tr) {
    pd <- pred_base |>
      dplyr::mutate(
        CRSwNP              = factor(CRSwNP, levels = c("No", "Yes")),
        Enrolled_Trial_name = factor(tr, levels = trial_lvls)
      )
    tryCatch(
      predict(fit_j, newdata = pd, type = "link", se.fit = TRUE),
      error = function(e) list(fit    = rep(NA_real_, n_pred),
                               se.fit = rep(NA_real_, n_pred))
    )
  })

  lp_mat <- sapply(pr_list, `[[`, "fit")    # n_pred × n_trials
  se_mat <- sapply(pr_list, `[[`, "se.fit")

  all_lp[, j] <- rowMeans(lp_mat, na.rm = TRUE)
  all_se[, j] <- rowMeans(se_mat, na.rm = TRUE)
  cat("  Imputation", imps[j], "done\n")
}

# ── 7.4  Rubin's rules across imputations ─────────────────────────────
Q_bar   <- rowMeans(all_lp, na.rm = TRUE)
U_bar   <- rowMeans(all_se^2, na.rm = TRUE)
B_var   <- apply(all_lp, 1, var, na.rm = TRUE)
T_var   <- U_bar + (1 + 1 / M) * B_var
SE_pool <- sqrt(T_var)

pred_pooled <- pred_base |>
  dplyr::mutate(
    pred   = exp(Q_bar),
    lower  = exp(Q_bar - 1.96 * SE_pool),
    upper  = exp(Q_bar + 1.96 * SE_pool),
    CRSwNP = factor(CRSwNP, levels = c("No", "Yes"),
                    labels = c("Without CRSwNP", "With CRSwNP"))
  )

# ── 7.5  Spline plot ──────────────────────────────────────────────────
ige_breaks_log <- log10(c(10, 30, 100, 300, 1000, 3000, 10000))
ige_labels     <- c("10", "30", "100", "300", "1,000", "3,000", "10,000")
x_range        <- range(x_seq)

p_spline <- ggplot(pred_pooled,
                   aes(x = IgE_log10, y = pred,
                       color = CRSwNP, fill = CRSwNP)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Without CRSwNP" = "grey50",
                                "With CRSwNP"    = "darkorange")) +
  scale_fill_manual(values  = c("Without CRSwNP" = "grey50",
                                "With CRSwNP"    = "darkorange")) +
  scale_x_continuous(breaks = ige_breaks_log, labels = ige_labels,
                     limits = x_range, name = NULL) +
  labs(y = "Predicted annual attack rate") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "top",
        legend.title     = element_blank(),
        axis.title       = element_text(face = "bold"),
        axis.title.x      = element_blank(),
        axis.text.x       = element_blank(),
        axis.ticks.x      = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank())

# ── 7.6  Density strip (distribution of IgE by CRSwNP, imp 1) ─────────
imp1_sp <- imp1 |>
  dplyr::mutate(
    CRSwNP = factor(CRSwNP, levels = c("No", "Yes"),
                    labels = c("Without CRSwNP", "With CRSwNP"))
  )

p_density <- ggplot(imp1_sp,
                    aes(x = IgE_log10, fill = CRSwNP, color = CRSwNP)) +
  geom_density(alpha = 0.25, linewidth = 0.5) +
  scale_color_manual(values = c("Without CRSwNP" = "grey50",
                                "With CRSwNP"    = "darkorange")) +
  scale_fill_manual(values  = c("Without CRSwNP" = "grey50",
                                "With CRSwNP"    = "darkorange")) +
  scale_x_continuous(breaks = ige_breaks_log, labels = ige_labels,
                     limits = x_range, name = "IgE (IU/mL)") +
  labs(y = "Density") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "none",
        axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# ── 7.7  Combine and save ─────────────────────────────────────────────
spline_combined <- cowplot::plot_grid(
  p_spline, p_density,
  ncol        = 1,
  align       = "v",
  rel_heights = c(3, 1)
)
spline_combined

ggsave(file.path(Figure_path, "Spline_IgE_CRSwNP.pdf"),
       plot = spline_combined, width = 6, height = 6,
       units = "in", device = "pdf")
