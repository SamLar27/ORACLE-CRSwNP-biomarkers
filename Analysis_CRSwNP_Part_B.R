########################################################################
# T2 biomarkers and CRSwNP - Part B - Discriminative value of Type 2 inflammatory biomarkers
########################################################################

########################################################################
# SECTION 1 — Packages
########################################################################

library(tidyverse)    # dplyr, tidyr, purrr, stringr, ggplot2, readr
library(readxl)
library(lubridate)

library(rspiro)       # pred_GLI, pred_NHANES3
library(gtsummary)
library(writexl)
library(broom)
library(patchwork)

user_name <- Sys.info()[["user"]]

########################################################################
# SECTION 2 — Paths, saving directory
########################################################################



# 2.1. Loading path
Loading_path <- file.path(
  "/Users", user_name,
  "Library/CloudStorage/OneDrive-USherbrooke/Recherche/Projets/Pneumologie/ORACLE/ORACLE_NP_IgE/Data/After_part_A"
)



Saving_path <- file.path(
  "/Users", user_name,
  "Library/CloudStorage/OneDrive-USherbrooke/Recherche/Projets/Pneumologie/ORACLE/ORACLE_NP_IgE/Data/After_part_B"
)

dir.create(Saving_path, showWarnings = FALSE, recursive = TRUE)

Figure_path <- file.path(
  "/Users", user_name,
  "Library/CloudStorage/OneDrive-USherbrooke/Recherche/Projets/Pneumologie/ORACLE/ORACLE_NP_IgE/Figures_and_Tables/Part_B"
)
dir.create(Figure_path, showWarnings = FALSE, recursive = TRUE)

########################################################################
# SECTION 3 — Loading the dataset
########################################################################

load(file.path(Loading_path, "DATA_Raw.RData"))      # loads DATA_Raw
load(file.path(Loading_path, "DATA_imputed.RData"))  # loads DATA_imputed

########################################################################
# SECTION 4 — Evaluation of the Discriminative value
#
# Evaluate the ability of FeNO, BEC, and IgE to detect CRSwNP.
# Strategy:
#   1. Compute ROC curve per biomarker in each of the 10 imputations
#   2. Pool AUC using Rubin's rules on the logit scale
#      (within-imputation variance via DeLong; between via var of logit-AUCs)
#   3. Plot: individual imputation curves (translucent) + mean ROC curve
########################################################################

library(pROC)
library(metafor)   # rma() for random-effects meta-analysis

# ── 4.1  Helper: ROC curves across MI + pooled AUC ───────────────────

compute_roc_mi <- function(biomarker, data = DATA_imputed) {
  imps <- sort(unique(data$.imp))

  # Per-imputation ROC
  roc_list <- purrr::map(imps, function(m) {
    d <- data |>
      dplyr::filter(.imp == m,
                    !is.na(.data[[biomarker]]),
                    !is.na(CRSwNP))
    pROC::roc(response  = d$CRSwNP,
              predictor = d[[biomarker]],
              levels    = c("No", "Yes"),
              direction = "<",
              quiet     = TRUE)
  })

  aucs    <- purrr::map_dbl(roc_list, ~ as.numeric(pROC::auc(.x)))
  m_imp   <- length(imps)

  # Logit-scale pooling (Rubin's rules)
  logit_aucs <- qlogis(aucs)
  Q_bar      <- mean(logit_aucs)
  B          <- var(logit_aucs)                              # between-imputation variance

  # Within-imputation variance (DeLong), converted to logit scale
  W_vec      <- purrr::map_dbl(roc_list, ~ as.numeric(pROC::var(.x)))
  W_logit    <- W_vec / (aucs * (1 - aucs))^2               # Delta method
  W_bar      <- mean(W_logit)

  T_var      <- W_bar + (1 + 1 / m_imp) * B

  list(
    roc_list   = roc_list,
    aucs       = aucs,
    pooled_auc = plogis(Q_bar),
    ci_low     = plogis(Q_bar - 1.96 * sqrt(T_var)),
    ci_high    = plogis(Q_bar + 1.96 * sqrt(T_var))
  )
}

# ── 4.2  Helper: mean ROC curve (average sensitivity on a fixed grid) ─

mean_roc_curve <- function(roc_list, n_grid = 200) {
  spec_grid <- seq(0, 1, length.out = n_grid)
  sens_mat  <- purrr::map(roc_list, function(r) {
    approx(x    = 1 - r$specificities,
           y    = r$sensitivities,
           xout = spec_grid,
           rule = 2)$y
  }) |> do.call(what = cbind)
  tibble::tibble(
    specificity = 1 - spec_grid,
    sensitivity = rowMeans(sens_mat)
  )
}

# ── 4.3  Plot helper ──────────────────────────────────────────────────

plot_roc_mi <- function(result, biomarker_label, colour) {
  # Individual imputation curves as data frames
  imp_df <- purrr::imap_dfr(result$roc_list, function(r, i) {
    tibble::tibble(
      specificity = r$specificities,
      sensitivity = r$sensitivities,
      imp         = i
    )
  })

  # Mean curve
  mean_df <- mean_roc_curve(result$roc_list)

  auc_label <- sprintf("AUC = %.3f (95%% CI: %.3f - %.3f)",
                       result$pooled_auc, result$ci_low, result$ci_high)

  ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    # Individual imputation curves
    geom_path(data = imp_df,
              aes(x = 1 - specificity, y = sensitivity, group = imp),
              colour = colour, alpha = 0.2, linewidth = 0.5) +
    # Mean curve
    geom_path(data = mean_df,
              aes(x = 1 - specificity, y = sensitivity),
              colour = colour, linewidth = 1.2) +
    annotate("text", x = 0.60, y = 0.08, hjust = 0,
             label = auc_label, size = 3.8, colour = colour, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                       labels = scales::label_number(drop0trailing = TRUE)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                       labels = scales::label_number(drop0trailing = TRUE)) +
    labs(
      x     = "1 - Specificity",
      y     = "Sensitivity",
      title = biomarker_label
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.title  = element_text(face = "bold"),
      plot.title  = element_text(face = "bold", size = 13)
    )
}

# NOTE: All analyses — AUC estimates, ROC curves, Youden thresholds,
# and biomarker comparisons — use the per-trial meta-analytic framework.
# compute_roc_mi() is retained only as an internal building block.

# ── 4.4  Meta-analytic helpers ────────────────────────────────────────

# Per-trial mean ROC curve (averaged across MI), single biomarker
trial_mean_roc_mi <- function(biomarker,
                               trials = unique(as.character(DATA_Raw$Enrolled_Trial_name)),
                               data   = DATA_imputed,
                               n_grid = 200) {
  spec_grid <- seq(0, 1, length.out = n_grid)
  purrr::map(trials, function(tr) {
    d_tr   <- data |> dplyr::filter(Enrolled_Trial_name == tr)
    imps   <- sort(unique(d_tr$.imp))
    n_subj <- nrow(d_tr |> dplyr::filter(.imp == imps[1],
                                          !is.na(.data[[biomarker]]),
                                          !is.na(CRSwNP)))
    roc_list <- purrr::map(imps, function(m) {
      d <- d_tr |>
        dplyr::filter(.imp == m, !is.na(.data[[biomarker]]), !is.na(CRSwNP))
      tryCatch(
        pROC::roc(response = d$CRSwNP, predictor = d[[biomarker]],
                  levels = c("No", "Yes"), direction = "<", quiet = TRUE),
        error = function(e) NULL
      )
    }) |> purrr::compact()
    if (length(roc_list) == 0) return(NULL)
    sens_mat <- purrr::map(roc_list, function(r) {
      approx(x = 1 - r$specificities, y = r$sensitivities,
             xout = spec_grid, rule = 2)$y
    }) |> do.call(what = cbind)
    list(trial = tr, n_subj = n_subj, spec_grid = spec_grid,
         mean_sens = rowMeans(sens_mat))
  }) |> purrr::compact()
}

# Per-trial mean ROC curve for a combined logistic model (two biomarkers)
trial_mean_roc_combined_mi <- function(var1, var2,
                                        trials = unique(as.character(DATA_Raw$Enrolled_Trial_name)),
                                        data   = DATA_imputed,
                                        n_grid = 200) {
  spec_grid <- seq(0, 1, length.out = n_grid)
  purrr::map(trials, function(tr) {
    d_tr   <- data |> dplyr::filter(Enrolled_Trial_name == tr)
    imps   <- sort(unique(d_tr$.imp))
    n_subj <- nrow(d_tr |> dplyr::filter(.imp == imps[1],
                                          !is.na(.data[[var1]]),
                                          !is.na(.data[[var2]]),
                                          !is.na(CRSwNP)))
    roc_list <- purrr::map(imps, function(m) {
      d <- d_tr |>
        dplyr::filter(.imp == m, !is.na(.data[[var1]]),
                      !is.na(.data[[var2]]), !is.na(CRSwNP))
      tryCatch({
        mod  <- glm(reformulate(c(var1, var2), "CRSwNP"), data = d, family = binomial)
        pred <- predict(mod, type = "response")
        pROC::roc(response = d$CRSwNP, predictor = pred,
                  levels = c("No", "Yes"), direction = "<", quiet = TRUE)
      }, error = function(e) NULL)
    }) |> purrr::compact()
    if (length(roc_list) == 0) return(NULL)
    sens_mat <- purrr::map(roc_list, function(r) {
      approx(x = 1 - r$specificities, y = r$sensitivities,
             xout = spec_grid, rule = 2)$y
    }) |> do.call(what = cbind)
    list(trial = tr, n_subj = n_subj, spec_grid = spec_grid,
         mean_sens = rowMeans(sens_mat))
  }) |> purrr::compact()
}

# Trial-size-weighted mean ROC curve across trials
weighted_mean_roc <- function(trial_roc_list) {
  weights  <- purrr::map_dbl(trial_roc_list, ~ .x$n_subj)
  weights  <- weights / sum(weights)
  sens_mat <- purrr::map(trial_roc_list, ~ .x$mean_sens) |> do.call(what = cbind)
  tibble::tibble(
    specificity = 1 - trial_roc_list[[1]]$spec_grid,
    sensitivity = as.vector(sens_mat %*% weights)
  )
}

# Meta-analytic ROC plot: trial-size-weighted pooled curve only
plot_meta_roc <- function(trial_roc_list, rma_obj, biomarker_label, colour) {
  summary_df <- weighted_mean_roc(trial_roc_list)
  auc_label  <- sprintf("AUC = %.3f [%.3f\u2013%.3f]",
                        plogis(rma_obj$b),
                        plogis(rma_obj$ci.lb),
                        plogis(rma_obj$ci.ub))
  ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_path(data = summary_df,
              aes(x = 1 - specificity, y = sensitivity),
              colour = colour, linewidth = 1.2) +
    annotate("text", x = 0.52, y = 0.08, hjust = 0,
             label = auc_label, size = 3.8, colour = colour, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                       labels = scales::label_number(drop0trailing = TRUE)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                       labels = scales::label_number(drop0trailing = TRUE)) +
    labs(x = "1 - Specificity",
         y = "Sensitivity",
         title = biomarker_label) +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = 13))
}

# Meta-analytic Youden threshold: per-trial mean (across MI), then REML across trials
youden_meta <- function(biomarker,
                         trials = unique(as.character(DATA_Raw$Enrolled_Trial_name)),
                         data   = DATA_imputed) {
  per_trial <- purrr::map_dfr(trials, function(tr) {
    tryCatch({
      d_tr <- data |> dplyr::filter(Enrolled_Trial_name == tr)
      imps <- sort(unique(d_tr$.imp))
      thrs <- purrr::map_dbl(imps, function(m) {
        d  <- d_tr |>
          dplyr::filter(.imp == m, !is.na(.data[[biomarker]]), !is.na(CRSwNP))
        r  <- pROC::roc(response = d$CRSwNP, predictor = d[[biomarker]],
                        levels = c("No", "Yes"), direction = "<", quiet = TRUE)
        co <- pROC::coords(r, "best", best.method = "youden",
                           ret = "threshold", transpose = FALSE)
        as.numeric(co[1, "threshold"])
      })
      tibble::tibble(trial = tr,
                     threshold = mean(thrs),
                     se        = sd(thrs) / sqrt(length(thrs)))
    }, error = function(e) {
      message("Skipping trial ", tr, ": ", conditionMessage(e)); NULL
    })
  })
  rma_thr <- metafor::rma(yi = threshold, sei = se,
                           data = per_trial, method = "REML")
  list(per_trial = per_trial,
       pooled    = as.numeric(rma_thr$b),
       ci_low    = as.numeric(rma_thr$ci.lb),
       ci_high   = as.numeric(rma_thr$ci.ub))
}


########################################################################
# SECTION 5 — Primary analysis: per-trial ROC → random-effects
#             meta-analysis
#
# Strategy:
#   1. Identify qualifying trials (≥1 CRSwNP case, so ROC is computable)
#   2. Within each trial: pool AUC across 10 imputations
#      (Rubin's rules on logit scale via compute_roc_mi)
#   3. Random-effects meta-analysis across trials (REML, logit scale)
#   4. Forest plot per biomarker
#   5. Pairwise biomarker comparison on meta-analytic logit-AUC estimates
########################################################################

# ── 5.1  All trials ───────────────────────────────────────────────────
# Part A already excluded trials with no CRSwNP data and trials where all
# subjects have CRSwNP = No, so all remaining trials are included.
qualifying_trials <- unique(as.character(DATA_Raw$Enrolled_Trial_name))

cat("Trials included (n =", length(qualifying_trials), "):",
    paste(qualifying_trials, collapse = ", "), "\n")

# ── 5.2  Helper: per-trial MI-pooled AUC → meta-analysis input ───────
# Returns logit(AUC) and its SE derived from the 95 % CI

trial_auc_mi <- function(biomarker, trials = unique(as.character(DATA_Raw$Enrolled_Trial_name)),
                          data = DATA_imputed) {
  purrr::map_dfr(trials, function(tr) {
    tryCatch({
      d_tr <- data |> dplyr::filter(Enrolled_Trial_name == tr)
      res  <- compute_roc_mi(biomarker, data = d_tr)

      logit_lo <- qlogis(max(res$ci_low,  1e-6))
      logit_hi <- qlogis(min(res$ci_high, 1 - 1e-6))
      se_logit <- (logit_hi - logit_lo) / (2 * 1.96)

      tibble::tibble(
        trial     = tr,
        auc       = res$pooled_auc,
        logit_auc = qlogis(res$pooled_auc),
        se_logit  = se_logit,
        ci_low    = res$ci_low,
        ci_high   = res$ci_high
      )
    }, error = function(e) {
      message("Skipping trial ", tr, ": ", conditionMessage(e))
      NULL
    })
  })
}

# ── 5.3  Run for each biomarker ───────────────────────────────────────
meta_FeNO <- trial_auc_mi("FeNO")
meta_BEC  <- trial_auc_mi("BEC")
meta_IgE  <- trial_auc_mi("IgE")

# ── 5.4  Random-effects meta-analysis (REML, logit scale) ────────────
run_rma <- function(df) {
  metafor::rma(yi = logit_auc, sei = se_logit,
               data = df, method = "REML")
}

rma_FeNO <- run_rma(meta_FeNO)
rma_BEC  <- run_rma(meta_BEC)
rma_IgE  <- run_rma(meta_IgE)

# Summary on AUC scale
summarise_rma <- function(rma_obj, label) {
  auc    <- plogis(rma_obj$b)
  ci_lo  <- plogis(rma_obj$ci.lb)
  ci_hi  <- plogis(rma_obj$ci.ub)
  cat(sprintf("%s — pooled AUC: %.3f (95%% CI: %.3f - %.3f) | I2 = %.1f%%\n",
              label, auc, ci_lo, ci_hi, rma_obj$I2))
}

summarise_rma(rma_FeNO, "FeNO")
summarise_rma(rma_BEC,  "BEC")
summarise_rma(rma_IgE,  "IgE")

# ── 5.5  Forest plot helper ───────────────────────────────────────────
plot_forest_rma <- function(df, rma_obj, biomarker_label, colour, i2_x = -0.3, right_expand = 0) {
  # Back-transform pooled estimate
  pooled_auc <- plogis(rma_obj$b[1])
  pooled_lo  <- plogis(rma_obj$ci.lb)
  pooled_hi  <- plogis(rma_obj$ci.ub)

  # Alphabetical order: A at top, Z at bottom
  # (ggplot draws level 1 at bottom, last level at top)
  df <- df |>
    dplyr::arrange(dplyr::desc(trial)) |>
    dplyr::mutate(trial = factor(trial, levels = unique(trial)))

  # Summary row
  summary_row <- tibble::tibble(
    trial   = factor("Pooled (RE)", levels = c(levels(df$trial), "Pooled (RE)")),
    auc     = pooled_auc,
    ci_low  = pooled_lo,
    ci_high = pooled_hi,
    is_pool = TRUE
  )
  n_trials <- nlevels(df$trial)
  plot_df <- df |>
    dplyr::mutate(is_pool = FALSE) |>
    dplyr::bind_rows(summary_row) |>
    dplyr::mutate(auc_text = sprintf("%.3f (%.3f - %.3f)", auc, ci_low, ci_high))

  # Diamond polygon for the pooled row
  y_pool     <- n_trials + 1
  half_h     <- 0.20
  diamond_df <- tibble::tibble(
    x = c(pooled_lo, pooled_auc, pooled_hi, pooled_auc),
    y = c(y_pool,    y_pool + half_h, y_pool, y_pool - half_h)
  )

  x_text <- 1.07

  ggplot(plot_df, aes(x = auc, y = trial)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey60") +
    # Trial rows: errorbar + point
    geom_errorbarh(
      data = plot_df |> dplyr::filter(!is_pool),
      aes(xmin = ci_low, xmax = ci_high),
      height = 0.25, linewidth = 0.8, colour = colour
    ) +
    geom_point(
      data = plot_df |> dplyr::filter(!is_pool),
      aes(x = auc), shape = 18, size = 3, colour = colour
    ) +
    # Pooled row: diamond spanning LL to UL
    geom_polygon(data = diamond_df, aes(x = x, y = y),
                 fill = "black", colour = "black", alpha = 0.9,
                 inherit.aes = FALSE) +
    # AUC text column
    geom_text(aes(x = x_text, y = trial, label = auc_text),
              hjust = 0, size = 3.2, colour = "black") +
    annotate("text", x = x_text, y = n_trials + 1.7,
             label = "AUC (95% CI)", hjust = 0, size = 3.4,
             fontface = "bold", colour = "black") +
    # I² annotation
    annotate("text", x = i2_x, y = 0.4,
             label = sprintf("I\u00b2 = %.1f%%", rma_obj$I2),
             hjust = 0, vjust = 1, size = 2.8,
             colour = "grey30", fontface = "italic") +
    scale_x_continuous(breaks = seq(0.0, 1.0, 0.2),
                       expand = expansion(add = c(0, right_expand))) +
    coord_cartesian(xlim = c(0.0, 1.0), clip = "off") +
    labs(x = "AUC", y = NULL, title = biomarker_label
        ) +
    theme_bw(base_size = 11) +
    theme(
      axis.title    = element_text(face = "bold"),
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      plot.margin   = margin(5, 110, 5, 5)
    )
}

# ── 5.6  Build and save ───────────────────────────────────────────────
fp_FeNO <- plot_forest_rma(meta_FeNO, rma_FeNO,
                            "FeNO",                  "steelblue")
fp_BEC  <- plot_forest_rma(meta_BEC,  rma_BEC,
                            "BEC","darkviolet")
fp_IgE  <- plot_forest_rma(meta_IgE,  rma_IgE,
                            "Total IgE",         "darkorange"   )

forest_panel <- fp_FeNO | fp_BEC | fp_IgE
forest_panel
fp_FeNO

ggsave(
  file.path(Figure_path, "Forest_AUC_primary.pdf"),
  plot   = forest_panel,
  width  = 14, height = 6, units = "in", device = "pdf"
)

# ── Pairwise biomarker comparison (meta-analytic logit-AUC estimates) ─
# z-test on difference of pooled logit-AUC ± pooled SE (Bonferroni × 3)
compare_rma_pair <- function(rma1, rma2, label) {
  diff <- as.numeric(rma1$b) - as.numeric(rma2$b)
  se   <- sqrt(rma1$se^2 + rma2$se^2)
  z    <- diff / se
  p    <- 2 * pnorm(-abs(z))
  tibble::tibble(
    Comparison  = label,
    AUC_1       = round(plogis(rma1$b), 3),
    AUC_2       = round(plogis(rma2$b), 3),
    AUC_diff    = round(plogis(rma1$b) - plogis(rma2$b), 4),
    z_statistic = round(z, 3),
    p_value     = round(p, 4),
    p_bonf      = round(min(p * 3, 1), 4)
  )
}

delong_meta <- dplyr::bind_rows(
  compare_rma_pair(rma_FeNO, rma_BEC,  "FeNO vs BEC"),
  compare_rma_pair(rma_FeNO, rma_IgE,  "FeNO vs IgE"),
  compare_rma_pair(rma_BEC,  rma_IgE,  "BEC vs IgE")
)
delong_meta



########################################################################
# SECTION 5b — Discriminative performance at selected thresholds
#              (restricted to qualifying trials, pooled across MI)
#
# Sens, Spec, PPV, NPV, LR+, LR− at clinically relevant thresholds.
# Pooled across 10 imputations: logit scale for proportions,
# log scale for likelihood ratios (Rubin's rules).
########################################################################

# ── Helper: pool proportion (logit + Rubin's) ─────────────────────────
pool_prop_ci <- function(x, digits_pct = 1) {
  m  <- length(x[!is.na(x)])
  x  <- x[!is.na(x)]
  lx <- qlogis(pmax(pmin(x, 1 - 1e-6), 1e-6))
  Q  <- mean(lx)
  B  <- var(lx)
  W  <- mean(x * (1 - x) / m)
  Tv <- W + (1 + 1/m) * B
  lo <- plogis(Q - 1.96 * sqrt(Tv)) * 100
  hi <- plogis(Q + 1.96 * sqrt(Tv)) * 100
  sprintf(paste0("%.", digits_pct, "f (%.", digits_pct, "f\u2013%.", digits_pct, "f)"),
          plogis(Q) * 100, lo, hi)
}

# ── Helper: pool likelihood ratio (log + Rubin's) ─────────────────────
pool_lr_ci <- function(x, digits = 2) {
  x  <- x[!is.na(x) & is.finite(x) & x > 0]
  m  <- length(x)
  lx <- log(x)
  Q  <- mean(lx)
  B  <- var(lx)
  Tv <- B * (1 + 1/m)
  sprintf(paste0("%.", digits, "f (%.", digits, "f\u2013%.", digits, "f)"),
          exp(Q), exp(Q - 1.96 * sqrt(Tv)), exp(Q + 1.96 * sqrt(Tv)))
}

# ── Core: metrics at one threshold, pooled across MI ──────────────────
diag_at_threshold <- function(biomarker, threshold,
                               data = DATA_imputed) {
  imps <- sort(unique(data$.imp))
  m <- purrr::map_dfr(imps, function(imp) {
    d <- data |>
      dplyr::filter(.imp == imp, !is.na(.data[[biomarker]]), !is.na(CRSwNP))
    pred_pos   <- d[[biomarker]] >= threshold
    actual_pos <- d$CRSwNP == "Yes"
    TP <- sum( pred_pos &  actual_pos)
    FP <- sum( pred_pos & !actual_pos)
    TN <- sum(!pred_pos & !actual_pos)
    FN <- sum(!pred_pos &  actual_pos)
    tibble::tibble(
      sens = TP / (TP + FN),
      spec = TN / (TN + FP),
      ppv  = if (TP + FP == 0) NA_real_ else TP / (TP + FP),
      npv  = TN / (TN + FN),
      lrp  = (TP / (TP + FN)) / (FP / (FP + TN)),
      lrm  = (FN / (TP + FN)) / (TN / (TN + FP))
    )
  })
  tibble::tibble(
    Threshold   = threshold,
    Sensitivity = pool_prop_ci(m$sens),
    Specificity = pool_prop_ci(m$spec),
    PPV         = pool_prop_ci(m$ppv),
    NPV         = pool_prop_ci(m$npv),
    `LR+`       = pool_lr_ci(m$lrp),
    `LR-`       = pool_lr_ci(m$lrm)
  )
}

make_diag_table <- function(biomarker, thresholds, unit_label) {
  purrr::map_dfr(thresholds, ~ diag_at_threshold(biomarker, .x)) |>
    dplyr::mutate(Threshold = paste0(Threshold, " ", unit_label))
}

tbl_FeNO <- make_diag_table("FeNO", c(25, 35, 50),       "ppb")
tbl_BEC  <- make_diag_table("BEC",  c(0.15, 0.30, 0.45), "x10\u2079/L")
tbl_IgE  <- make_diag_table("IgE",  c(150, 300, 600),    "kU/L")

# ── Export to Word ─────────────────────────────────────────────────────
footnote_txt <- paste(
  "Values are pooled across 10 multiple imputations using Rubin's rules",
  "(logit scale for proportions, log scale for likelihood ratios).",
  "PPV = Positive Predictive Value; NPV = Negative Predictive Value;",
  "LR+ = Positive Likelihood Ratio; LR- = Negative Likelihood Ratio.",
  "Prevalence of CRSwNP in this cohort: ~12.5%."
)

make_word_table <- function(tbl, title, footnote, path) {
  landscape <- officer::prop_section(
    page_size = officer::page_size(orient = "landscape")
  )
  flextable::flextable(tbl) |>
    flextable::set_header_labels(
      Threshold   = "Threshold",
      Sensitivity = "Sensitivity\n% (95% CI)",
      Specificity = "Specificity\n% (95% CI)",
      PPV         = "PPV\n% (95% CI)",
      NPV         = "NPV\n% (95% CI)",
      `LR+`       = "LR+\n(95% CI)",
      `LR-`       = "LR-\n(95% CI)"
    ) |>
    flextable::add_header_lines(values = title) |>
    flextable::bold(part = "header") |>
    flextable::align(align = "center", part = "header") |>
    flextable::align(j = "Threshold", align = "left",   part = "body") |>
    flextable::align(j = 2:7,         align = "center", part = "body") |>
    flextable::add_footer_lines(values = footnote) |>
    flextable::fontsize(size = 10, part = "all") |>
    flextable::fontsize(size = 8,  part = "footer") |>
    flextable::autofit() |>
    flextable::save_as_docx(path = path, pr_section = landscape)
}


# ── Combined table: all three biomarkers in one landscape Word document ─
tbl_combined <- dplyr::bind_rows(
  tibble::tibble(Threshold = "FeNO",              Sensitivity = "", Specificity = "",
                 PPV = "", NPV = "", `LR+` = "", `LR-` = "", .is_header = TRUE),
  tbl_FeNO |> dplyr::mutate(.is_header = FALSE),
  tibble::tibble(Threshold = "Blood eosinophils", Sensitivity = "", Specificity = "",
                 PPV = "", NPV = "", `LR+` = "", `LR-` = "", .is_header = TRUE),
  tbl_BEC  |> dplyr::mutate(.is_header = FALSE),
  tibble::tibble(Threshold = "Total IgE",         Sensitivity = "", Specificity = "",
                 PPV = "", NPV = "", `LR+` = "", `LR-` = "", .is_header = TRUE),
  tbl_IgE  |> dplyr::mutate(.is_header = FALSE)
)

header_rows <- which( tbl_combined$.is_header)
data_rows   <- which(!tbl_combined$.is_header)
tbl_display <- tbl_combined |> dplyr::select(-.is_header)

landscape <- officer::prop_section(
  page_size = officer::page_size(orient = "landscape")
)

ft_combined <- flextable::flextable(tbl_display) |>
  flextable::set_header_labels(
    Threshold   = "Threshold",
    Sensitivity = "Sensitivity\n% (95% CI)",
    Specificity = "Specificity\n% (95% CI)",
    PPV         = "PPV\n% (95% CI)",
    NPV         = "NPV\n% (95% CI)",
    `LR+`       = "LR+\n(95% CI)",
    `LR-`       = "LR-\n(95% CI)"
  ) |>
  flextable::add_header_lines(
    values = "Discriminative performance of T2 biomarkers for CRSwNP detection"
  ) |>
  flextable::bold(i = header_rows, part = "body") |>
  flextable::bg(i = header_rows, bg = "#D9D9D9", part = "body") |>
  flextable::merge_h(i = header_rows, part = "body") |>
  flextable::padding(i = data_rows, j = "Threshold",
                     padding.left = 15, part = "body") |>
  flextable::bold(part = "header") |>
  flextable::align(align = "center", part = "header") |>
  flextable::align(j = "Threshold", align = "left",   part = "body") |>
  flextable::align(j = 2:7,         align = "center", part = "body") |>
  flextable::add_footer_lines(values = footnote_txt) |>
  flextable::fontsize(size = 10, part = "all") |>
  flextable::fontsize(size = 8,  part = "footer") |>
  flextable::autofit()

flextable::save_as_docx(ft_combined,
  path       = file.path(Figure_path, "DiagTable_Combined.docx"),
  pr_section = landscape)

# ── Youden threshold (meta-analytic: per-trial REML) ──────────────────
# Per trial: mean Youden threshold across 10 imputations.
# Across trials: REML random-effects meta-analysis.

youden_FeNO <- youden_meta("FeNO")
youden_BEC  <- youden_meta("BEC")
youden_IgE  <- youden_meta("IgE")

thr_FeNO <- youden_FeNO$pooled
thr_BEC  <- youden_BEC$pooled
thr_IgE  <- youden_IgE$pooled

cat(sprintf("Youden threshold (meta-analytic) — FeNO: %.2f ppb | BEC: %.2f x10\u2079/L | IgE: %.2f kU/L\n",
            thr_FeNO, thr_BEC, thr_IgE))

youden_row <- function(biomarker, threshold, unit) {
  row <- diag_at_threshold(biomarker, threshold)
  row$Threshold <- sprintf("Youden: %.2f %s", threshold, unit)
  row
}

tbl_FeNO_y <- dplyr::bind_rows(tbl_FeNO, youden_row("FeNO", thr_FeNO, "ppb"))
tbl_BEC_y  <- dplyr::bind_rows(tbl_BEC,  youden_row("BEC",  thr_BEC,  "x10\u2079/L"))
tbl_IgE_y  <- dplyr::bind_rows(tbl_IgE,  youden_row("IgE",  thr_IgE,  "kU/L"))

# Rebuild combined table with Youden rows highlighted
tbl_combined2 <- dplyr::bind_rows(
  tibble::tibble(Threshold = "FeNO",              Sensitivity = "", Specificity = "",
                 PPV = "", NPV = "", `LR+` = "", `LR-` = "",
                 .is_header = TRUE, .is_youden = FALSE),
  tbl_FeNO_y |> dplyr::mutate(.is_header = FALSE, .is_youden = grepl("Youden", Threshold)),
  tibble::tibble(Threshold = "Blood eosinophils", Sensitivity = "", Specificity = "",
                 PPV = "", NPV = "", `LR+` = "", `LR-` = "",
                 .is_header = TRUE, .is_youden = FALSE),
  tbl_BEC_y  |> dplyr::mutate(.is_header = FALSE, .is_youden = grepl("Youden", Threshold)),
  tibble::tibble(Threshold = "Total IgE",         Sensitivity = "", Specificity = "",
                 PPV = "", NPV = "", `LR+` = "", `LR-` = "",
                 .is_header = TRUE, .is_youden = FALSE),
  tbl_IgE_y  |> dplyr::mutate(.is_header = FALSE, .is_youden = grepl("Youden", Threshold))
)

header_rows2 <- which( tbl_combined2$.is_header)
youden_rows2 <- which( tbl_combined2$.is_youden)
data_rows2   <- which(!tbl_combined2$.is_header)
tbl_display2 <- tbl_combined2 |> dplyr::select(-.is_header, -.is_youden)

ft2 <- flextable::flextable(tbl_display2) |>
  flextable::set_header_labels(
    Threshold   = "Threshold",
    Sensitivity = "Sensitivity\n% (95% CI)",
    Specificity = "Specificity\n% (95% CI)",
    PPV         = "PPV\n% (95% CI)",
    NPV         = "NPV\n% (95% CI)",
    `LR+`       = "LR+\n(95% CI)",
    `LR-`       = "LR-\n(95% CI)"
  ) |>
  flextable::add_header_lines(
    values = "Discriminative performance of T2 biomarkers for CRSwNP detection"
  ) |>
  flextable::bold(i = header_rows2, part = "body") |>
  flextable::bg(i = header_rows2, bg = "#D9D9D9", part = "body") |>
  flextable::merge_h(i = header_rows2, part = "body") |>
  flextable::italic(i = youden_rows2, part = "body") |>
  flextable::bg(i = youden_rows2, bg = "#FFF9C4", part = "body") |>
  flextable::padding(i = data_rows2, j = "Threshold",
                     padding.left = 15, part = "body") |>
  flextable::bold(part = "header") |>
  flextable::align(align = "center", part = "header") |>
  flextable::align(j = "Threshold", align = "left",   part = "body") |>
  flextable::align(j = 2:7,         align = "center", part = "body") |>
  flextable::add_footer_lines(values = paste(
    footnote_txt,
    "Youden threshold (italic, highlighted) = threshold maximising Sensitivity + Specificity - 1,",
    "derived from a per-trial REML random-effects meta-analysis of Youden thresholds."
  )) |>
  flextable::fontsize(size = 10, part = "all") |>
  flextable::fontsize(size = 8,  part = "footer") |>
  flextable::autofit()

flextable::save_as_docx(ft2,
  path       = file.path(Figure_path, "DiagTable_Combined.docx"),
  pr_section = landscape)


########################################################################
# SECTION 6 — Combined FeNO + BEC: does combining improve the AUC?
#
# Strategy: per-trial logistic model → pool AUC across MI within each
# trial → REML meta-analysis across trials.  Same pipeline as Section 5.
# No direct analysis on the full pooled dataset.
########################################################################

# compute_roc_combined is kept as a building block called inside
# trial_auc_combined_mi() — NOT called on the full population directly.
compute_roc_combined <- function(var1, var2, data = DATA_imputed) {
  imps <- sort(unique(data$.imp))
  roc_list <- purrr::map(imps, function(m) {
    d <- data |>
      dplyr::filter(.imp == m,
                    !is.na(.data[[var1]]), !is.na(.data[[var2]]),
                    !is.na(CRSwNP))
    mod  <- glm(reformulate(c(var1, var2), "CRSwNP"), data = d, family = binomial)
    pred <- predict(mod, type = "response")
    pROC::roc(response = d$CRSwNP, predictor = pred,
              levels = c("No", "Yes"), direction = "<", quiet = TRUE)
  })
  aucs       <- purrr::map_dbl(roc_list, ~ as.numeric(pROC::auc(.x)))
  m_imp      <- length(imps)
  logit_aucs <- qlogis(aucs)
  Q_bar      <- mean(logit_aucs)
  B          <- var(logit_aucs)
  W_vec      <- purrr::map_dbl(roc_list, ~ as.numeric(pROC::var(.x)))
  W_logit    <- W_vec / (aucs * (1 - aucs))^2
  T_var      <- mean(W_logit) + (1 + 1/m_imp) * B
  list(roc_list   = roc_list, aucs = aucs,
       pooled_auc = plogis(Q_bar),
       ci_low     = plogis(Q_bar - 1.96 * sqrt(T_var)),
       ci_high    = plogis(Q_bar + 1.96 * sqrt(T_var)))
}

# ── Per-trial meta-analytic approach ──────────────────────────────────
# Fit combined model within each trial × imputation, pool across MI,
# then meta-analyze across trials (REML) — same pipeline as individual biomarkers.

trial_auc_combined_mi <- function(var1, var2,
                                   trials = unique(as.character(DATA_Raw$Enrolled_Trial_name)),
                                   data   = DATA_imputed) {
  purrr::map_dfr(trials, function(tr) {
    tryCatch({
      d_tr <- data |> dplyr::filter(Enrolled_Trial_name == tr)
      imps <- sort(unique(d_tr$.imp))
      roc_list <- purrr::map(imps, function(m) {
        d <- d_tr |>
          dplyr::filter(.imp == m,
                        !is.na(.data[[var1]]), !is.na(.data[[var2]]),
                        !is.na(CRSwNP))
        mod  <- glm(reformulate(c(var1, var2), "CRSwNP"),
                    data = d, family = binomial)
        pred <- predict(mod, type = "response")
        pROC::roc(response = d$CRSwNP, predictor = pred,
                  levels = c("No", "Yes"), direction = "<", quiet = TRUE)
      })
      aucs       <- purrr::map_dbl(roc_list, ~ as.numeric(pROC::auc(.x)))
      m_imp      <- length(imps)
      logit_aucs <- qlogis(aucs)
      Q_bar      <- mean(logit_aucs)
      B          <- var(logit_aucs)
      W_vec      <- purrr::map_dbl(roc_list, ~ as.numeric(pROC::var(.x)))
      W_logit    <- W_vec / (aucs * (1 - aucs))^2
      T_var      <- mean(W_logit) + (1 + 1/m_imp) * B
      logit_lo   <- qlogis(max(plogis(Q_bar - 1.96*sqrt(T_var)), 1e-6))
      logit_hi   <- qlogis(min(plogis(Q_bar + 1.96*sqrt(T_var)), 1-1e-6))
      tibble::tibble(
        trial     = tr,
        auc       = plogis(Q_bar),
        logit_auc = Q_bar,
        se_logit  = (logit_hi - logit_lo) / (2 * 1.96),
        ci_low    = plogis(Q_bar - 1.96*sqrt(T_var)),
        ci_high   = plogis(Q_bar + 1.96*sqrt(T_var))
      )
    }, error = function(e) {
      message("Skipping trial ", tr, ": ", conditionMessage(e)); NULL
    })
  })
}

meta_combined <- trial_auc_combined_mi("FeNO", "BEC")
rma_combined  <- metafor::rma(yi = logit_auc, sei = se_logit,
                               data = meta_combined, method = "REML")

cat(sprintf("Combined FeNO+BEC (per-trial meta): AUC = %.3f (%.3f - %.3f) | I2 = %.1f%%\n",
            plogis(rma_combined$b), plogis(rma_combined$ci.lb),
            plogis(rma_combined$ci.ub), rma_combined$I2))

# Formal comparison via z-test on pooled logit-AUC estimates
compare_rma <- function(rma_comb, rma_single, label) {
  diff_logit <- as.numeric(rma_comb$b) - as.numeric(rma_single$b)
  se_diff    <- sqrt(rma_comb$se^2 + rma_single$se^2)
  z          <- diff_logit / se_diff
  p          <- 2 * pnorm(-abs(z))
  tibble::tibble(
    Comparison  = label,
    AUC_comb    = round(plogis(rma_comb$b),   3),
    AUC_single  = round(plogis(rma_single$b), 3),
    AUC_diff    = round(plogis(rma_comb$b) - plogis(rma_single$b), 4),
    z_statistic = round(z, 3),
    p_value     = round(p, 4)
  )
}

comparison_meta <- dplyr::bind_rows(
  compare_rma(rma_combined, rma_FeNO, "Combined vs FeNO (per-trial meta)"),
  compare_rma(rma_combined, rma_BEC,  "Combined vs BEC  (per-trial meta)")
)
comparison_meta

# Forest plot adding combined model column
fp_combined <- plot_forest_rma(meta_combined, rma_combined,
                                "FeNO + BEC (combined)", "black", i2_x = -0.15, right_expand = 0.03)
forest_panel_comb <- fp_combined
forest_panel_comb

ggsave(file.path(Figure_path, "Forest_AUC_combined.pdf"),
       plot = forest_panel_comb, width = 7, height = 6,
       units = "in", device = "pdf")

# Meta-analytic comparison: combined vs individual (uses comparison_meta computed above)
comparison_meta

# Meta-analytic ROC comparison plot
# Per-trial mean ROC curves (averaged across MI), then trial-size-weighted summary curve
trial_roc_FeNO <- trial_mean_roc_mi("FeNO")
trial_roc_BEC  <- trial_mean_roc_mi("BEC")
trial_roc_comb <- trial_mean_roc_combined_mi("FeNO", "BEC")

meta_feno <- weighted_mean_roc(trial_roc_FeNO)
meta_bec  <- weighted_mean_roc(trial_roc_BEC)
meta_comb <- weighted_mean_roc(trial_roc_comb)

p_combined_comparison <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_path(data = meta_feno, aes(x = 1 - specificity, y = sensitivity),
            colour = "steelblue",  linewidth = 1) +
  geom_path(data = meta_bec,  aes(x = 1 - specificity, y = sensitivity),
            colour = "darkviolet", linewidth = 1) +
  geom_path(data = meta_comb, aes(x = 1 - specificity, y = sensitivity),
            colour = "black",      linewidth = 1.3) +
  annotate("text", x = 0.48, y = 0.16, hjust = 0, size = 3.5, colour = "steelblue",
           fontface = "bold",
           label = sprintf("FeNO: AUC = %.3f [%.3f\u2013%.3f]",
                           plogis(rma_FeNO$b), plogis(rma_FeNO$ci.lb), plogis(rma_FeNO$ci.ub))) +
  annotate("text", x = 0.48, y = 0.10, hjust = 0, size = 3.5, colour = "darkviolet",
           fontface = "bold",
           label = sprintf("BEC: AUC = %.3f [%.3f\u2013%.3f]",
                           plogis(rma_BEC$b), plogis(rma_BEC$ci.lb), plogis(rma_BEC$ci.ub))) +
  annotate("text", x = 0.48, y = 0.04, hjust = 0, size = 3.5, colour = "black",
           fontface = "bold",
           label = sprintf("FeNO + BEC: AUC = %.3f [%.3f\u2013%.3f]",
                           plogis(rma_combined$b), plogis(rma_combined$ci.lb), plogis(rma_combined$ci.ub))) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                     labels = scales::label_number(drop0trailing = TRUE)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                     labels = scales::label_number(drop0trailing = TRUE)) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 12) +
  theme(axis.title = element_text(face = "bold"))

p_combined_comparison

ggsave(file.path(Figure_path, "ROC_combined_FeNO_BEC.pdf"),
       plot = p_combined_comparison, width = 6, height = 5,
       units = "in", device = "pdf")

########################################################################
# SECTION 7 — Individual biomarker ROC curves (meta-analytic)
#
# Per-trial mean ROC curves + trial-size-weighted summary curve for
# each of the three T2 biomarkers, displayed as a 3-panel figure.
########################################################################

# Compute per-trial ROC curves for IgE (FeNO and BEC already computed above)
trial_roc_IgE <- trial_mean_roc_mi("IgE")

# Individual meta-analytic ROC plots
roc_plot_FeNO <- plot_meta_roc(trial_roc_FeNO, rma_FeNO,
                                "Fractional exhaled nitric oxide", "steelblue") +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
roc_plot_BEC  <- plot_meta_roc(trial_roc_BEC,  rma_BEC,
                                "Blood eosinophil count",          "darkviolet") +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
roc_plot_IgE  <- plot_meta_roc(trial_roc_IgE,  rma_IgE,
                                "Total immunoglobulin E",          "darkorange") +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))

# Multipanel: 3 biomarkers side by side
roc_panel <- (roc_plot_FeNO | roc_plot_BEC | roc_plot_IgE) +
  plot_annotation(
    tag_levels = "A", tag_prefix = "(", tag_suffix = ")",
    theme = theme(plot.tag = element_text(size = 16, face = "bold"))
  )
roc_panel

ggsave(file.path(Figure_path, "ROC_individual_biomarkers.pdf"),
       plot = roc_panel, width = 15, height = 5,
       units = "in", device = "pdf")

