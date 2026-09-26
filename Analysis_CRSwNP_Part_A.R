# T2 biomarkers and CRSwNP - Part A - Preparing the dataframe

#######################################################################################################################################
#######################################################################################################################################
# Remove everything in the environment 
rm(list = ls())

# After removing large objects, prompt R to release memory back to the OS
gc()
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 1.	Loading the packages
#--------------------------------------------------------------------------------------------------------------------------------------
# 1.1.	From library

# Data manipulation
library(tidyverse)   # includes dplyr, tidyr, purrr, stringr, forcats, tibble, readr, ggplot2, lubridate
library(magrittr)
library(readxl)

# Statistical analysis
library(MASS)
library(mice)
library(mitools)
library(rms)
library(MuMIn)
library(boot)
library(AICcmodavg)
library(caret)
library(gam)
library(rstatix)
library(metafor)
library(lme4)
library(glmmTMB)
library(leaps)

# Survival analysis
library(survival)
library(survivalAnalysis)
library(ggsurvfit)
library(tidycmprsk)

# Visualization
library(ggpubr)
library(ggvenn)
library(cowplot)
library(pheatmap)
library(circlize)
library(forplo)
library(patchwork)
library(gridExtra)
library(grid)
library(DiagrammeR)
library(rsvg)

# Reporting & tables
library(gtsummary)
library(tableone)
library(table1)
library(Gmisc, quietly = TRUE)
library(glue)
library(knitr)
library(flextable)

# Utilities
library(writexl)
library(rspiro)
library(remotes)
library(devtools)
#--------------------------------------------------------------------------------------------------------------------------------------
# 1.2.	From local path

# loading the MIanalysis package
user_name <- Sys.info()[["user"]]  # Automatically get the username
MIAnalysis_path <- file.path("/Users", user_name,
                             "Library/CloudStorage/OneDrive-USherbrooke/Recherche/R_packages/MIAnalysis/R")
devtools::load_all(MIAnalysis_path)
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 2.	Defining the path

#--------------------------------------------------------------------------------------------------------------------------------------
# 2.1. Loading path
Load_path <- file.path(
  "/Users", user_name,
  "Library/CloudStorage/OneDrive-USherbrooke/Recherche/Projets/Pneumologie/ORACLE/ORACLE2_Imputation/Data/6_After_part_D_multiple_imputation/Main_FCS_single_pmm_logreg_polyreg_pmm"
)

#--------------------------------------------------------------------------------------------------------------------------------------
# 2.2. Saving data path
Saving_path <- file.path(
  "/Users", user_name,
  "Library/CloudStorage/OneDrive-USherbrooke/Recherche/Projets/Pneumologie/ORACLE/ORACLE_NP_IgE/Data/After_part_A"
)

#--------------------------------------------------------------------------------------------------------------------------------------
# 2.3. Saving figure path
Figure_path <- file.path(
  "/Users", user_name,
  "Library/CloudStorage/OneDrive-USherbrooke/Recherche/Projets/Pneumologie/ORACLE/ORACLE_NP_IgE/Figures_and_Tables/Part_A"
)


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 3.	Importing the dataset 

#--------------------------------------------------------------------------------------------------------------------------------------
# 3.1.	Read the dataframe
Load_file <- "ORACLE_MI_FULL_long_imp0_to_m.rds"
full_path  <- file.path(Load_path, Load_file)
file.exists(full_path)

data_long <- readRDS(full_path)

#--------------------------------------------------------------------------------------------------------------------------------------
# 3.2.	Convert to data.frame
ORACLE_after_imputation <- as.data.frame(data_long)
colnames(ORACLE_after_imputation)



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 4. Selection of the subjects to include in the analysis

# 4.1 Work on the non-imputed dataset (imp == 0) for subject counts
imp0 <- ORACLE_after_imputation |> dplyr::filter(.imp == 0)

# Counts before any selection
total_n    <- nrow(imp0)
total_rcts <- dplyr::n_distinct(imp0$Enrolled_Trial_name)

# --- Step 1: Exclude trials with no CRSwNP data collected (all NA) ---
trials_step1 <- imp0 |>
  dplyr::group_by(Enrolled_Trial_name) |>
  dplyr::summarise(any_nonNA = any(!is.na(CRSwNP)), .groups = "drop") |>
  dplyr::filter(any_nonNA) |>
  dplyr::pull(Enrolled_Trial_name)

excl1_trials <- imp0 |>
  dplyr::filter(!Enrolled_Trial_name %in% trials_step1) |>
  dplyr::pull(Enrolled_Trial_name) |>
  unique() |> as.character() |> sort()
excl1_n    <- nrow(imp0 |> dplyr::filter(!Enrolled_Trial_name %in% trials_step1))
excl1_rcts <- length(excl1_trials)

# Intermediate counts after step 1
step1_n    <- nrow(imp0 |> dplyr::filter(Enrolled_Trial_name %in% trials_step1))
step1_rcts <- length(trials_step1)

# --- Step 2: Exclude trials where all subjects have CRSwNP = "No" (no cases) ---
trials_step2 <- imp0 |>
  dplyr::filter(Enrolled_Trial_name %in% trials_step1) |>
  dplyr::group_by(Enrolled_Trial_name) |>
  dplyr::summarise(any_yes = any(CRSwNP == "Yes", na.rm = TRUE), .groups = "drop") |>
  dplyr::filter(any_yes) |>
  dplyr::pull(Enrolled_Trial_name)

excl2_trials <- imp0 |>
  dplyr::filter(Enrolled_Trial_name %in% trials_step1,
                !Enrolled_Trial_name %in% trials_step2) |>
  dplyr::pull(Enrolled_Trial_name) |>
  unique() |> as.character() |> sort()
excl2_n    <- nrow(imp0 |> dplyr::filter(Enrolled_Trial_name %in% trials_step1,
                                          !Enrolled_Trial_name %in% trials_step2))
excl2_rcts <- length(excl2_trials)

# Filter all imputations to keep only trials passing both steps
ORACLE_adult <- ORACLE_after_imputation |>
  dplyr::filter(Enrolled_Trial_name %in% trials_step2)
ORACLE_adult$Enrolled_Trial_name <- droplevels(ORACLE_adult$Enrolled_Trial_name)

# Final counts
final_n    <- nrow(imp0 |> dplyr::filter(Enrolled_Trial_name %in% trials_step2))
final_rcts <- dplyr::n_distinct(ORACLE_adult$Enrolled_Trial_name[ORACLE_adult$.imp == 0])

cat("Step 1 - Excluded trials (no CRSwNP data collected):", paste(excl1_trials, collapse = ", "), "\n")
cat("Step 1 - Excluded subjects:", excl1_n, "from", excl1_rcts, "RCTs\n")
cat("Step 2 - Excluded trials (all CRSwNP = No):", paste(excl2_trials, collapse = ", "), "\n")
cat("Step 2 - Excluded subjects:", excl2_n, "from", excl2_rcts, "RCTs\n")
cat("Final dataset:", final_n, "subjects from", final_rcts, "RCTs\n")

# 4.2 Generation of the flowchart
flowchart <- DiagrammeR::grViz(sprintf("
digraph consort {
  graph [layout = dot, rankdir = TB, fontname = 'Helvetica', splines = ortho, nodesep = 0.6]
  node [shape = rectangle, fontname = 'Helvetica', fontsize = 13,
        margin = '0.3,0.2', width = 4.5, height = 0.8]

  A [label = 'ORACLE2 Dataset\\n%d subjects from %d RCTs',
     fillcolor = '#DDEEFF', style = filled]
  C [label = 'Analysed Dataset\\n%d subjects from %d RCTs',
     fillcolor = '#DDEEFF', style = filled]
  E1 [label = 'Excluded (n = %d from %d RCTs)\\nNo CRSwNP data collected',
      fillcolor = '#FFDDDD', style = filled, width = 4.2]
  E2 [label = 'Excluded (n = %d from %d RCTs)\\nAll subjects have no CRSwNP',
      fillcolor = '#FFDDDD', style = filled, width = 4.2]

  node [shape = point, width = 0.01, style = invis, label = '']
  mid1
  mid2

  A    -> mid1 [arrowhead = none]
  mid1 -> mid2 [arrowhead = none]
  mid2 -> C
  mid1 -> E1
  mid2 -> E2

  { rank = same; mid1; E1 }
  { rank = same; mid2; E2 }
}
",
total_n, total_rcts,
final_n, final_rcts,
excl1_n, excl1_rcts,
excl2_n, excl2_rcts
))
flowchart

# 4.3 Save flowchart as PDF
svg_tmp <- tempfile(fileext = ".svg")
writeLines(DiagrammeRsvg::export_svg(flowchart), svg_tmp)
rsvg::rsvg_pdf(svg_tmp, file = file.path(Figure_path, "Flowchart_selection.pdf"))
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 5.	Variables formatting

## 5.1.	Creating of categorical variable

# 5.1.1.	Age_by_group
ORACLE_adult$Age_by_group <- cut(
  ORACLE_adult$Age,
  breaks = c(-10, 40, 50, 60, 100000),
  labels = c("<40", "≥40 to 50", "≥50 to <60", "≥60")
)

# 5.1.2.	BMI_by_group
ORACLE_adult$BMI_by_group <- cut(
  ORACLE_adult$BMI,
  breaks = c(-10, 25, 30, 35, 100000),
  labels = c("<25", "≥25 to 30", "≥30 to <35", "≥35")
)

# 5.1.3.	BEC_by_group
ORACLE_adult$BEC_by_group <- cut(
  ORACLE_adult$BEC,
  breaks = c(0, 0.15, 0.3, Inf),
  labels = c("<1.5", "≥1.5 to <0.3", "≥0.3"),
  right  = FALSE
)

# 5.1.4.	FeNO_by_group
ORACLE_adult$FeNO_by_group <- cut(
  ORACLE_adult$FeNO,
  breaks = c(0, 25, 50, Inf),
  labels = c("<25", "≥25 to <50", "≥50"),
  right  = TRUE
)

# 5.1.5.	IgE_by_group
ORACLE_adult$IgE_by_group <- cut(
  ORACLE_adult$IgE,
  breaks = c(0, 150, 600, 100000),
  labels = c("<150", "≥150 to <600", "≥600")
)

# 5.1.6.	ACQ_score_0W_by_group
ORACLE_adult$ACQ_score_0W_by_group <- cut(
  ORACLE_adult$ACQ_score_0W,
  breaks = c(-10, 1.5, 3, 100000),
  labels = c("<1.5", "≥1.5 to <3", "≥3")
)

# 5.1.7.	Attack_history_last12mo_yes_no
ORACLE_adult <- ORACLE_adult |>
  dplyr::mutate(
    Attack_history_last12mo_yes_no = factor(
      Attack_history_last12mo_yes_no,
      levels = c("No", "Yes")
    )
  )

# 5.1.8.	Attack_history_cat
ORACLE_adult <- ORACLE_adult |>
  dplyr::mutate(
    Attack_history_cat = dplyr::case_when(
      Attack_history_num == 0   ~ "0",
      Attack_history_num == 1   ~ "1",
      Attack_history_num >= 2   ~ "≥2",
      is.na(Attack_history_num) ~ NA_character_,
      TRUE                      ~ NA_character_
    ),
    Attack_history_cat = factor(Attack_history_cat, levels = c("0", "1", "≥2"))
  )

# 5.1.9.	Attack_number_during_followup_yes_no
ORACLE_adult <- ORACLE_adult |>
  dplyr::mutate(
    Attack_number_during_followup_yes_no = dplyr::case_when(
      Attack_number_during_followup == 0   ~ "No",
      Attack_number_during_followup >= 1   ~ "Yes",
      is.na(Attack_number_during_followup) ~ NA_character_,
      TRUE                                 ~ NA_character_
    ),
    Attack_number_during_followup_yes_no = factor(
      Attack_number_during_followup_yes_no,
      levels = c("No", "Yes")
    )
  )

## 5.2.	Creating of variable combining others

# 5.2.1.	LTRA_or_LAMA_or_Theophylline
ORACLE_adult$LTRA_or_LAMA_or_Theophylline <- dplyr::case_when(
  ORACLE_adult$LTRA          == "Yes" |
    ORACLE_adult$LAMA        == "Yes" |
    ORACLE_adult$Theophylline == "Yes" ~ "Yes",
  TRUE ~ "No"
) |>
  factor(levels = c("No", "Yes"))

# 5.2.2.	mOCS_dose_supraphysiologic
ORACLE_adult <- ORACLE_adult |>
  dplyr::mutate(
    mOCS_dose_supraphysiologic = dplyr::case_when(
      is.na(mOCS_dose) ~ NA_character_,
      mOCS_dose > 5    ~ "Yes",
      TRUE             ~ "No"
    ),
    mOCS_dose_supraphysiologic = factor(
      mOCS_dose_supraphysiologic,
      levels = c("No", "Yes")
    )
  )

# 5.2.3.	GINA_step
ORACLE_adult <- ORACLE_adult |>
  dplyr::mutate(
    GINA_step = dplyr::case_when(
      # Step 1: No ICS and No LTRA
      ICS_Dose_Category == "No" & LTRA == "No"                                                ~ "Step 1",
      # Step 2: Low ICS + No LTRA + No LABA OR No ICS + LTRA + No LABA
      (ICS_Dose_Category == "Low"  & LABA == "No") |
        (ICS_Dose_Category == "No" & LTRA == "Yes" & LABA == "No")                           ~ "Step 2",
      # Step 3: Low ICS + Additional controller
      ICS_Dose_Category == "Low" & (LABA == "Yes" | LTRA_or_LAMA_or_Theophylline == "Yes")  ~ "Step 3",
      # Step 4: Medium ICS
      ICS_Dose_Category == "Medium"                                                           ~ "Step 4",
      # Step 5: High ICS or supraphysiologic mOCS
      ICS_Dose_Category == "High" | mOCS_dose_supraphysiologic == "Yes"                      ~ "Step 5",
      TRUE                                                                                    ~ NA_character_
    ),
    GINA_step = factor(
      GINA_step,
      levels  = c("Step 1", "Step 2", "Step 3", "Step 4", "Step 5"),
      ordered = TRUE
    )
  )

# 5.2.4.	GINA_step_numeric
ORACLE_adult$GINA_step_numeric <- as.numeric(ORACLE_adult$GINA_step)


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 6.	Generating the final dataset

# 6.1.	Raw dataset : Not imputed
DATA_Raw <- ORACLE_adult |>
  dplyr::filter(.imp == 0)

# 6.2.	Imputed dataset
DATA_imputed <- ORACLE_adult |>
  dplyr::filter(.imp %in% 1:10)
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 7.	Creation of Tables


# 8.1.	Creation of Table 1
Table_1_CRSwNP <- tbl_summary(
  DATA_Raw |>
    dplyr::mutate(
      GINA_step = factor(GINA_step, levels = c("Step 3", "Step 4", "Step 5"))
    ),
  by = CRSwNP,
  include = c(
    #Enrolled_Trial_name,
    Age,
    Sex,
    BMI,
    Smoking_history_yes_no,
    Psychiatric_disease,
    Eczema,
    Allergic_Rhinitis,
    Airborne_allergen_sensitisation,
    CRSsNP,
    GINA_step,
    LTRA,
    LAMA,
    Theophylline,
    mOCS,
    Intranasal_steroid,
    Attack_history_cat,
    ICU_or_ETI_history,
    ACQ_score_0W,
    FEV1_preBD_PCT_0W,
    Tif_preBD_0W,
    Tif_preBD_PCT_0W,
    BEC,
    FeNO,
    IgE,
    Follow_up_duration_days,
    Attack_number_during_followup_yes_no
  ),
  type  = list(Sex ~ "dichotomous"),
  value = list(Sex ~ "Female"),
  statistic = list(
    all_continuous()  ~ "{median} ({p25}, {p75})",
    all_categorical() ~ "{n} / {N} ({p}%)"
  ),
  missing_text = "(Missing)",
  label = list(
    Age                              ~ "Age (y)",
    Sex                              ~ "Female sex",
    BMI                              ~ "BMI (kg/m²)",
    Smoking_history_yes_no           ~ "Former smokers",
    Psychiatric_disease              ~ "Psychiatric disease",
    Eczema                           ~ "Eczema",
    Allergic_Rhinitis                ~ "Allergic Rhinitis",
    Airborne_allergen_sensitisation  ~ "Airborne allergen sensitisation",
    CRSsNP                           ~ "CRSsNP",
    GINA_step                        ~ "GINA treatment step",
    LTRA                             ~ "LTRA",
    LAMA                             ~ "LAMA",
    Theophylline                     ~ "Theophylline",
    mOCS                             ~ "mOCS",
    Intranasal_steroid               ~ "Intranasal corticosteroids",
    Attack_history_cat               ~ "Severe Attack in the past 12 mo",
    ICU_or_ETI_history               ~ "History of ICU admission or intubation",
    ACQ_score_0W                     ~ "ACQ-5 Score",
    FEV1_preBD_PCT_0W                ~ "FEV1 (% of predicted)",
    Tif_preBD_0W                     ~ "FEV1/FVC",
    Tif_preBD_PCT_0W                 ~ "FEV1/FVC (% of predicted)",
    BEC                              ~ "BEC (×10⁹ cells per L)",
    FeNO                             ~ "FeNO (ppb)",
    IgE                              ~ "Total IgE (kU/L)",
    Follow_up_duration_days          ~ "Follow-up duration (Day)",
    Attack_number_during_followup_yes_no ~ "Severe Attack"
  )
) |>
  gtsummary::add_overall() |>
  gtsummary::add_p() |>
  gtsummary::bold_p() |>
  gtsummary::modify_spanning_header(c(stat_1, stat_2) ~ "**CRSwNP**") |>
  gtsummary::modify_table_body(function(tbl) {
    
    # Remove missing rows for categorical/dichotomous variables
    tbl <- tbl |>
      dplyr::filter(!(row_type == "missing" & var_type %in% c("categorical", "dichotomous")))
    
    # Define section headers: label and the variable they should appear before
    headers <- list(
      list(before = "Age",                    label = "Demographic"),
      list(before = "Psychiatric_disease",    label = "Comorbidities"),
      list(before = "GINA_step",                label = "Medication"),
      list(before = "Attack_history_cat",     label = "Exacerbation history"),
      list(before = "ACQ_score_0W",           label = "Symptoms"),
      list(before = "FEV1_preBD_PCT_0W",      label = "Lung function pre-BD"),
      list(before = "BEC",                    label = "Inflammatory biomarkers"),
      list(before = "Follow_up_duration_days",label = "Follow-up")
    )
    
    # Insert headers in reverse order so row indices remain valid
    for (h in rev(headers)) {
      idx <- min(which(tbl$variable == h$before))
      new_row <- tibble::tibble(
        variable  = paste0("header_", gsub("[^a-zA-Z]", "_", h$label)),
        var_label = h$label,
        label     = h$label,
        row_type  = "label",
        var_type  = NA_character_
      )
      tbl <- dplyr::bind_rows(tbl[seq_len(idx - 1), ], new_row, tbl[idx:nrow(tbl), ])
    }
    
    # Indent rows: variable labels get 4 spaces, subcategory levels get 8 spaces
    tbl <- tbl |>
      dplyr::mutate(label = dplyr::case_when(
        startsWith(variable, "header_") ~ label,
        row_type %in% c("level", "missing") ~ paste0("\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0", label),
        TRUE                            ~ paste0("\u00a0\u00a0\u00a0\u00a0", label)
      ))
    tbl
  }) |>
  gtsummary::modify_table_styling(
    columns = label,
    rows = startsWith(variable, "header_"),
    text_format = "bold"
  ) |>
  gtsummary::modify_header(
    label  ~ "",
    stat_0 ~ "**Overall**  \nN = {n}",
    stat_1 ~ "**No**  \nN = {n}",
    stat_2 ~ "**Yes**  \nN = {n}"
  )

# Save as Word with abbreviation legend (no reference number)
abbreviations <- paste(
  "Abbreviations:",
  "ACQ, Asthma Control Questionnaire;",
  "BEC, Blood Eosinophil Count;",
  "BMI, Body Mass Index;",
  "CRSsNP, Chronic Rhinosinusitis without Nasal Polyps;",
  "CRSwNP, Chronic Rhinosinusitis with Nasal Polyps;",
  "ETI, Endotracheal Intubation;",
  "FeNO, Fractional exhaled Nitric Oxide;",
  "FEV1, Forced Expiratory Volume in 1 second;",
  "FVC, Forced Vital Capacity;",
  "GINA, Global Initiative for Asthma;",
  "ICU, Intensive Care Unit;",
  "IgE, Immunoglobulin E;",
  "LAMA, Long-Acting Muscarinic Antagonist;",
  "LTRA, Leukotriene Receptor Antagonist;",
  "mOCS, Maintenance Oral Corticosteroids;",
  "ppb, parts per billion;",
  "preBD, pre-Bronchodilator."
)

Table_1_CRSwNP |>
  gtsummary::remove_footnote_header(everything()) |>
  gtsummary::as_flex_table() |>
  flextable::footnote(
    i = 1, j = "stat_1",   # † on the "CRSwNP" spanning header
    value = flextable::as_paragraph(
      "8 subjects with not available status of CRSwNP were not included in the current table but included in the analysis with their imputed value."
    ),
    ref_symbols = "\u2020",   # †
    part = "header"
  ) |>
  flextable::footnote(
    i = 2, j = c("stat_0", "stat_1", "stat_2"),   # ‡ on stat column headers
    value = flextable::as_paragraph(
      "Median (Q1, Q3) for continuous variables; n/N (%) for categorical variables."
    ),
    ref_symbols = "\u2021",   # ‡
    part = "header"
  ) |>
  flextable::footnote(
    i = 2, j = "p.value",   # § on p-value column
    value = flextable::as_paragraph(
      "Wilcoxon rank sum test; Pearson\u2019s Chi-squared test; Fisher\u2019s exact test."
    ),
    ref_symbols = "\u00a7",   # §
    part = "header"
  ) |>
  flextable::fontsize(size = 10, part = "all") |>
  flextable::padding(padding = 3, part = "all") |>
  flextable::add_footer_lines(values = abbreviations) |>
  flextable::fontsize(size = 9, part = "footer") |>
  flextable::set_table_properties(layout = "autofit", width = 1) |>
  flextable::save_as_docx(path = file.path(Figure_path, "Table1_CRSwNP.docx"))

summary(DATA_Raw$Enrolled_Trial_name)



# 8.2.	Creation of Table for the trials, ethnicity and regions

Table_CRSwNP_trials_ethnicity_region <- tbl_summary(
  DATA_Raw |>
    dplyr::mutate(
      Ethnicity = dplyr::recode(Ethnicity,
                                "American_Indian_or_Alaska_Native"          = "American Indian or Alaska Native",
                                "Black_or_African_American"                 = "Black or African American",
                                "Native_Hawaiian_or_other_Pacific_Islander" = "Native Hawaiian or other Pacific Islander",
                                "Multiple"                                  = "Multiple ethnicities"
      ),
      Region = dplyr::recode(Region,
                             "North_America" = "North America",
                             "South_America" = "South America",
                             "South_Africa"  = "South Africa"
      )
    ),
  by = CRSwNP,
  include = c(
    Enrolled_Trial_name,
    Ethnicity,
    Region
  ),
  statistic = list(
    all_continuous()  ~ "{median} ({p25}, {p75})",
    all_categorical() ~ "{n} / {N} ({p}%)"
  ),
  missing_text = "(Missing)",
  label = list(
    Enrolled_Trial_name ~ "Trial name",
    Ethnicity           ~ "Ethnicity",
    Region              ~ "Region"
  )
) |>
  gtsummary::add_overall() |>
  gtsummary::add_p(test    = list(all_categorical() ~ "chisq.test"),
                   include = c(Ethnicity, Region)) |>
  gtsummary::bold_p() |>
  gtsummary::modify_spanning_header(c(stat_1, stat_2) ~ "**CRSwNP**") |>
  gtsummary::modify_table_body(function(tbl) {
    tbl |>
      dplyr::mutate(label = dplyr::case_when(
        row_type %in% c("level", "missing") ~ paste0("\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0", label),
        TRUE                                ~ label
      ))
  }) |>
  gtsummary::modify_header(label ~ "")
Table_CRSwNP_trials_ethnicity_region

Table_CRSwNP_trials_ethnicity_region |>
  gtsummary::as_flex_table() |>
  flextable::fontsize(size = 10, part = "all") |>
  flextable::padding(padding = 3, part = "all") |>
  flextable::set_table_properties(layout = "autofit", width = 1) |>
  flextable::save_as_docx(path = file.path(Figure_path, "Table2_CRSwNP_trials_ethnicity_region.docx"))

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 9.	Saving the datasets

save(DATA_Raw,     file = file.path(Saving_path, "DATA_Raw.RData"))
save(DATA_imputed, file = file.path(Saving_path, "DATA_imputed.RData"))

