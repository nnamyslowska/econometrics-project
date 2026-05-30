# =============================================================================
# Kinship Intensity and Political Freedom: An Ordered Logit / Probit Analysis
# =============================================================================
#
# DEPENDENT VARIABLE: fh_status (Freedom House Freedom Status)
#   3 ordered categories. Built here as an ORDERED FACTOR:
#       Not Free  <  Partly Free  <  Free      (higher = MORE free)
#   -> ordered logit / ordered probit (proportional-odds models)
#
# KEY INDEPENDENT VARIABLE: kinship (ekne kinship intensity index) [0, 1]
#   0 = nuclear-family society ; 1 = tight kinship / cousin-marriage networks
#   Theory (Schulz et al. 2019; Enke 2019): tight kinship fosters in-group
#   loyalty and personalised exchange, undermining the impersonal institutions
#   that sustain political freedom. Expected sign: NEGATIVE (less free).
#
# CONTROLS (general model):
#   log_gdppc  (CENTRED)  log_pop   urban_pct   unemp   trade   internet
# INTERACTION (requirement c): kinship x log_gdppc_c
#   -> does wealth offset the kinship penalty on freedom?
#
# DROPPED from earlier draft (bad / circular controls):
#   cpi, govt_eff  -> same construct as the outcome (corruption / state quality)
#   dem_duration   -> near-circular with current freedom status
#
# DATA:
#   QoG Standard Dataset Jan 2026 (Teorell et al.) -- qog_std_cs_jan26.csv
#   Kinship intensity index -- kinship_df.csv
# =============================================================================


# -----------------------------------------------------------------------------
# 0. PACKAGES
# -----------------------------------------------------------------------------
# install.packages(c("MASS","brant","pscl","generalhoslem","erer","DescTools",
#                    "stargazer","ggplot2","dplyr","lmtest","car","sandwich",
#                    "nortest"))
# WNE::linktest is a Faculty-of-Economics package; install from your lab source
# if available. A manual fallback (linktest_ordered) is provided below.

library(MASS)          # polr() -- ordered logit / probit
library(brant)         # brant() -- proportional-odds (parallel regression) test
library(pscl)          # pR2()   -- pseudo-R2 statistics
library(generalhoslem) # logitgof(), lipsitz.test(), pulkrob.chisq()
library(erer)          # ocME()  -- ordered-model marginal effects per category
library(DescTools)     # Desc(), PseudoR2()
library(stargazer)     # publication-quality table
library(ggplot2)
library(dplyr)
library(lmtest)        # coeftest, lrtest
library(car)           # linearHypothesis()
library(sandwich)

Sys.setenv(LANG = "en")
options(scipen = 100)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# =============================================================================
# SECTION 1: LOAD & MERGE
# =============================================================================

qog     <- read.csv("data/qog_std_cs_jan26.csv", stringsAsFactors = FALSE)
kinship <- read.csv("data/kinship_df.csv",        stringsAsFactors = FALSE)

cat("QoG dimensions:     ", dim(qog),     "\n")
cat("Kinship dimensions: ", dim(kinship), "\n")

# Inner join on ISO-3 country codes (ccodealp in QoG; isocode in kinship)
df <- merge(qog, kinship, by.x = "ccodealp", by.y = "isocode", all = FALSE)
cat("Countries after inner join:", nrow(df), "\n")


qog_lookup <- qog[, c("ccodealp", "cname")]
names(qog_lookup) <- c("iso", "qog_country")

kinship_lookup <- kinship
names(kinship_lookup)[names(kinship_lookup) == "isocode"] <- "iso"

merge_audit <- merge(qog_lookup, kinship_lookup, by = "iso", all = TRUE)

merge_audit$source_status <- ifelse(
  !is.na(merge_audit$qog_country) & !is.na(merge_audit$kinship_score),
  "Matched in both",
  ifelse(
    !is.na(merge_audit$qog_country) & is.na(merge_audit$kinship_score),
    "QoG only",
    "Kinship only"
  )
)

cat("\n--- Merge audit summary ---\n")
print(table(merge_audit$source_status))

cat("\n--- Countries not matched in both datasets ---\n")
print(merge_audit[merge_audit$source_status != "Matched in both", ],
      row.names = FALSE)

write.csv(merge_audit, "merge_audit_qog_kinship.csv", row.names = FALSE)

# =============================================================================
# SECTION 2: BUILD THE DEPENDENT VARIABLE  (VERIFY CODING FIRST!)
# =============================================================================
# Inspect how fh_status is stored in YOUR merged file before recoding.
cat("\n--- Raw fh_status values ---\n")
print(table(df$fh_status, useNA = "ifany"))

# QoG convention is typically NUMERIC: 1 = Free, 2 = Partly Free, 3 = Not Free.
# >>> CONFIRM against the codebook entry for fh_status. <<<
# We build an ordered factor running LOW->HIGH freedom so a POSITIVE coefficient
# means "more free":  Not Free (1) < Partly Free (2) < Free (3).

if (is.numeric(df$fh_status)) {
  df$freedom <- factor(df$fh_status,
                       levels = c(3, 2, 1),                 # 3=NotFree ... 1=Free
                       labels = c("Not Free", "Partly Free", "Free"),
                       ordered = TRUE)
} else {
  # If stored as text labels, just order them explicitly:
  df$freedom <- factor(df$fh_status,
                       levels = c("Not Free", "Partly Free", "Free"),
                       ordered = TRUE)
}

cat("\n--- Ordered DV (freedom) ---\n")
print(table(df$freedom, useNA = "ifany"))


# =============================================================================
# SECTION 3: MISSING-VALUE DIAGNOSTICS, COVERAGE SCREEN & CANDIDATE SELECTION
# =============================================================================
# Self-contained variable-selection record:
#   (a) missingness of EVERY variable in the merged dataset
#   (b) a coverage screen that removes variables too sparse to use
#   (c) how many variables survive as VIABLE candidates
#   (d) the THEORY-DRIVEN shortlist actually taken forward, each justified
#   (e) the candidate dataset + its maximum available sample size

# --- 3a. Missingness of every column in the merged data ----------------------
missing_all <- data.frame(
  variable    = names(df),
  n_missing   = sapply(df, function(x) sum(is.na(x))),
  pct_missing = round(sapply(df, function(x) mean(is.na(x))) * 100, 1),
  row.names   = NULL
)
missing_all <- missing_all[order(-missing_all$pct_missing), ]

cat("Total variables in merged dataset:", nrow(missing_all), "\n")
print(missing_all, row.names = FALSE)

# --- 3b. Coverage screen: how many variables are VIABLE? ---------------------
# A variable is viable only if its missingness is low enough not to collapse a
# global cross-section. Threshold set at 20% (adjustable).
threshold <- 20   # max % missing allowed

cat("\n--- Variables by missingness band (answers 'how many are viable?') ---\n")
bands <- cut(missing_all$pct_missing,
             breaks = c(-Inf, 5, 10, 20, 50, 80, Inf),
             labels = c("<=5%", "5-10%", "10-20%", "20-50%", "50-80%", ">80%"))
print(table(bands))

viable <- missing_all[missing_all$pct_missing <= threshold, ]
cat("\nVariables passing the", threshold, "% coverage screen:",
    nrow(viable), "of", nrow(missing_all), "\n")
# Everything above the threshold (mostly EU-only, survey, or historical series)
# is removed HERE, before any theoretical selection.

cat("\n--- First 100 variables passing the coverage screen ---\n")
print(head(viable[, c("variable", "n_missing", "pct_missing")], 100), row.names = FALSE)

write.csv(viable, "viable_variables_under_20pct_missing.csv", row.names = FALSE)
cat("\nFull viable-variable list saved to: viable_variables_under_20pct_missing.csv\n") # 534 variables left

# --- 3c. THEORY-DRIVEN candidate shortlist among viable variables ------------
# The coverage screen only tells us which variables have enough observations.
# It does NOT decide which variables belong in the model.
#
# From the viable pool, we manually keep variables that are defensible controls
# in a cross-country model of political freedom. The selection follows
# determinants-of-democracy theory, especially modernization theory, structural
# country characteristics, information access, economic openness, and rentier-
# state theory.
#
# Candidate variables:
#
#   fh_status
#     DEPENDENT VARIABLE -- Freedom House status, later recoded as an ordered
#     factor: Not Free < Partly Free < Free.
#
#   kinship_score
#     KEY EXPLANATORY VARIABLE -- kinship intensity. The expected sign is
#     negative: tighter kinship structures are expected to be associated with
#     lower political freedom.
#
#   wdi_gdpcapcon2015
#     ECONOMIC DEVELOPMENT -- GDP per capita. Included as the main
#     modernization-theory control and as the base variable for the interaction
#     with kinship. Because GDP and democracy may be jointly determined, the GDP
#     coefficient is treated as a control association, not as a causal estimate.
#
#   wdi_pop
#     COUNTRY SIZE -- population size. Included as a scale/governability control.
#
#   wdi_popurb
#     URBANISATION -- captures structural modernization and social mobilisation.
#
#   wdi_trade
#     TRADE OPENNESS -- captures exposure to international markets, norms, and
#     external pressure. It has relatively higher missingness, so it should be
#     checked in sensitivity analysis.
#
#   wdi_internet OR dr_ig
#     INFORMATION ACCESS -- internet use may support information flows and civic
#     mobilisation. It may also be partly endogenous to development and regime
#     type, so it is interpreted as a control, not as a causal effect.
#
#   wdi_unempilo
#     ECONOMIC STRESS -- unemployment may proxy economic grievance. This is the
#     weakest theoretical control and can be removed by the general-to-specific
#     procedure if insignificant.
#
#   wdi_oilrent
#     RENTIER-STATE CONTROL -- oil rents as % of GDP. Included because natural
#     resource rents can allow states to remain politically unfree despite high
#     income, urbanisation, or trade openness.
#
# Deliberately excluded:
#
#   ti_cpi, wbgi_gee, wbgi_rle, wbgi_vae, vdem_*, ciri_*, fh_*
#     BAD CONTROLS / OUTCOME-LIKE VARIABLES -- these are corruption,
#     governance, rights, or democracy measures. They are too close to the
#     dependent variable or may lie on the causal pathway:
#     kinship -> institutions/state quality/rights -> freedom.
#
#   bmr_demdur
#     CIRCULAR CONTROL -- democratic duration mechanically overlaps with current
#     freedom/democracy status.
#
#   wdi_lifexp, wdi_fertility and similar development outcomes
#     OVER-CONTROL RISK -- these are likely consequences of broader development
#     and may be correlated with kinship, so they are not used as baseline
#     controls.

key_vars <- c("fh_status", "kinship_score", "wdi_gdpcapcon2015",
              "wdi_popurb", "wdi_pop", "wdi_trade", "wdi_internet", "wdi_unempilo", "wdi_oilrent")

# --- 3d. Verify shortlisted variables exist and pass coverage screen ----------

missing_from_data <- setdiff(key_vars, names(df))

if (length(missing_from_data) > 0) {
  stop("These key_vars are not present in df: ",
       paste(missing_from_data, collapse = ", "))
}

key_check <- missing_all[missing_all$variable %in% key_vars, ]
key_check$viable <- key_check$pct_missing <= threshold

cat("\n--- Coverage check on the theory-driven candidate shortlist ---\n")
print(key_check[order(key_check$pct_missing), ], row.names = FALSE)

stopifnot(all(key_check$viable))

cat("\nAll", length(key_vars), "candidate variables pass the",
    threshold, "% screen.\n")


# --- 3e. Candidate dataset and maximum complete-case sample ------------------

id_vars <- c("cname", "ccodealp")
df_candidates <- df[, c(id_vars, key_vars)]

cat("\nCandidate dataset dimensions:", dim(df_candidates), "\n")

cat("Complete cases across ALL candidate model variables:",
    sum(complete.cases(df[, key_vars])), "of", nrow(df), "\n")


candidate_coverage <- data.frame(
  variable = key_vars,
  n_missing = sapply(df[key_vars], function(x) sum(is.na(x))),
  pct_missing = round(sapply(df[key_vars], function(x) mean(is.na(x))) * 100, 1),
  viable = sapply(df[key_vars], function(x) mean(is.na(x)) * 100 <= threshold),
  row.names = NULL
)

candidate_coverage <- candidate_coverage[order(candidate_coverage$pct_missing), ]

print(candidate_coverage, row.names = FALSE)


candidate_sample_summary <- data.frame(
  total_countries_after_merge = nrow(df),
  complete_cases_all_candidate_variables = sum(complete.cases(df[, key_vars])),
  countries_lost_due_to_missingness = nrow(df) - sum(complete.cases(df[, key_vars])),
  pct_complete = round(mean(complete.cases(df[, key_vars])) * 100, 1)
)

print(candidate_sample_summary, row.names = FALSE)

# =============================================================================
# SECTION 4: VARIABLE CONSTRUCTION AND COMPLETE-CASE MODEL SAMPLE
# =============================================================================
# In Section 3 we only diagnosed missingness and selected the candidate variables.
# Here we actually construct the variables used in the models and remove countries
# with missing values in any model variable.

# --- Check raw distributions before transformations --------------------------

dist_raw <- df %>%
  transmute(
    gdppc    = wdi_gdpcapcon2015,
    pop      = wdi_pop,
    urban_pct = wdi_popurb,
    trade    = wdi_trade,
    internet = wdi_internet,
    unemp    = wdi_unempilo,
    oilrent  = wdi_oilrent,
    kinship  = kinship_score
  ) %>%
  tidyr::pivot_longer(everything(),
                      names_to = "variable",
                      values_to = "value")

ggplot(dist_raw, aes(x = value)) +
  geom_histogram(bins = 25, fill = "grey70", colour = "white") +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(title = "Raw Distributions of Candidate Variables",
       x = NULL,
       y = "Number of countries") +
  theme_minimal(base_size = 12)


# --- 4a. Construct model variables before removing missing values -------------

df_model_raw <- df %>%
  transmute(
    country   = cname,
    iso       = ccodealp,
    freedom   = freedom,                  # ordered DV: Not Free < Partly Free < Free
    kinship   = kinship_score,            # key explanatory variable [0,1]
    log_gdppc = log(wdi_gdpcapcon2015),   # log GDP per capita
    log_pop   = log(wdi_pop),             # log population
    urban_pct = wdi_popurb,               # urbanisation (%)
    trade     = wdi_trade,                # trade openness (% of GDP)
    internet  = wdi_internet,             # internet users (%)
    unemp     = wdi_unempilo,             # unemployment (%)
    log_oilrent = log1p(wdi_oilrent)      # oil rents (% of GDP)
  )

cat("\nModel dataset before complete-case filtering:",
    nrow(df_model_raw), "countries\n")


# --- 4b. Check whether logs created invalid values ----------------------------
# This should usually be zero, but it is safer to verify.

cat("\n--- Invalid log values check ---\n")
cat("Missing/invalid log GDP per capita:",
    sum(is.na(df_model_raw$log_gdppc) | is.infinite(df_model_raw$log_gdppc)), "\n")
cat("Missing/invalid log population:",
    sum(is.na(df_model_raw$log_pop) | is.infinite(df_model_raw$log_pop)), "\n")


# --- 4c. Identify countries excluded due to missing model variables -----------

model_vars <- c("freedom", "kinship", "log_gdppc", "log_pop", "urban_pct",
                "trade", "internet", "unemp", "log_oilrent")

excluded_missing <- df_model_raw[!complete.cases(df_model_raw[, model_vars]), ]

cat("\nCountries excluded due to missing model variables:",
    nrow(excluded_missing), "\n")

print(excluded_missing[, c("iso", "country", model_vars)], row.names = FALSE)


# Optional: clearer missingness map by country
missingness_by_country <- excluded_missing[, c("iso", "country")]

for (v in model_vars) {
  missingness_by_country[[paste0("missing_", v)]] <- is.na(excluded_missing[[v]])
}

cat("\n--- Missingness pattern among excluded countries ---\n")
print(missingness_by_country, row.names = FALSE)


# --- 4d. Create final complete-case modelling sample --------------------------

df_model <- df_model_raw %>%
  filter(complete.cases(across(all_of(model_vars))))

cat("\nFinal complete-case model sample:",
    nrow(df_model), "countries\n")


# --- 4e. Centre GDP after defining the final sample ---------------------------
# This makes the kinship coefficient interpretable at average income.

df_model$log_gdppc_c <- df_model$log_gdppc - mean(df_model$log_gdppc)

cat("\nMean log GDP per capita used for centering:",
    round(mean(df_model$log_gdppc), 3), "\n")

cat("\nOutcome distribution in final model sample:\n")
print(table(df_model$freedom))


# --- 4f. Sample-loss summary --------------------------------------------------

sample_loss_summary <- data.frame(
  stage = c("Merged dataset", "Complete-case model sample", "Excluded due to missing model variables"),
  n_countries = c(nrow(df_model_raw), nrow(df_model), nrow(excluded_missing)),
  pct_of_merged = round(c(
    100,
    nrow(df_model) / nrow(df_model_raw) * 100,
    nrow(excluded_missing) / nrow(df_model_raw) * 100
  ), 1)
)

cat("\n--- Sample-loss summary ---\n")
print(sample_loss_summary, row.names = FALSE)

write.csv(df_model, "data/df_model_clean_complete_cases.csv", row.names = FALSE)
df_model

# =============================================================================
# SECTION 5: DESCRIPTIVE STATISTICS AND EXPLORATORY DATA ANALYSIS
# =============================================================================

# --- 5a. Descriptive statistics ----------------------------------------------

Desc(df_model$freedom, main = "Freedom Status (ordered DV)")
Desc(df_model$kinship, main = "Kinship intensity score")

continuous_vars <- c("kinship", "log_gdppc", "log_pop", "urban_pct",
                     "trade", "internet", "unemp", "log_oilrent")

cat("\n--- Summary statistics: continuous variables ---\n")
print(summary(df_model[, continuous_vars]))


# --- 5b. Outcome distribution ------------------------------------------------

ggplot(df_model, aes(x = freedom)) +
  geom_bar(fill = "grey70", colour = "white") +
  labs(title = "Distribution of Freedom Status",
       x = "Freedom status",
       y = "Number of countries") +
  theme_minimal(base_size = 13)


# --- 5c. Kinship distribution ------------------------------------------------

ggplot(df_model, aes(x = kinship)) +
  geom_histogram(bins = 20, fill = "grey70", colour = "white") +
  labs(title = "Distribution of Kinship Intensity",
       x = "Kinship score",
       y = "Number of countries") +
  theme_minimal(base_size = 13)


# --- 5d. Kinship by freedom category -----------------------------------------

ggplot(df_model, aes(x = freedom, y = kinship)) +
  geom_boxplot(fill = "grey70", colour = "grey30") +
  labs(title = "Kinship Intensity by Freedom Status",
       x = "Freedom status",
       y = "Kinship score") +
  theme_minimal(base_size = 13)


# --- 5e. GDP by freedom category ---------------------------------------------

ggplot(df_model, aes(x = freedom, y = log_gdppc)) +
  geom_boxplot(fill = "grey70", colour = "grey30") +
  labs(title = "GDP per Capita by Freedom Status",
       x = "Freedom status",
       y = "Log GDP per capita") +
  theme_minimal(base_size = 13)


# --- 5f. Correlation matrix for continuous regressors -------------------------

cor_matrix <- round(cor(df_model[, continuous_vars], use = "complete.obs"), 2)

cor_df <- as.data.frame(as.table(cor_matrix))
names(cor_df) <- c("Variable_1", "Variable_2", "Correlation")

ggplot(cor_df, aes(x = Variable_1, y = Variable_2, fill = Correlation)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = Correlation), size = 3) +
  scale_fill_gradient2(midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Correlation Matrix of Continuous Regressors",
       x = NULL,
       y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# =============================================================================
# SECTION 6: GENERAL MODELS -- LPM, ORDERED LOGIT, ORDERED PROBIT (req. a, d)
# =============================================================================
# Following the lab approach, we first estimate a general ordered model.
# The LPM is reported only as a benchmark; ordered logit/probit are preferred
# because the dependent variable is ordinal: Not Free < Partly Free < Free.

# Quick checks before estimation
stopifnot(is.ordered(df_model$freedom))
stopifnot(all(complete.cases(df_model)))

cat("\nModel sample used in Section 6:", nrow(df_model), "countries\n")
cat("Outcome coding:\n")
print(table(df_model$freedom))

# General specification
form_general <- freedom ~ kinship + log_gdppc_c + log_pop + urban_pct +
  trade + internet + unemp + log_oilrent + kinship:log_gdppc_c

# LPM: treat ordered outcome as numeric 1/2/3, benchmark only
df_model$freedom_num <- as.numeric(df_model$freedom)   # Not Free = 1 ... Free = 3
LPM_general <- lm(update(form_general, freedom_num ~ .), data = df_model)

# Ordered logit and ordered probit
ologit_general  <- polr(form_general, data = df_model,
                        method = "logistic", Hess = TRUE)

oprobit_general <- polr(form_general, data = df_model,
                        method = "probit", Hess = TRUE)

cat("\n=== GENERAL LPM benchmark ===\n")
print(summary(LPM_general))

cat("\n=== GENERAL ordered logit ===\n")
print(summary(ologit_general))

cat("\n=== GENERAL ordered probit ===\n")
print(summary(oprobit_general))


library(texreg)

screenreg(
  list(LPM_general, ologit_general, oprobit_general),
  custom.model.names = c("LPM", "Ordered logit", "Ordered probit"),
  digits = 3,
  stars = c(0.001, 0.01, 0.05, 0.1),
  custom.note = "*** p < 0.001; ** p < 0.01; * p < 0.05; . p < 0.1"
)

htmlreg(
  list(LPM_general, ologit_general, oprobit_general),
  file = "section6_general_models_table.html",
  custom.model.names = c("LPM", "Ordered logit", "Ordered probit"),
  digits = 3,
  stars = c(0.001, 0.01, 0.05, 0.1),
  custom.note = "*** p < 0.001; ** p < 0.01; * p < 0.05; . p < 0.1",
  caption = "General model comparison: LPM, ordered logit, and ordered probit",
  caption.above = TRUE
)


# Helper: p-values for polr slope coefficients
# polr() reports t-values, interpreted here using the standard normal approximation.
polr_pvals <- function(model) {
  ct <- coef(summary(model))
  slopes <- ct[!rownames(ct) %in% names(model$zeta), , drop = FALSE]
  pv <- 2 * pnorm(abs(slopes[, "t value"]), lower.tail = FALSE)
  setNames(pv, rownames(slopes))
}

cat("\n--- General ordered logit: coefficient p-values ---\n")
print(round(polr_pvals(ologit_general), 4))

cat("\n--- General ordered probit: coefficient p-values ---\n")
print(round(polr_pvals(oprobit_general), 4))


# Compare logit and probit by information criteria
cat("\n--- Information criteria: logit vs probit ---\n")
cat("AIC  logit/probit:", AIC(ologit_general),  AIC(oprobit_general), "\n")
cat("BIC  logit/probit:", BIC(ologit_general),  BIC(oprobit_general), "\n")

# We report both ordered logit and ordered probit.
# The ordered logit is carried forward for the main interpretation because its
# coefficients can be discussed in odds terms and because the Brant test applies
# directly to the proportional-odds logit specification.

cat("\n--- LR test: ordered logit general model vs null ---\n")
ologit_null <- polr(freedom ~ 1, data = df_model,
                    method = "logistic", Hess = TRUE)
print(lrtest(ologit_general, ologit_null))

cat("\n--- LR test: ordered probit general model vs null ---\n")
oprobit_null <- polr(freedom ~ 1, data = df_model,
                     method = "probit", Hess = TRUE)
print(lrtest(oprobit_general, oprobit_null))

# =============================================================================
# SECTION 7: GENERAL-TO-SPECIFIC SELECTION (req. b)
# =============================================================================
# Lab-style approach:
#   1. Start from the general ordered logit model.
#   2. Identify the most insignificant non-protected variable.
#   3. Drop one variable at a time.
#   4. After each step, compare the restricted model with the ORIGINAL general model.
#   5. If the LR test p-value >= 0.05, the restriction is accepted.
#   6. If the LR test p-value < 0.05, stop and keep the previous accepted model.
#
# Ordered logit is used for model selection because it had slightly better AIC/BIC
# and allows proportional-odds diagnostics. Ordered probit will be estimated later
# using the same final formula as a robustness check.

# --- Helper: coefficient table with p-values for polr models -----------------

polr_table <- function(model) {
  ct <- coef(summary(model))
  slopes <- ct[!rownames(ct) %in% names(model$zeta), , drop = FALSE]
  
  p_values <- 2 * pnorm(abs(slopes[, "t value"]), lower.tail = FALSE)
  
  out <- data.frame(
    Estimate  = slopes[, "Value"],
    Std_Error = slopes[, "Std. Error"],
    t_value   = slopes[, "t value"],
    p_value   = p_values,
    Signif = cut(
      p_values,
      breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
      labels = c("***", "**", "*", ".", "")
    ),
    row.names = rownames(slopes)
  )
  
  # Round only numeric columns
  numeric_cols <- sapply(out, is.numeric)
  out[numeric_cols] <- round(out[numeric_cols], 4)
  
  return(out)
}


# --- Step 0: Original general ordered logit model -----------------------------
# This is the unrestricted model from Section 6.

gts_general <- ologit_general

cat("\n=== Step 0: GENERAL ordered logit model ===\n")
print(summary(gts_general))

cat("\n--- Step 0 coefficient table ---\n")
print(polr_table(gts_general))

ologit_drop_all_insig <- polr(
  freedom ~ kinship + log_gdppc_c + log_pop +
    trade + log_oilrent + kinship:log_gdppc_c,
  data = df_model,
  method = "logistic",
  Hess = TRUE
)

anova(gts_general, ologit_drop_all_insig)
# can drop all three insignificant variables together, but we still follow the general-to-specific one by one

# --- Step 1: Drop unemp ------------------------------------------------------

ologit_gts_1 <- polr(
  freedom ~ kinship + log_gdppc_c + log_pop + urban_pct +
    trade + internet + log_oilrent + kinship:log_gdppc_c,
  data = df_model,
  method = "logistic",
  Hess = TRUE
)

cat("\n=== Step 1: Ordered logit without unemp ===\n")
print(summary(ologit_gts_1))

cat("\n--- Step 1 coefficient table ---\n")
print(polr_table(ologit_gts_1))

cat("\n--- LR test: Can we drop unemp? ---\n")
print(anova(gts_general, ologit_gts_1))
# we can remove unemp

# --- Step 2: Drop urban_pct as well ------------------------------------------
# We now test whether unemp and urban_pct can be jointly removed.

ologit_gts_2 <- polr(
  freedom ~ kinship + log_gdppc_c + log_pop +
    trade + internet + log_oilrent + kinship:log_gdppc_c,
  data = df_model,
  method = "logistic",
  Hess = TRUE
)

cat("\n=== Step 2: Ordered logit without unemp and urban_pct ===\n")
print(summary(ologit_gts_2))

cat("\n--- Step 2 coefficient table ---\n")
print(polr_table(ologit_gts_2))

cat("\n--- LR test: Can we jointly drop unemp and urban_pct? ---\n")
print(anova(gts_general, ologit_gts_2))
# we can remove urban_pct

# --- Step 3: Drop internet as well -------------------------------------------
# We now test whether unemp, urban_pct, and internet can be jointly removed.

ologit_gts_3 <- polr(
  freedom ~ kinship + log_gdppc_c + log_pop +
    trade + log_oilrent + kinship:log_gdppc_c,
  data = df_model,
  method = "logistic",
  Hess = TRUE
)

cat("\n=== Step 3: Ordered logit without unemp, urban_pct, and internet ===\n")
print(summary(ologit_gts_3))

cat("\n--- Step 3 coefficient table ---\n")
print(polr_table(ologit_gts_3))

cat("\n--- LR test: Can we jointly drop unemp, urban_pct, and internet? ---\n")
print(anova(gts_general, ologit_gts_3))
# we can remove internet

# --- Final model if Step 3 was accepted --------------------------------------

ologit_final <- ologit_gts_3
final_form <- formula(ologit_final)

cat("\n=== FINAL selected ordered logit model ===\n")
print(summary(ologit_final))

cat("\n--- Final ordered logit coefficient table ---\n")
print(polr_table(ologit_final))

# --- Matching final probit and LPM using the same final formula ---------------

oprobit_final <- polr(
  final_form,
  data = df_model,
  method = "probit",
  Hess = TRUE
)

LPM_final <- lm(update(final_form, freedom_num ~ .), data = df_model)

cat("\n=== FINAL matching ordered probit model ===\n")
print(summary(oprobit_final))

cat("\n--- Final ordered probit coefficient table ---\n")
print(polr_table(oprobit_final))

cat("\n=== FINAL matching LPM benchmark ===\n")
print(summary(LPM_final))

# --- Compare general and final models ----------------------------------------

cat("\n--- Information criteria: general vs final ordered logit ---\n")
cat("AIC general/final:", AIC(ologit_general), AIC(ologit_final), "\n")
cat("BIC general/final:", BIC(ologit_general), BIC(ologit_final), "\n")

cat("\n--- LR test: final ordered logit vs null ---\n")
ologit_null <- polr(freedom ~ 1, data = df_model,
                    method = "logistic", Hess = TRUE)

print(lrtest(ologit_final, ologit_null))

# --- Publication-style comparison table --------------------------------------

screenreg(
  list(LPM_general, ologit_general, oprobit_general,
       LPM_final, ologit_final, oprobit_final),
  custom.model.names = c("LPM general", "Logit general", "Probit general",
                         "LPM final", "Logit final", "Probit final"),
  digits = 3,
  stars = c(0.001, 0.01, 0.05, 0.1),
  custom.note = "*** p < 0.001; ** p < 0.01; * p < 0.05; . p < 0.1"
)

htmlreg(
  list(LPM_general, ologit_general, oprobit_general,
       LPM_final, ologit_final, oprobit_final),
  file = "section7_general_vs_final_models.html",
  custom.model.names = c("LPM general", "Logit general", "Probit general",
                         "LPM final", "Logit final", "Probit final"),
  digits = 3,
  stars = c(0.001, 0.01, 0.05, 0.1),
  custom.note = "*** p < 0.001; ** p < 0.01; * p < 0.05; . p < 0.1",
  caption = "General-to-specific model selection: general and final specifications",
  caption.above = TRUE
)

# --- Multicollinearity check using final LPM as proxy -------------------------

cat("\n--- VIF check on final LPM specification ---\n")
print(vif(LPM_final))

cat("\n--- VIF check with interaction-aware option ---\n")
print(vif(LPM_final, type = "predictor"))

# save
htmlreg(
  list(LPM_general, ologit_general, oprobit_general,
       LPM_final, ologit_final, oprobit_final),
  file = "section7_general_vs_final_models.html",
  custom.model.names = c("LPM general", "Logit general", "Probit general",
                         "LPM final", "Logit final", "Probit final"),
  digits = 3,
  stars = c(0.001, 0.01, 0.05, 0.1),
  custom.note = "*** p < 0.001; ** p < 0.01; * p < 0.05; . p < 0.1",
  caption = "General-to-specific model selection: general and final specifications",
  caption.above = TRUE
)

# =============================================================================
# SECTION 8: MODEL DIAGNOSTICS FOR FINAL ORDERED LOGIT
# =============================================================================

# --- 8a. Brant test: proportional odds / parallel regression assumption -------
# H0: proportional odds / parallel regression assumption holds
# H1: proportional odds assumption is violated
cat("\n--- Brant test: proportional odds assumption ---\n")
brant(ologit_final)

# --- 8b. Lipsitz goodness-of-fit test ----------------------------------------

cat("\n--- Hosmer-Lemeshow for ordered models (logitgof) ---\n")
print(tryCatch(
  logitgof(df_model$freedom, fitted(ologit_final), g = 10, ord = TRUE),
  error = function(e) paste("HL note:", e$message)))

cat("\n--- Lipsitz test ---\n")
print(tryCatch(lipsitz.test(ologit_final),
               error = function(e) paste("Lipsitz note:", e$message)))

# add that we cannot use Pulkstenis-Robinson test because we don't have the categorical predictor

# --- 8c. Pseudo-R2 measures --------------------------------------------------

cat("\n--- Pseudo-R2 for final ordered logit ---\n")
pR2(ologit_final)

cat("\n--- Pseudo-R2 for final ordered probit ---\n")
pR2(oprobit_final)


# =============================================================================
# SECTION 9: HYPOTHESIS TESTING THROUGH PREDICTED PROBABILITIES
# =============================================================================
# We use predicted probabilities from the final ordered logit model.
# Ordered probit is used as robustness, but the interpretation is based mainly
# on the final ordered logit specification selected in Section 7.

# Use ordered probit for predicted probabilities because the Brant test
# indicated problems with the proportional-odds assumption in ordered logit.
pred_model <- oprobit_final

# --- 9a. Helper function: predicted probabilities ----------------------------

predict_probs <- function(model, newdata) {
  probs <- as.data.frame(predict(model, newdata = newdata, type = "probs"))
  cbind(newdata, probs)
}


# --- 9b. Define baseline values for controls ---------------------------------
# Continuous controls are held at their sample means.
# GDP is centred, so log_gdppc_c = 0 means average GDP per capita.

base_log_pop     <- mean(df_model$log_pop)
base_trade       <- mean(df_model$trade)
base_log_oilrent <- mean(df_model$log_oilrent)


# =============================================================================
# H1 and H2: Kinship effect and kinship-GDP interaction
# =============================================================================
# H1: Higher kinship intensity is expected to reduce political freedom.
# H2: The kinship effect is expected to depend on national wealth.

gdp_levels <- data.frame(
  gdp_label = c("Low GDP", "Average GDP", "High GDP"),
  log_gdppc_c = c(
    quantile(df_model$log_gdppc_c, 0.25),
    0,
    quantile(df_model$log_gdppc_c, 0.75)
  )
)

kinship_grid <- expand.grid(
  kinship = seq(0, 1, length.out = 100),
  gdp_label = gdp_levels$gdp_label
)

kinship_grid$log_gdppc_c <- gdp_levels$log_gdppc_c[
  match(kinship_grid$gdp_label, gdp_levels$gdp_label)
]

kinship_grid$log_pop <- base_log_pop
kinship_grid$trade <- base_trade
kinship_grid$log_oilrent <- base_log_oilrent

pred_kinship <- predict_probs(pred_model, kinship_grid)

pred_kinship_long <- pred_kinship %>%
  tidyr::pivot_longer(
    cols = c("Not Free", "Partly Free", "Free"),
    names_to = "freedom_status",
    values_to = "probability"
  )

ggplot(pred_kinship_long,
       aes(x = kinship, y = probability, linetype = gdp_label)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ freedom_status) +
  labs(title = "Predicted Freedom Probabilities by Kinship and GDP Level",
       x = "Kinship intensity",
       y = "Predicted probability",
       linetype = "GDP level") +
  theme_minimal(base_size = 12)


# --- 9c. Compact prediction table for H1/H2 ----------------------------------

kinship_table_data <- expand.grid(
  kinship = c(0.1, 0.9),
  gdp_label = gdp_levels$gdp_label
)

kinship_table_data$log_gdppc_c <- gdp_levels$log_gdppc_c[
  match(kinship_table_data$gdp_label, gdp_levels$gdp_label)
]

kinship_table_data$log_pop <- base_log_pop
kinship_table_data$trade <- base_trade
kinship_table_data$log_oilrent <- base_log_oilrent

kinship_pred_table <- predict_probs(pred_model, kinship_table_data)

cat("\n--- Predicted probabilities: low vs high kinship at different GDP levels ---\n")

kinship_pred_table_print <- kinship_pred_table

num_cols <- sapply(kinship_pred_table_print, is.numeric)
kinship_pred_table_print[num_cols] <- round(kinship_pred_table_print[num_cols], 3)

print(kinship_pred_table_print)

# =============================================================================
# H3: Oil-rent effect
# =============================================================================
# H3: Higher oil wealth is expected to reduce political freedom.

oil_grid <- data.frame(
  kinship = mean(df_model$kinship),
  log_gdppc_c = 0,
  log_pop = base_log_pop,
  trade = base_trade,
  log_oilrent = seq(
    quantile(df_model$log_oilrent, 0.05),
    quantile(df_model$log_oilrent, 0.95),
    length.out = 100
  )
)

pred_oil <- predict_probs(pred_model, oil_grid)

pred_oil_long <- pred_oil %>%
  tidyr::pivot_longer(
    cols = c("Not Free", "Partly Free", "Free"),
    names_to = "freedom_status",
    values_to = "probability"
  )

ggplot(pred_oil_long,
       aes(x = log_oilrent, y = probability)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ freedom_status) +
  labs(title = "Predicted Freedom Probabilities by Oil-Rent Dependence",
       x = "Log(1 + oil rents)",
       y = "Predicted probability") +
  theme_minimal(base_size = 12)


# --- 9d. Compact prediction table for H3 -------------------------------------

oil_table_data <- data.frame(
  kinship = mean(df_model$kinship),
  log_gdppc_c = 0,
  log_pop = base_log_pop,
  trade = base_trade,
  log_oilrent = quantile(df_model$log_oilrent, c(0.10, 0.50, 0.90))
)

oil_table_data$oil_level <- c("Low oil rents", "Median oil rents", "High oil rents")

oil_pred_table <- predict_probs(pred_model, oil_table_data)

cat("\n--- Predicted probabilities: oil-rent scenarios ---\n")

oil_pred_table_print <- oil_pred_table

num_cols <- sapply(oil_pred_table_print, is.numeric)
oil_pred_table_print[num_cols] <- round(oil_pred_table_print[num_cols], 3)

print(oil_pred_table_print)


