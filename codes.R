# =============================================================================
# Kinship Intensity and Human Development: A Tobit Analysis
# =============================================================================
# Dependent variable  : undp_hdi  (UNDP Human Development Index, bounded [0,1])
# Key independent var : kinship_score (ekne, 0 = nuclear family, 1 = tight kinship)
# Model type          : Tobit (Limited Dependent Variable) — upper censoring at 1
# Data sources        : QoG Standard Dataset (Jan 2026) + kinship_df (ekne index)
# =============================================================================


# --- 0. Packages --------------------------------------------------------------

# Install if not already present
# install.packages(c("AER", "VGAM", "lmtest", "sandwich", "stargazer",
#                    "ggplot2", "dplyr", "tidyr", "car", "nortest"))

library(AER)        # tobit() function
library(VGAM)       # alternative tobit / marginal effects helpers
library(lmtest)     # resettest (Ramsey-RESET), bptest (Breusch-Pagan)
library(sandwich)   # heteroskedasticity-robust standard errors
library(stargazer)  # publication-quality regression tables
library(ggplot2)    # plots
library(dplyr)      # data wrangling
library(tidyr)      # pivot helpers
library(car)        # vif(), linearHypothesis()
library(nortest)    # Anderson-Darling normality test


# --- 1. Load data -------------------------------------------------------------

qog     <- read.csv("data/qog_std_cs_jan26.csv", stringsAsFactors = FALSE)
kinship <- read.csv("data/kinship_df.csv", stringsAsFactors = FALSE)

# Quick sanity check
dim(qog)      # should be ~194 rows, 1300+ columns
dim(kinship)  # should be 217 rows, 3 columns


# --- 2. Merge datasets --------------------------------------------------------

# Both datasets share ISO-3 country codes (ccodealp in QoG, isocode in kinship)
df <- merge(qog, kinship,
            by.x = "ccodealp",
            by.y = "isocode",
            all   = FALSE)   # inner join — keep only matched countries

cat("Countries after merge:", nrow(df), "\n")


# --- 3. Missing value diagnostics --------------------------------------------
# Check missingness across ALL columns in the merged dataset

missing_summary <- data.frame(
  variable = names(df),
  n_missing = sapply(df, function(x) sum(is.na(x))),
  pct_missing = round(sapply(df, function(x) mean(is.na(x))) * 100, 1)
)

# Print only columns that actually have missing values, sorted descending
missing_nonzero <- missing_summary[missing_summary$n_missing > 0, ]
missing_nonzero <- missing_nonzero[order(-missing_nonzero$n_missing), ]

cat("\n--- Missing values (columns with at least 1 NA) ---\n")
print(missing_nonzero, row.names = FALSE)

# Print missingness for the key variables we will use
key_vars <- c("undp_hdi", "kinship_score",
              "wdi_gdpcapcon2015", "wdi_pop", "wdi_unempilo",
              "wbgi_cce", "wbgi_rle", "vdem_libdem", "ti_cpi")

cat("\n--- Missing values for key model variables ---\n")
print(missing_summary[missing_summary$variable %in% key_vars, ],
      row.names = FALSE)


# --- 4. Select and transform variables ---------------------------------------

# Note on Tobit rationale:
# HDI is bounded [0,1]. Several high-income countries cluster near the ceiling
# (Norway 0.967, Switzerland 0.966, Iceland 0.964, etc.). Standard OLS would
# underestimate the true effect at the upper boundary. Upper-censored Tobit
# corrects for this. We set the upper limit at 1.

df_model <- df %>%
  transmute(
    country          = cname,
    iso              = ccodealp,
    
    # Dependent variable
    hdi              = undp_hdi,
    
    # Key independent variable
    kinship          = kinship_score,
    
    # Controls (log-transform skewed variables for linearity)
    log_gdppc        = log(wdi_gdpcapcon2015),   # World Bank GDP per capita (const. 2015 USD)
    log_pop          = log(wdi_pop),              # Population size
    unemp            = wdi_unempilo,              # ILO unemployment rate (%)
    
    # Additional governance context (for robustness checks — not in main model)
    rule_of_law      = wbgi_rle,
    ctrl_corruption  = wbgi_cce,
    libdem           = vdem_libdem
  ) %>%
  filter(!is.na(hdi),
         !is.na(kinship),
         !is.na(log_gdppc),
         !is.na(log_pop),
         !is.na(unemp))

cat("\nObservations in working sample:", nrow(df_model), "\n")

# Descriptive summary of working sample
cat("\n--- Descriptive statistics (working sample) ---\n")
summary(df_model[ , c("hdi","kinship","log_gdppc","log_pop","unemp")])


# --- 5. Exploratory visualisation -------------------------------------------

# Distribution of HDI — check for ceiling clustering
ggplot(df_model, aes(x = hdi)) +
  geom_histogram(bins = 30, fill = "#378ADD", colour = "white", alpha = 0.85) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#D85A30") +
  labs(title = "Distribution of HDI (working sample)",
       subtitle = "Dashed line = upper censoring point (1.0)",
       x = "Human Development Index", y = "Count") +
  theme_minimal(base_size = 13)

# Bivariate relationship: kinship vs HDI
ggplot(df_model, aes(x = kinship, y = hdi)) +
  geom_point(alpha = 0.6, colour = "#378ADD") +
  geom_smooth(method = "lm", se = TRUE, colour = "#D85A30") +
  labs(title = "Kinship intensity vs Human Development Index",
       x = "Kinship score (0 = nuclear, 1 = tight kinship)",
       y = "HDI") +
  theme_minimal(base_size = 13)


# --- 6. General-to-specific model selection ----------------------------------
# Start with all theoretically motivated variables, then drop non-significant ones.
# AIC/BIC guide model comparison; Ramsey-RESET checks functional form.

# 6a. GENERAL model (all controls + interaction term)
# Interaction: kinship may matter more in poor countries (tight kinship networks
# substitute for weak formal institutions when GDP is low)
tobit_general <- tobit(
  hdi ~ kinship + log_gdppc + log_pop + unemp + kinship:log_gdppc,
  left  = -Inf,   # no left censoring (min HDI in sample is 0.385)
  right = 1,      # upper censoring at 1
  data  = df_model
)
summary(tobit_general)

# 6b. INTERMEDIATE model — drop interaction if not significant
tobit_intermediate <- tobit(
  hdi ~ kinship + log_gdppc + log_pop + unemp,
  left  = -Inf,
  right = 1,
  data  = df_model
)
summary(tobit_intermediate)

# 6c. SPECIFIC (FINAL) model — drop further non-significant controls
# (Adjust based on your results — example keeps kinship + gdppc + pop)
tobit_final <- tobit(
  hdi ~ kinship + log_gdppc + log_pop,
  left  = -Inf,
  right = 1,
  data  = df_model
)
summary(tobit_final)

# Model comparison via AIC and BIC
AIC(tobit_general, tobit_intermediate, tobit_final)
BIC(tobit_general, tobit_intermediate, tobit_final)


# --- 7. Publication-quality table (all three models) -------------------------

stargazer(
  tobit_general, tobit_intermediate, tobit_final,
  type          = "text",                   # change to "latex" for paper
  title         = "Tobit Models: Effect of Kinship Intensity on HDI",
  dep.var.label = "Human Development Index (upper-censored at 1)",
  covariate.labels = c(
    "Kinship intensity",
    "Log GDP per capita",
    "Log population",
    "Unemployment rate",
    "Kinship × Log GDP p.c."
  ),
  column.labels   = c("General", "Intermediate", "Final"),
  omit.stat       = c("f","ser"),
  star.cutoffs    = c(0.10, 0.05, 0.01),
  notes           = "† p<0.10, * p<0.05, ** p<0.01. Upper censoring at HDI = 1.",
  notes.append    = FALSE,
  out             = "table_tobit_models.txt"
)


# --- 8. Marginal effects (three kinds required for Tobit) --------------------
# For a Tobit model E[y | x, y < c] we compute three marginal effects:
#   (1) Unconditional (population-average): dE[y*]/dx — latent index
#   (2) Conditional on being uncensored:    dE[y | y < 1]/dx
#   (3) Probability of being uncensored:    dP(y < 1)/dx

beta   <- coef(tobit_final)
sigma  <- tobit_final$scale
x_mean <- colMeans(model.matrix(tobit_final))   # evaluate at means

# Upper limit
upper <- 1
z     <- (upper - x_mean %*% beta) / sigma
phi_z <- dnorm(z)
Phi_z <- pnorm(z)

# (1) Unconditional marginal effects = beta_j (latent index slope)
me_unconditional <- beta[-length(beta)]   # drop log(sigma)
cat("\n--- (1) Unconditional marginal effects (on latent HDI) ---\n")
print(round(me_unconditional, 5))

# (2) Conditional marginal effects = beta_j * Phi(z)
me_conditional <- me_unconditional * as.numeric(Phi_z)
cat("\n--- (2) Conditional marginal effects (given HDI < 1) ---\n")
print(round(me_conditional, 5))

# (3) Marginal effect on probability of being uncensored = -phi(z)/sigma * beta_j
me_prob <- -(phi_z / sigma) * me_unconditional
cat("\n--- (3) Marginal effect on P(HDI reaches ceiling) ---\n")
print(round(me_prob, 5))


# --- 9. Hypotheses verification -----------------------------------------------

cat("\n--- Hypothesis tests (final model) ---\n")

# H1: kinship has a negative effect on HDI
# Coefficient on kinship should be significantly negative
# Read from summary(tobit_final) — t-statistic and p-value

# H2: GDP per capita positively predicts HDI (sanity / secondary check)
# Coefficient on log_gdppc should be significantly positive

# Formal Wald test for H1: kinship coefficient = 0
linearHypothesis(tobit_final, "kinship = 0")


# --- 10. Diagnostic tests -----------------------------------------------------

# We extract residuals from the underlying OLS on uncensored observations
# (AER::tobit residuals) for the diagnostic tests below.

resid_tobit <- residuals(tobit_final, type = "response")   # raw residuals
fitted_vals <- fitted(tobit_final)

# 10a. Ramsey-RESET test — checks for omitted nonlinearities / functional form
#      We run RESET on an auxiliary OLS (same spec) as a proxy;
#      Tobit RESET is typically done this way in practice.
ols_proxy <- lm(hdi ~ kinship + log_gdppc + log_pop, data = df_model)
resettest(ols_proxy, power = 2:3, type = "fitted")

# 10b. Breusch-Pagan test — homoscedasticity
bptest(ols_proxy)

# 10c. White's test — heteroscedasticity (no cross-terms assumption)
bptest(ols_proxy, ~ fitted(ols_proxy) + I(fitted(ols_proxy)^2))

# 10d. Breusch-Godfrey test — no autocorrelation (less critical in cross-section,
#      but run to satisfy checklist; residuals ordered by HDI)
bgtest(ols_proxy, order = 2)

# 10e. Normality of residuals (Anderson-Darling)
ad.test(resid_tobit)

# 10f. Linktest (Pregibon) — model specification
#      Regress DV on yhat and yhat^2; yhat^2 should be non-significant
yhat  <- fitted(ols_proxy)
yhat2 <- yhat^2
linktest_model <- lm(hdi ~ yhat + yhat2, data = df_model)
summary(linktest_model)
# Interpretation: if yhat2 is non-significant → model is correctly specified


# --- 11. Pseudo-R² statistics ------------------------------------------------

# McFadden's R²
loglik_full  <- logLik(tobit_final)
loglik_null  <- logLik(update(tobit_final, . ~ 1))
R2_mcfadden  <- 1 - as.numeric(loglik_full) / as.numeric(loglik_null)
cat("\nMcFadden R²:", round(R2_mcfadden, 4), "\n")

# Note: McKelvey-Zavoina, count R², and adjusted count R² are primarily for
# binary/ordered models. For Tobit, the conventional pseudo-R² is McFadden's,
# or report the variance explained in the latent variable:
R2_latent <- var(fitted_vals) / (var(fitted_vals) + sigma^2)
cat("R² (latent variable):", round(R2_latent, 4), "\n")


# --- 12. Robustness checks (optional, mention in appendix) --------------------

# 12a. OLS baseline — compare coefficient direction and magnitude
ols_baseline <- lm(hdi ~ kinship + log_gdppc + log_pop, data = df_model)
summary(ols_baseline)

# 12b. Add unemployment back (check sensitivity)
tobit_robustness <- tobit(
  hdi ~ kinship + log_gdppc + log_pop + unemp,
  left  = -Inf,
  right = 1,
  data  = df_model
)
summary(tobit_robustness)

# 12c. Compare Tobit vs OLS in one table (Appendix)
stargazer(
  ols_baseline, tobit_final,
  type          = "text",
  title         = "Appendix: OLS vs Tobit comparison",
  column.labels = c("OLS baseline", "Tobit (final)"),
  out           = "table_ols_vs_tobit.txt"
)

# =============================================================================
# End of script
# =============================================================================