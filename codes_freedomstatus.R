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

cat("Countries with missing Freedom House status:", sum(is.na(df$freedom)), "\n")

# =============================================================================
# SECTION 3: MISSING-VALUE DIAGNOSTICS (key variables)
# =============================================================================
key_vars <- c("fh_status", "kinship_score", "wdi_gdpcapcon2015", "wdi_pop",
              "wdi_popurb", "wdi_unempilo", "wdi_trade", "dr_ig")

missing_key <- data.frame(
  variable    = key_vars,
  n_missing   = sapply(df[key_vars], function(x) sum(is.na(x))),
  pct_missing = round(sapply(df[key_vars], function(x) mean(is.na(x))) * 100, 1),
  row.names   = NULL
)
cat("\n--- Missing values: key model variables ---\n")
print(missing_key, row.names = FALSE)


# =============================================================================
# SECTION 4: VARIABLE CONSTRUCTION  (complete-case working sample)
# =============================================================================
df_model <- df %>%
  transmute(
    country   = cname,
    iso       = ccodealp,
    freedom   = freedom,                  # ordered DV
    kinship   = kinship_score,            # key X  [0,1]
    log_gdppc = log(wdi_gdpcapcon2015),   # log GDP per capita
    log_pop   = log(wdi_pop),             # log population
    urban_pct = wdi_popurb,               # urbanisation (%)
    unemp     = wdi_unempilo,             # unemployment (%)
    trade     = wdi_trade,                # trade openness (% GDP)
    internet  = dr_ig                     # internet users (%)
  ) %>%
  filter(complete.cases(.))

# CENTRE log_gdppc so the kinship main effect is interpreted at MEAN income
df_model$log_gdppc_c <- df_model$log_gdppc - mean(df_model$log_gdppc)

cat("\nWorking sample (complete cases):", nrow(df_model), "countries\n")
cat("Outcome distribution:\n"); print(table(df_model$freedom))
# NOTE: if 'trade' (highest missingness) costs too many countries, drop it and
#       re-run; report the larger N in the paper.


# =============================================================================
# SECTION 5: DESCRIPTIVE STATISTICS (lab-standard Desc())
# =============================================================================
Desc(df_model$freedom, main = "Freedom Status (ordered DV)")
Desc(df_model$kinship, main = "Kinship intensity score (ekne)")

cat("\n--- Summary statistics (continuous regressors) ---\n")
summary(df_model[, c("kinship","log_gdppc","log_pop","urban_pct",
                     "unemp","trade","internet")])


# =============================================================================
# SECTION 6: VISUAL INSPECTION
# =============================================================================
# Plot 1 -- DV category counts
ggplot(df_model, aes(x = freedom, fill = freedom)) +
  geom_bar(alpha = 0.85, colour = "white") +
  labs(title = "Distribution of Freedom Status", x = NULL, y = "Count") +
  theme_minimal(base_size = 13) + theme(legend.position = "none")

# Plot 2 -- kinship distribution by freedom category (core relationship)
ggplot(df_model, aes(x = freedom, y = kinship, fill = freedom)) +
  geom_boxplot(alpha = 0.8) +
  labs(title = "Kinship intensity across Freedom Status",
       subtitle = "Expected: less free countries have higher kinship",
       x = NULL, y = "Kinship score (ekne)") +
  theme_minimal(base_size = 13) + theme(legend.position = "none")

# Plot 3 -- kinship vs income, coloured by freedom (motivates interaction)
ggplot(df_model, aes(x = kinship, y = log_gdppc, colour = freedom)) +
  geom_point(alpha = 0.8, size = 2.2) +
  labs(title = "Kinship vs income by freedom status",
       subtitle = "Motivates the kinship x log GDP interaction",
       x = "Kinship score (ekne)", y = "log GDP per capita") +
  theme_minimal(base_size = 13)

# Plot 4 -- correlation heatmap (continuous variables)
cor_vars   <- c("kinship","log_gdppc","log_pop","urban_pct","unemp","trade","internet")
cor_matrix <- round(cor(df_model[, cor_vars], use = "complete.obs"), 2)
cor_df     <- as.data.frame(as.table(cor_matrix)); names(cor_df) <- c("V1","V2","r")
ggplot(cor_df, aes(V1, V2, fill = r)) +
  geom_tile(colour = "white") + geom_text(aes(label = r), size = 3) +
  scale_fill_gradient2(low = "#C94040", mid = "white", high = "#3A7DC9",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Correlation matrix -- regressors", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# =============================================================================
# SECTION 7: GENERAL MODELS -- LPM, ORDERED LOGIT, ORDERED PROBIT (req. a, d)
# =============================================================================
form_general <- freedom ~ kinship + log_gdppc_c + log_pop + urban_pct +
  unemp + trade + internet + kinship:log_gdppc_c

# LPM: treat ordered outcome as numeric 1/2/3 (benchmark only; not preferred)
df_model$freedom_num <- as.numeric(df_model$freedom)   # NotFree=1 ... Free=3
LPM_general <- lm(update(form_general, freedom_num ~ .), data = df_model)

# Ordered logit and ordered probit
ologit_general  <- polr(form_general, data = df_model, method = "logistic", Hess = TRUE)
oprobit_general <- polr(form_general, data = df_model, method = "probit",   Hess = TRUE)

cat("\n=== GENERAL ordered logit ===\n");  print(summary(ologit_general))
cat("\n=== GENERAL ordered probit ===\n"); print(summary(oprobit_general))

# Helper: p-values for polr slope coefficients (polr gives no p-values directly)
polr_pvals <- function(model) {
  ct <- coef(summary(model))
  slopes <- ct[!rownames(ct) %in% names(model$zeta), , drop = FALSE]  # drop thresholds
  pv <- pnorm(abs(slopes[, "t value"]), lower.tail = FALSE) * 2
  setNames(pv, rownames(slopes))
}
cat("\n--- General ordered logit: p-values ---\n"); print(round(polr_pvals(ologit_general), 4))

# Choose logit vs probit by information criteria (lower = better)
cat("\nAIC  logit/probit:", AIC(ologit_general),  AIC(oprobit_general), "\n")
cat("BIC  logit/probit:", BIC(ologit_general),  BIC(oprobit_general), "\n")
# We carry the ordered LOGIT forward (brant + odds interpretation); the probit
# is reported alongside for comparison.


# =============================================================================
# SECTION 8: GENERAL-TO-SPECIFIC SELECTION (req. b)
# =============================================================================
# Rule (lab method): at each step drop the single most-insignificant control,
# then verify with anova() against the ORIGINAL GENERAL MODEL that ALL dropped
# variables are jointly = 0 (p >= 0.05 -> safe to drop). Protected from removal:
# kinship, log_gdppc_c, and the interaction (hierarchy principle).

protected <- c("kinship", "log_gdppc_c", "kinship:log_gdppc_c")

gts_polr <- function(general_model, data, protected, method = "probit",
                     alpha = 0.05) {
  full <- general_model
  current_terms <- attr(terms(formula(general_model)), "term.labels")
  repeat {
    mod <- polr(reformulate(current_terms, response = "freedom"),
                data = data, method = method, Hess = TRUE)
    pv  <- polr_pvals(mod)
    # candidate controls = current terms that are NOT protected and ARE signific…?
    cand <- setdiff(current_terms, protected)
    # map term -> its p-value (interaction term name may differ in coef table)
    cand_p <- pv[intersect(names(pv), cand)]
    cand_p <- cand_p[cand_p >= alpha]          # only insignificant ones
    if (length(cand_p) == 0) {
      cat("\nGTS STOP: all remaining controls significant.\n")
      return(mod)
    }
    drop_var <- names(which.max(cand_p))       # most insignificant
    reduced_terms <- setdiff(current_terms, drop_var)
    reduced <- polr(reformulate(reduced_terms, response = "freedom"),
                    data = data, method = method, Hess = TRUE)
    lr <- anova(full, reduced)                 # joint test vs GENERAL model
    p_joint <- lr$"Pr(Chi)"[2]
    cat(sprintf("\nStep: drop '%s' (p=%.3f) | joint test vs general p=%.3f -> %s\n",
                drop_var, max(cand_p), p_joint,
                ifelse(p_joint >= alpha, "DROP", "KEEP & STOP")))
    if (is.na(p_joint) || p_joint < alpha) {
      cat("Cannot jointly drop -> keep current model as final.\n")
      return(mod)
    }
    current_terms <- reduced_terms
  }
}

oprobit_final <- gts_polr(ologit_general, df_model, protected, method = "probit")
cat("\n=== FINAL ordered logit ===\n"); print(summary(oprobit_final))
cat("\n--- Final model p-values ---\n"); print(round(polr_pvals(oprobit_final), 4))

# Re-estimate the matching final probit and final LPM on the SAME final formula
final_form    <- formula(oprobit_final)
oprobit_final <- polr(final_form, data = df_model, method = "probit", Hess = TRUE)
LPM_final     <- lm(update(final_form, freedom_num ~ .), data = df_model)

cat("\n--- VIF (final LPM, as proxy for multicollinearity) ---\n")
print(vif(LPM_final))

# =============================================================================
# SECTION 9 (FIXED): PUBLICATION TABLE using texreg
# =============================================================================
install.packages("texreg")
library(texreg)

# Console display (equivalent to stargazer type="text")
screenreg(list(LPM_general, ologit_general, oprobit_general,
               LPM_final,   oprobit_final,   oprobit_final),
          custom.model.names = c("LPM-gen","oLogit-gen","oProbit-gen",
                                 "LPM-fin","oLogit-fin","oProbit-fin"),
          digits = 3)

# HTML output for pasting into Word (File > Open in Word)
htmlreg(list(LPM_general, ologit_general, oprobit_general,
             LPM_final,   oprobit_final,   oprobit_final),
        file = "results_table.html",
        custom.model.names = c("LPM-gen","oLogit-gen","oProbit-gen",
                               "LPM-fin","oLogit-fin","oProbit-fin"),
        digits = 3,
        caption = "Kinship and Freedom Status: general vs final models",
        caption.above = TRUE)

# If you want LaTeX instead:
# texreg(list(...), file = "results_table.tex", digits = 3)


# =============================================================================
# SECTION 10: HYPOTHESIS VERIFICATION (req. -- joint & single significance)
# =============================================================================
# (a) Joint significance of all regressors: final vs intercept-only
null_mod <- polr(freedom ~ 1, data = df_model, method = "logistic", Hess = TRUE)
cat("\n--- LR test: final model vs null (joint significance) ---\n")
print(lrtest(oprobit_final, null_mod))

# (b) Key hypothesis H1: kinship reduces freedom (single-coefficient test)
cat("\n--- Coefficient on kinship (H1) ---\n")
print(round(polr_pvals(oprobit_final)["kinship"], 4))

# (c) Secondary hypothesis H2: interaction kinship x income
cat("\n--- Interaction term (H2) ---\n")
print(round(polr_pvals(oprobit_final)[grep("kinship:", names(polr_pvals(oprobit_final)))], 4))


# =============================================================================
# SECTION 11: MARGINAL EFFECTS for the FINAL model, per category (req. e)
# =============================================================================
# Ordered marginal effects differ per outcome category and SUM TO ZERO across
# categories for each variable. Reported in percentage-point terms.
cat("\n--- Marginal effects (ordered logit, at means) ---\n")
me_ologit <- ocME(oprobit_final)          # erer::ocME -> ME per category
print(me_ologit)
# me_ologit$out holds the ME matrices; interpret e.g.:
# "A one-unit rise in kinship lowers P(Free) by X pp and raises P(Not Free) by Y pp."


# =============================================================================
# SECTION 12: PSEUDO-R2 STATISTICS (req. f)
# =============================================================================
cat("\n--- pscl::pR2 (McFadden etc.) ---\n")
print(pR2(oprobit_final))

cat("\n--- McKelvey-Zavoina & others (DescTools) ---\n")
print(tryCatch(
  PseudoR2(oprobit_final, which = c("McFadden","McKelveyZavoina","Nagelkerke","CoxSnell")),
  error = function(e) paste("PseudoR2 note:", e$message)))

# Count R2 and adjusted Count R2 (classification-based)
count_r2 <- function(model) {
  obs  <- model$model[[1]]
  pred <- predict(model, type = "class")
  n    <- length(obs)
  ncorrect <- sum(pred == obs)
  nmode    <- max(table(obs))                       # modal-category count
  c(CountR2     = ncorrect / n,
    AdjCountR2  = (ncorrect - nmode) / (n - nmode))
}
cat("\n--- Count R2 / Adjusted Count R2 ---\n")
print(round(count_r2(oprobit_final), 3))
# "The model correctly classifies about XX% of countries."


# =============================================================================
# SECTION 13: LINKTEST -- specification (req. g)
# =============================================================================
# We want _hat significant and _hatsq INSIGNIFICANT (no misspecification).
# Lab-standard: WNE::linktest(oprobit_final)  -- if the WNE package is installed.
# Manual fallback for ordered models:
linktest_ordered <- function(model) {
  beta <- coef(model)
  X    <- model.matrix(model)[, names(beta), drop = FALSE]
  yhat <- as.vector(X %*% beta)
  dd   <- data.frame(y = model$model[[1]], yhat = yhat, yhat2 = yhat^2)
  lt   <- polr(y ~ yhat + yhat2, data = dd, method = model$method, Hess = TRUE)
  ct   <- coef(summary(lt))
  ct   <- ct[c("yhat","yhat2"), , drop = FALSE]
  p    <- pnorm(abs(ct[, "t value"]), lower.tail = FALSE) * 2
  cbind(round(ct, 4), p.value = round(p, 4))
}
cat("\n--- Linktest (ordered, manual) ---\n")
print(tryCatch(linktest_ordered(oprobit_final),
               error = function(e) paste("linktest note:", e$message)))
# If installed:  WNE::linktest(oprobit_final)


# =============================================================================
# SECTION 14: GOODNESS-OF-FIT -- Hosmer-Lemeshow, Lipsitz, Pulkstenis-Robinson (req. h)
# =============================================================================
# All share H0: the model fits well (p >= 0.05 -> no evidence of poor fit).
cat("\n--- Hosmer-Lemeshow for ordered models (logitgof) ---\n")
print(tryCatch(
  logitgof(df_model$freedom, fitted(oprobit_final), g = 10, ord = TRUE),
  error = function(e) paste("HL note:", e$message)))

cat("\n--- Lipsitz test ---\n")
print(tryCatch(lipsitz.test(oprobit_final),
               error = function(e) paste("Lipsitz note:", e$message)))

cat("\n--- Pulkstenis-Robinson test ---\n")
# NOTE: pulkrob.chisq requires at least one CATEGORICAL predictor in the model.
# Our regressors are continuous, so this test may not run as-is. To enable it,
# add a categorical control (e.g. income tercile) and pass its name below.
# Example:
#   df_model$inc_grp <- cut(df_model$log_gdppc, quantile(df_model$log_gdppc,
#                           c(0,1/3,2/3,1)), include.lowest = TRUE,
#                           labels = c("low","mid","high"))
#   m_cat <- polr(update(final_form, . ~ . + inc_grp), data = df_model)
#   pulkrob.chisq(m_cat, c("inc_grp"))
print(tryCatch(pulkrob.chisq(oprobit_final, character(0)),
               error = function(e) paste("PR note: needs a categorical predictor -", e$message)))


# =============================================================================
# SECTION 15: PROPORTIONAL-ODDS ASSUMPTION -- Brant test (req. i)
# =============================================================================
# H0: parallel-regression (proportional-odds) assumption holds.
# Omnibus p >= 0.05 -> no evidence of violation.
cat("\n--- Brant test (proportional odds) ---\n")
print(tryCatch(brant(oprobit_final),
               error = function(e) paste("Brant note:", e$message)))
# If violated: consider generalized ordered logit (VGAM::vglm, cumulative,
# parallel = FALSE) or a partial-proportional-odds model.


# =============================================================================
# SECTION 16: ODDS RATIOS (optional, aids interpretation)
# =============================================================================
cat("\n--- Odds ratios (exp(beta)) for the final ordered logit ---\n")
print(round(exp(coef(oprobit_final)), 3))

# =============================================================================
# END
# =============================================================================