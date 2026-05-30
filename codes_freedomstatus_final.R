# =============================================================================
# Kinship intensity and political freedom: an ordered logit / probit analysis
# =============================================================================
# Dependent variable: Freedom House status, recoded as an ordered factor
#   Not Free < Partly Free < Free   (higher category = more free)
# Key explanatory variable: kinship intensity index (Enke 2019)
# Data: QoG Standard Dataset (Jan 2026) + kinship intensity index
# =============================================================================

library(MASS)          # polr(): ordered logit / probit
library(brant)         # brant(): proportional odds test
library(pscl)          # pR2(): pseudo R2
library(generalhoslem) # logitgof(), lipsitz.test(), pulkrob.chisq()
library(erer)          # ocME(): marginal effects for ordered models
library(DescTools)     # Desc(), PseudoR2()
library(texreg)        # screenreg(), htmlreg(): result tables
library(ggplot2)
library(dplyr)
library(tidyr)
library(lmtest)        # lrtest()
library(car)           # vif()

Sys.setenv(LANG = "en")
options(scipen = 100)

# Set the working directory to the script location (RStudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# =============================================================================
# SECTION 1: load and merge
# =============================================================================

qog     <- read.csv("data/qog_std_cs_jan26.csv", stringsAsFactors = FALSE)
kinship <- read.csv("data/kinship_df.csv",        stringsAsFactors = FALSE)

cat("QoG dimensions:    ", dim(qog),     "\n")
cat("Kinship dimensions:", dim(kinship), "\n")

# Inner join on ISO-3 codes: keep only countries present in both datasets,
# because the kinship score is the key explanatory variable.
df <- merge(qog, kinship, by.x = "ccodealp", by.y = "isocode", all = FALSE)
cat("Countries after inner join:", nrow(df), "\n")

# Merge audit: record which countries matched and which were dropped.
qog_lookup <- qog[, c("ccodealp", "cname")]
names(qog_lookup) <- c("iso", "qog_country")
kinship_lookup <- kinship
names(kinship_lookup)[names(kinship_lookup) == "isocode"] <- "iso"

merge_audit <- merge(qog_lookup, kinship_lookup, by = "iso", all = TRUE)
merge_audit$source_status <- ifelse(
  !is.na(merge_audit$qog_country) & !is.na(merge_audit$kinship_score), "Matched in both",
  ifelse(!is.na(merge_audit$qog_country) & is.na(merge_audit$kinship_score),
         "QoG only", "Kinship only"))

cat("\nMerge audit summary:\n")
print(table(merge_audit$source_status))

# The relevant exclusions are the QoG countries with no kinship score.
cat("\nQoG countries without a kinship score:\n")
print(merge_audit[merge_audit$source_status == "QoG only", c("iso", "qog_country")],
      row.names = FALSE)

write.csv(merge_audit, "merge_audit_qog_kinship.csv", row.names = FALSE)


# =============================================================================
# SECTION 2: build the dependent variable
# =============================================================================
# QoG codes fh_status numerically: 1 = Free, 2 = Partly Free, 3 = Not Free.
# We reverse it into an ordered factor so that a higher category means more free.

cat("\nRaw fh_status values:\n")
print(table(df$fh_status, useNA = "ifany"))

df$freedom <- factor(df$fh_status,
                     levels = c(3, 2, 1),
                     labels = c("Not Free", "Partly Free", "Free"),
                     ordered = TRUE)

cat("\nOrdered dependent variable (freedom):\n")
print(table(df$freedom, useNA = "ifany"))


# =============================================================================
# SECTION 3: missing value screen and candidate selection
# =============================================================================
# We first compute missingness for every variable, then keep only variables
# with acceptable coverage, and finally select a theory-driven candidate set.

missing_all <- data.frame(
  variable    = names(df),
  n_missing   = sapply(df, function(x) sum(is.na(x))),
  pct_missing = round(sapply(df, function(x) mean(is.na(x))) * 100, 1),
  row.names   = NULL
)
missing_all <- missing_all[order(-missing_all$pct_missing), ]
cat("\nTotal variables in merged dataset:", nrow(missing_all), "\n")
cat("Ten variables with the most missing values:\n")
print(head(missing_all, 10), row.names = FALSE)

# Coverage screen: a variable is viable if missingness is at most 20 percent.
threshold <- 20
bands <- cut(missing_all$pct_missing,
             breaks = c(-Inf, 5, 10, 20, 50, 80, Inf),
             labels = c("<=5%", "5-10%", "10-20%", "20-50%", "50-80%", ">80%"))
cat("\nVariables by missingness band:\n")
print(table(bands))

viable <- missing_all[missing_all$pct_missing <= threshold, ]
cat("\nVariables passing the", threshold, "percent screen:",
    nrow(viable), "of", nrow(missing_all), "\n")
write.csv(viable, "viable_variables_under_20pct_missing.csv", row.names = FALSE)

# Theory-driven candidate set. Passing the coverage screen is not enough:
# governance, democracy, rights and corruption indicators (vdem_*, wbgi_*,
# ciri_*, fh_*, bmr_*, ti_cpi, ...) are excluded because they are either
# alternative measures of the outcome or lie on the pathway from kinship to
# freedom, which would create bad-control bias. Democratic duration is
# excluded because it is mechanically tied to current freedom status, and
# development outcomes such as life expectancy or fertility are excluded to
# avoid over-controlling.
key_vars <- c("fh_status", "kinship_score", "wdi_gdpcapcon2015",
              "wdi_popurb", "wdi_pop", "wdi_trade", "wdi_internet",
              "wdi_unempilo", "wdi_oilrent")

missing_from_data <- setdiff(key_vars, names(df))
if (length(missing_from_data) > 0) {
  stop("These key_vars are not present in df: ",
       paste(missing_from_data, collapse = ", "))
}

key_check <- missing_all[missing_all$variable %in% key_vars, ]
key_check$viable <- key_check$pct_missing <= threshold
cat("\nCoverage check on the candidate set:\n")
print(key_check[order(key_check$pct_missing), ], row.names = FALSE)
stopifnot(all(key_check$viable))

# Maximum sample available across all candidate variables.
cat("\nComplete cases across all candidate variables:",
    sum(complete.cases(df[, key_vars])), "of", nrow(df), "\n")


# =============================================================================
# SECTION 4: variable construction and complete-case sample
# =============================================================================
# GDP per capita, population and oil rents are right skewed, so we log them.
# Oil rents contain zeros, so we use log(1 + oil rents) to keep those countries.

df_model_raw <- df %>%
  transmute(
    country     = cname,
    iso         = ccodealp,
    freedom     = freedom,
    kinship     = kinship_score,
    log_gdppc   = log(wdi_gdpcapcon2015),
    log_pop     = log(wdi_pop),
    urban_pct   = wdi_popurb,
    trade       = wdi_trade,
    internet    = wdi_internet,
    unemp       = wdi_unempilo,
    log_oilrent = log1p(wdi_oilrent)
  )

model_vars <- c("freedom", "kinship", "log_gdppc", "log_pop", "urban_pct",
                "trade", "internet", "unemp", "log_oilrent")

# Complete-case sample: a country is kept only if no model variable is missing.
df_model <- df_model_raw %>%
  filter(complete.cases(across(all_of(model_vars))))

cat("\nFinal complete-case sample:", nrow(df_model), "countries\n")

# Centre log GDP so the kinship main effect is read at average income.
df_model$log_gdppc_c <- df_model$log_gdppc - mean(df_model$log_gdppc)

# Sample-loss summary used in the paper.
sample_loss <- data.frame(
  stage = c("Merged dataset", "Complete-case model sample",
            "Excluded for missing model variables"),
  n_countries = c(nrow(df_model_raw), nrow(df_model),
                  nrow(df_model_raw) - nrow(df_model)),
  pct_of_merged = round(c(100,
                          nrow(df_model) / nrow(df_model_raw) * 100,
                          (nrow(df_model_raw) - nrow(df_model)) / nrow(df_model_raw) * 100), 1)
)
cat("\nSample-loss summary:\n")
print(sample_loss, row.names = FALSE)

cat("\nOutcome distribution in the final sample:\n")
print(table(df_model$freedom))

write.csv(df_model, "data/df_model_clean_complete_cases.csv", row.names = FALSE)


# =============================================================================
# SECTION 5: descriptive statistics and exploratory graphs
# =============================================================================

Desc(df_model$freedom, main = "Freedom status (ordered DV)")
Desc(df_model$kinship,  main = "Kinship intensity score")

continuous_vars <- c("kinship", "log_gdppc", "log_pop", "urban_pct",
                     "trade", "internet", "unemp", "log_oilrent")
cat("\nSummary statistics, continuous variables:\n")
print(summary(df_model[, continuous_vars]))

# Graph 1: distribution of the dependent variable.
g_freedom <- ggplot(df_model, aes(x = freedom)) +
  geom_bar(fill = "grey70", colour = "white") +
  labs(title = "Distribution of freedom status",
       x = "Freedom status", y = "Number of countries") +
  theme_minimal(base_size = 13)
print(g_freedom)
ggsave("fig_freedom_distribution.png", g_freedom, width = 7, height = 4.5, dpi = 300)

# Graph 2: kinship by freedom status (core descriptive relationship).
g_kin_box <- ggplot(df_model, aes(x = freedom, y = kinship)) +
  geom_boxplot(fill = "grey70", colour = "grey30") +
  labs(title = "Kinship intensity by freedom status",
       x = "Freedom status", y = "Kinship score") +
  theme_minimal(base_size = 13)
print(g_kin_box)
ggsave("fig_kinship_by_freedom.png", g_kin_box, width = 7, height = 4.5, dpi = 300)

# Graph 3: GDP per capita by freedom status.
g_gdp_box <- ggplot(df_model, aes(x = freedom, y = log_gdppc)) +
  geom_boxplot(fill = "grey70", colour = "grey30") +
  labs(title = "GDP per capita by freedom status",
       x = "Freedom status", y = "Log GDP per capita") +
  theme_minimal(base_size = 13)
print(g_gdp_box)
ggsave("fig_gdp_by_freedom.png", g_gdp_box, width = 7, height = 4.5, dpi = 300)

# Graph 4: correlation heatmap of continuous regressors.
cor_matrix <- round(cor(df_model[, continuous_vars], use = "complete.obs"), 2)
cor_df <- as.data.frame(as.table(cor_matrix))
names(cor_df) <- c("Variable_1", "Variable_2", "Correlation")
g_corr <- ggplot(cor_df, aes(x = Variable_1, y = Variable_2, fill = Correlation)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = Correlation), size = 3) +
  scale_fill_gradient2(midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Correlation matrix of continuous regressors", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(g_corr)
ggsave("fig_correlation_matrix.png", g_corr, width = 7, height = 6, dpi = 300)


# =============================================================================
# SECTION 6: general models (LPM, ordered logit, ordered probit)
# =============================================================================
stopifnot(is.ordered(df_model$freedom))

form_general <- freedom ~ kinship + log_gdppc_c + log_pop + urban_pct +
  trade + internet + unemp + log_oilrent + kinship:log_gdppc_c

# LPM benchmark only: treats the ordered outcome as numeric 1/2/3.
df_model$freedom_num <- as.numeric(df_model$freedom)
LPM_general     <- lm(update(form_general, freedom_num ~ .), data = df_model)
ologit_general  <- polr(form_general, data = df_model, method = "logistic", Hess = TRUE)
oprobit_general <- polr(form_general, data = df_model, method = "probit",   Hess = TRUE)

# Helper: coefficient table with p-values for polr models (polr reports only t).
polr_table <- function(model) {
  ct <- coef(summary(model))
  slopes <- ct[!rownames(ct) %in% names(model$zeta), , drop = FALSE]
  p_values <- 2 * pnorm(abs(slopes[, "t value"]), lower.tail = FALSE)
  out <- data.frame(
    Estimate  = round(slopes[, "Value"], 4),
    Std_Error = round(slopes[, "Std. Error"], 4),
    t_value   = round(slopes[, "t value"], 4),
    p_value   = round(p_values, 4),
    Signif = cut(p_values, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                 labels = c("***", "**", "*", ".", "")),
    row.names = rownames(slopes)
  )
  out
}

cat("\nGeneral ordered logit:\n");  print(polr_table(ologit_general))
cat("\nGeneral ordered probit:\n"); print(polr_table(oprobit_general))

cat("\nInformation criteria (logit / probit):\n")
cat("AIC:", AIC(ologit_general), "/", AIC(oprobit_general), "\n")
cat("BIC:", BIC(ologit_general), "/", BIC(oprobit_general), "\n")

# Joint significance against an intercept-only model.
ologit_null  <- polr(freedom ~ 1, data = df_model, method = "logistic", Hess = TRUE)
oprobit_null <- polr(freedom ~ 1, data = df_model, method = "probit",   Hess = TRUE)
cat("\nLR test, general ordered logit vs null:\n");  print(lrtest(ologit_general, ologit_null))
cat("\nLR test, general ordered probit vs null:\n"); print(lrtest(oprobit_general, oprobit_null))


# =============================================================================
# SECTION 7: general-to-specific selection
# =============================================================================
# We start from the general ordered logit, drop weak controls one at a time,
# and check each restriction with a likelihood ratio test against the original
# general model. kinship and log_gdppc_c are kept because they form the
# interaction term (hierarchy principle).

cat("\nStep 0, general ordered logit:\n")
print(polr_table(ologit_general))

# Step 1: drop unemp (least significant control).
ologit_gts_1 <- polr(freedom ~ kinship + log_gdppc_c + log_pop + urban_pct +
                       trade + internet + log_oilrent + kinship:log_gdppc_c,
                     data = df_model, method = "logistic", Hess = TRUE)
cat("\nStep 1, can we drop unemp?\n")
print(anova(ologit_general, ologit_gts_1))

# Step 2: also drop urban_pct.
ologit_gts_2 <- polr(freedom ~ kinship + log_gdppc_c + log_pop +
                       trade + internet + log_oilrent + kinship:log_gdppc_c,
                     data = df_model, method = "logistic", Hess = TRUE)
cat("\nStep 2, can we jointly drop unemp and urban_pct?\n")
print(anova(ologit_general, ologit_gts_2))

# Step 3: also drop internet.
ologit_gts_3 <- polr(freedom ~ kinship + log_gdppc_c + log_pop +
                       trade + log_oilrent + kinship:log_gdppc_c,
                     data = df_model, method = "logistic", Hess = TRUE)
cat("\nStep 3, can we jointly drop unemp, urban_pct and internet?\n")
print(anova(ologit_general, ologit_gts_3))

# All three restrictions are accepted (p well above 0.05), so this is the final model.
ologit_final <- ologit_gts_3
final_form   <- formula(ologit_final)

oprobit_final <- polr(final_form, data = df_model, method = "probit", Hess = TRUE)
LPM_final     <- lm(update(final_form, freedom_num ~ .), data = df_model)

cat("\nFinal ordered logit:\n");  print(polr_table(ologit_final))
cat("\nFinal ordered probit:\n"); print(polr_table(oprobit_final))

cat("\nInformation criteria, general vs final ordered logit:\n")
cat("AIC:", AIC(ologit_general), "/", AIC(ologit_final), "\n")
cat("BIC:", BIC(ologit_general), "/", BIC(ologit_final), "\n")
cat("\nLR test, final ordered logit vs null:\n")
print(lrtest(ologit_final, ologit_null))

# Publication table: general and final LPM, ordered logit, ordered probit.
screenreg(list(LPM_general, ologit_general, oprobit_general,
               LPM_final, ologit_final, oprobit_final),
          custom.model.names = c("LPM general", "Logit general", "Probit general",
                                 "LPM final", "Logit final", "Probit final"),
          digits = 3, stars = c(0.001, 0.01, 0.05, 0.1),
          custom.note = "*** p<0.001; ** p<0.01; * p<0.05; . p<0.1")

htmlreg(list(LPM_general, ologit_general, oprobit_general,
             LPM_final, ologit_final, oprobit_final),
        file = "table_general_vs_final.html",
        custom.model.names = c("LPM general", "Logit general", "Probit general",
                               "LPM final", "Logit final", "Probit final"),
        digits = 3, stars = c(0.001, 0.01, 0.05, 0.1),
        custom.note = "*** p<0.001; ** p<0.01; * p<0.05; . p<0.1",
        caption = "General and final models", caption.above = TRUE)

# Multicollinearity check on the final LPM. Ordinary VIF is inflated by the
# interaction, so the interaction-aware version is the one to read.
cat("\nVIF, final LPM (interaction-aware):\n")
print(vif(LPM_final, type = "predictor"))


# =============================================================================
# SECTION 8: diagnostics for the final model
# =============================================================================

# 8a. Brant test: proportional odds assumption (ordered logit).
cat("\nBrant test (proportional odds):\n")
print(brant(ologit_final))

# 8b. Goodness of fit: ordinal Hosmer-Lemeshow and Lipsitz.
cat("\nOrdinal Hosmer-Lemeshow test:\n")
print(logitgof(df_model$freedom, fitted(ologit_final), g = 10, ord = TRUE))

cat("\nLipsitz test:\n")
print(lipsitz.test(ologit_final))

# The Pulkstenis-Robinson test cannot be applied here: it requires at least one
# categorical predictor to form covariate-pattern strata, and the final model
# uses only continuous regressors. We report this rather than force the test.
cat("\nPulkstenis-Robinson test:\n")
print(tryCatch(pulkrob.chisq(ologit_final, character(0)),
               error = function(e) "Not applicable: the model has no categorical predictor."))

# 8c. Pseudo R2: McFadden, McKelvey-Zavoina, plus count R2 and adjusted count R2.
cat("\nMcFadden pseudo R2 (logit / probit):\n")
print(pR2(ologit_final)["McFadden"])
print(pR2(oprobit_final)["McFadden"])

cat("\nMcKelvey-Zavoina pseudo R2:\n")
print(tryCatch(PseudoR2(ologit_final, which = "McKelveyZavoina"),
               error = function(e) paste("note:", e$message)))
print(tryCatch(PseudoR2(oprobit_final, which = "McKelveyZavoina"),
               error = function(e) paste("note:", e$message)))

count_r2 <- function(model) {
  obs  <- as.character(model.frame(model)[[1]])
  pred <- as.character(predict(model, type = "class"))
  
  n <- length(obs)
  ncorrect <- sum(pred == obs)
  nmode <- max(table(obs))
  
  c(
    CountR2 = round(ncorrect / n, 3),
    AdjCountR2 = round((ncorrect - nmode) / (n - nmode), 3)
  )
}

cat("\nCount R2 and adjusted count R2 (logit):\n")
print(count_r2(ologit_final))

cat("\nCount R2 and adjusted count R2 (probit):\n")
print(count_r2(oprobit_final))

# 8d. Linktest for specification (ordered version). We want the linear term
# significant and the squared term insignificant.
linktest_ordered <- function(model) {
  beta <- coef(model)
  X    <- model.matrix(model)[, names(beta), drop = FALSE]
  yhat <- as.vector(X %*% beta)
  dd   <- data.frame(y = model.frame(model)[[1]], yhat = yhat, yhat2 = yhat^2)
  lt   <- polr(y ~ yhat + yhat2, data = dd, method = model$method, Hess = TRUE)
  polr_table(lt)
}
cat("\nLinktest, final ordered logit:\n")
print(linktest_ordered(ologit_final))
cat("\nLinktest, final ordered probit:\n")
print(linktest_ordered(oprobit_final))

# 8e. Marginal effects per category for the final model. Ordered marginal
# effects differ by outcome and sum to zero across categories for each variable.
# We use the ordered probit, which is the preferred specification for
# interpretation (see the predicted probabilities below).
cat("\nMarginal effects, final ordered probit (ocME):\n")
print(tryCatch(ocME(oprobit_final), error = function(e) paste("note:", e$message)))
cat("\nMarginal effects, final ordered logit (ocME):\n")
print(tryCatch(ocME(ologit_final), error = function(e) paste("note:", e$message)))


# =============================================================================
# SECTION 9: predicted probabilities and hypothesis testing
# =============================================================================
# The Brant test rejects proportional odds for the ordered logit, mostly because
# of oil rents, so predicted probabilities are based on the ordered probit and
# the ordered logit is treated as a robustness model.
pred_model <- oprobit_final

predict_probs <- function(model, newdata) {
  probs <- as.data.frame(predict(model, newdata = newdata, type = "probs"))
  cbind(newdata, probs)
}

base_log_pop     <- mean(df_model$log_pop)
base_trade       <- mean(df_model$trade)
base_log_oilrent <- mean(df_model$log_oilrent)

# H1 and H2: kinship and the kinship-by-GDP interaction.
gdp_levels <- data.frame(
  gdp_label  = c("Low GDP", "Average GDP", "High GDP"),
  log_gdppc_c = c(quantile(df_model$log_gdppc_c, 0.25), 0,
                  quantile(df_model$log_gdppc_c, 0.75))
)

kinship_grid <- expand.grid(kinship = seq(0, 1, length.out = 100),
                            gdp_label = gdp_levels$gdp_label)
kinship_grid$log_gdppc_c <- gdp_levels$log_gdppc_c[match(kinship_grid$gdp_label,
                                                         gdp_levels$gdp_label)]
kinship_grid$log_pop     <- base_log_pop
kinship_grid$trade       <- base_trade
kinship_grid$log_oilrent <- base_log_oilrent

pred_kinship_long <- predict_probs(pred_model, kinship_grid) %>%
  pivot_longer(cols = c("Not Free", "Partly Free", "Free"),
               names_to = "freedom_status", values_to = "probability")

g_pp_kin <- ggplot(pred_kinship_long,
                   aes(x = kinship, y = probability, linetype = gdp_label)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ freedom_status) +
  labs(title = "Predicted freedom probabilities by kinship and GDP level",
       x = "Kinship intensity", y = "Predicted probability", linetype = "GDP level") +
  theme_minimal(base_size = 12)
print(g_pp_kin)
ggsave("fig_pp_kinship_gdp.png", g_pp_kin, width = 9, height = 5, dpi = 300)

kinship_table_data <- expand.grid(kinship = c(0.1, 0.9),
                                  gdp_label = gdp_levels$gdp_label)
kinship_table_data$log_gdppc_c <- gdp_levels$log_gdppc_c[match(kinship_table_data$gdp_label,
                                                               gdp_levels$gdp_label)]
kinship_table_data$log_pop     <- base_log_pop
kinship_table_data$trade       <- base_trade
kinship_table_data$log_oilrent <- base_log_oilrent
kinship_pred_table <- predict_probs(pred_model, kinship_table_data)
cat("\nPredicted probabilities, low vs high kinship by GDP level:\n")
print(kinship_pred_table, row.names = FALSE)
write.csv(kinship_pred_table, "pp_kinship_gdp.csv", row.names = FALSE)

# H3: oil-rent dependence.
oil_grid <- data.frame(
  kinship = mean(df_model$kinship), log_gdppc_c = 0,
  log_pop = base_log_pop, trade = base_trade,
  log_oilrent = seq(quantile(df_model$log_oilrent, 0.05),
                    quantile(df_model$log_oilrent, 0.95), length.out = 100)
)
pred_oil_long <- predict_probs(pred_model, oil_grid) %>%
  pivot_longer(cols = c("Not Free", "Partly Free", "Free"),
               names_to = "freedom_status", values_to = "probability")

g_pp_oil <- ggplot(pred_oil_long, aes(x = log_oilrent, y = probability)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ freedom_status) +
  labs(title = "Predicted freedom probabilities by oil-rent dependence",
       x = "Log(1 + oil rents)", y = "Predicted probability") +
  theme_minimal(base_size = 12)
print(g_pp_oil)
ggsave("fig_pp_oilrent.png", g_pp_oil, width = 9, height = 5, dpi = 300)

oil_table_data <- data.frame(
  kinship = mean(df_model$kinship), log_gdppc_c = 0,
  log_pop = base_log_pop, trade = base_trade,
  log_oilrent = quantile(df_model$log_oilrent, c(0.10, 0.50, 0.90)),
  oil_level = c("Low oil rents", "Median oil rents", "High oil rents")
)
oil_pred_table <- predict_probs(pred_model, oil_table_data)
cat("\nPredicted probabilities, oil-rent scenarios:\n")
print(oil_pred_table, row.names = FALSE)
write.csv(oil_pred_table, "pp_oilrent.csv", row.names = FALSE)
