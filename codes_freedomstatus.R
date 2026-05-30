# =============================================================================
# Kinship Intensity, Governance, and Democracy: Two Tobit Models
# =============================================================================
#
# MODEL 1 — Dependent variable: wbgi_gee (Government Effectiveness, World Bank WGI)
#   Theoretical range [-2.5, 2.5]; observed [-2.01, 2.20]
#   Tobit rationale: GEE is derived from a latent governance measurement model
#   whose underlying scale is censored at both tails (±2.5) by the estimation
#   method. 18 observations cluster near the upper bound (>1.5) and 10 near the
#   lower bound (<-1.5). Two-sided censoring: left = -2.5, right = 2.5.
#   Interaction: kinship × log_gdppc
#   Rationale: the dampening effect of kinship on state capacity may be weaker
#   in wealthier countries that have resources to build impersonal institutions
#   regardless of cultural legacy.
#
# MODEL 2 — Dependent variable: vdem_libdem (Liberal Democracy Index, V-Dem)
#   Theoretical range [0, 1]; observed [0.012, 0.890]
#   Tobit rationale: 14 observations pile up near the floor (< 0.05) — these are
#   hard autocracies where liberal democracy is effectively zero and cannot fall
#   further. This is a left-censored Tobit (left = 0, right = Inf).
#   Interaction: kinship × dem_duration
#   Rationale: in long-established democracies, formal institutions may have
#   displaced kin-based political organisation; the kinship–democracy link
#   should therefore weaken with accumulated democratic experience.
#
# WHY NOT TAUTOLOGICAL?
#   Kinship score (ekne, Enke 2019 / Schulz et al. 2019) is built from
#   pre-industrial anthropological data: cousin-marriage rates, clan structures,
#   co-residence patterns, and lineage organisation — all sourced from the
#   Ethnographic Atlas (Murdock 1967) and historical church records.
#   GEE and libdem are contemporary institutional outcomes measured by entirely
#   independent surveys and expert codings. There is no shared input.
#
# CONTROLS (8 variables, shared across both models):
#   log_gdppc    — economic development (log GDP p.c., const. 2015 USD, WB)
#   log_pop      — country size (log total population, WB)
#   urban_pct    — urbanisation % (WB): urban settings erode kin-network reliance
#   unemp        — unemployment rate % (ILO): labour market stress
#   dem_duration — years of continuous democracy (BMR): institutional memory
#   cpi          — Corruption Perceptions Index 0–100 (TI): related but distinct
#   trade        — trade openness % of GDP (WB): external institutional pressure
#   internet     — internet users % of population: civic information access
#
# DATA:
#   QoG Standard Dataset Jan 2026 (Teorell et al.) — qog_std_cs_jan26.csv
#   Kinship intensity index — kinship_df.csv
#   Working sample: N = 144 (complete cases across both models)
# =============================================================================


# --- 0. Packages --------------------------------------------------------------

# install.packages(c("censReg", "AER", "DescTools",
#                    "stargazer", "ggplot2", "dplyr", "lmtest", "nortest", "car"))

library("censReg")    # censReg() — Tobit estimator (lab-standard)
library("AER")        # tobit()   — alternative; gives $x, $scale for manual ME
library("DescTools")  # Desc()    — descriptive statistics (lab-standard)
library("stargazer")  # publication-quality regression tables
library("ggplot2")    # visualisation
library("dplyr")      # data wrangling
library("lmtest")     # resettest, bptest, bgtest
library("nortest")    # ad.test() — Anderson-Darling normality test
library("car")        # linearHypothesis(), vif()

Sys.setenv(LANG = "en")
options(scipen = 100)


# =============================================================================
# SECTION 1: LOAD & MERGE
# =============================================================================

qog     <- read.csv("data/qog_std_cs_jan26.csv",  stringsAsFactors = FALSE)
kinship <- read.csv("data/kinship_df.csv",         stringsAsFactors = FALSE)

cat("QoG dimensions:     ", dim(qog),     "\n")
cat("Kinship dimensions: ", dim(kinship), "\n")

# Inner join on ISO-3 country codes (ccodealp in QoG; isocode in kinship)
df <- merge(qog, kinship,
            by.x = "ccodealp",
            by.y = "isocode",
            all  = FALSE)

cat("Countries after inner join:", nrow(df), "\n")


# =============================================================================
# SECTION 2: MISSING VALUE DIAGNOSTICS
# =============================================================================

# --- 2a. Full merged dataset: every column with at least one NA ---
missing_all <- data.frame(
  variable    = names(df),
  n_missing   = sapply(df, function(x) sum(is.na(x))),
  pct_missing = round(sapply(df, function(x) mean(is.na(x))) * 100, 1),
  row.names   = NULL
)

missing_nonzero <- missing_all[missing_all$n_missing > 0, ]
missing_nonzero <- missing_nonzero[order(-missing_nonzero$n_missing), ]

cat("\n--- Missing values: all columns with NA > 0 (sorted descending) ---\n")
print(missing_nonzero, row.names = FALSE)

# --- 2b. Key model variables only ---
key_vars <- c("wbgi_gee", "vdem_libdem", "kinship_score",
              "wdi_gdpcapcon2015", "wdi_pop", "wdi_popurb",
              "wdi_unempilo", "bmr_demdur", "ti_cpi",
              "wdi_trade", "dr_ig")

cat("\n--- Missing values: key model variables ---\n")
print(missing_all[missing_all$variable %in% key_vars, ], row.names = FALSE)


# =============================================================================
# SECTION 3: BUILD WORKING DATASET
# =============================================================================
# Use the intersection of complete cases across BOTH models so all comparisons
# are made on the same set of countries.

df_model <- df %>%
  transmute(
    country      = cname,
    iso          = ccodealp,
    
    # Dependent variables
    govt_eff     = wbgi_gee,                  # Gov. Effectiveness [-2.5, 2.5]
    libdem       = vdem_libdem,               # Liberal Democracy  [0, 1]
    
    # Key independent variable
    kinship      = kinship_score,             # ekne [0,1]
    
    # Controls
    log_gdppc    = log(wdi_gdpcapcon2015),    # log GDP per capita
    log_pop      = log(wdi_pop),              # log population
    urban_pct    = wdi_popurb,                # urbanisation (%)
    unemp        = wdi_unempilo,              # unemployment rate (%)
    dem_duration = bmr_demdur,               # years of democracy
    cpi          = ti_cpi,                   # corruption index [0–100]
    trade        = wdi_trade,                # trade openness (% GDP)
    internet     = dr_ig                     # internet users (%)
  ) %>%
  filter(complete.cases(.))

cat("\nWorking sample (complete cases, both models):", nrow(df_model), "countries\n")


# =============================================================================
# SECTION 4: DESCRIPTIVE STATISTICS
# =============================================================================

# Desc() — lab-standard deep descriptive for the two DVs
Desc(df_model$govt_eff, main = "Government Effectiveness (wbgi_gee)")
Desc(df_model$libdem,   main = "Liberal Democracy Index (vdem_libdem)")
Desc(df_model$kinship,  main = "Kinship intensity score (ekne)")

# Compact summary for all model variables
cat("\n--- Summary statistics (N =", nrow(df_model), ") ---\n")
summary(df_model[ , c("govt_eff", "libdem", "kinship",
                      "log_gdppc", "log_pop", "urban_pct", "unemp",
                      "dem_duration", "cpi", "trade", "internet")])


# =============================================================================
# SECTION 5: VISUAL INSPECTION
# =============================================================================

# --- Plot 1: Distribution of Government Effectiveness ---
# Look for clustering near ±2.5 (Tobit censoring justification)
ggplot(df_model, aes(x = govt_eff)) +
  geom_histogram(bins = 30, fill = "#3A7DC9", colour = "white", alpha = 0.85) +
  geom_vline(xintercept =  2.5, linetype = "dashed",
             colour = "#C94040", linewidth = 0.8) +
  geom_vline(xintercept = -2.5, linetype = "dashed",
             colour = "#C94040", linewidth = 0.8) +
  annotate("text", x =  2.35, y = Inf, label = "Upper bound +2.5",
           colour = "#C94040", angle = 90, vjust = 1.3, hjust = 1.2, size = 3.2) +
  annotate("text", x = -2.35, y = Inf, label = "Lower bound −2.5",
           colour = "#C94040", angle = 90, vjust = 1.3, hjust = 1.2, size = 3.2) +
  labs(title    = "Distribution of Government Effectiveness (Model 1 DV)",
       subtitle = paste0("N = ", nrow(df_model),
                         "  |  18 obs > 1.5  |  10 obs < −1.5",
                         "  →  two-sided Tobit [-2.5, 2.5]"),
       x = "Government Effectiveness (WGI)", y = "Count") +
  theme_minimal(base_size = 13)

# --- Plot 2: Distribution of Liberal Democracy ---
# Look for floor clustering near 0 (left-censored Tobit justification)
ggplot(df_model, aes(x = libdem)) +
  geom_histogram(bins = 30, fill = "#5AA65A", colour = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "#C94040", linewidth = 0.8) +
  annotate("text", x = 0.015, y = Inf, label = "Floor = 0",
           colour = "#C94040", angle = 90, vjust = 1.3, hjust = 1.2, size = 3.2) +
  labs(title    = "Distribution of Liberal Democracy Index (Model 2 DV)",
       subtitle = paste0("N = ", nrow(df_model),
                         "  |  14 obs < 0.05  |  max = 0.890",
                         "  →  left-censored Tobit [0, +Inf]"),
       x = "Liberal Democracy Index (V-Dem)", y = "Count") +
  theme_minimal(base_size = 13)

# --- Plot 3: Kinship vs Government Effectiveness ---
ggplot(df_model, aes(x = kinship, y = govt_eff)) +
  geom_point(alpha = 0.55, colour = "#3A7DC9", size = 2) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "#C94040", linewidth = 0.9) +
  labs(title    = "Kinship intensity vs Government Effectiveness",
       subtitle = "Expected negative relationship (H1)",
       x = "Kinship score (ekne) — 0: nuclear family, 1: tight kinship",
       y = "Government Effectiveness (WGI)") +
  theme_minimal(base_size = 13)

# --- Plot 4: Kinship vs Liberal Democracy ---
ggplot(df_model, aes(x = kinship, y = libdem)) +
  geom_point(alpha = 0.55, colour = "#5AA65A", size = 2) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "#C94040", linewidth = 0.9) +
  labs(title    = "Kinship intensity vs Liberal Democracy",
       subtitle = "Expected negative relationship (H2)",
       x = "Kinship score (ekne) — 0: nuclear family, 1: tight kinship",
       y = "Liberal Democracy Index (V-Dem)") +
  theme_minimal(base_size = 13)

# --- Plot 5: Interaction motivation — kinship × log_gdppc on govt_eff ---
# Split by income tercile to visualise whether the slope varies by wealth
df_model$income_group <- cut(
  df_model$log_gdppc,
  breaks         = quantile(df_model$log_gdppc, probs = c(0, 1/3, 2/3, 1)),
  labels         = c("Low income", "Middle income", "High income"),
  include.lowest = TRUE
)

ggplot(df_model, aes(x = kinship, y = govt_eff)) +
  geom_point(alpha = 0.5, colour = "#3A7DC9", size = 1.8) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "#C94040", linewidth = 0.8) +
  facet_wrap(~income_group) +
  labs(title    = "Kinship vs Government Effectiveness — by income tercile",
       subtitle = "Motivates kinship × log GDP p.c. interaction in Model 1",
       x = "Kinship score (ekne)", y = "Government Effectiveness") +
  theme_minimal(base_size = 12)

# --- Plot 6: Interaction motivation — kinship × dem_duration on libdem ---
# Split by democracy duration tercile
df_model$dem_group <- cut(
  df_model$dem_duration,
  breaks         = quantile(df_model$dem_duration, probs = c(0, 1/3, 2/3, 1)),
  labels         = c("Short democracy", "Medium democracy", "Long democracy"),
  include.lowest = TRUE
)

ggplot(df_model, aes(x = kinship, y = libdem)) +
  geom_point(alpha = 0.5, colour = "#5AA65A", size = 1.8) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "#C94040", linewidth = 0.8) +
  facet_wrap(~dem_group) +
  labs(title    = "Kinship vs Liberal Democracy — by democratic experience",
       subtitle = "Motivates kinship × democracy duration interaction in Model 2",
       x = "Kinship score (ekne)", y = "Liberal Democracy Index (V-Dem)") +
  theme_minimal(base_size = 12)

# --- Plot 7: Correlation heatmap — all model variables ---
cor_vars <- c("govt_eff", "libdem", "kinship", "log_gdppc", "log_pop",
              "urban_pct", "unemp", "dem_duration", "cpi", "trade", "internet")

cor_matrix <- round(cor(df_model[, cor_vars], use = "complete.obs"), 2)
cor_df     <- as.data.frame(as.table(cor_matrix))
names(cor_df) <- c("Var1", "Var2", "Corr")

ggplot(cor_df, aes(x = Var1, y = Var2, fill = Corr)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = Corr), size = 2.8, colour = "black") +
  scale_fill_gradient2(low  = "#C94040", mid = "white", high = "#3A7DC9",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  labs(title = "Correlation matrix — all model variables",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# =============================================================================
# SECTION 6: OLS BASELINES
# =============================================================================
# Estimate OLS for both models first, mirroring the lab script approach.
# These serve as: (a) benchmark for Tobit comparison, (b) basis for
# diagnostic tests (RESET, BP, BG, linktest) later in the analysis.

OLS_gee <- lm(
  govt_eff ~ kinship + log_gdppc + log_pop + urban_pct +
    unemp + dem_duration + cpi + trade + internet +
    kinship:log_gdppc,
  data = df_model
)
cat("\n--- OLS baseline: Government Effectiveness ---\n")
summary(OLS_gee)

OLS_libdem <- lm(
  libdem ~ kinship + log_gdppc + log_pop + urban_pct +
    unemp + dem_duration + cpi + trade + internet +
    kinship:dem_duration,
  data = df_model
)
cat("\n--- OLS baseline: Liberal Democracy ---\n")
summary(OLS_libdem)


# =============================================================================
# SECTION 7: GENERAL TOBIT MODELS
# =============================================================================

# --- Model 1: Government Effectiveness — two-sided Tobit [-2.5, +2.5] --------
# Rationale: WGI scores are derived from a latent measurement model whose scale
# is anchored at ±2.5; observations at the extremes are measurement-censored.

tobit_gee <- censReg(
  govt_eff ~ kinship + log_gdppc + log_pop + urban_pct +
    unemp + dem_duration + cpi + trade + internet +
    kinship:log_gdppc,
  left  = -2.5,
  right =  2.5,
  data  = df_model
)

cat("\n=== GENERAL MODEL 1: Government Effectiveness (Tobit, ±2.5) ===\n")
summary(tobit_gee)

cat("\n--- Conditional marginal effects E(y | uncensored) — Model 1 ---\n")
summary(margEff(tobit_gee))


# --- Model 2: Liberal Democracy — left-censored Tobit [0, +Inf] --------------
# Rationale: 14 countries pile up at the floor of the [0,1] scale (hard
# autocracies where measured liberal democracy is effectively zero). No
# observations cluster near the ceiling (max = 0.890), so no right censoring.

tobit_libdem <- censReg(
  libdem ~ kinship + log_gdppc + log_pop + urban_pct +
    unemp + dem_duration + cpi + trade + internet +
    kinship:dem_duration,
  left  = 0,
  right = Inf,
  data  = df_model
)

cat("\n=== GENERAL MODEL 2: Liberal Democracy (Tobit, left-censored at 0) ===\n")
summary(tobit_libdem)

cat("\n--- Conditional marginal effects E(y | uncensored) — Model 2 ---\n")
summary(margEff(tobit_libdem))


# =============================================================================
# NOTE: General-to-specific selection, full marginal effects (3 kinds),
#       diagnostic tests, hypothesis verification, and publication table
#       continue in the next section of the analysis.
# =============================================================================