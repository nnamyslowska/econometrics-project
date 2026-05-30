# =============================================================================
# Kinship Intensity and Human Development: A Tobit Analysis
# =============================================================================
#
# DEPENDENT VARIABLE: undp_hdi — UNDP Human Development Index [0, 1]
#
# WHY TOBIT (AND WHY HDI)?
#   HDI is a composite index bounded strictly between 0 and 1 by construction.
#   The upper boundary is a genuine constraint: highly developed countries
#   (Norway, Switzerland, Iceland, etc.) cluster near the ceiling of 1,
#   meaning OLS would systematically underestimate effects at the top of the
#   distribution. This is precisely the data structure Tobit is designed for —
#   a latent variable (true human development potential) that is censored above
#   at 1. We apply upper censoring at right = 1 and no left censoring (the
#   observed minimum of 0.385 is well away from the theoretical floor of 0).
#
#   HDI is also the conceptually cleanest link to kinship: it captures health,
#   education, and living standards simultaneously — all three dimensions that
#   tight kinship networks are theorised to affect (Schulz et al. 2019, Science).
#   Societies with high kinship intensity tend to concentrate resources within
#   the clan, limiting investment in broad public goods and universal education.
#
# KEY INDEPENDENT VARIABLE: kinship_score (ekne) [0, 1]
#   0 = nuclear-family society; 1 = tight kinship / cousin-marriage networks
#
# CONTROLS (8 variables + 1 interaction):
#   log_gdppc     — economic development (log World Bank GDP p.c., const. 2015 USD)
#   log_pop       — country size (log total population)
#   urban_pct     — urbanisation (% of population in urban areas)
#   unemp         — unemployment rate (ILO, %)
#   dem_duration  — years of continuous democracy (BMR)
#   cpi           — Corruption Perceptions Index [0–100] (Transparency Intl.)
#   govt_eff      — government effectiveness (World Bank WGI)
#   internet      — internet users as % of population
#   INTERACTION   — kinship × log_gdppc
#     Rationale: kinship networks may substitute for weak formal institutions
#     in poor countries, making the kinship effect contingent on wealth level.
#
# DATA SOURCES:
#   QoG Standard Dataset, Jan 2026 (Teorell et al.)
#   kinship_df — ekne kinship intensity index
# =============================================================================


# --- 0. Packages --------------------------------------------------------------

# install.packages(c("censReg", "AER", "DescTools",
#                    "stargazer", "ggplot2", "dplyr", "lmtest", "nortest"))

library("censReg")    # censReg() — Tobit estimator (lab-standard)
library("AER")        # tobit()   — alternative; gives $x and $scale for ME
library("DescTools")  # Desc()    — descriptive statistics (lab-standard)
library("stargazer")  # publication-quality regression tables
library("ggplot2")    # visualisation
library("dplyr")      # data wrangling
library("lmtest")     # resettest, bptest, bgtest
library("nortest")    # ad.test() — Anderson-Darling normality test
library("car")        # linearHypothesis()

Sys.setenv(LANG = "en")
options(scipen = 100)


# =============================================================================
# SECTION 1: DATA LOADING & MERGING
# =============================================================================

qog     <- read.csv("data/qog_std_cs_jan26.csv",  stringsAsFactors = FALSE)
kinship <- read.csv("data/kinship_df.csv",         stringsAsFactors = FALSE)

cat("QoG dimensions:     ", dim(qog),     "\n")
cat("Kinship dimensions: ", dim(kinship), "\n")
cat("Kinship columns:    ", names(kinship), "\n")

# Inner join on ISO-3 country codes
# QoG uses 'ccodealp'; kinship uses 'isocode'
df <- merge(qog, kinship,
            by.x  = "ccodealp",
            by.y  = "isocode",
            all   = FALSE)

cat("\nCountries after inner join:", nrow(df), "\n")


# =============================================================================
# SECTION 2: MISSING VALUE DIAGNOSTICS
# =============================================================================

# --- 2a. Full dataset: all columns with at least one NA ---
missing_all <- data.frame(
  variable    = names(df),
  n_missing   = sapply(df, function(x) sum(is.na(x))),
  pct_missing = round(sapply(df, function(x) mean(is.na(x))) * 100, 1),
  row.names   = NULL
)

missing_nonzero <- missing_all[missing_all$n_missing > 0, ]
missing_nonzero <- missing_nonzero[order(-missing_nonzero$n_missing), ]

cat("\n--- Missing values: all columns with NA > 0 (sorted) ---\n")
print(missing_nonzero, row.names = FALSE)

# --- 2b. Key model variables only ---
key_vars <- c("undp_hdi", "kinship_score",
              "wdi_gdpcapcon2015", "wdi_pop", "wdi_popurb",
              "wdi_unempilo", "bmr_demdur", "ti_cpi",
              "wdi_trade", "wbgi_gee", "dr_ig")

cat("\n--- Missing values: key model variables ---\n")
print(missing_all[missing_all$variable %in% key_vars, ], row.names = FALSE)


# =============================================================================
# SECTION 3: VARIABLE CONSTRUCTION
# =============================================================================

df_model <- df %>%
  transmute(
    country      = cname,
    iso          = ccodealp,
    
    # Dependent variable
    hdi          = undp_hdi,                    # HDI [0,1], upper-censored at 1
    
    # Key independent variable
    kinship      = kinship_score,               # ekne [0,1]
    
    # Controls
    log_gdppc    = log(wdi_gdpcapcon2015),      # log GDP per capita
    log_pop      = log(wdi_pop),                # log population
    urban_pct    = wdi_popurb,                  # urbanisation (%)
    unemp        = wdi_unempilo,                # unemployment rate (%)
    dem_duration = bmr_demdur,                  # years of democracy
    cpi          = ti_cpi,                      # corruption perception [0–100]
    govt_eff     = wbgi_gee,                    # government effectiveness
    internet     = dr_ig                        # internet users (%)
  ) %>%
  filter(complete.cases(.))

cat("\nWorking sample (complete cases):", nrow(df_model), "countries\n")


# =============================================================================
# SECTION 4: DESCRIPTIVE STATISTICS
# =============================================================================

# Desc() mirrors the lab-script approach
Desc(df_model$hdi,     main = "HDI — Dependent Variable")
Desc(df_model$kinship, main = "Kinship intensity score")

# Summary table for all model variables
cat("\n--- Summary statistics (working sample, N =", nrow(df_model), ") ---\n")
summary(df_model[ , c("hdi", "kinship", "log_gdppc", "log_pop",
                      "urban_pct", "unemp", "dem_duration",
                      "cpi", "govt_eff", "internet")])


# =============================================================================
# SECTION 5: VISUAL INSPECTION
# =============================================================================

# --- Plot 1: HDI distribution — check ceiling clustering ---
ggplot(df_model, aes(x = hdi)) +
  geom_histogram(bins = 30, fill = "#3A7DC9", colour = "white", alpha = 0.85) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#C94040", linewidth = 0.8) +
  annotate("text", x = 0.985, y = Inf, label = "Ceiling = 1",
           colour = "#C94040", angle = 90, vjust = 1.4, hjust = 1.2, size = 3.5) +
  labs(title    = "Distribution of HDI (working sample)",
       subtitle = paste0("N = ", nrow(df_model),
                         "  |  ", sum(df_model$hdi >= 0.95),
                         " observations above 0.95 — upper censoring justification"),
       x = "Human Development Index", y = "Count") +
  theme_minimal(base_size = 13)

# --- Plot 2: Kinship distribution ---
ggplot(df_model, aes(x = kinship)) +
  geom_histogram(bins = 25, fill = "#5AA65A", colour = "white", alpha = 0.85) +
  labs(title    = "Distribution of kinship intensity score",
       subtitle = "0 = nuclear-family society   |   1 = tight kinship / cousin-marriage",
       x = "Kinship score", y = "Count") +
  theme_minimal(base_size = 13)

# --- Plot 3: Kinship vs HDI — core bivariate relationship ---
ggplot(df_model, aes(x = kinship, y = hdi)) +
  geom_point(alpha = 0.55, colour = "#3A7DC9", size = 2) +
  geom_smooth(method = "lm", se = TRUE, colour = "#C94040", linewidth = 0.9) +
  labs(title    = "Kinship intensity vs Human Development Index",
       subtitle = "OLS fit shown; expected negative slope",
       x = "Kinship score (ekne)",
       y = "HDI") +
  theme_minimal(base_size = 13)

# --- Plot 4: Kinship vs HDI, coloured by log GDP per capita ---
# Motivates the interaction term: does the kinship–HDI relationship
# differ across income levels?
ggplot(df_model, aes(x = kinship, y = hdi, colour = log_gdppc)) +
  geom_point(alpha = 0.7, size = 2.2) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.7) +
  scale_colour_gradient(low = "#FDE68A", high = "#1E3A5F",
                        name = "Log GDP p.c.") +
  labs(title    = "Kinship vs HDI — coloured by income level",
       subtitle = "Motivates the kinship × log GDP interaction term",
       x = "Kinship score (ekne)", y = "HDI") +
  theme_minimal(base_size = 13)

# --- Plot 5: Kinship vs HDI split by income tercile ---
df_model$income_group <- cut(df_model$log_gdppc,
                             breaks = quantile(df_model$log_gdppc,
                                               probs = c(0, 1/3, 2/3, 1)),
                             labels = c("Low income", "Middle income", "High income"),
                             include.lowest = TRUE)

ggplot(df_model, aes(x = kinship, y = hdi)) +
  geom_point(alpha = 0.5, colour = "#3A7DC9", size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "#C94040", linewidth = 0.8) +
  facet_wrap(~income_group) +
  labs(title    = "Kinship vs HDI by income tercile",
       subtitle = "Interaction effect visualised across subgroups",
       x = "Kinship score (ekne)", y = "HDI") +
  theme_minimal(base_size = 12)

# --- Plot 6: Correlation matrix heatmap (key numeric variables) ---
cor_vars <- c("hdi","kinship","log_gdppc","log_pop","urban_pct",
              "unemp","dem_duration","cpi","govt_eff","internet")
cor_matrix <- round(cor(df_model[, cor_vars], use = "complete.obs"), 2)

cor_df <- as.data.frame(as.table(cor_matrix))
names(cor_df) <- c("Var1","Var2","Corr")

ggplot(cor_df, aes(x = Var1, y = Var2, fill = Corr)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = Corr), size = 3, colour = "black") +
  scale_fill_gradient2(low = "#C94040", mid = "white", high = "#3A7DC9",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  labs(title = "Correlation matrix — model variables",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# =============================================================================
# SECTION 6: OLS BASELINE (as in lab script)
# =============================================================================
# Estimate OLS first — serves as a benchmark and as the basis for diagnostic
# tests later (Ramsey-RESET, Breusch-Pagan, Breusch-Godfrey, linktest)

OLS_general <- lm(
  hdi ~ kinship + log_gdppc + log_pop + urban_pct +
    unemp + dem_duration + cpi + govt_eff + internet +
    kinship:log_gdppc,
  data = df_model
)
summary(OLS_general)


# =============================================================================
# SECTION 7: GENERAL TOBIT MODEL
# =============================================================================
# Upper censoring at 1 (HDI ceiling); no left censoring (min observed = 0.385)
# Using censReg() — the estimator used throughout the lab script

tobit_general <- censReg(
  hdi ~ kinship + log_gdppc + log_pop + urban_pct +
    unemp + dem_duration + cpi + govt_eff + internet +
    kinship:log_gdppc,
  left  = 0,   # no left censoring
  right = 1,      # upper censoring at HDI = 1
  data  = df_model
)

summary(tobit_general)

# --- Quick marginal effects preview (conditional, i.e. E[y | y < 1]) ---
# Using margEff() from censReg, as in the lab
cat("\n--- Conditional marginal effects E(y | y < 1) — general model ---\n")
summary(margEff(tobit_general))

# =============================================================================
# NOTE: General-to-specific variable selection, diagnostic tests,
#       three types of marginal effects, hypothesis testing, and the
#       publication-quality table PUT HERE.
# =============================================================================