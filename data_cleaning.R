# ============================================================
# PROJECT: Kinship & Unemployment Benefits (Tobit Model)
# SCRIPT:  01_data_import_explore.R
# ============================================================

# ---- 0. PACKAGES -------------------------------------------
# Install any missing packages first (run once)
install.packages(c("tidyverse", "haven", "readxl", "skimr", "remotes"))
remotes::install_github("xmarquez/democracyData")

library(tidyverse)   # data wrangling + ggplot2
library(haven)       # read .dta (Stata) files
library(skimr)       # rich summary statistics
library(democracyData) # Polity V download


# ============================================================
# 1. KINSHIP INDEX
#    Source: partner's .dta file
# ============================================================

kinship_raw <- read_dta("data/kinship_data/Data/CountryData.dta")

# Explore
glimpse(kinship_raw)
skim(kinship_raw)
head(kinship_raw)

# Keep only what we need
kinship_df <- kinship_raw %>%
  select(isocode, country, kinship_score)

cat("\n--- Kinship: countries available ---\n")
cat("N countries:", nrow(kinship_df), "\n")
summary(kinship_df$kinship_score)


# ============================================================
# 2. ILO — UNEMPLOYMENT BENEFITS EXPENDITURE (% GDP)
#    File: SDG_0131_SEX_SOC_RT_A-filtered-2026-05-25.csv
#    Note: SDG 1.3.1 — social protection coverage rates
# ============================================================

ilo_raw <- read_csv(
  "data/SDG_0131_SEX_SOC_RT_A-filtered-2026-05-25.csv",
  show_col_types = FALSE
)

# Explore structure first — ILO files vary in layout
glimpse(ilo_raw)
head(ilo_raw, 20)
names(ilo_raw)

# ILO CSVs typically have columns like:
#   ref_area, ref_area.label, sex, sex.label, classif1, classif1.label,
#   time, obs_value, obs_status, ...
# We want: unemployment benefits, both sexes, all years
# Check what classif1 categories exist:
ilo_raw %>%
  distinct(classif1, classif1.label) %>%
  print(n = 50)

# Check sex categories:
ilo_raw %>% distinct(sex, sex.label)

# Filter to: Total sex + Unemployment function + keep obs_value
# (adjust column names below if they differ in your file)
ilo_benefits <- ilo_raw %>%
  filter(
    sex.label  == "Sex: Total",
    str_detect(classif1.label, regex("unemploy", ignore_case = TRUE))
  ) %>%
  select(
    iso3c       = ref_area,
    country_ilo = ref_area.label,
    year        = time,
    unemp_benefits_pct = obs_value
  )

cat("\n--- ILO: years available ---\n")
print(table(ilo_benefits$year))

cat("\n--- ILO: share of zeros ---\n")
mean(ilo_benefits$unemp_benefits_pct == 0, na.rm = TRUE)


# ============================================================
# 3. WORLD BANK — GDP PER CAPITA PPP (constant 2021 USD)
#    File: API_NY.GDP.PCAP.PP.KD_DS2_en_csv_v2_1700.csv
#    IMPORTANT: World Bank CSVs have 4 header rows → skip = 4
# ============================================================

gdp_raw <- read_csv(
  "data/API_NY.GDP.PCAP.PP.KD_DS2_en_csv_v2_1700.csv",
  skip = 4,
  show_col_types = FALSE
)

glimpse(gdp_raw)
head(gdp_raw)

# World Bank format: wide — one column per year
# Columns: "Country Name", "Country Code", "Indicator Name",
#          "Indicator Code", "1960", "1961", ..., "2024"

# Pivot to long format
gdp_long <- gdp_raw %>%
  select(-`Indicator Name`, -`Indicator Code`) %>%
  rename(country_wb = `Country Name`, iso3c = `Country Code`) %>%
  pivot_longer(
    cols      = matches("^[0-9]{4}$"),
    names_to  = "year",
    values_to = "gdp_pc_ppp"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(!is.na(gdp_pc_ppp))

cat("\n--- GDP: years available ---\n")
range(gdp_long$year)
cat("N country-years:", nrow(gdp_long), "\n")


# ============================================================
# 4. WORLD BANK — UNEMPLOYMENT RATE (% of labour force)
#    File: API_SL.UEM.TOTL.ZS_DS2_en_csv_v2_115692.csv
# ============================================================

unemp_raw <- read_csv(
  "data/API_SL.UEM.TOTL.ZS_DS2_en_csv_v2_115692.csv",
  skip = 4,
  show_col_types = FALSE
)

unemp_long <- unemp_raw %>%
  select(-`Indicator Name`, -`Indicator Code`) %>%
  rename(country_wb = `Country Name`, iso3c = `Country Code`) %>%
  pivot_longer(
    cols      = matches("^[0-9]{4}$"),
    names_to  = "year",
    values_to = "unemp_rate"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(!is.na(unemp_rate))

glimpse(unemp_long)
cat("\n--- Unemployment rate: years ---\n")
range(unemp_long$year)


# ============================================================
# 5. WORLD BANK — POPULATION 65+ (% of total)
#    File: API_SP.POP.65UP.TO.ZS_DS2_en_csv_v2_115532.csv
# ============================================================

pop65_raw <- read_csv(
  "data/API_SP.POP.65UP.TO.ZS_DS2_en_csv_v2_115532.csv",
  skip = 4,
  show_col_types = FALSE
)

pop65_long <- pop65_raw %>%
  select(-`Indicator Name`, -`Indicator Code`) %>%
  rename(country_wb = `Country Name`, iso3c = `Country Code`) %>%
  pivot_longer(
    cols      = matches("^[0-9]{4}$"),
    names_to  = "year",
    values_to = "pop_over_65"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(!is.na(pop_over_65))

glimpse(pop65_long)
cat("\n--- Pop 65+: years ---\n")
range(pop65_long$year)


# ============================================================
# 6. POLITY V — DEMOCRACY INDEX
#    Downloaded automatically via democracyData package
# ============================================================

cat("\nDownloading Polity 5 data...\n")
polity_raw <- download_polity_annual()

glimpse(polity_raw)

# polity2 is the cleaned version (handles -66/-77/-88 transitions)
# Range: -10 (full autocracy) to +10 (full democracy)
# We use polity2 and will filter to our target year

polity_df <- polity_raw %>%
  select(
    iso3c   = scode,           # 3-letter country code
    country_polity = polityIV_country,
    year,
    polity2                    # cleaned annual polity score
  ) %>%
  filter(!is.na(polity2))

cat("\n--- Polity: years available ---\n")
range(polity_df$year, na.rm = TRUE)
cat("Polity2 range:", range(polity_df$polity2, na.rm = TRUE), "\n")


# ============================================================
# 7. QUICK OVERLAP CHECK — which year to use?
#    Goal: maximise countries with data in ALL datasets
# ============================================================

# The kinship index is cross-sectional (one value per country)
# All other datasets are annual — we need to pick ONE year
# Candidate years: 2015–2019 (pre-COVID, good ILO coverage)

cat("\n=== COVERAGE CHECK BY YEAR ===\n")

for (yr in 2015:2019) {
  ilo_yr    <- ilo_benefits %>% filter(year == yr) %>% pull(iso3c)
  gdp_yr    <- gdp_long     %>% filter(year == yr) %>% pull(iso3c)
  unemp_yr  <- unemp_long   %>% filter(year == yr) %>% pull(iso3c)
  pop65_yr  <- pop65_long   %>% filter(year == yr) %>% pull(iso3c)
  pol_yr    <- polity_df    %>% filter(year == yr) %>% pull(iso3c)
  kinship_c <- kinship_df$isocode
  
  overlap <- Reduce(intersect,
                    list(ilo_yr, gdp_yr, unemp_yr, pop65_yr, pol_yr, kinship_c))
  
  cat("Year", yr, "→", length(overlap), "countries in all 6 datasets\n")
}

# ============================================================
# NEXT STEP: run 02_data_merge.R once you pick the year
# ============================================================