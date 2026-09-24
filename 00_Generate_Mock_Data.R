# ============================================================================
# 00_Generate_Mock_Data.R
# Purpose: Generate a synthetic cohort for demonstrating the core calculations.
# The data are entirely simulated and do not reproduce the study data.
# ============================================================================

set.seed(20260924)

N <- 10000L

# The example uses pan-cancer plus 12 site-specific endpoints.
OUTCOMES <- c("pan", "162", "150", "151", "153", "154", "155",
              "156", "157", "189", "188", "14", "193")

RISK_FACTORS <- c(
  "eversmoke", "drink_2", "adiposity_2", "TPA_2",
  "proc_meat_2", "redmeat_2", "veg_2", "fruit_2",
  "HT2DM", "HLD", "Hgallstone"
)

sample_binary <- function(n, p) rbinom(n, size = 1L, prob = p)

sex <- factor(sample(c("Female", "Male"), N, replace = TRUE,
                     prob = c(0.54, 0.46)),
              levels = c("Female", "Male"))
country <- factor(sample(c("England", "Scotland", "Wales"), N,
                         replace = TRUE, prob = c(0.84, 0.10, 0.06)))
startage <- pmin(pmax(rnorm(N, mean = 57, sd = 7), 40.1), 70)
base_date <- as.Date("2010-01-01")
DOB <- base_date - round(startage * 365.25)
startage <- as.numeric(base_date - DOB) / 365.25
birth_year <- as.integer(format(DOB, "%Y"))
birth_cohort <- cut(
  birth_year,
  breaks = c(1910, 1930, 1940, 1950, 1960, 1970, 1980, 1990),
  right = FALSE,
  include.lowest = TRUE,
  labels = c("<1930", "1930-1939", "1940-1949", "1950-1959",
             "1960-1969", "1970-1979", "1980+")
)

# Education and income are categorical covariates with an additional "prefer not to answer" category.
edu_2 <- factor(sample(c("Low", "High", "Prefer not to answer"), N,
                        replace = TRUE, prob = c(0.18, 0.80, 0.02)))
inc_2 <- factor(sample(c("Low", "High", "Prefer not to answer"), N,
                        replace = TRUE, prob = c(0.22, 0.74, 0.04)))

# Sex-dependent risk-factor prevalences.
p_male <- sex == "Male"

eversmoke <- sample_binary(N, ifelse(p_male, 0.51, 0.40))
drink_2 <- sample_binary(N, ifelse(p_male, 0.25, 0.16))
adiposity_2 <- sample_binary(N, ifelse(p_male, 0.35, 0.38))
TPA_2 <- sample_binary(N, 0.25)
proc_meat_2 <- sample_binary(N, ifelse(p_male, 0.43, 0.21))
redmeat_2 <- sample_binary(N, ifelse(p_male, 0.29, 0.25))
veg_2 <- sample_binary(N, ifelse(p_male, 0.21, 0.13))
fruit_2 <- sample_binary(N, ifelse(p_male, 0.09, 0.18))
HT2DM <- sample_binary(N, ifelse(p_male, 0.14, 0.10))
HLD <- sample_binary(N, ifelse(p_male, 0.12, 0.10))
Hgallstone <- sample_binary(N, ifelse(p_male, 0.12, 0.15))

# Add a small amount of dependence between selected health conditions while
# retaining simple 0/1 variables for the counterfactual PAF calculations.
HT2DM[adiposity_2 == 1 & runif(N) < 0.08] <- 1L
HLD[drink_2 == 1 & runif(N) < 0.03] <- 1L
Hgallstone[adiposity_2 == 1 & runif(N) < 0.05] <- 1L

# Exposures and age affect incidence across the simulated endpoints.
lp_common <-
  0.018 * (startage - 55) +
  0.10 * (sex == "Male") +
  0.06 * (edu_2 == "Low") +
  0.04 * (inc_2 == "Low")

mock_data <- data.frame(
  id = seq_len(N),
  sex = sex,
  country = country,
  DOB = DOB,
  startage = startage,
  birth_cohort = birth_cohort,
  edu_2 = edu_2,
  inc_2 = inc_2,
  eversmoke = eversmoke,
  drink_2 = drink_2,
  adiposity_2 = adiposity_2,
  TPA_2 = TPA_2,
  proc_meat_2 = proc_meat_2,
  redmeat_2 = redmeat_2,
  veg_2 = veg_2,
  fruit_2 = fruit_2,
  HT2DM = HT2DM,
  HLD = HLD,
  Hgallstone = Hgallstone
)

effect_sizes <- c(
  eversmoke = 0.28, drink_2 = 0.18, adiposity_2 = 0.16,
  TPA_2 = 0.10, proc_meat_2 = 0.10, redmeat_2 = 0.08,
  veg_2 = 0.08, fruit_2 = 0.07, HT2DM = 0.18,
  HLD = 0.28, Hgallstone = 0.14
)

# Each site has a randomly chosen illustrative baseline hazard. 
site_hazards <- setNames(runif(length(OUTCOMES) - 1L, 0.002, 0.003),
                         setdiff(OUTCOMES, "pan"))
risk_lp <- lp_common
for (rf in RISK_FACTORS) {
  risk_lp <- risk_lp + effect_sizes[[rf]] *
    ifelse(is.na(mock_data[[rf]]), 0, mock_data[[rf]])
}
risk_multiplier <- exp(risk_lp - mean(risk_lp))
censor_time <- runif(N, min = 5, max = 15)
first_cancer_time <- rep(Inf, N)

for (icd in setdiff(OUTCOMES, "pan")) {
  event_time <- rexp(N, rate = site_hazards[[icd]] * risk_multiplier)
  first_cancer_time <- pmin(first_cancer_time, event_time)
  followup_time <- pmin(event_time, censor_time)
  mock_data[[paste0("CA_", icd)]] <- as.integer(event_time <= censor_time)
  mock_data[[paste0("CA_", icd, "_py")]] <- followup_time
  mock_data[[paste0("CA_", icd, "_endage")]] <- startage + followup_time
}

# Other qualifying cancers fill the gap between the listed sites and the
# combined non-sex-specific cancer endpoint; they are not separate outputs.
other_time <- rexp(N, rate = 0.002 * risk_multiplier)
first_cancer_time <- pmin(first_cancer_time, other_time)
pan_time <- pmin(first_cancer_time, censor_time)
mock_data$CA_pan <- as.integer(first_cancer_time <= censor_time)
mock_data$CA_pan_py <- pan_time
mock_data$CA_pan_endage <- startage + pan_time

saveRDS(mock_data, "mock_data.rds")

rm()
