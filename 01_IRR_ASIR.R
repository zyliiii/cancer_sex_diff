# ============================================================================
# 01_IRR_ASIR.R
# Purpose: Sex-specific ASIRs, male-to-female Poisson IRRs, and ASIR gaps.
# ============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(splines)
  library(epitools)
})

DATA <- readRDS("mock_data.rds")

OUTCOMES <- c("pan", "162", "150", "151", "153", "154", "155",
              "156", "157", "189", "188", "14", "193")

AGE_BREAKS <- seq(40, 85, 5)
AGE_LABELS <- paste(seq(40, 80, 5), seq(44, 84, 5), sep = "-")
SEGI_DOLL_WEIGHTS <- c(6, 6, 5, 4, 4, 3, 2, 1, 0.5)

split_by_age <- function(data, event_var, endage_var) {
  d <- data[complete.cases(data[, c("startage", event_var, endage_var)]), ,
            drop = FALSE]
  tmp <- data.frame(
    startage = d$startage,
    endage = d[[endage_var]],
    status = d[[event_var]]
  )
  sp <- survival::survSplit(
    Surv(startage, endage, status) ~ ., data = tmp,
    cut = AGE_BREAKS
  )
  sp <- sp[sp$startage >= 40 & sp$startage < 85, , drop = FALSE]
  sp$age_group <- cut(
    sp$startage, breaks = AGE_BREAKS, labels = AGE_LABELS,
    right = FALSE
  )
  sp$person_years <- sp$endage - sp$startage
  sp
}

age_specific_rates <- function(data, event_var, endage_var) {
  sp <- split_by_age(data, event_var, endage_var)
  agg <- aggregate(
    cbind(events = status, person_years = person_years) ~ age_group,
    data = sp, FUN = sum
  )
  all_ages <- data.frame(
    age_group = factor(AGE_LABELS, levels = AGE_LABELS),
    stringsAsFactors = FALSE
  )
  all_ages <- merge(all_ages, agg, by = "age_group", all.x = TRUE,
                    sort = FALSE)
  all_ages$events[is.na(all_ages$events)] <- 0
  all_ages$person_years[is.na(all_ages$person_years)] <- 0
  all_ages
}

compute_asir <- function(data, event_var, endage_var) {
  age_tab <- age_specific_rates(data, event_var, endage_var)
  gamma_result <- epitools::ageadjust.direct(
    count = age_tab$events,
    pop = age_tab$person_years,
    stdpop = SEGI_DOLL_WEIGHTS
  )

  w <- SEGI_DOLL_WEIGHTS
  W <- sum(w)
  var_rate <- sum(
    (w / W)^2 * ifelse(
      age_tab$person_years > 0,
      age_tab$events / age_tab$person_years^2,
      0
    )
  )

  data.frame(
    asir = unname(gamma_result["adj.rate"]) * 1e5,
    asir_lower = unname(gamma_result["lci"]) * 1e5,
    asir_upper = unname(gamma_result["uci"]) * 1e5,
    asir_variance = var_rate * 1e10,
    events = sum(age_tab$events),
    person_years = sum(age_tab$person_years)
  )
}

compute_irr <- function(data, event_var, endage_var) {
  cohort_breaks <- seq(1920, 1980, by = 10)
  cohort_labels <- paste0("[", head(cohort_breaks, -1), ", ", tail(cohort_breaks, -1), ")")
  d <- data.frame(
    id = data$id,
    startage = data$startage,
    endage = data[[endage_var]],
    status = data[[event_var]],
    sex = factor(data$sex, levels = c("Female", "Male")),
    birth_cohort = cut(
      as.numeric(substr(as.character(data$DOB), 1, 4)),
      breaks = cohort_breaks, labels = cohort_labels,
      include.lowest = TRUE, right = FALSE
    )
  )
  sp <- survSplit(
    Surv(startage, endage, status) ~ ., data = d,
    cut = seq(1, 100, by = 1),
    start = "agestart", end = "agestop", event = "event",
    na.action = na.pass
  )
  sp$pyears <- sp$agestop - sp$agestart
  cases_male <- sum(sp$event[sp$sex == "Male"], na.rm = TRUE)
  cases_female <- sum(sp$event[sp$sex == "Female"], na.rm = TRUE)
  if (sum(sp$event, na.rm = TRUE) == 0) stop("No events found.")

  fit <- glm(
    event ~ sex + ns(agestart, df = 4) + birth_cohort + offset(log(pyears)),
    family = poisson, data = sp
  )
  sex_male <- summary(fit)$coefficients["sexMale", ]
  beta <- sex_male["Estimate"]
  se <- sex_male["Std. Error"]
  data.frame(
    irr = unname(exp(beta)),
    irr_lower = unname(exp(beta - 1.96 * se)),
    irr_upper = unname(exp(beta + 1.96 * se)),
    irr_p = unname(sex_male["Pr(>|z|)"]),
    male_events = cases_male,
    female_events = cases_female
  )
}

compute_outcome <- function(data, icd) {
  event_var <- paste0("CA_", icd)
  endage_var <- paste0("CA_", icd, "_endage")
  asir_m <- compute_asir(data[data$sex == "Male", , drop = FALSE],
                         event_var, endage_var)
  asir_f <- compute_asir(data[data$sex == "Female", , drop = FALSE],
                         event_var, endage_var)
  irr <- compute_irr(data, event_var, endage_var)

  difference <- asir_m$asir - asir_f$asir
  difference_se <- sqrt(asir_m$asir_variance + asir_f$asir_variance)
  z <- qnorm(0.975)

  data.frame(
    outcome = icd,
    irr = irr$irr,
    irr_lower = irr$irr_lower,
    irr_upper = irr$irr_upper,
    irr_p = irr$irr_p,
    asir_male = asir_m$asir,
    asir_male_lower = asir_m$asir_lower,
    asir_male_upper = asir_m$asir_upper,
    asir_female = asir_f$asir,
    asir_female_lower = asir_f$asir_lower,
    asir_female_upper = asir_f$asir_upper,
    asir_difference = difference,
    asir_difference_lower = difference - z * difference_se,
    asir_difference_upper = difference + z * difference_se,
    male_events = irr$male_events,
    female_events = irr$female_events
  )
}

IRR_ASIR_results <- do.call(rbind, lapply(OUTCOMES, function(icd) {
  compute_outcome(DATA, icd)
}))

write.csv(IRR_ASIR_results, "IRR_ASIR_results.csv", row.names = FALSE)
