# ============================================================================
# 02_PAF_ARD_AP.R
# Purpose: Sex-specific 5-year PAFs, joint PAFs, ARDs, and APs.
# ============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(AF)
})

source("AFcoxph_joint.R")

DATA <- readRDS("mock_data.rds")

PAF_input <- list(
  pan = c("eversmoke", "drink_2", "adiposity_2", "TPA_2",
          "proc_meat_2", "redmeat_2", "veg_2", "fruit_2", "HT2DM", "HLD", "Hgallstone"),
  `162` = c("eversmoke", "drink_2", "TPA_2",
            "proc_meat_2", "redmeat_2", "veg_2", "fruit_2"),
  `150` = c("eversmoke", "drink_2", "adiposity_2", "TPA_2",
            "proc_meat_2", "veg_2", "fruit_2", "HT2DM", "HLD"),
  `151` = c("eversmoke", "drink_2", "adiposity_2",
            "proc_meat_2", "fruit_2", "HT2DM", "HLD", "Hgallstone"),
  `153` = c("eversmoke", "drink_2", "adiposity_2", "TPA_2",
            "proc_meat_2", "redmeat_2", "veg_2", "fruit_2", "HT2DM", "HLD", "Hgallstone"),
  `154` = c("eversmoke", "drink_2", "adiposity_2", "TPA_2",
            "proc_meat_2", "redmeat_2", "veg_2", "fruit_2", "HT2DM", "HLD", "Hgallstone"),
  `155` = c("eversmoke", "drink_2", "adiposity_2", "TPA_2", "HT2DM", "HLD", "Hgallstone"),
  `156` = c("adiposity_2", "TPA_2", "HT2DM", "HLD", "Hgallstone"),
  `157` = c("eversmoke", "drink_2", "adiposity_2",
            "proc_meat_2", "redmeat_2", "fruit_2", "HT2DM", "HLD", "Hgallstone"),
  `189` = c("eversmoke", "adiposity_2", "HT2DM", "HLD", "Hgallstone"),
  `188` = c("eversmoke", "TPA_2", "veg_2", "fruit_2", "HT2DM", "HLD", "Hgallstone"),
  `14` = c("eversmoke", "drink_2", "adiposity_2",
           "proc_meat_2", "redmeat_2", "veg_2", "fruit_2", "HLD"),
  `193` = c("eversmoke", "drink_2", "adiposity_2", "TPA_2", "proc_meat_2", "HT2DM")
)

OUTCOMES <- names(PAF_input)

T_POINT <- 5

fit_cox <- function(data, icd, rfs) {
  event_var <- paste0("CA_", icd)
  time_var <- paste0("CA_", icd, "_py")
  needed <- c("startage", "birth_cohort", "edu_2", "inc_2", "country",
              event_var, time_var, rfs)
  d <- data[complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
  for (rf in rfs) {
    observed <- unique(as.character(d[[rf]]))
    if (!setequal(observed, c("0", "1"))) {
      stop(sprintf("%s / %s: %s must contain both 0 and 1 among Cox complete cases; observed: %s",
                   icd, as.character(unique(d$sex)), rf,
                   paste(observed, collapse = ", ")))
    }
    d[[rf]] <- as.numeric(as.character(d[[rf]]))
  }

  formula <- as.formula(paste0(
    "Surv(", time_var, ", ", event_var, ") ~ ",
    "startage + ridge(birth_cohort, edu_2, inc_2) + country + ",
    paste(rfs, collapse = " + ")
  ))
  fit <- coxph(formula, data = d, ties = "breslow")
  list(fit = fit, data = d)
}

paf_one <- function(fit, data, rf) {
  obj <- AF::AFcoxph(
    object = fit,
    data = data,
    exposure = rf,
    times = T_POINT
  )
  paf <- as.numeric(obj$AF.est[1])
  variance <- as.numeric(obj$AF.var[1])
  p_value <- 2 * stats::pnorm(-abs(paf / sqrt(variance)))
  data.frame(
    risk_factor = rf,
    paf = paf,
    paf_p = p_value,
    stringsAsFactors = FALSE
  )
}

paf_joint <- function(fit, data, rfs) {
  obj <- AFcoxph_joint(
    object = fit,
    data = data,
    exposures = rfs,
    times = T_POINT
  )
  paf <- as.numeric(obj$AF.est[1])
  variance <- as.numeric(obj$AF.var[1])
  p_value <- 2 * stats::pnorm(-abs(paf / sqrt(variance)))
  data.frame(
    risk_factor = "__JOINT__",
    paf = paf,
    paf_p = p_value,
    stringsAsFactors = FALSE
  )
}

compute_paf_for_sex <- function(data, icd, sex_value, rfs) {
  d <- subset(data, sex == sex_value)
  fitted <- fit_cox(d, icd, rfs)
  single <- do.call(rbind, lapply(rfs, function(rf) {
    paf_one(fitted$fit, fitted$data, rf)
  }))
  joint <- paf_joint(fitted$fit, fitted$data, rfs)
  out <- rbind(single, joint)
  out$sex <- sex_value
  out$outcome <- icd
  out$n_cox <- fitted$fit$n
  out$n_events <- fitted$fit$nevent
  out
}

# ASIR point estimates are reused for the ARD calculation.
if (!file.exists("IRR_ASIR_results.csv")) {
  stop("Run 01_IRR_ASIR.R before running this script.")
}
asir_lookup <- read.csv("IRR_ASIR_results.csv", stringsAsFactors = FALSE)[,
                        c("outcome", "asir_male", "asir_female")]

paf_results <- do.call(rbind, lapply(OUTCOMES, function(icd) {
  rfs <- PAF_input[[icd]]
  rbind(
    compute_paf_for_sex(DATA, icd, "Male", rfs),
    compute_paf_for_sex(DATA, icd, "Female", rfs)
  )
}))

paf_wide <- reshape(
  paf_results[, c("outcome", "risk_factor", "sex", "paf", "paf_p")],
  idvar = c("outcome", "risk_factor"),
  timevar = "sex",
  direction = "wide"
)
names(paf_wide) <- sub("^paf\\.", "paf_", names(paf_wide))
names(paf_wide) <- sub("^paf_p\\.", "paf_p_", names(paf_wide))

paf_wide <- merge(paf_wide, asir_lookup, by = "outcome", all.x = TRUE)
paf_wide$ard <- with(paf_wide,
                     asir_male * paf_Male - asir_female * paf_Female)
paf_wide$ap <- with(paf_wide,
                    ard / (asir_male - asir_female))
paf_wide$ap_percent <- paf_wide$ap * 100


write.csv(paf_wide, "PAF_ARD_AP_results.csv", row.names = FALSE)
