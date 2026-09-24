# ============================================================================
# 03_Bootstrap_ARD_CI.R
# Purpose: Obtain percentile bootstrap CIs for ARD using B = 100 replicates.
# ============================================================================

suppressPackageStartupMessages({
  library(survival)
})

B <- 100L
T_POINT <- 5

if (!file.exists("PAF_ARD_AP_results.csv")) {
  stop("Run 02_PAF_ARD_AP.R before running this bootstrap script.")
}
point_estimates <- read.csv("PAF_ARD_AP_results.csv",
                            stringsAsFactors = FALSE)[,
                            c("outcome", "risk_factor", "ard", "ap_percent")]
names(point_estimates)[names(point_estimates) == "ard"] <- "ard_estimate"

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

PAF_point <- function(fit, data, exposures, time_point) {
  bh <- survival::basehaz(fit, centered = TRUE)
  index <- findInterval(time_point, bh$time)
  H0 <- if (index == 0L) 0 else bh$hazard[index]
  lp_observed <- stats::predict(fit, newdata = data, type = "lp")
  data_counterfactual <- data
  data_counterfactual[exposures] <- lapply(
    data_counterfactual[exposures], function(x) 0
  )
  lp_counterfactual <- stats::predict(
    fit, newdata = data_counterfactual, type = "lp"
  )
  survival_observed <- mean(exp(-H0 * exp(lp_observed)))
  survival_counterfactual <- mean(exp(-H0 * exp(lp_counterfactual)))
  1 - (1 - survival_counterfactual) / (1 - survival_observed)
}

compute_asir_point <- function(data, event_var, endage_var) {
  age_breaks <- seq(40, 85, 5)
  age_labels <- paste(seq(40, 80, 5), seq(44, 84, 5), sep = "-")
  weights <- c(6, 6, 5, 4, 4, 3, 2, 1, 0.5)
  d <- data[complete.cases(data[, c("startage", event_var, endage_var)]), ,
            drop = FALSE]
  tmp <- data.frame(
    startage = d$startage,
    endage = d[[endage_var]],
    status = d[[event_var]]
  )
  sp <- survival::survSplit(
    Surv(startage, endage, status) ~ ., data = tmp,
    cut = age_breaks
  )
  sp <- sp[sp$startage >= 40 & sp$startage < 85, , drop = FALSE]
  sp$age_group <- cut(sp$startage, age_breaks, labels = age_labels,
                      right = FALSE)
  sp$py <- sp$endage - sp$startage
  age_tab <- aggregate(cbind(events = status, py = py) ~ age_group,
                       data = sp, FUN = sum)
  all_ages <- data.frame(age_group = factor(age_labels, levels = age_labels))
  age_tab <- merge(all_ages, age_tab, by = "age_group", all.x = TRUE,
                   sort = FALSE)
  age_tab$events[is.na(age_tab$events)] <- 0
  age_tab$py[is.na(age_tab$py)] <- 0
  sum(weights * ifelse(age_tab$py > 0, age_tab$events / age_tab$py, 0)) /
    sum(weights) * 1e5
}

one_bootstrap <- function(b, data, outcomes, paf_input) {
  set.seed(20260924 + b)
  male <- data[data$sex == "Male", , drop = FALSE]
  female <- data[data$sex == "Female", , drop = FALSE]
  male <- male[sample.int(nrow(male), nrow(male), replace = TRUE), , drop = FALSE]
  female <- female[sample.int(nrow(female), nrow(female), replace = TRUE), , drop = FALSE]

  result <- list()
  for (icd in outcomes) {
    rfs <- paf_input[[icd]]
    event_var <- paste0("CA_", icd)
    endage_var <- paste0("CA_", icd, "_endage")
    asir_m <- compute_asir_point(male, event_var, endage_var)
    asir_f <- compute_asir_point(female, event_var, endage_var)

    paf_m <- fit_cox(male, icd, rfs)
    paf_f <- fit_cox(female, icd, rfs)
    paf_rows <- c(
      setNames(lapply(rfs, function(rf) {
        PAF_point(paf_m$fit, paf_m$data, rf, T_POINT)
      }), rfs),
      list(`__JOINT__` = PAF_point(paf_m$fit, paf_m$data, rfs, T_POINT))
    )
    paf_rows_f <- c(
      setNames(lapply(rfs, function(rf) {
        PAF_point(paf_f$fit, paf_f$data, rf, T_POINT)
      }), rfs),
      list(`__JOINT__` = PAF_point(paf_f$fit, paf_f$data, rfs, T_POINT))
    )

    for (rf in c(rfs, "__JOINT__")) {
      ard <- asir_m * paf_rows[[rf]] - asir_f * paf_rows_f[[rf]]
      result[[length(result) + 1L]] <- data.frame(
        b = b, outcome = icd, risk_factor = rf, ard = ard
      )
    }
  }
  do.call(rbind, result)
}

bootstrap_list <- lapply(seq_len(B), one_bootstrap, data = DATA,
                         outcomes = OUTCOMES, paf_input = PAF_input)
bootstrap_results <- do.call(rbind, bootstrap_list)

bootstrap_ci <- do.call(rbind, lapply(split(
  bootstrap_results,
  list(bootstrap_results$outcome, bootstrap_results$risk_factor)
), function(d) {
  q <- stats::quantile(d$ard, probs = c(0.025, 0.975), na.rm = TRUE)
  data.frame(
    outcome = d$outcome[1],
    risk_factor = d$risk_factor[1],
    ard_lower = unname(q[1]),
    ard_upper = unname(q[2])
  )
}))

bootstrap_ci <- merge(point_estimates, bootstrap_ci,
                      by = c("outcome", "risk_factor"), sort = FALSE)
bootstrap_ci <- bootstrap_ci[, c("outcome", "risk_factor", "ard_estimate",
                                "ard_lower", "ard_upper", "ap_percent")]

write.csv(bootstrap_results, "ARD_bootstrap_replicates.csv", row.names = FALSE)
write.csv(bootstrap_ci, "ARD_bootstrap_CI.csv", row.names = FALSE)

