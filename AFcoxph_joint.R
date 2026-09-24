# ============================================================================
# Joint attributable fraction for a Cox model.
#
# This is based on the AF::AFcoxph implementation used in the analysis.  The
# counterfactual data set differs by setting all supplied binary exposures to
# their low-risk reference level (0) simultaneously.
# ============================================================================

AFcoxph_joint <- function(object, data, exposures, times, clusterid) {
  if (!identical(object$method, "breslow")) {
    stop("AFcoxph_joint requires a Cox model fitted with ties = 'breslow'.")
  }
  if (length(exposures) < 1L) {
    stop("At least one exposure must be supplied.")
  }
  if (!all(exposures %in% names(data))) {
    stop("All exposures must be columns in data.")
  }
  if (!all(vapply(data[exposures], function(x) {
    all(na.omit(x) %in% c(0, 1))
  }, logical(1)))) {
    stop("All joint exposures must be coded 0/1.")
  }

  formula <- object$formula
  data <- as.data.frame(data)
  rownames(data) <- seq_len(nrow(data))
  model_matrix <- model.matrix(object = formula, data = data)
  complete <- as.integer(rownames(model_matrix))
  data <- data[complete, , drop = FALSE]

  object_detail <- survival::coxph.detail(object)
  if (missing(times)) times <- object_detail$time

  data0 <- data
  data0[exposures] <- lapply(data0[exposures], function(x) 0)

  response_terms <- stats::delete.response(stats::terms(object))
  design <- as.matrix(stats::model.matrix(response_terms, data = data))[, -1,
                                                                         drop = FALSE]
  design0 <- as.matrix(stats::model.matrix(response_terms, data = data0))[, -1,
                                                                            drop = FALSE]

  detail_time <- object_detail$time
  dH0 <- object_detail$hazard
  H0 <- cumsum(dH0)
  H0step <- stats::stepfun(detail_time, c(0, H0))

  response_text <- deparse(formula[[2]])
  response_vars <- all.vars(formula[[2]])
  if (length(response_vars) < 2L) {
    stop("The model must use a two-column Surv(time, event) response.")
  }
  endvar <- response_vars[1]
  eventvar <- response_vars[2]

  ord <- order(data[[endvar]])
  data <- data[ord, , drop = FALSE]
  design <- design[ord, , drop = FALSE]
  design0 <- design0[ord, , drop = FALSE]
  data0 <- data0[ord, , drop = FALSE]

  n <- nrow(data)
  npar <- length(object$coef)
  n_cases <- sum(data[[eventvar]])

  epred <- stats::predict(object, newdata = data, type = "risk")
  epred0 <- stats::predict(object, newdata = data0, type = "risk")
  score_beta <- as.matrix(stats::residuals(object, type = "score"))[ord, , drop = FALSE]

  H0res <- rep(0, n)
  dH0_untied <- rep(dH0, object_detail$nevent) /
    rep(object_detail$nevent, object_detail$nevent)
  H0res[data[[eventvar]] == 1] <- dH0_untied * n

  E <- matrix(0, nrow = n, ncol = npar)
  means <- as.matrix(object_detail$means)
  means <- means[rep(seq_len(nrow(means)), object_detail$nevent), , drop = FALSE]
  E[data[[eventvar]] == 1, ] <- means

  has_cluster <- !missing(clusterid)
  if (has_cluster) {
    clusters <- data[[clusterid]]
    n_cluster <- length(unique(clusters))
  }

  S_est <- numeric(length(times))
  S0_est <- numeric(length(times))
  AF_var <- numeric(length(times))
  S_var <- numeric(length(times))
  S0_var <- numeric(length(times))

  for (i in seq_along(times)) {
    time_i <- times[i]
    score_S <- exp(-H0step(time_i) * epred)
    score_S0 <- exp(-H0step(time_i) * epred0)
    score_H0 <- H0res * (data[[endvar]] <= time_i)

    score_equations <- cbind(score_S, score_S0, score_H0, score_beta)
    if (has_cluster) {
      split_scores <- split(seq_len(nrow(score_equations)), clusters)
      score_equations <- do.call(rbind, lapply(split_scores, function(ii) {
        colSums(score_equations[ii, , drop = FALSE])
      }))
    }

    meat <- stats::var(score_equations, na.rm = TRUE)
    hessian_S <- c(
      -1, 0,
      mean(epred * score_S),
      colMeans(design * H0step(time_i) * epred * score_S)
    )
    hessian_S0 <- c(
      0, -1,
      mean(epred0 * score_S0),
      colMeans(design0 * H0step(time_i) * epred0 * score_S0)
    )
    hessian_H0 <- c(
      0, 0, -1,
      -colMeans(E * score_H0, na.rm = TRUE)
    )
    hessian_beta <- cbind(matrix(0, nrow = npar, ncol = 3),
                          -solve(stats::vcov(object)) / n)
    bread <- rbind(hessian_S, hessian_S0, hessian_H0, hessian_beta)
    bread_inv <- solve(bread)

    sandwich <- bread_inv %*% meat %*% t(bread_inv)
    if (has_cluster) {
      sandwich <- sandwich * n_cluster / n^2
    } else {
      sandwich <- sandwich / n
    }
    sandwich <- sandwich[1:2, 1:2, drop = FALSE]

    S_est[i] <- mean(score_S, na.rm = TRUE)
    S0_est[i] <- mean(score_S0, na.rm = TRUE)
    gradient <- matrix(c(
      -(1 - S0_est[i]) / (1 - S_est[i])^2,
      1 / (1 - S_est[i])
    ), nrow = 2)
    AF_var[i] <- as.numeric(t(gradient) %*% sandwich %*% gradient)
    S_var[i] <- sandwich[1, 1]
    S0_var[i] <- sandwich[2, 2]
  }

  AF_est <- 1 - (1 - S0_est) / (1 - S_est)
  out <- list(
    AF.est = AF_est,
    AF.var = AF_var,
    S.est = S_est,
    S0.est = S0_est,
    S.var = S_var,
    S0.var = S0_var,
    objectcall = object$call,
    exposures = exposures,
    outcome = response_text,
    object = object,
    formula = formula,
    n = n,
    n.cases = n_cases,
    n.cluster = if (has_cluster) n_cluster else 0L,
    times = times
  )
  class(out) <- "jointAF"
  out
}
