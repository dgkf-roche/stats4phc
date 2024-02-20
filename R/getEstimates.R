
# Helper function to check the user input - bins and quantiles arguments
checkBinsQuantiles <- function(bins, quantiles, score) {

  checkmate::assert(
    checkmate::check_numeric(bins, any.missing = FALSE, sorted = TRUE),
    checkmate::check_null(bins)
  )
  checkmate::assert(
    checkmate::check_integerish(quantiles, lower = 1),
    checkmate::check_null(quantiles)
  )

  # Check combination of quantiles and bins
  if (!is.null(quantiles) && !is.null(bins)) {
    stop("bins and quantiles cannot be specified together, choose one and set the other to NULL")
  }

  if (is.null(quantiles) && is.null(bins)) {
    quantiles <- 10
  }

  # Further check for bins
  if (length(bins) > 1) {
    if (bins[1] > min(score)) {
      stop(paste("The first element of bins must be <= min(score), not", bins[1]))
    }
    if (bins[length(bins)] < max(score)) {
      stop(paste("The last element of bins must be >= max(score), not", bins[length(bins)]))
    }
  }

  # Scalar bins -> number of intervals
  if (length(bins) == 1 && bins <= 1) {
    stop("bins must be > 1 when provided as a scalar (i.e. number of bins)")
  }

  # Check discrete score
  lu <- length(unique(score))
  if (length(bins) == 1 && bins > lu) {
    warning(
      paste0(
        "The number of `bins` (", bins, ") ",
        "is > the number of unique score values (", lu, "). ",
        "The results may be unreliable."
      )
    )
  }
  if (length(bins) > 1 && length(bins) - 1 > lu) {
    warning(
      paste0(
        "The number of `bins` (", length(bins) - 1, ") ",
        "is > the number of unique score values (", lu, "). ",
        "The results may be unreliable."
      )
    )
  }
  if (!is.null(quantiles) && lu <= 10) {
    warning(
      "Using the quantile method for non-continuous score. ",
      "The results may be unreliable."
    )
  }

  return(list(bins = bins, quantiles = quantiles))
}


# Helper function to check the user input - errorbar.sem argument
checkErrorbarSem <- function(errorbar.sem) {
  checkmate::assert(
    checkmate::check_number(errorbar.sem, lower = 0), # this checks >= 0
    checkmate::check_null(errorbar.sem)
  )
  if (!is.null(errorbar.sem)) {
    stopifnot("`errorbar.sem` must be > 0" = errorbar.sem > 0)
  }
  return(errorbar.sem)
}


#' Binned Risk Estimates
#'
#' Calculates bins based on number of evenly spaced bins or n-tiles.
#' Determines average risk within bins, used for risk estimates.
#'
#' @inheritParams riskProfile
#' @param quantiles Numeric; quantiles to split bins.
#' @param bins Numeric; number of evenly spaced bins or bin locations.
#' @param right Logical indicating right closed interval. Defaults to `TRUE`.
#' @param errorbar.sem Scalar numeric representing the number of standard error from the means
#' (SEM) used to calculate risk error bar.
#'
#' @return A data frame with 4 columns
#' (score, score percentile, outcome, estimate).
#' Additionally, there is an attribute "errorbar" holding the error-bar data if
#' `errorbar.sem` was specified.
#'
#' @seealso [getASISest()] [getCGAMest()] [getGAMest()] [getMSPLINEest()] [getPAVAest()]
#'
#' @export
#'
#' @examples
#' # Read in example data
#' auroc <- read.csv(system.file("extdata", "sample.csv", package = "stats4phc"))
#' rscore <- auroc$predicted
#' truth <- as.numeric(auroc$actual)
#'
#' getBINNEDest(outcome = truth, score = rscore)
#'
getBINNEDest <- function(outcome,
                         score,
                         quantiles = NULL,
                         bins = NULL,
                         right = TRUE,
                         errorbar.sem = NULL) {

  # Argument checks
  checkmate::assert_numeric(outcome)
  checkmate::assert_numeric(score, len = length(outcome))
  checkmate::assert_flag(right)
  errorbar.sem <- checkErrorbarSem(errorbar.sem)

  # Check bins and quantiles
  bqp <- checkBinsQuantiles(bins = bins, quantiles = quantiles, score = score)

  # Retrieve data summaries
  df <- getSummaries(
    outcome = outcome, score = score,
    quantiles = bqp$quantiles, bins = bqp$bins,
    right = right
  )

  if (!is.null(errorbar.sem)) {
    errorbar <- getERRORest(binlvl = df[["binlvl"]], z = errorbar.sem) %>%
      mutate(
        percentile = df[["binlvl"]][["riskpercentile"]],
        bin.mid = df[["binlvl"]][["avg.outcome"]]
      )
  } else {
    errorbar <- NULL
  }

  # Create binned.data
  binned.data <- data.frame(
    score = df[["binlvl"]][["avg.risk"]],
    percentile = df[["binlvl"]][["riskpercentile"]],
    outcome = NA,
    estimate = df[["binlvl"]][["avg.outcome"]]
  )
  attr(binned.data, "errorbar") <- errorbar

  return(binned.data)
}



#' Get Summaries: observation level & bin level dataframes
#'
#' Given vectors of outcomes and risk scores, the function will return a list of
#' observation level and bin level summary dataframes for the data.
#'
#' Observation level dataframe has columns for outcome, riskscore, risk percentile, bin number,
#' and corresponding minimum and maximum score for that bin.
#'
#' Bin level dataframe has columns indicating bin number and the observation count,
#' number of events, average outcome, average risk, and standard deviation of risk,
#' within each of the bins. Risk percentile and bin intervals are also provided.
#'
#' @return List of observation level and bin level dataframes.
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#' auroc <- read.csv(system.file("extdata", "sample.csv", package = "stats4phc"))
#' truth <- as.numeric(auroc$actual)
#' rscore <- auroc$predicted
#'
#' # Bin by quantiles
#' getSummaries(
#'   outcome = truth,
#'   score = rscore,
#'   quantiles = 10
#' )
#'
#' # Bin by specific percentiles
#' getSummaries(
#'   outcome = truth,
#'   score = rscore,
#'   quantiles = 0,
#'   bin = c(0, 0.25, 0.5, 0.8, 1)
#' )
#'
getSummaries <- function(outcome,
                         score,
                         quantiles = NULL,
                         bins = NULL,
                         right = TRUE) {

  stopifnot(
    !is.null(quantiles) | !is.null(bins),
    quantiles != 0 | bins != 0,
    is.logical(right)
  )

  cdf.fit <- ecdf(score)

  if (!is.null(quantiles)) {
    bin.int <- Hmisc::cut2(score, g = quantiles)

  } else if (length(bins) > 1) {
    bin.int <- cut(cdf.fit(score), breaks = bins, include.lowest = TRUE, right = right)

  } else if (length(bins) == 1) {
    bin.int <- cut(score, breaks = bins, include.lowest = TRUE, right = right)

  } else {
    stop("Unrecognized option")
  }
  
  # get numeric label of interval
  bin.num <- as.numeric(bin.int)

  # get interval borders (min, max)
  min_max <- strsplit(
    gsub("(?![-,.])[[:punct:]]", "", trimws(as.character(bin.int)), perl = TRUE),
    ","
  )
  min_max <- lapply(
    min_max, \(x) `if`(length(x) < 2, rep(x, 2), x)
  )

  # Observation level data: outcome, rscore, observation percentile, and interval
  obslvl <- data.frame(
    outcome = outcome,
    score = score,
    riskpercentile = cdf.fit(score),
    bin = bin.num,
    interval = as.character(bin.int),
    min = vapply(min_max, "[[", character(1), 1),
    max = vapply(min_max, "[[", character(1), 2)
  ) %>%
    arrange(.data$score)

  # Bin level data: within each bin: n, avg risk, sd risk, quantile, error bar
  binlvl <- obslvl %>%
    group_by(.data$bin, .data$interval) %>%
    summarise(
      n = dplyr::n(),
      events = sum(.data$outcome),
      avg.outcome = mean(.data$outcome),
      sd.outcome = sd(.data$outcome, na.rm = TRUE),
      avg.risk = mean(.data$score, na.rm = TRUE),
      sd.risk = sd(.data$score),
      riskpercentile = max(.data$riskpercentile),
      .groups = "drop"
    )

  # Replace NAs with 0
  binlvl[is.na(binlvl)] <- 0

  return(list(obslvl = obslvl, binlvl = binlvl))
}


# ERROR BAR Estimates
getERRORest <- function(binlvl, z) {

  stopifnot(
    is.data.frame(binlvl),
    is.numeric(z)
  )

  binlvl %>%
    mutate(
      bin.low = .data$avg.outcome - (z * (.data$sd.outcome / sqrt(.data$n))),
      bin.low = ifelse(.data$bin.low < 0, 0, .data$bin.low),
      bin.high = .data$avg.outcome + (z * (.data$sd.outcome / sqrt(.data$n))),
      midquantile = .data$riskpercentile - (diff(c(0, .data$riskpercentile)) / 2)
    ) %>%
    select(all_of(c("midquantile", "bin.high", "bin.low"))) %>%
    tidyr::replace_na(list(bin.high = 0, bin.low = 0))
}


#' PAVA Risk Estimates
#'
#' Determines isotonic regression estimates via pava, given a vector of binary outcomes,
#' and a vector of scores.
#'
#' @inheritParams riskProfile
#' @param weights Vector of numerics to specify PAVA observation weighting.
#' @param ties String to specify how ties should be handled for PAVA.
#' @param low_events Numeric, specifying number of events in the lowest bin.
#' @param low_nonevents Numeric, specifying number of nonevents in the lowest bin.
#' @param high_events Numeric, specifying number of events in the highest bin.
#' @param high_nonevents Numeric, specifying number of nonevents in the highest bin.
#' @param hilo_obs Numeric, specifying number of observations in the highest and lowest bins.
#'
#' @return A data frame with 4 columns
#' (score, score percentile, outcome, estimate).
#'
#' @seealso [getASISest()] [getBINNEDest()] [getCGAMest()] [getGAMest()] [getMSPLINEest()]
#'
#' @export
#'
#' @examples
#' # Read in example data
#' auroc <- read.csv(system.file("extdata", "sample.csv", package = "stats4phc"))
#' rscore <- auroc$predicted
#' truth <- as.numeric(auroc$actual)
#'
#' tail(getPAVAest(outcome = truth, score = rscore), 10)
#'
getPAVAest <- function(outcome,
                       score,
                       weights = rep(1, length(outcome)),
                       ties = "primary",
                       low_events = NULL,
                       low_nonevents = NULL,
                       high_events = NULL,
                       high_nonevents = NULL,
                       hilo_obs = NULL) {

  checkmate::assert_numeric(outcome)
  checkmate::assert_numeric(score, len = length(outcome))
  checkmate::assert_numeric(weights, any.missing = FALSE, len = length(outcome))
  checkmate::assert_character(ties, any.missing = FALSE, len = 1)
  checkmate::assert(
    checkmate::check_integerish(low_events, lower = 1, any.missing = FALSE, len = 1),
    checkmate::check_null(low_events)
  )
  checkmate::assert(
    checkmate::check_integerish(low_nonevents, lower = 1, any.missing = FALSE, len = 1),
    checkmate::check_null(low_nonevents)
  )
  checkmate::assert(
    checkmate::check_integerish(high_events, lower = 1, any.missing = FALSE, len = 1),
    checkmate::check_null(high_events)
  )
  checkmate::assert(
    checkmate::check_integerish(high_nonevents, lower = 1, any.missing = FALSE, len = 1),
    checkmate::check_null(high_nonevents)
  )
  checkmate::assert(
    checkmate::check_integerish(hilo_obs, lower = 1, any.missing = FALSE, len = 1),
    checkmate::check_null(hilo_obs)
  )

  pava.est <- isotone::gpava(z = score, y = outcome, weights = weights, ties = ties)$x

  # if constrained, then replace percentiles....
  check <- any(
    is.numeric(low_events), is.numeric(low_nonevents), is.numeric(high_events),
    is.numeric(high_nonevents), is.numeric(hilo_obs)
  )
  if (check) {
    percentile <- getConstraints(
      outcome = outcome, rscore = score,
      low_events = low_events, low_nonevents = low_nonevents,
      high_events = high_events, high_nonevents = high_nonevents,
      hilo_obs = hilo_obs
    )
  } else {
    percentile <- ecdf(score)(score)
  }

  return(data.frame(score, percentile, outcome, estimate = pava.est))
}


#' Constrained Risk Percentile Estimates
#'
#' Adjusts PAVA risk percentile estimates for the first and last bin, to meet criteria
#' for events, non-events, or total observation count.
#'
#' @inheritParams riskProfile
#'
#' @return A vector of constrained risk percentiles.
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#' auroc <- read.csv(system.file("extdata", "sample.csv", package = "stats4phc"))
#' truth <- as.numeric(auroc$actual)
#' rscore <- auroc$predicted
#'
#' getConstraints(outcome = truth, rscore = rscore, low_events = 3, high_nonevents = 3)
#'
getConstraints <- function(outcome,
                           rscore,
                           low_events = NULL, # min events in lower bin (useful for a PC)
                           low_nonevents = NULL, # min non-events in lower bin
                           high_events = NULL, # min events in upper bin
                           high_nonevents = NULL, # min non-events in upper bin (useful for a PC)
                           hilo_obs = NULL) { # min total obs in upper AND lower bin

  n <- length(rscore)
  rscore <- seq_along(rscore) / length(rscore)
  # rscore_cons: scores with binning constraint (same as rscore if no constraint specified)
  rscore_cons <- rscore

  if (!is.null(low_events) && low_events == 0) low_events <- NULL
  if (!is.null(low_nonevents) && low_nonevents == 0) low_nonevents <- NULL
  if (!is.null(high_events) && high_events == 0) high_events <- NULL
  if (!is.null(high_nonevents) && high_nonevents == 0) high_nonevents <- NULL
  if (!is.null(hilo_obs) && hilo_obs == 0) hilo_obs <- NULL

  if (is.numeric(low_events) && is.numeric(low_nonevents)) {
    warning(
      paste(
        "Specified both a minimum number of events and non-events for the lower bin.",
        "Combining for total observations instead."
      )
    )
    hilo_obs <- round(low_events + low_nonevents)
    rscore_cons[1:hilo_obs] <- min(rscore)

  } else if (is.numeric(high_events) && is.numeric(high_nonevents)) {
    warning(
      paste(
        "Specified both a minimum number of events and non-events for the upper bin.",
        "Combining for total observations instead."
      )
    )
    hilo_obs <- round(high_events + high_nonevents)
    rscore_cons[(n + 1 - hilo_obs):n] <- max(rscore)

  } else {
    # Apply upper / lower constraints
    if (is.numeric(low_events)) {
      low_events <- round(low_events)
      indlo <- match(TRUE, cumsum(outcome) == low_events)
      rscore_cons[1:indlo] <- min(rscore)
    }
    if (is.numeric(low_nonevents)) {
      low_nonevents <- round(low_nonevents)
      indlo <- match(TRUE, cumsum(1 - outcome) == low_nonevents)
      rscore_cons[1:indlo] <- min(rscore)
    }
    if (is.numeric(high_nonevents)) {
      high_nonevents <- round(high_nonevents)
      indhi <- match(TRUE, cumsum(1 - outcome[n:1]) == high_nonevents)
      rscore_cons[(n + 1 - indhi):n] <- max(rscore)
    }

    if (is.numeric(high_events)) {
      high_events <- round(high_events)
      indhi <- match(TRUE, cumsum(outcome[n:1]) == high_events)
      rscore_cons[(n + 1 - indhi):n] <- max(rscore)
    }

    if (is.numeric(hilo_obs)) {
      hilo_obs <- round(hilo_obs)
      rscore_cons[1:hilo_obs] <- min(rscore)
      rscore_cons[(n + 1 - hilo_obs):n] <- max(rscore)
    }
  }

  return(as.vector(rscore_cons))
}


#' GAM Risk Estimates
#'
#' Fits a Generalized Additive Model to estimate risk, given a vector of binary outcome,
#' and a vector of scores.
#'
#' @inheritParams riskProfile
#' @param k Numeric to specify the upper limit of basis functions to fit for GAM.
#' See [mgcv::s()] for more details. Defaults to -1.
#' @param bs Character string to specify spline type.
#' See [mgcv::s()] for more details. Defaults to `"tp"`.
#' @param method Character string to specify method type.
#' See [mgcv::s()] for more details. Defaults to "REML".
#' @param logscores Logical; if `TRUE`, fit gam on log scores. Defaults to `FALSE`.
#' @param fitonPerc Logical; if `TRUE`, fit gam on risk percentiles. Defaults to `TRUE`.
#'
#' @return A data frame with 4 columns
#' (score, score percentile, outcome, estimate).
#'
#' @seealso [getASISest()] [getBINNEDest()] [getCGAMest()] [getMSPLINEest()] [getPAVAest()]
#'
#' @export
#'
#' @examples
#' # Read in example data
#' auroc <- read.csv(system.file("extdata", "sample.csv", package = "stats4phc"))
#' rscore <- auroc$predicted
#' truth <- as.numeric(auroc$actual)
#'
#' tail(getGAMest(outcome = truth, score = rscore), 10)
#'
getGAMest <- function(outcome,
                      score,
                      k = -1,
                      bs = "tp",
                      method = "REML",
                      logscores = FALSE,
                      fitonPerc = TRUE) {

  checkmate::assert_numeric(outcome)
  checkmate::assert_numeric(score, len = length(outcome))
  checkmate::assert_number(k)
  checkmate::assert_character(bs, any.missing = FALSE, len = 1)
  checkmate::assert_character(method, any.missing = FALSE, len = 1)
  checkmate::assert_flag(logscores)
  checkmate::assert_flag(fitonPerc)

  mygrid <- ecdf(score)(score)
  df <- data.frame(outcome = outcome, score = score, perc = mygrid)

  # mgcv::s does not work in formula, need to define it here
  s <- mgcv::s

  if (fitonPerc && logscores) {
    formula <- outcome ~ s(log(perc), k = k, bs = bs)
  } else if (fitonPerc && !logscores) {
    formula <- outcome ~ s(perc, k = k, bs = bs)
  } else if (!fitonPerc && logscores) {
    formula <- outcome ~ s(log(score), k = k, bs = bs)
  } else {
    formula <- outcome ~ s(score, k = k, bs = bs)
  }

  gam.fit <- mgcv::gam(formula, data = df, family = "binomial", method = method)

  gam.est <- mgcv::predict.gam(gam.fit, type = "response")

  return(data.frame(score, percentile = mygrid, outcome, estimate = gam.est))
}



#' Constrained GAM (cgam) Risk Estimates
#'
#' Fits a Constrained Generalized Additive Model to estimate risk,
#' given a vector of binary outcomes and a vector of scores.
#'
#' @inheritParams riskProfile
#' @param numknots Numeric to specify the number of knots.
#' Passed to the `smoother` function. Defaults to 3.
#' @param smoother Character string to specify the smoother (from cgam package).
#' Defaults to "s.incr".
#' @param logscores Logical; if `TRUE`, fit gam on log scores. Defaults to `FALSE`.
#' @param fitonPerc Logical; if `TRUE`, fit gam on risk percentiles. Defaults to `TRUE`.
#'
#' @return A data frame with 4 columns
#' (score, score percentile, outcome, estimate).
#'
#' @seealso [getASISest()] [getBINNEDest()] [getGAMest()] [getMSPLINEest()] [getPAVAest()]
#'
#' @export
#'
#' @examples
#' # Read in example data
#' auroc <- read.csv(system.file("extdata", "sample.csv", package = "stats4phc"))
#' rscore <- auroc$predicted
#' truth <- as.numeric(auroc$actual)
#'
#' tail(getCGAMest(outcome = truth, score = rscore), 10)
#'
getCGAMest <- function(outcome,
                       score,
                       numknots = 0,
                       smoother = "s.incr",
                       logscores = FALSE,
                       fitonPerc = TRUE) {

  checkmate::assert_numeric(outcome)
  checkmate::assert_numeric(score, len = length(outcome))
  checkmate::assert_number(numknots)
  checkmate::assert_character(smoother, any.missing = FALSE, len = 1)
  checkmate::assert_flag(logscores)
  checkmate::assert_flag(fitonPerc)

  mygrid <- ecdf(score)(score)
  df <- data.frame(outcome = outcome, score = score, perc = mygrid)

  # cgam::s does not work in formula, need to define it here
  assign(smoother, do.call(`::`, list(pkg = "cgam", name = smoother)))

  formula <- as.formula(
    paste0(
      "outcome ~ ",
      smoother, "(",
      `if`(logscores, "log("),
      `if`(fitonPerc, "perc", "score"),
      `if`(logscores, ")"),
      ", numknots = ", numknots, ")"
    )
  )

  cgam.fit <- cgam::cgam(formula, data = df, family = "binomial")

  cgam.est <- cgam::predict.cgam(cgam.fit, type = "response")$fit

  return(data.frame(score, percentile = mygrid, outcome, estimate = cgam.est))
}


### MSPLINE ESTIMATES ###

# Function written by willtownes, joseph paulson.
mspline <- function(x, y, k = 10, lower = NA, upper = NA) {
  # fits a monotonic spline to data
  # small values of k= more smoothing (flatter curves)
  # large values of k= more flexible (wiggly curves)
  # k is related to effective degrees of freedom and number of knots
  # use unconstrained gam to get rough parameter estimates
  # lower, upper optional bounds on the function
  # basically a slight modifimessageion of an example in the mgcv::pcls documentation
  dat <- data.frame(x = x, y = y)
  s <- mgcv::s # mgcv::s does not work in formula, need to define it here
  init_gam <- mgcv::gam(y ~ s(x, k = k, bs = "cr"))
  # Create Design matrix, constraints etc. for monotonic spline....
  sm <- mgcv::smoothCon(s(x, k = k, bs = "cr"), dat, knots = NULL)[[1]]
  mc <- mgcv::mono.con(sm$xp, lower = lower, upper = upper) # monotonicity constraints
  M <- list(
    X = sm$X, y = y, # design matrix, outcome
    C = matrix(0, 0, 0), # equality constraints (none)
    Ain = mc$A, bin = mc$b, # inequality constraints
    sp = init_gam$sp, p = sm$xp, # initial guesses for param estimates
    S = sm$S, # smoothness penalty matrix
    w = y * 0 + 1, off = 0 # weights, offset
  )
  # fit spine using penalized constrained least squares
  p <- mgcv::pcls(M)
  return(list(sm = sm, p = p))
}

# Function written by joseph paulson
predict.mspline <- function(msp, x) {
  # using the monotone spline msp, predict values for the vector x
  as.vector(mgcv::Predict.matrix(msp$sm, data.frame(x = x)) %*% msp$p)
}



#' Monotone Spline Risk Estimates
#'
#' Fits a Monotone constrained Generalized Additive Model (GAM) to estimate risk,
#' given a vector of binary outcomes and a vector of scores.
#'
#' @inheritParams getGAMest
#'
#' @return A data frame with 4 columns
#' (score, score percentile, outcome, estimate).
#'
#' @seealso [getASISest()] [getBINNEDest()] [getCGAMest()] [getGAMest()] [getPAVAest()]
#'
#' @export
#'
#' @examples
#' # Read in example data
#' auroc <- read.csv(system.file("extdata", "sample.csv", package = "stats4phc"))
#' rscore <- auroc$predicted
#' truth <- as.numeric(auroc$actual)
#'
#' tail(getMSPLINEest(outcome = truth, score = rscore), 10)
#'
getMSPLINEest <- function(outcome,
                          score,
                          k = 10,
                          fitonPerc = TRUE) {

  checkmate::assert_numeric(outcome)
  checkmate::assert_numeric(score, len = length(outcome))
  checkmate::assert_integerish(k, len = 1, any.missing = FALSE)
  checkmate::assert_flag(fitonPerc)

  stopifnot(!is.null(k), is.logical(fitonPerc))

  scorefit <- ecdf(score)
  mygrid <- scorefit(score)

  if (!fitonPerc) {
    fitspl <- mspline(x = score, y = outcome, k = k)
    mspline.est <- predict.mspline(fitspl, score)
  } else {
    fitspl <- mspline(x = mygrid, y = outcome, k = k)
    mspline.est <- predict.mspline(fitspl, mygrid)
  }

  mspline.est[mspline.est < 0] <- 0
  mspline.est[mspline.est > 1] <- 1

  return(data.frame(score, percentile = mygrid, outcome, estimate = mspline.est))
}


#' "As is" estimates
#'
#' This function does no estimation, but uses the score as it is
#' (it works like an identity function).
#'
#' @inheritParams getGAMest
#'
#' @return A data frame with 4 columns
#' (score, score percentile, outcome, estimate).
#'
#' @seealso [getBINNEDest()] [getCGAMest()] [getGAMest()] [getMSPLINEest()] [getPAVAest()]
#'
#' @export
#'
#' @examples
#' # Read in example data
#' auroc <- read.csv(system.file("extdata", "sample.csv", package = "stats4phc"))
#' rscore <- auroc$predicted
#' truth <- as.numeric(auroc$actual)
#'
#' tail(getASISest(outcome = truth, score = rscore), 10)
#'
getASISest <- function(outcome, score) {
  # Argument checks
  checkmate::assert_numeric(outcome)
  checkmate::assert_numeric(score, len = length(outcome))
  return(
    data.frame(
      score = score,
      percentile = ecdf(score)(score),
      outcome,
      estimate = score
    )
  )
}


# Returns Risk Estimates.
# This function calls all the other getXXXest functions
# Used for Predictiveness Curve Data
getEsts <- function(methods, outcome, score) {

  # check methods
  stopifnot(is.list(methods))

  # Get function names (and perform checks)
  fun.names <- getEstMethods(methods, with.names = TRUE)

  # Run estimations
  m.est <- lapply(methods, \(x) getEst(x, outcome = outcome, score = score))

  # Special case: get errorbar data if existing
  idx.binned <- fun.names == "binned"
  m.er <- lapply(
    which(idx.binned),
    \(i) attr(m.est[[i]], "errorbar")
  )

  # Bind together and convert to long data frame
  m.est <- bind_rows(m.est, .id = "method")
  rownames(m.est) <- NULL

  # Bind together and convert to long data frame; otherwise return NULL
  m.error <- bind_rows(m.er, .id = "method")
  if (nrow(m.error) == 0) {
    m.error <- NULL
  }

  # Return indexes of step methods and risk methods
  idx.step <- fun.names %in% c("binned", "pava")
  names(idx.step) <- names(fun.names)
  idx.asis <- fun.names == "asis"

  # Check if asis was called multiple times
  if (sum(idx.asis) > 1) {
    stop("Please use 'asis' just once (as it does not have any further arguments).")
  }

  return(
    list(
      plotdata = m.est, errorbardata = m.error,
      idx.step = idx.step, idx.asis = idx.asis,
      idx.binned = idx.binned, idx.pava = fun.names == "pava"
    )
  )
}

# Define estimation functions
est.funs <- function() {
  list(
    gam = getGAMest,
    cgam = getCGAMest,
    mspline = getMSPLINEest,
    binned = getBINNEDest,
    pava = getPAVAest,
    asis = getASISest
  )
}

# Generic and S3 methods for estimations
getEst <- function(x, outcome, score) {
  UseMethod("getEst")
}


getEst.list <- function(x, outcome, score) {
  do.call(
    est.funs()[[x[["method"]]]],
    append(list(outcome = outcome, score = score), x[names(x) != "method"])
  )
}


getEst.character <- function(x, outcome, score) {
  do.call(
    est.funs()[[x]],
    list(outcome = outcome, score = score)
  )
}


getEst.function <- function(x, outcome, score) {

  # Run estimation
  out <- x(outcome = outcome, score = score)

  # Check output
  check1 <- is.data.frame(out) &&
    identical(c("estimate", "outcome", "percentile", "score"), sort(colnames(out)))
  if (!check1) {
    stop(
      paste(
        "User defined estimation functions must return a data.frame of 4 columns:",
        "score - the predictions,",
        "percentile - the percentile of score,",
        "outcome - the original outcome,",
        "and estimate - the estimated value"
      )
    )
  }
  
  check2 <- vapply(out, is.numeric, FUN.VALUE = logical(1), USE.NAMES = TRUE)
  if (!all(check2)) {
    stop(
      paste(
        "All columns of the returned data.frame in user defined estimation function",
        "must be numeric.",
        paste0("`", names(which(!check2)), "`", collapse = ", "), "is/are not numeric."
      )
    )
  }

  return(out)
}

# Generic and S3 methods for estimation method as string (actual estimation function used)
getEstMethod <- function(x) {
  UseMethod("getEstMethod")
}

#' @export
getEstMethod.character <- function(x) {
  x
}

#' @export
getEstMethod.list <- function(x) {
  x[["method"]]
}

#' @export
getEstMethod.function <- function(x) {
  "udf"
}

getEstMethods <- function(x, with.names) {
  vapply(x, getEstMethod, FUN.VALUE = character(1), USE.NAMES = with.names)
}
