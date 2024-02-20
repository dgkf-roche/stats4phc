#' Calibration plot
#'
#' Calibration curve risk estimates
#'
#' @inheritParams riskProfile
#' @param methods Character vector of method names (case-insensitive) for plotting curves or
#' a named list where elements are method function and its arguments.
#' Default is set to `list(gam = list(method = "gam", fitonPerc = FALSE))`.
#'
#' Full options are: `c("binned", "pava", "mspline", "gam", "cgam")`.
#'
#' To specify arguments per method, use lists. For example:
#' ```
#' list(
#'   pava = list(method = "pava", ties = "primary"),
#'   mspline = list(method = "mspline", fitonPerc = TRUE),
#'   gam = list(method = "gam", bs = "tp", logscores = FALSE),
#'   bin = list(method = "binned", bins = 10),
#' )
#' ```
#' See section "Estimation" for more details.
#' @param include Character vector (case-insensitive, partial matching) or `NULL` specifying
#' what quantities to include in the plot.
#'
#' Default is: `c("loess", "citl")`.
#'
#' Full options are: `c("loess", "citl", "rug", "datapoints")` or `NULL`.
#' "loess" adds a Loess fit, "citl" stands for "Calibration in the large",
#' "rug" adds rug ticks of `score` by `outcome` (top x-axis: `score` for `outcome == 1`,
#' bottom x-axis: `score` for `outcome == 0`),
#' "datapoints" adds jittered `score` by `outcome` (slightly shifted away from 0 / 1 y-values),
#' "`NULL`" stands for no extra information.
#' @param plot.raw Logical to show percentiles or raw values.
#' Defaults to `TRUE` (i.e. raw `score`).
#' @param rev.order Logical to reverse ordering of scores. Defaults to `FALSE`.
#' @param margin.type Type of additional margin plot, can be one of
#' `c("density", "histogram", "boxplot", "violin", "densigram")`.
#' See [ggExtra::ggMarginal()] for more details.
#' @param ... Additional arguments passed to [ggExtra::ggMarginal()].
#'
#' @inheritSection riskProfile Estimation
#'
#' @return A list containing the plot and data, plus `citl` data if they were requested.
#'
#' @export
#'
#' @seealso [riskProfile()] [sensSpec()]
#'
#' [getPAVAest()] [getBINNEDest()] [getGAMest()] [getCGAMest()] [getMSPLINEest()]
#' [getASISest()]
#'
#' @examples
#' # Read in example data
#' auroc <- read.csv(system.file("extdata", "sample.csv", package = "stats4phc"))
#' rscore <- auroc$predicted_calibrated
#' truth <- as.numeric(auroc$actual)
#'
#' # Default calibration plot
#' p1 <- calibrationProfile(outcome = truth, score = rscore)
#' p1$plot
#'
#' # Specifying multiple estimation methods
#' # By default, all the methods fit on percentiles
#' calibrationProfile(
#'   outcome = truth,
#'   score = rscore,
#'   methods = c("gam", "mspline", "binned")
#' )$plot
#'
#' # Specifying multiple estimation methods with parameters
#' calibrationProfile(
#'   outcome = truth,
#'   score = rscore,
#'   methods = list(
#'     gam = list(method = "gam", fitonPerc = FALSE, k = 3),
#'     mspline = list(method = "mspline"),
#'     bin = list(method = "binned", quantiles = 5)
#'   )
#' )$plot
#'
#' # Additional quantities and marginal histogram with specified number of bins
#' calibrationProfile(
#'   outcome = truth,
#'   score = rscore,
#'   include = c("rug", "datapoints", "citl"),
#'   # or use partial matching: include = c("r", "d", "c"),
#'   margin.type = "histogram",
#'   bins = 100 # passed to ggExtra::ggMarginal
#' )$plot
#'
calibrationProfile <- function(outcome,
                               score,
                               methods = list(gam = list(method = "gam", fitonPerc = FALSE)),
                               include = c("loess", "citl"),
                               plot.raw = TRUE,
                               rev.order = FALSE,
                               margin.type = NULL,
                               ...) {

  # Argument checks (except outcome, score, methods - below)
  checkmate::assert(
    checkmate::check_character(include),
    checkmate::check_null(include)
  )
  if (is.character(include)) {
    include <- matchArgSubset(tolower(include), choices = c("loess", "citl", "rug", "datapoints"))
  }
  checkmate::assert_flag(plot.raw)
  checkmate::assert_flag(rev.order)
  checkmate::assert(
    checkmate::check_character(margin.type, len = 1, any.missing = FALSE),
    checkmate::check_null(margin.type)
  )

  if (plot.raw) {
    xvar <- "score"
  } else {
    xvar <- "percentile"
  }

  # Standardize/Check outcome, scores
  op <- inputCheck(outcome = outcome, score = score)

  # Order Data by scores
  tempdf <- orderInputs(outcome = op$outcome, score = op$score, rev.order = rev.order)
  score <- tempdf$score
  outcome <- tempdf$outcome

  # Check methods
  methods <- methodCheck(methods = methods)
  method.names <- names(methods)
  if ("asis" %in% getEstMethods(methods, with.names = FALSE)) {
    stop('"asis" method is not suitable for this plot. Please remove it.')
  }

  # Get estimates
  pc.ests <- getEsts(methods = methods, outcome = outcome, score = score)
  PC.data <- pc.ests$plotdata
  step.methods <- method.names[pc.ests$idx.step]

  # Calculate percentiles
  ecdf.score <- ecdf(score)
  percentile <- ecdf.score(score)

  # Calculate Calibration in the large
  citl.data <- data.frame(
    outcome = mean(outcome),
    score = mean(score),
    percentile = ecdf.score(mean(score)),
    method = "Calibration In The Large"
  )

  # Subset PC data for plotting
  smoothPC <- PC.data[!PC.data$method %in% step.methods, , drop = FALSE]
  stepPC <- PC.data[PC.data$method %in% step.methods, , drop = FALSE]

  # data.frame with user inputs
  ddf <- data.frame(score, percentile, outcome)

  # Shape type storage
  shapes <- c()

  # Add empty scatterplot layer for ggMarginal
  if (!is.null(margin.type)) {
    p <- ggplot() +
      geom_point(
        data = ddf,
        aes(x = .data[[xvar]], y = .data$outcome), shape = NA, na.rm = TRUE
      )
  } else {
    p <- ggplot()
  }

  # Build plot
  p <- p +
    geom_line(
      data = smoothPC,
      aes(x = .data[[xvar]], y = .data$estimate, linetype = .data$method, colour = .data$method),
      alpha = 0.8,
      linewidth = 0.5
    ) +
    geom_step(
      data = stepPC,
      aes(x = .data[[xvar]], y = .data$estimate, linetype = .data$method, colour = .data$method),
      direction = "vh",
      alpha = 0.8,
      linewidth = 0.5
    ) +
    geom_abline(
      aes(slope = 1, intercept = 0, linewidth = "Identity line"),
      colour = "gray50",
      linetype = "solid"
    )

  # Add loess if requested
  if ("loess" %in% include) {
    p <- p + geom_smooth(
      data = ddf,
      aes(x = .data[[xvar]], y = .data$outcome, linetype = "loess", colour = "loess"),
      method = "loess", formula = y ~ x, se = FALSE,
      linewidth = 0.5
    )
  }

  # Add calibration in the large if requested
  if ("citl" %in% include) {
    p <- p + geom_point(
      data = citl.data,
      aes(x = .data[[xvar]], y = .data$outcome, shape = .data$method),
      colour = "red",
      size = 3,
      stroke = 1
    )
    shapes <- c(shapes, c("Calibration In The Large" = 4))
  }

  # Add datapoints if requested
  if ("datapoints" %in% include) {
    p <- p + geom_jitter(
      data = data.frame(
        score,
        percentile,
        outcome = ifelse(outcome == 0, -0.1, 1.1),
        method = "Data points"
      ),
      aes(x = .data[[xvar]], y = .data$outcome, shape = .data$method),
      colour = "black",
      size = 1.5,
      alpha = 0.4,
      position = position_jitter(seed = 5, height = 0.03)
    )
    shapes <- c(shapes, c("Data points" = 16))
  }

  # Add rug if requested
  if ("rug" %in% include) {
    p <- p + geom_rug(
      data = ddf[ddf$outcome == 0, ],
      aes(x = .data[[xvar]]),
      sides = "b",
      show.legend = FALSE
    ) + geom_rug(
      data = ddf[ddf$outcome == 1, ],
      aes(x = .data[[xvar]]),
      sides = "t",
      show.legend = FALSE
    )
  }

  # Fix shape legend ...
  if (all(c("citl", "datapoints") %in% include)) {
    shape_guide <- guide_legend(
      override.aes = list(
        colour = c("Calibration In The Large" = "red", "Data points" = "black"),
        alpha = 1, size = 2.5
      ),
      order = 3
    )
  } else if (any(c("citl", "datapoints") %in% include)) {
    shape_guide <- guide_legend(
      override.aes = list(alpha = 1, size = 2.5),
      order = 3
    )
  } else {
    shape_guide <- NULL
  }

  # Finalize graph
  p <- p +
    scale_linewidth_manual(values = c("Identity line" = 0.5)) +
    scale_shape_manual(values = shapes) +
    scale_colour_hue(l = 45) +
    scale_x_continuous(n.breaks = 6) +
    scale_y_continuous(n.breaks = 6) +
    labs(
      title = "Calibration Plot",
      x = "Predicted Probability",
      y = "Observed",
      linetype = "Estimation Method",
      linewidth = NULL,
      colour = "Estimation Method",
      shape = `if`(all(c("citl", "datapoints") %in% include), "Points", NULL)
    ) +
    theme_bw() +
    theme(legend.key.width = unit(2, "line")) +
    guides(
      linetype = guide_legend(order = 1),
      colour = guide_legend(order = 1),
      linewidth = guide_legend(order = 2),
      shape = shape_guide
    )

  # Add margin plot
  if (!is.null(margin.type)) {
    p <- ggExtra::ggMarginal(p, type = margin.type, margins = "x", ...)
  }

  return(list(plot = p, data = dplyr::as_tibble(PC.data), citl = citl.data))
}
