#' Risk profile plot
#'
#' Predictiveness curve, PPV, NPV and 1-NPV risk estimates
#'
#' @param outcome Vector of binary outcome for each observation.
#' @param score Numeric vector of continuous predicted risk score.
#' @param methods Character vector of method names (case-insensitive) for plotting curves or
#' a named list where elements are method function and its arguments.
#' Default is set to `"asis"`.
#'
#' Full options are: `c("asis", "binned", "pava", "mspline", "gam", "cgam")`.
#'
#' To specify arguments per method, use lists. For example:
#' ```
#' list(
#'   pava = list(method = "pava", ties = "primary"),
#'   mspline = list(method = "mspline", fitonPerc = TRUE),
#'   gam = list(method = "gam", bs = "tp", logscores = FALSE),
#'   bin = list(method = "binned", bins = 10),
#'   risk = list(method = "asis")
#' )
#' ```
#' See section "Estimation" for more details.
#' @param prev.adj `NULL` (default) or scalar numeric between 0 and 1 for prevalence adjustment.
#' @param show.prev Logical, show prevalence value in the graph. Defaults to `TRUE`.
#' @param show.nonparam.pv Logical, show non-parametric calculation of PVs. Defaults to `TRUE`.
#' @param show.best.pv Logical, show best possible PVs. Defaults to `TRUE`.
#' @param include Character vector (case-insensitive, partial matching) specifying what quantities
#' to include in the plot.
#'
#' Default is: `c("PC", "PPV", "1-NPV")`.
#'
#' Full options are: `c("NPV", "PC", "PPV", "1-NPV")`.
#' @param plot.raw Logical to show percentiles or raw values.
#' Defaults to `FALSE` (i.e. percentiles).
#' @param rev.order Logical, reverse ordering of scores. Defaults to `FALSE`.
#'
#' @section Estimation:
#' The `methods` argument specifies the estimation method.
#' You can provide either a vector of strings, any of
#' ```
#' c("asis", "binned", "pava", "mspline", "gam", "cgam")
#' ```
#' (`"asis"` is not available for `calibrationProfile`),
#' or a named list of lists.
#' In the latter case, the inner list must have an element "method",
#' which specifies the estimation function (one of those above),
#' and optionally other elements, which are passed to the estimation function.
#' For example:
#' ```
#' list(
#'   gam = list(method = "gam", k = 3),
#'   c_gam = list(method = "cgam", numknots = 3)
#' )
#' ```
#'
#' To see what arguments are available for each estimation method,
#' see the documentation of that function.
#' The naming convention is `getXest`,
#' where `X` stands for the estimation method, for example [getGAMest()].
#'
#' "gam", "cgam", and "mspline" always fit on percentiles by default.
#' To change this, use `fitonPerc = FALSE`, for example 
#' ```
#' list(gam = list(method = "gam", fitonPerc = FALSE))
#' ```
#'
#' "gam" and "cgam" methods are wrappers of [mgcv::gam()] and [cgam::cgam()], respectively.
#' The default values of function arguments (like `k`, the number of knots in [mgcv::s()])
#' mirror the package defaults.
#'
#' @return A list containing the plot and data, plus `errorbar` data if they were requested 
#' (through `"binned"` estimation method with a parameter `errorbar.sem`).
#'
#' @export
#'
#' @seealso [calibrationProfile()] [sensSpec()]
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
#' # Default plot includes 1-NPV, PPV, and a predictiveness curve (PC) based on risk-cutoff
#' p1 <- riskProfile(outcome = truth, score = rscore)
#' p1$plot
#' p1$data
#'
#' # Show also NPV
#' p2 <- riskProfile(
#'   outcome = truth,
#'   score = rscore,
#'   include = c("PC", "NPV", "PPV", "1-NPV")
#'   # or use partial matching: include = c("PC", "N", "PPV", "1")
#' )
#' p2$plot
#' p2$data
#'
#' # All estimates of prediction curve
#' p3 <- riskProfile(
#'   outcome = truth,
#'   score = rscore,
#'   methods = c("mspline", "gam", "cgam", "binned", "pava", "asis"),
#'   include = c("PC", "PPV", "1-NPV")
#' )
#' p3$plot
#'
#' # Specifying method arguments (note each list has a "method" element)
#' p4 <- riskProfile(
#'   outcome = truth,
#'   score = rscore,
#'   methods = list(
#'     "gam" = list(method = "gam", bs = "tp", logscores = FALSE, fitonPerc = TRUE),
#'     "risk" = list(method = "asis"), # no available arguments for this method
#'     "bin" = list(method = "binned", quantiles = 10, errorbar.sem = 1.2)
#'   )
#' )
#' p4$plot
#'
#' # Compare multiple GAMs in terms of Predictiveness Curves
#' p5 <- riskProfile(
#'   outcome = truth,
#'   score = rscore,
#'   methods = list(
#'     "gam_3" = list(method = "gam", k = 3),
#'     "gam_4" = list(method = "gam", k = 4),
#'     "gam_7" = list(method = "gam", k = 7)
#'   ),
#'   include = "PC"
#' )
#' p5$plot
#'
#' # Using logistic regression as user-defined estimation function, fitting on percentiles
#' # Function needs to take exactly these two arguments
#' my_est <- function(outcome, score) {
#'   # Calculate percentiles
#'   perc <- ecdf(score)(score)
#'   # Fit
#'   m <- glm(outcome ~ perc, family = "binomial")
#'   # Generate predictions
#'   preds <- predict(m, type = "response")
#'   # Return a data.frame with exactly these columns
#'   return(
#'     data.frame(
#'       score = score,
#'       percentile = perc,
#'       outcome = outcome,
#'       estimate = preds
#'     )
#'   )
#' }
#' p6 <- riskProfile(
#'   outcome = truth,
#'   score = rscore,
#'   methods = list(my_lr = my_est)
#' )
#' p6$plot
#'
#' # Using cgam as user-defined estimation function
#' # Note that you can also use the predefined cgam using methods = "cgam"
#' # Attach needed library
#' # Watch out for masking of mgcv::s and cgam::s if both are attached
#' library(cgam, quietly = TRUE) 
#' # Function needs to take exactly these two arguments
#' my_est <- function(outcome, score) {
#'   # Fit on raw predictions with space = "E"
#'   m <- cgam(
#'     outcome ~ s.incr(score, numknots = 5, space = "E"),
#'     family = "binomial"
#'   )
#'   # Generate predictions and convert to vector
#'   preds <- predict(m, type = "response")$fit
#'   # Return a data.frame with exactly these columns
#'   out <- data.frame(
#'     score = score,
#'     percentile = ecdf(score)(score),
#'     outcome = outcome,
#'     estimate = preds
#'   )
#'   return(out)
#' }
#'
#' p7 <- riskProfile(
#'   outcome = truth,
#'   score = rscore,
#'   methods = list(my_cgam = my_est)
#' )
#' p7$plot
#' 
#' # Prevalence adjustment to 0.1
#' p8 <- riskProfile(outcome = truth, score = rscore, prev.adj = 0.1)
#' p8$plot
#'
riskProfile <- function(outcome,
                        score,
                        methods = "asis",
                        prev.adj = NULL,
                        show.prev = TRUE,
                        show.nonparam.pv = TRUE,
                        show.best.pv = TRUE,
                        include = c("PC", "PPV", "1-NPV"),
                        plot.raw = FALSE,
                        rev.order = FALSE) {

  # Argument checks (except outcome, score, methods - below)
  checkmate::assert_number(prev.adj, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_flag(show.nonparam.pv)
  checkmate::assert_flag(show.best.pv)
  checkmate::assert_flag(show.prev)
  include <- matchArgSubset(toupper(include), choices = c("PC", "PPV", "NPV", "1-NPV"))
  checkmate::assert_flag(plot.raw)
  checkmate::assert_flag(rev.order)

  # Standardize/Check outcome, scores
  op <- inputCheck(outcome = outcome, score = score)

  # Order Data by scores
  tempdf <- orderInputs(outcome = op$outcome, score = op$score, rev.order = rev.order)
  score <- tempdf$score
  outcome <- tempdf$outcome

  # Check methods
  methods <- methodCheck(methods = methods)
  method.names <- names(methods)

  # Calculate prevalence and percentiles
  prev <- mean(outcome, na.rm = TRUE)

  # Get the plot settings
  show.pc <- "PC" %in% include
  show.pv <- any(c("PPV", "NPV", "1-NPV") %in% include)
  show.one.only <- sum(c("PC", "PPV", "NPV", "1-NPV") %in% include) == 1

  if (plot.raw) {
    xvar <- "score"
  } else {
    xvar <- "percentile"
  }

  # Prediction Curve Data
  if (show.pc) {
    pc.ests <- getEsts(methods = methods, outcome = outcome, score = score)
    PC.data <- mutate(pc.ests$plotdata, pv = "PC")
    errorbar.data <- pc.ests$errorbardata
    step.methods <- method.names[pc.ests$idx.step]
  } else {
    PC.data <- data.frame(
      method = character(0), pv = character(0),
      percentile = numeric(0), score = numeric(0), estimate = numeric(0)
    )
    pc.ests <- errorbar.data <- NULL
    step.methods <- character(0)
  }

  # Predictive Value Data
  if (show.pv) {
    PV.data <- getPVdata(outcome = outcome, score = score, methods = methods, pc.ests = pc.ests)
  } else {
    PV.data <- data.frame(
      method = character(0), score = numeric(0), percentile = numeric(0),
      outcome = numeric(0), estimate = numeric(0),
      MNPV = numeric(0), NPV = numeric(0), PPV = numeric(0)
    )
  }

  # show.nonparam.pv
  if (show.nonparam.pv) {
    tmp <- nonParametricPV(outcome = outcome, score = score) %>%
      mutate(method = "non-parametric", estimate = NA)
    PV.data <- bind_rows(PV.data, tmp)
  }

  # Dataset of inputs
  df.in <- data.frame(outcome, score, percentile = ecdf(score)(score))

  # Adjust based on user defined prevalence
  if (!is.null(prev.adj)) {

    cdf.cases <- ecdf(score[outcome == 1])
    cdf.controls <- ecdf(score[outcome == 0])

    df.in$percentile <- adjPrevPerc(
      perc = df.in$score, prev.new = prev.adj,
      cdf.case = cdf.cases, cdf.control = cdf.controls
    )

    if (show.pc) {
      PC.data <- adjPrevPC(
        dat = PC.data, prev = prev, prev.new = prev.adj,
        cdf.case = cdf.cases, cdf.control = cdf.controls
      )
    }

    if (show.pv) {
      PV.data <- adjPrevPV(
        dat = PV.data, prev = prev, prev.new = prev.adj,
        cdf.case = cdf.cases, cdf.control = cdf.controls
      )
    }

    prev <- prev.adj
  }

  # pivot PV data
  if (show.pv) {
    PV.data <- PV.data %>%
      rename(`1-NPV` = "MNPV") %>%
      tidyr::pivot_longer(
        cols = all_of(c("1-NPV", "NPV", "PPV")), names_to = "pv", values_to = "pvValue"
      ) %>%
      mutate(pv = factor(.data$pv, levels = c("NPV", "1-NPV", "PPV"))) %>%
      arrange(.data$method, .data$pv, .data$percentile)
  } else {
    PV.data <- data.frame(
      method = character(0), score = numeric(0), percentile = numeric(0),
      outcome = numeric(0), estimate = numeric(0),
      pv = character(0), pvValue = numeric(0)
    )
  }

  # If showing percentiles, add row with 0th percentile
  if (!plot.raw) {
    if (show.pc) {
      PC.data <- add0thPercPC(PC.data)
    }
    if (show.pv) {
      PV.data <- add0thPercPV(PV.data)
    }
  }

  # Subset PC data
  smoothPC <- PC.data[!PC.data$method %in% step.methods, , drop = FALSE]
  stepPC <- PC.data[PC.data$method %in% step.methods, , drop = FALSE]

  # Subset PV data
  smoothPV <- PV.data[
    PV.data$pv %in% include & !PV.data$method %in% step.methods, ,
    drop = FALSE
  ]
  stepPV <- PV.data[
    PV.data$pv %in% include & PV.data$method %in% step.methods, ,
    drop = FALSE
  ]

  # Different aes based on what is to be shown
  # (if one kind of PV value, use both coloring and linetype for distinguishing estimation methods)
  if (show.one.only) {
    aes.pc <- aes(
      x = .data[[xvar]], y = .data$estimate,
      colour = .data$method, linetype = .data$method
    )
    aes.pv <- aes(
      x = .data[[xvar]], y = .data$pvValue,
      colour = .data$method, linetype = .data$method
    )
  } else {
    aes.pc <- aes(
      x = .data[[xvar]], y = .data$estimate,
      colour = .data$pv, linetype = .data$method
    )
    aes.pv <- aes(
      x = .data[[xvar]], y = .data$pvValue,
      colour = .data$pv, linetype = .data$method
    )
  }

  # Build plot
  p <- ggplot() +
    geom_line(aes.pc, data = smoothPC, alpha = 0.8) +
    geom_step(aes.pc, data = stepPC, alpha = 0.8, direction = "vh") +
    geom_line(aes.pv, data = smoothPV, alpha = 0.8) +
    geom_step(aes.pv, data = stepPV, alpha = 0.8, direction = "vh") +
    geom_hline(yintercept = prev, alpha = 0.8, col = "black", linetype = "dashed")

  # Best PVs
  if (show.best.pv) {
    if (!plot.raw) {
      df.in <- bind_rows(dplyr::tibble(percentile = 0, score = NA), df.in)
    }
    best <- df.in %>%
      select(all_of(c("percentile", "score"))) %>%
      distinct() %>%
      arrange(.data$percentile, .data$score) %>%
      mutate(
        PPV = bestPPV(perc = .data$percentile, prev = prev),
        `1-NPV` = bestMNPV(perc = .data$percentile, prev = prev),
        NPV = 1 - .data$`1-NPV`
      ) %>%
      tidyr::pivot_longer(
        cols = c("PPV", "1-NPV", "NPV"), names_to = "pv", values_to = "pvValue"
      ) %>%
      filter(.data$pv %in% include) %>%
      mutate(
        method = "Best PVs",
        pv = paste("Best", .data$pv)
      )

    if (show.one.only) {
      p <- p +
        geom_line(
          data = best,
          aes(x = .data[[xvar]], y = .data$pvValue, linewidth = .data$pv),
          colour = "gray70"
        ) +
        scale_linewidth_manual(
          values = c("Best 1-NPV" = 0.3, "Best PPV" = 0.3, "Best NPV" = 0.3),
          name = "Best PVs"
        )
    } else {
      p <- p +
        geom_line(data = best, aes(x = .data[[xvar]], y = .data$pvValue, colour = .data$pv))
    }
  }

  # Add errorbars
  if (show.pc && !is.null(errorbar.data)) {
    p <- p +
      geom_point(
        data = errorbar.data,
        aes(x = .data$midquantile, y = .data$bin.mid),
        alpha = 0.8,
        size = 0.2
      ) +
      geom_errorbar(
        data = errorbar.data,
        aes(
          x = .data$midquantile,
          ymin = .data$bin.low,
          ymax = .data$bin.high,
          width = .02
        ),
        alpha = 0.7,
        linewidth = 0.2,
        inherit.aes = FALSE
      )
  }

  # Set always the same colours for PVs
  if (!show.one.only) {
    clrs <- predictionColours(include, show.best = show.best.pv)
    p <- p + scale_colour_manual(values = clrs, breaks = names(clrs))
  } else {
    p <- p + scale_colour_hue(l = 45)
  }

  # Finalize plot
  p <- p +
    labs(
      title = "Predictiveness Plot",
      x = ifelse(plot.raw, "Prediction Score", "Risk Percentile"),
      y = "Predicted Risk / Predictive Value",
      linetype = "Estimation Method",
      colour = ifelse(show.one.only, "Estimation Method", "Predictive Quantity")
    ) +
    scale_x_continuous(n.breaks = 6) +
    scale_y_continuous(n.breaks = 6) +
    theme_bw() +
    theme(legend.key.width = unit(2, "line"))
  
  # Add prevalence annotation if requested
  if (show.prev) {
    # x-value for plotting prevalence label
    prev_x <- ifelse(plot.raw, min(score), 0)
    prev_nudge_x <- ifelse(plot.raw, (max(score) - min(score)) / 10, 0.1)
    prev_nudge_y <- ggplot2::layer_scales(p)$y$get_limits()[2] / 10
    
    p <- p + annotate(
      geom = "text",
      x = prev_x + prev_nudge_x,
      y = ifelse(prev > 0.8, prev - prev_nudge_y, prev + prev_nudge_y),
      label = paste0("Prevalence: ", "\n", round(prev, 3)),
      colour = "black",
      alpha = 0.8,
      size = 3.5
    )
  }

  if (show.one.only) {
    p <- p + guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 1))
  } else {
    p <- p + guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2))
  }

  if (show.best.pv) {
    PV.data <- bind_rows(PV.data, mutate(best, pv = gsub("Best ", "", .data$pv)))
  }

  return(
    list(
      plot = p,
      data = bind_rows(
        dplyr::as_tibble(PC.data),
        dplyr::as_tibble(PV.data)
      ),
      errorbar = errorbar.data
    )
  )
}
