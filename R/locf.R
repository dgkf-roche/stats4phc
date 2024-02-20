
#' zipFastener for two dataframes of unequal length
#'
#' The following function acts like a “zip fastener” for combining two dataframes.
#' It takes the first row of the first data frame and places it above of
#' the first row of the second data frame and so on.
#'
#' @param df1 dataframe 1
#' @param df2 dataframe 2
#'
#' @return Zipped data.frame
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#' df1 <- data.frame(a = 1:3, b = 1:3, c = 1:3)
#' df2 <- data.frame(a = letters[1:3], b = letters[1:3], c = letters[1:3])
#'
#' zipFastener(df1, df2)
#'
zipFastener <- function(df1, df2) {

  if (ncol(df1) != ncol(df2)) {
    stop("the no. of columns has to be equal to merge them by zip feeding")
  }

  d1 <- nrow(df1)
  d2 <- nrow(df2)

  if (d1 != d2) {
    stop("the no. of rows has to be equal to merge them by zip feeding")
  }

  # zip fastener preperations
  i1 <- 1:d1 # index vector 1
  i2 <- 1:d1 + d1 # index vector 2

  # zip fastener operations
  index <- as.vector(matrix(c(i1, i2), ncol = d1, byrow = TRUE))
  index <- index[!is.na(index)] # remove NAs

  colnames(df2) <- colnames(df1) # keep 1st colnames
  res <- rbind(df1, df2)[index, ] # reorder data frame

  return(res)
}


#' last observation carried forward: for binned
#'
#' Function assumes columns in data frame are: percentile, estimate, method
#'
#' @param df dataframe containing predcurve estimates
#'
#' @return dataframe with adjusted points
#'
#' @keywords internal
#' @noRd
#'
locf.binned <- function(df) {

  df.copy <- df
  # keep same x, slide y back one spot
  df.copy$estimate <- df.copy$estimate[c(2:length(df.copy$estimate), NA)]
  combined.df <- zipFastener(df, df.copy)
  combined.df <- combined.df[-nrow(combined.df), ] # removes a redundant point

  return(combined.df)
}

#' last observation carried forward: for pava
#'
#' Function assumes columns in data frame are: percentile, estimate, method
#'
#' @param df dataframe containing predcurve estimates
#'
#' @return dataframe with adjusted points
#'
#' @keywords internal
#' @noRd
#'
locf.pava <- function(df) {

  # save first and last points
  first <- df[1, ]
  last <- df[nrow(df) - 1, ]

  # adjust data to remove extra points
  diff_y <- diff(c(0, df$estimate))
  pava.int <- df[diff_y != 0, ]

  # re-add first and last points
  df <- rbind(first, pava.int, last)
  df <- df[order(df$percentile), ] # order on percentile

  df.copy <- df
  # keep same x, slide y back one spot
  df.copy$percentile <- df.copy$percentile[
    c(2:(length(df.copy$percentile)), NA)
  ]
  combined.df <- zipFastener(df, df.copy)

  # Combine smooth estimate to fixed df
  combined.df <- na.omit(combined.df)
  combined.df <- combined.df[-nrow(combined.df), ] # removes a redundant point

  return(combined.df)
}

#' Last observation carried forward for specific estimation methods
#'
#' Iterates locf for multiple methods: calls locf.binned and locf.pava
#'
#' @param df dataframe containing predcurve estimates
#' @param method.binned vector of "binned" methods
#' @param method.pava vector of "pava" methods
#'
#' @return dataframe with adjusted points
#'
#' @keywords internal
#' @noRd
#'
locf <- function(df, method.binned, method.pava) {

  complete.data <- df[!df$method %in% c(method.binned, method.pava), ]
  binned.data <- df[df$method %in% method.binned, ]
  pava.data <- df[df$method %in% method.pava, ]

  # pava locf
  if (nrow(pava.data) >= 1) {
    pava.plotdata <- locf.pava(pava.data)
    complete.data <- bind_rows(complete.data, pava.plotdata)
  }

  # binned locf
  if (nrow(binned.data) >= 1) {
    binned.plotdata <- locf.binned(binned.data)
    complete.data <- bind_rows(complete.data, binned.plotdata)
  }
  rownames(complete.data) <- NULL

  return(complete.data)
}
