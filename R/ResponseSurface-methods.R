#' Method for plotting response surface objects
#'
#' @param x Output of \code{\link{fitSurface}}
#' @param color Character indicating on what values surface coloring will be
#'   based.
#'
#'   If \code{color = "z-score"}, surface coloring will be based on median of
#'   standartized off-axis Z-scores. Median function can be replaced by other
#'   function using an optional \code{colorfun} argument which will be passed to
#'   \code{plotResponseSurface}. Color breaks are determined here by standard
#'   deviation of off-axis Z-scores. For \code{color = "maxR"}, coloring will be
#'   based on values of maxR statistic and the quantile of its distribution
#'   (bootstrapped or not). If \code{color = "occupancy"}, coloring will be
#'   based on calculated occupancy rate for the respective dose combination.
#'
#' @param ... Further parameters passed to \code{\link{plotResponseSurface}}.
#'   \code{colorBy} argument in this method is computed automatically and thus
#'   cannot be passed to \code{\link{plotResponseSurface}}.
#' @export
plot.ResponseSurface <- function(x, color = c("z-score", "maxR", "occupancy"), ...) {

  color <- match.arg(color)
  inputs <- as.list(substitute(list(...)))[-1L]

  ## Green is synergy, red is antagonism
  if(!exists("colorPalette", inputs)) {
    inputs$colorPalette <- c("red", rep("grey70", 2), "blue")
    if (x$fitResult$coef["b"] > x$fitResult$coef["m1"]) {
      inputs$colorPalette <- rev(inputs$colorPalette)
    }
  }

  if (color == "z-score") {
    boundary <- sd(x$offAxisTable[["z.score"]])
    inputs$colorBy <- x$offAxisTable[, c("d1", "d2", "z.score")]
    if (!exists("breaks", inputs)) inputs$breaks <- c(-Inf, -boundary, 0, boundary, Inf)
    if (!exists("main", inputs)) inputs$main <- "Z-scores"
  } else if (color == "maxR") {
    inputs$colorBy <- x$maxR$Ymean[, c("d1", "d2", "R")]
    q <- attr(x$maxR$Ymean, "q")
    if (!exists("breaks", inputs)) inputs$breaks <- c(-Inf, -q, 0, q, Inf)
    if (!exists("main", inputs)) inputs$main <- "maxR"
  } else if (color == "occupancy") {
    inputs$colorBy <- x$occupancy
    inputs$colorPalette <- c("#EFF3FF", "#BDD7E7", "#6BAED6", "#2171B5")
    if (!exists("breaks", inputs)) inputs$breaks <- c(0, 0.25, 0.5, 0.75, 1)
    if (!exists("main", inputs)) inputs$main <- "Occupancy rate"
  }

  inputs$data <- x$data
  inputs$fitResult <- x$fitResult
  inputs$transforms <- x$transforms
  inputs$null_model <- x$null_model

  do.call(plotResponseSurface, inputs)

}

#' Method for plotting of contours based on maxR statistics
#'
#' @param x Output of \code{\link{fitSurface}}
#' @param ... Further parameters passed to \code{\link{plot.maxR}}
#' @export
contour.ResponseSurface <- function(x, ...) {

  if (!exists("maxR", x))
    stop("maxR statistics were not found.")

  plot(x$maxR, ...)

}

#' Method for contour plots with isoboles based on null model predictions
#'
#' @param x Output of \code{\link{fitSurface}}
#' @param grid.len Number of concentrations to plot for each compound in the
#'   contour plot. An evenly spaced grid of doses will be generated for each
#'   compound given its respective observed minimum and maximum doses. Note that
#'   \code{grid.len^2} computations will be needed later so this number should
#'   stay reasonably low.
#' @param logScale If \code{logScale = TRUE}, then grid of doses is evenly
#'   spaced in the logarithmic scale.
#' @param ... Further parameters passed to ...
#' @import ggplot2
#' @export
image.ResponseSurface <- function(x, grid.len = 100, logScale = TRUE, ...) {

  ## Generate evenly spaced grid either on a linear or a log-scale
  genSeq <- function(doses) {

    if (logScale) {
      ## Log-scale removed zero dose
      doses <- setdiff(doses, 0)
      seq.range <- log(range(doses))
      c(0, exp(seq(seq.range[1], seq.range[2], length.out = grid.len - 1)))
    } else {
      ## Linear scale
      seq.range <- range(doses)
      seq(seq.range[1], seq.range[2], length.out = grid.len)
    }

  }

  coefs <- coef(x$fitResult)

  ## Generate a grid of doses for Compound 1 and predict the response
  doses1 <- genSeq(unique(x$data$d1))
  resp1 <- L4(doses1, coefs["h1"], coefs["b"], coefs["m1"], coefs["e1"])
  ## Generate a grid of doses for Compound 2 and predict the response
  doses2 <- genSeq(unique(x$data$d2))
  resp2 <- L4(doses2, coefs["h2"], coefs["b"], coefs["m2"], coefs["e2"])
  ## Combine both compounds and their marginal predictions
  data <- rbind(data.frame("d1" = doses1, "d2" = 0, "effect" = resp1),
                data.frame("d1" = 0, "d2" = doses2, "effect" = resp2))

  ## Based on marginal data, generate null model predictions
  predSurface <- predictOffAxis(data, x$fitResult,
                                null_model = x$null_model)$predSurface

  melt.surface <- data.frame("d1" = rep(doses1, length(doses2)),
                             "d2" = rep(doses2, each = length(doses1)),
                             "effect" = as.numeric(predSurface))

  labnames <- c("Response", "Compound 1", "Compound 2")
  if (!is.null(attr(x$data, "orig.colnames"))) {
    labnames <- unlist(attr(x$data, "orig.colnames"))
  }

  p <- ggplot(melt.surface,
              aes_string(x = "d1", y = "d2",
                         z = "effect", fill = "effect")) +
    theme_bw() +
    geom_tile() +
    labs(x = labnames[2], y = labnames[3]) +
    scale_fill_gradientn(labnames[1],
                         colours = c("steelblue", "lightsteelblue", "lightblue",
                                     "floralwhite", "beige", "khaki",
                                     "orange1", "tomato3", "red")) +
    geom_contour(bins = 7, col = "black", size = 0.2)

  if (logScale) {
    p <- p +
      scale_x_log10(breaks = unique(x$data$d1)) +
      scale_y_log10(breaks = unique(x$data$d1))
  } else {
    p <- p +
      scale_x_continuous(breaks = unique(x$data$d1)) +
      scale_y_continuous(breaks = unique(x$data$d1))
  }

  p

}


#' Summary of \code{ResponseSurface} object
#'
#' @param object Output of \code{\link{fitSurface}}
#' @param ... Further parameters
#' @export
summary.ResponseSurface <- function(object, ...) {

  ans <- list()
  ans$marginalFit <- summary(object$fitResult)
  ans$null_model <- object$null_model
  ans$shared_asymptote <- object$fitResult$shared_asymptote

  if (!is.null(object$meanR)) ans$meanR <- summary(object$meanR)
  if (!is.null(object$maxR)) ans$maxR <- summary(object$maxR)

  ans$occup <- mean(object$occupancy$occupancy)

  class(ans) <- "summary.ResponseSurface"
  ans
}

#' Print method for the summary function of \code{ResponseSurface} object
#'
#' @param x Summary of \code{ResponseSurface} object
#' @param ... Further parameters
#' @export
print.summary.ResponseSurface <- function(x, ...) {

  cat("Null model: ")
  if (x$null_model == "loewe" & x$shared_asymptote == TRUE)
    cat("Standard Loewe Additivity")
  else if (x$null_model == "loewe" & x$shared_asymptote == FALSE)
    cat("Generalized Loewe Additivity")
  else if (x$null_model == "hsa" & x$shared_asymptote == TRUE)
    cat("Highest Single Agent with shared maximal response")
  else if (x$null_model == "hsa" & x$shared_asymptote == FALSE)
    cat("Highest Single Agent with differing maximal response")
  else
    cat(x$null_model)

  cat("\n")
  cat("Mean occupancy rate:", x$occup)
  cat("\n\n")
  print(x$marginalFit)
  cat("\n")

  if (!is.null(x$meanR)) print(x$meanR)
  if (!is.null(x$maxR)) print(x$maxR)

  if (is.null(x$meanR) & is.null(x$maxR)) {
    cat("\n\n")
    cat("No test statistics were computed.")
  }

  cat("\n")
}


#' Predicted values of the response surface according to the given null model
#'
#' @param object Output of \code{\link{fitSurface}}
#' @param ... Further parameters
#' @export
fitted.ResponseSurface <- function(object, ...) {

  doseInput <- object$data[, c("d1", "d2")]
  parmInput <- coef(object$fitResult)

  switch(object$null_model,
         "loewe" = generalizedLoewe(doseInput, parmInput)$response,
         "hsa" = hsa(doseInput, parmInput))

}
