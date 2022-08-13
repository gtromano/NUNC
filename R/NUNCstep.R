#' NUNC Step
#' 
#' @description Performs a single step of the NUNC algorithm
#' 
#' @details This function tests a single window of data for a change in distribution using the 
#' NUNC Local algorithm. This is performed by comparing the data inside the sliding window for a 
#' change in distribution using a cost function based on the empirical CDFs for the data.
#' If the difference between the CDFs is significant enough, then a change is detected.
#' The comparison is performed at K quantiles of the ECDF, and the information aggregated to enhance
#' power.
#' 
#' This function allows a user to specify the quantiles of the distribution, and so incorporates
#' both NUNC Local and NUNC Semi-Parametric
#' 
#' The different methods require different arguments. These arguments can be specified
#' through the \code{params} argument:
#' \describe{
#' \item{if \code{method = "local"}}{An extra parameter to set the maximum amount of comparisons 
#'  to check can be specified via extra parameter \code{max_checks}. If not specified this 
#'  defaults to \code{w}, the size of the window.}
#' }
#'
#' @param win A window of data to test for a change in distribution.
#' @param threshold A positive numeric that gives the threshold for the test.
#' @param quantiles A vector of length K giving the quantiles at which the ECDFs are compared.
#' @param method The method used by NUNC.
#' @param params Additional arguments passed to NUNC - see details.
#'
#' @return Returns an object of class \code{NUNCout}, which is a \code{\link[=list]{list}}
#' of items containing
#' \describe{
#' \item{\code{$t}}{Time change is detected - indicated by the right edge of the window.}
#' \item{\code{$changepoint}}{Location of the changepoint in the data stream.}
#' \item{\code{$method}}{Method that has been specified by the user.}
#' \item{\code{$zVals}}{If the method is set to \code{"global"} then this returns the
#' value of the empirical CDF for the historic data at the K quantiles.} 
#' \item{\code{$data}}{The data inputted into the function.}
#' }
#' 
#' 
#' @seealso \link{NUNC} and \link{NUNCoffline}

#' 
#' @export
#'
#' @examples
#' 
#' window_data <- c(rnorm(100), rnorm(100, 3, 9))
#' NUNCstep(window_data, threshold = 3, method = "local")
#' 

NUNCstep <- function(win,
                     threshold,
                     quantiles = quantile(w, probs = seq(from = 1 / 15, by = 1 / 15, length = 14)),
                     method = c("global", "local"),
                     params = list()) {
  
  
  # setting up the method
  method <- match.arg(method)
  
  # checks on the data
  if(!is.numeric(win)) stop("Please provide a window of observations.")
  
  w <- length(win)
  
  if(!is.numeric(threshold) | threshold < 0)
    stop("threshold must be a positive numeric")
  
  if(!is.numeric(quantiles)){
    stop("quantiles must be a vector of numerics")
  }
  
  K <- length(quantiles)
  
  # first arguments checks
  if (w < 2) stop("The size of the window must be greater than 1.")
  
  # extra parameters check
  if (method == "local" && !("max_checks" %in% names(params))) {
    if (!is.null(names(params)))
      warning("the step for NUNC Local requires 'max_checks' argument: defaulting to w. Check documentation for more details.")
    params["max_checks"] <- length(win)
  } else if (method == "local" && params$max_checks > w) {
    stop("Checks to perform on local algorithm must be smaller than the window size.")
  } else if (method == "global" &&  !("past_info" %in% names(params))) {
    stop("the step for NUNC global requires a list of previous information. Check documentation for more details.")
  }
  
  if (method == "local") {
    J <- params$max_checks
    if (J == w)
      grid <- 1:(w-1)
    else {
      grid <- unique(.local_grid(w, beta, J))
      # if (J > w/2)
      #   grid <- unique(c(grid, .local_grid(w, beta, J, "left")))
    }
  }
  
  # running the actual function
  out <- switch (method, 
                 global = .stepGlobal(win, params$past_info, threshold, quantiles),
                 local = .NUNCLocalOffline(win, w, threshold, K, grid)
  )
  class(out) <- "NUNCout"
  out$data <- win
  out
}