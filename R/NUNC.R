#' NUNC Online
#' 
#' @description Online detection of changes in the distribution of streamed data.
#' 
#' @details This is the online version of the NUNC algorithm. For an offline version that
#' takes an dataset that has been observed in full, rather than a stream of data, see 
#' \link{NUNCoffline}.
#' 
#' Three different variants of NUNC exist: "NUNC Local", "NUNC Global", and "NUNC Semi-Parametric".
#' Each of these three variants tests for changes in distribution through use of a cost function
#' that makes a comparison between the pre and post change empirical CDFs for the data. 
#' This comparison is aggregated over K quantiles, as this enhances the power of the test by 
#' comparing both the centre, and tails, of the estimated distributions.
#' 
#' The three different methods can be described as follows.
#' \describe{
#' \item{NUNC Local}{This method searches for a change in distribution inside the points of 
#' data contained in the sliding window. An approximation for this algorithm can also be specified,
#' that only searches a subset of the points in the sliding window for a change in order to 
#' enhance computational efficiency.}
#' \item{NUNC Global}{This method tests if the data in the sliding window is drawn from a 
#' different distribution to the historic data.}
#' \item{NUNC Semi-Parametric}{This method is identical to NUNC Local, except the quantiles are
#' pre-specified rather than being estimated from the windowed data.}
#' }
#' 
#' 
#' The different methods require different arguments. These arguments can be specified
#' through the \code{params} argument:
#' \describe{
#' \item{if \code{method = "local"}}{An extra parameter to set the maximum amount of comparisons 
#'  to check can be specified via extra parameter \code{max_checks}. If not specified this 
#'  defaults to \code{w}, the size of the window.}
#' \item{if \code{method = "semi-param"}}{then parameter \code{quantiles} is required to perform 
#'  the semi-parametric estimation. This fixes the K quantiles the CDFs are compared at, rather than
#'  estimating them from the data.}
#' }
#'
#' @param dataFUN A function to pull the data. Needs to take no arguments and return
#' a single numeric value at each time point.
#' @param w  A positive integer for the size of the sliding window.
#' @param beta A positive numeric for the thresholds
#' @param K The positive integer for the number of quantiles to use.
#' @param method The NUNC method: one of "local", "global", or "semi-parametric".
#' @param params Additional arguments passed to NUNC - see details.
#' 
#'   
#' @return Returns an object of class \code{NUNCout}, which is a \code{\link[=list]{list}}
#' of items containing
#' \describe{
#' \item{\code{$t}}{Time change is detected - indicated by the right edge of the window.}
#' \item{\code{$changepoint}}{Location of the changepoint in the data stream.}
#' \item{\code{$zVals}}{If the method is set to \code{"global"} then this returns the
#' value of the empirical CDF for the historic data at the K quantiles.}
#' \item{\code{$method}}{Method that has been specified by the user.} 
#' }
#' 
#' 
#' @export
#'
#' @examples
#' set.seed(42)
#' onlinedata <- c(rnorm(300, 1, 10), rcauchy(1200, 12, .1))
#' 
#' f <- function() {
#'   out <- onlinedata[i]
#'   i <<- i + 1
#'   # uncomment and run with x11() or windows() to see real time plotting
#'   # plot(tail(onlinedata[1:i], 80), ylab = "test data", type = "l")
#'   out
#' }
#' 
#' i <- 1; NUNC(f, 50, beta = 10, K = 15, params = list(max_checks = 3), method = "local")
#' i <- 1; NUNC(f, 50, beta = 10, K = 15)



NUNC <- function(dataFUN,
           w,
           beta = find_beta(0.05, 1e3, w, K, method),
           K = 15,
           method = c("global", "local", "semi-param"),
           params = list()) {
  
  method <- match.arg(method)
  
  if( !is.numeric(w) | w <= 0 | round(w) != w)
    stop("w must be a positive integer")
  
  if( !is.numeric(beta) | beta < 0)
    stop("beta must be a positive numeric")
  
  
  
  # checks on the data generating function
  if(is.function(match.fun(dataFUN)))
    out_check <- dataFUN()
  else
    stop("Please provvide a data generating function.")
  if(!is.numeric(out_check)) stop("Data generating function provvided does not return a numeric output.")
  if(length(out_check) != 1) {
    warning("Length of the output from data generating function is greater than 1. Using only the first element.")
    f <- function() return(dataFUN()[1])
  } else {
    f <- match.fun(dataFUN)
  }
  
  
  # extra parameters check
  if (method == "local" && !("max_checks" %in% names(params))) {
    if (!is.null(names(params)))
      warning("NUNC Local requires 'max_checks' argument: defaulting to w. Check documentation for more details.")
    params["max_checks"] <- w
  } else if (method == "local" && params$max_checks > w) {
    stop("Checks to perform on local algorithm must be smaller than the window size.")
  }  else if (method == "semi-param" && !is.numeric(params$quantiles)) {
    stop("NUNC semi parametric needs a vector of quantiles of lenght K.")
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
                 global = .NUNCGlobal(f, w, beta, K),
                 local = .NUNCLocal(f, w, beta, K, grid),
                 "semi-param" = .NUNCsemiParam(f, w, beta, params$quantiles)
  )
  class(out) <- "NUNCout"
  out
}
