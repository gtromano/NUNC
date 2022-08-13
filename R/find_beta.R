#' Find Beta
#' 
#' @description Find the threshold value for a given probability of false alarm.
#' 
#' @details This function provides a threshold that offers a theoretical bound on the probability
#' of false alarm by time t at some level alpha. In some cases, however, this bound may be 
#' conservative and cause a loss of test power.
#'
#' @param alpha A numeric between zero and one specifying the probability of false alarm.
#' @param t The time by which there will be this probability of a false alarm.
#' @param w A positive integer for the window size used by NUNC.
#' @param K A positive integer for the number of quantiles used by NUNC
#' @param method The variant of NUNC used, one of "global", "local", or "semi-parametric".
#'
#' @return A numeric for the penalty value
#' 
#' @seealso \link{NUNC} and \link{NUNCoffline}
#' 
#' @export
#'
#' @examples
#' 
#' #Return the threshold corresponding to a 5% probability of false alarm by time 1000 for NUNC
#' #local using a window of size 50 and 20 quantiles.
#' 
#' beta <- find_beta(alpha = 0.05, t = 1000, w = 50, K = 20, method = "local")
#' 

find_beta <- function(alpha, t, w, K, method = c("global", "local", "semi-param")) {
  
  if( !is.numeric(alpha) | alpha >= 1 | alpha <= 0){
    stop("alpha is the probability of false alarm between zero and one")
  }
  
  if( !is.numeric(t) | round(t) != t | t < 1){
    stop("t is a positive integer")
  }
  
  if( !is.numeric(w) | w <= 0 | round(w) != w)
    stop("w must be a positive integer")
  
  if( !is.numeric(K) | K <= 0 | round(K) != K)
    stop("K must be a positive integer")
  
  method <- match.arg(method)
  
  if(method != "global"){   #number of multiple tests
    bonferroni <- t*(t - w + 1)
  }else{
    bonferroni <- t -w + 1
  }
  
  beta_one <- 1 - (8/K) * log(alpha / bonferroni)
  
  beta_two <- 1 + 2*sqrt(2*log(bonferroni/alpha))
  
  return(max(beta_one, beta_two))
  
  
}
