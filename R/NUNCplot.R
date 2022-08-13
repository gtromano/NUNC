#'NUNC Plot
#'
#'@description Plotting function for NUNC output
#'
#'@details This function takes an output from NUNC and produces a plot of the data and, where 
#'applicable, the changepoint and time it is detected.
#'
#'@param NUNCobj The output from \code{NUNCoffline}.
#'
#'@param xlab The x axis label
#'
#'@param ylab The y axis label
#'
#'@param main The main plot title
#'
#'@return A plot of the data with vertical lines indicating the changepoint and detection time.
#'
#'@export
#'
#'@examples
#'
#'#Example detection of change in distribution and subsequent plot
#'
#'x <- c(rnorm(1000), rnorm(1000, 5, 3))
#'
#'NUNCrun <- NUNCoffline(x, 100, 12, 15, "local")
#'
#'NUNCplot(NUNCrun)
#'
#'

NUNCplot <- function(NUNCobj, xlab = NA, ylab = NA, main = NA){
  
  if( class(NUNCobj) != "NUNCout"){
    stop("Must input an object of class NUNCout")
  }
  
  
  if( !is.na(xlab)){
    
    if(!is.character(xlab)){
      stop("xlab must be a vector of characters for x axis label")
    }
    
  }
  
  if( !is.na(ylab)){
    
    if(!is.character(ylab)){
      stop("ylab must be a vector of characters for y axis label")
    }
    
  }
  
  if( !is.na(main)){
    
    if(!is.character(main)){
      stop("main must be a vector of characters for plot title")
    }
    
  }
  
  time_vec <- 1:length(NUNCobj$data)
  
  plot(NUNCobj$data, xlab = xlab, ylab = ylab, main = main)
  
  if(NUNCobj$changepoint != -1){
    abline(v = NUNCobj$changepoint, col = "blue")
    abline(v = NUNCobj$t, col = "green")
    legend(x = "topleft",c("Changepoint", "Detection"), 
           col = c("blue", "green"), lty = 1)
  }
  
} 
