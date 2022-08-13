# finds the grid of points for which to search over
# W is the window size
# beta is the penalty value
# J is the number of checks to be performed

.local_grid <- function(W, beta, J, side = "right"){
  
  if(J >= 2){
    
    # now find the maximum value for which detection is possible
    if( .local_bound(floor(W/2), W) >= beta){
      # uniroot is way faster than sapply for large windows
      root <- floor(uniroot(function(x) .local_bound(x, W) - beta,
                                interval = c((floor(W/2)), (W-1)))$root)
      tau_star <- which((floor(W/2)) : (W-1) == root)
      
    }else{
      tau_star <- W / 2
    }
    
    #note now we only need to search over tau_star + 1 points at most
    # if(J > (tau_star + 1 ))
    #   J <- tau_star + 1
    
    # if (side == "right")
    #   trimmed_grid <- floor(W/2) : ( floor(W/2) - 1 + tau_star)
    # else
    #   trimmed_grid <- ( floor(W/2) - tau_star) : floor(W/2)
    
    trimmed_grid <- ( floor(W/2) + 1 - tau_star) : ( floor(W/2) - 1 + tau_star)
    
     if(J > length(trimmed_grid))
       J <- length(trimmed_grid)
    
    
    c <- log(2*W-1)
    # probs <- 1 / (1 +  (2 * ((1:(J-1)) + 10) * exp( (-c / J) * (2 * (1:(J-1)) - 1) ) ) )
    probs <- seq(1/J, by = 1/J, length.out = J - 1) # equally spaced grid
    spaced_grid <- round(quantile(trimmed_grid, probs = probs))
    spaced_grid <- c(floor(W/2), spaced_grid)  #always test the centre
    
  } else {
    spaced_grid <- floor(W/2)   #only test centre, if only one point to be tested
  }
  
  return(spaced_grid)
}


# calculates the bound from proposition one
.local_bound <- function(tau, W){
  return(-(tau / W) * log(tau / W) ) - ( (1 - (tau / W)) * log(1 - (tau / W)) )
}