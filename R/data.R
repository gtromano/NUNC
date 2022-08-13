#' Data from a Dualshock4 accelometer
#'
#' The accelerometer data is a recording of 100 instances of movement of a Dualshock4 controller.
#' The movement is reflected in the y axis value recorded by the motion sensor in the controller;
#' when the controller moves, the sensor detects this motion as a change in the y-axis value.
#' 
#' In each scenario a change has been inserted at a known time indicating one of four different
#' actions on the controller:
#' \describe{
#'   \item{Picking}{Picking the controller up off a surface.}
#'   \item{Sitting}{Sitting with the controller and starting to move it.}
#'   \item{Sliding}{Sliding the controller along a surface.}
#'   \item{Shaking}{Shaking the controller.}
#' }
#' 
#'
#' @format Accelerometer is a list of 100 data frames. Each data frame contains 3 variables:
#' \describe{
#'  \item{y}{A vector of length 2000 indicating the y axis position recorded by the controller.}
#'  \item{changepoint}{The time the action occurs.}
#'  \item{action}{The action that is performed.}
#' }
#' 
"accelerometer"
