#' @title Logistic function
#'
#' @description Computes the logistic function.
#'
#' @param x Points where the logistic function must be computed.
#' @param x50 x-value of the sigmoid's midpoint of the logistic function.
#' @param xd Minus the inverse of the logistic growth rate (steepness of the curve).
#' @details The function implements the logistic function which equation is \deqn{L(x)=1/(1+exp((x-x50)/xd))} where x50 is the x-value of the sigmoid's midpoint of the logistic function, and xd is minus the inverse of the logistic growth rate (steepness of the curve).
#'
#'
#' @return {Vector containing the values of the logistic function in each of the x-values.}
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' }
#' @export

Logistic<-function(x,x50,xd){
  if(x[1]==0){x[1]<-0.0001}
  Lo<-1/(1+exp((x-x50)/xd))
     return(Lo)
}
