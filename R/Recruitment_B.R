#' @title Beverton-Holt Recruitment Model
#'
#' @description Computes the recruitment using Beverton-Holt Recruitment Model.
#'
#'
#' @param SSB Maturity biomass for each year (spawning stock)
#' @param a Maximum number of recruitments produced.
#' @param b Spawning stock needed to produce recruitment equal to half maximum.
#' @details The formulation of the Beverton-Holt Recruitment Model equation is: \deqn{R=a*SSB/(b+SSB)} where SSB is the Maturity biomass (spawning stock), a is the maximum number of recruitments produced and b is the spawning stock needed to produce recruitment equal to half maximum.
#'
#' @return {Vector containing the recruitment for each year.}
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' }
#' @export

RBH<-function(SSB,a,b){
  niter<-length(SSB)
  av<-rep(a,niter)
  bv<-rep(b,niter)
  R<-(av*SSB)/(bv+SSB)
  return(R)
}
