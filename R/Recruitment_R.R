#' @title Ricker Recruitment Model
#'
#' @description Computes the recruitment using Ricker Model.
#'
#'
#' @param SSB Maturity biomass for each year (spawning stock)
#' @param a Recruits-per-spawner at low stock levels.
#' @param b Relates to the rate of decreasing of recruits-per-spawner as SSB increases.
#' @details The formulation of the Ricker equation is: \deqn{R=a*SSB*exp(-b*SSB)} where SSB is the Maturity biomass (spawning stock), a is the recruits-per-spawner at low stock levels and b relates to the rate of decreasing of recruits-per-spawner as SSB increases.
#'
#' @return {Vector containing the recruitment for each year.}
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' }
#' @export

RRK<-function(SSB,a,b){
  niter<-length(SSB)
  av<-rep(a,niter)
  bv<-rep(b,niter)
  R<-(a*SSB)*exp(-b*SSB)
  return(R)
}
