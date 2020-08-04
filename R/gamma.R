#' @title Gamma Selectivity Function
#'
#' @description Computes Gamma selectivity function.
#'
#' @param x Points where Gamma function must be computed (see Details).
#' @param gamma the size of the mesh used by the fleet (according to the units of the remaining parameters, for example, cm). Negative values or zero are NOT valid.
#' @param beta the scale parameter. Negative values or zero are NOT valid.
#' @param alpha the shape parameter. Negative values are NOT valid.
#' @details The function implements Gamma selectivity function which equation is \deqn{S(x)<-((x/((alpha-1) beta gamma))^(alpha-1))exp(alpha-1-(1/beta gamma)),} where gamma is the size of the mesh, alpha is the shape parameter and beta is the scale parameter.
#'
#'
#' @return {SG: }{Vector containing the values of Gamma function in each of the x-values.}
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' }
#' @export
gamma_SEL<-function(x,alpha,beta,gamma){
  if(x[1]==0){x[1]<-0.0001}
  sel_gamma<-((x/((alpha-1)*beta*gamma))^(alpha-1))*exp(alpha-1-(x/(beta*gamma)))
  return(SG=sel_gamma)}
