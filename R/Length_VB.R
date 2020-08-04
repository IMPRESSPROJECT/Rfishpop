#' @title Von Bertalanffy Growth Model (Length)
#'
#' @description Computes the length using Von Bertalanffy Growth Model.
#'
#'
#' @param ages Vector of ages.
#' @param L_inf Asymptotic average maximum body size.
#' @param t0 Hypothetical age at which the species has zero length.
#' @param k Growth rate coefficient that determines how quickly the maximum is attained.
#' @details The function implements Von Bertalanffy Growth Model which equation is
#' \deqn{L(x)=L_inf*(1-exp(-k*(x-t0)))} where L_inf is the asymptotic average maximum body size, t0 is hypothetical age at which the species has zero length, and k is growth rate coefficient and x is the vector of ages where the function must be computed.
#' Note if you want the length at some time of the year you only need to sum the corresponding value t to the ages before to introduce them into the function.
#' For example, in order to obtain the length at March you need to introduce ages+t where t=0.3. The year is a unit and hence t is between 0 and 1.
#'
#' @return {L: }{Vector containing the length for each age.}
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' }
#'
#' @export

Length_VB<-function(L_inf,k,ages,t0){
  L<-L_inf*(1-exp(-k*(ages-t0)))
  return(L)
}

