#' @title Length-Weight relationship
#'
#' @description Computes the weight from the length using the relationship explained in Details.
#'
#'
#' @param L Vector of lengths by age.
#' @param a Allometric growth parameter.
#' @param b Scaling constant.
#' @details The function implements the following Length-Weight relationship
#'  \deqn{W=a*L^b} where L is the vector of lengths by age, a is the allometric growth parameter, b scaling constant, and W is the age weight vector.
#'
#' @return {Vector containing the weight for each age.}
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' }
#' @export

Weight<-function(L,a,b){
  number_ages<-length(L)
  av<-rep(a,number_ages)
  bv<-rep(b,number_ages)
  W<-av*L^bv
  return(W)
}
