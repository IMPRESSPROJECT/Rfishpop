#' @title Andersen Selectivity Function
#'
#' @description Computes Andersen selectivity function.
#'
#' @param x Points where Andersen function (see Details) must be computed.
#' @param p1 Define the value of \eqn{x} at which the transition between the two gaussian functions happens; that is, \eqn{x=p5/exp(p1)}. Negative values are NOT valid.
#' @param p3 Ascending slope parameter. Negative values or zero are NOT valid.
#' @param p4 Descending slope parameter.  Negative values or zero are NOT valid.
#' @param p5 Define the value of \eqn{x} at which the transition between the two gaussian functions happens; that is, \eqn{x=p5/exp(p1)}.  Negative values or zero are NOT valid.
#' @details The function implements Andersen selectivity function whose equation is \deqn{S(x)<-p0+p2exp(-(ln(p5/x)-p1)^2/p4) if ln(p5/x)<=p1} and \deqn{S(x)<-p0+p2exp(-(ln(p5/x)-p1)^2/p3) if ln(p5/x)>p1.} We fixed \eqn{p0=0}, which is the beginning size for the plateau, and \eqn{p2=1}, which is the maximum value attainable by the function. Remember that \eqn{p3} and \eqn{p4} are the ascending and descending slope parameters, respectively, whereas \eqn{p1} and \eqn{p5} define the value of \eqn{x} at which the transition between the two gaussian functions happens; that is, \eqn{x=p5/exp(p1)}.
#'
#'
#' @return {SA: }{Vector containing the values of Andersen function in each of the x-values.}
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' }
#' @export



andersen<-function(x,p1,p3,p4,p5){
  if(x[1]==0){x[1]<-0.0001}
  p0<-0
  p2<-1

  number_ages<-length(x);sel_andersen<-1:number_ages

  for (i in 1:number_ages){
    if(log(p5/x[i])<=p1){
      sel_andersen[i]<-p0+p2*exp(-((log(p5/x[i])-p1)^2/p4))
    }


    if(log(p5/x[i])>p1){
      sel_andersen[i]<-p0+p2*exp(-((log(p5/x[i])-p1)^2/p3))
    }

  }
  return(SA=sel_andersen)
}
