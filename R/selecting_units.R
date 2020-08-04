#' @title Selecting Units
#'
#' @description Returns some of the matrices provided by Sum.Pop.Mod function but in tonnes instead of Kg.
#'
#'
#' @param Pop.Mod A list containing the components returned by Population.Modeling function (main function).
#' @param Elements A vector specifing which of the following elements must be reported by the function.\itemize{
#' \item{"C":}{Weight of the catches for each year and iteration.}
#' \item{"BIO":}{Total biomass for each year and iteration.}
#' \item{"SSB":}{Maturity biomass for each year (spawning stock biomass) and iteration.}}
#' @return A list containing the objects specified before using argument "Elements" in tonnes.
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' }
#' @examples
#'
#' #selecting_units(Pop.Mod,c("C","BIO","SSB"))
#' @export

selecting_units<-function(Pop.Mod,Elements){

  E<-Sum.Pop.Mod(Pop.Mod,Elements)


  E$C=E$C/1000
  E$BIO=E$BIO/1000
  E$SSB=E$SSB/1000

  return(E)
}
