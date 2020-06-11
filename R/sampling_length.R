#' @title Sampling length (stock or captures)
#'
#' @description Return a sample of the stock or capture length from the corresponding distribution computed using Distribution.length function.
#'
#'
#' @param L.D The distribution length returned by Distribution.length function.
#' @param sample.size The sample size of the desired sample.
#' @return An array containing in each column the corresponding length sample for each iteration (third dimension).
#' @details The function returns a length sample (stock or capture) generating random values from the computed length distribution function.
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' }
#' @examples
#' # First we introduce the basic parameters to define the population.
#' # Note that N0 is equal to 10000 individuals, and hence below we are
#' # consistent with this unit when we introduce the biological and
#' # stock-recruitment parameters.
#' ctrPop<-list(years=seq(1980,2020,by=1),niter=2,N0=10000,ages=0:15,minFage=4,
#' maxFage=7,ts=0,tc=0.5,tseed=NULL)
#'
#' # Now, we introduce the biological parameters of the population.
#' # Note that L_inf is in cm, and a and b parameters allow us to relate
#' # the length in cm with the weight in Kg.
#' number_ages<-length(ctrPop$ages);number_years<-length(ctrPop$years)
#' M<-matrix(rep(0.4,number_ages*number_years),ncol = number_years)
#' colnames(M)<-ctrPop$years
#' rownames(M)<-ctrPop$ages
#' ctrBio<-list(M=M,CV_M=0.2, L_inf=124.5, t0=0, k=0.164, CV_L=0, CV_LC=0, a=4.5*10^(-6), b=3.1049,
#'            a50_Mat=3, ad_Mat=-0.5,CV_Mat=0)
#'
#' # We continue introducing the fishing parameters.
#' # Below, we have different objects ctrSEL depending on which selectivity function is used.
#' # Logistic selectivity
#' ctrSEL<-list(type="Logistic", par=list(a50_Sel=1.5, ad_Sel=-1),CV_SEL=0)
#'
#' f=matrix(rep(0.5,number_years),ncol=number_years,nrow=2,byrow=TRUE)
#' ctrFish<-list(f=f,ctrSEL=ctrSEL)
#'
#' # Finally, we show below the three possible stock recruitment relationship.
#' # The values of the parameters of Beverton-Holt Recruitment Model and Ricker
#' # Recruitment Model are ones suitables when the biomass is measured in Kg and
#' # the recruitment is measured as number of individuals.
#'
#' a_BH=10000; b_BH=400; CV_REC_BH=0
#' # If the spawning stock recruitment relationship is Beverton-Holt Recruitment Model:
#' SR<-list(type="BH",par=c(a_BH,b_BH,CV_REC_BH))
#'
#' # The following lines allow us to use the described function.
#' Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)
#'
#' # We need to compute the distribution length.
#'
#' L.D<-Distribution.length(Pop.Mod,CV=0.2,Type="LengthC")
#'
#' # Now, from the distribution function we generated a sample.
#'
#' our.sample<-Sampling_length(L.D,sample.size=100)
#' @export

Sampling_length<-function(L.D, sample.size){
  column.names <- colnames(L.D)
  number_years<-dim(L.D)[2]
  niter<-dim(L.D)[3]
  Our.sample<-array(0, dim=c(sample.size, number_years, niter),dimnames=list(1:sample.size,column.names,1:niter))


  for (i in 1:niter){
    Our.sample[,,i]<-apply(L.D[,,i],2,f<-function(v,n=sample.size){random.generator(v,n)})
  }

return(Our.sample)




}
