#' @title Sampling Biomass at each year
#'
#' @description Returns a determined number of biomass samples.
#'
#'
#' @param Pop.Mod A list containing the components returned by Population.Modeling function (main function).
#' @param CV The biomass coefficient of variation. Default value 0, which means that the function returns the biomass computed in the main function of the package, NOT a sample.
#' @param niter_sampling The number of samples to be computed if the Pop.Mod object refers to a deterministic framework. If such object is stochastic (niter>1) for each iteration one sample is computed and hence a value of this parameter is not required.
#' @details A log-normal distribution is used to compute the biomass samples. More precisely, for each year and iteration the value of the biomass in the sample comes from a log-normal distribution centered in the corresponding value of biomass and variability determined by the CV.
#'
#' @return An array containing the samples of the total biomass for each year. The number of samples is equal to niter_sampling in the deterministic framework and to niter in the stochastic one.
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
#' ctrBio<-list(M=M,CV_M=0.2, L_inf=124.5, t0=0, k=0.164, CV_L=0.2, CV_LC=0.2, a=4.5*10^(-6), b=3.1049,
#'            a50_Mat=3, ad_Mat=-0.5,CV_Mat=0.2)
#'
#' # We continue introducing the fishing parameters.
#' # Below, we have different objects ctrSEL depending on which selectivity function is used.
#' # Constant selectivity
#' ctrSEL<-list(type="cte", par=list(cte=0.5),CV_SEL=0.2)
#'
#' # Logistic selectivity
#' ctrSEL<-list(type="Logistic", par=list(a50_Sel=1.5, ad_Sel=-1),CV_SEL=0.2)
#'
#' # Gamma selectivity
#' ctrSEL<-list(type="Gamma", par=list(gamma=10,alpha=15, beta=0.03),CV_SEL=0.05)
#'
#' # Andersen selectivity
#' ctrSEL<-list(type="Andersen", par=list(p1=2,p3=0.2,p4=0.2,p5=40),CV_SEL=0.05)
#'
#' f=matrix(rep(0.5,number_years),ncol=number_years,nrow=2,byrow=TRUE)
#' ctrFish<-list(f=f,ctrSEL=ctrSEL)
#'
#' # Finally, we show below the three possible stock recruitment relationship.
#' # The values of the parameters of Beverton-Holt Recruitment Model and Ricker
#' # Recruitment Model are ones suitables when the biomass is measured in Kg and
#' # the recruitment is measured as number of individuals.
#'
#' a_BH=10000; b_BH=400; CV_REC_BH=0.2; a_RK=10; b_RK=0.0002; CV_REC_RK=0.2
#' # If the spawning stock recruiment relationship is constant:
#' SR<-list(type="cte",par=NULL)
#' # If the spawning stock recruitment relationship is Beverton-Holt Recruitment Model:
#' SR<-list(type="BH",par=c(a_BH,b_BH,CV_REC_BH))
#' # If the spawning stock recruitment relationship is Ricker Recruitment Model:
#' SR<-list(type="RK",par=c(a_RK,b_RK,CV_REC_RK))
#'
#' # The following lines allow us to use the described function.
#' Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)
#' # Now,we can extract the biomass.
#' # Deterministic biomass:
#' B<-Sampling_Biomass(Pop.Mod,CV=0)
#' # Stochastic biomass:
#' B<-Sampling_Biomass(Pop.Mod,CV=0.2)
#'
#' # If niter=1, then the following lines allow us to obtain the biomass samples.
#'
#' # B<-Sampling_Biomass(Pop.Mod,CV=0)
#' # B<-Sampling_Biomass(Pop.Mod,CV=0.2,niter_sampling=1000)
#' @export

Sampling_Biomass<-function(Pop.Mod,CV,niter_sampling){
  if(is.null(CV)){CV=0}

  if(is.numeric(Pop.Mod$Info$seed)){set.seed(Pop.Mod$Info$seed)}
  niter<-dim(Pop.Mod$Matrices$N)[3]
  biomass<-(Sum.Pop.Mod(Pop.Mod,c("BIO"))$BIO)
  N<-Pop.Mod$Matrices$N
  years<-colnames(N)
  number_years<-dim(N)[2]
  if(niter==1){



  if(CV>0){
    column.names <- years
    row.names <- ""
    matrix.names <- 1:niter_sampling
    B_S<-array(rep(0,number_years), dim=c(1, number_years, niter_sampling),dimnames=list("",column.names,matrix.names))
    for (j in 1:number_years){
    v<-(CV*(biomass))^2
    m<-(biomass)
    mu<-log(m^2/sqrt(m^2+v))
    sigma<-sqrt(log(v/m^2+1))




    for (ind in 1:niter_sampling){
      B_S[,j,ind]<-(stats::rlnorm(1, meanlog =mu[j], sdlog =sigma[j]))
    }

    }

  } else {return(biomass)}
}


  if(niter>1){

    column.names <- years
    row.names <- ""
    matrix.names <- 1:niter
    B_S<-array(rep(0,number_years), dim=c(1, number_years, niter),dimnames=list("",column.names,matrix.names))


    if(CV>0){
      for (ind in 1:niter){
      for (j in 1:number_years){
        v<-(CV*(biomass[,,ind]))^2
        m<-biomass[,,ind]
        mu<-log(m^2/sqrt(m^2+v))
        sigma<-sqrt(log(v/m^2+1))





          B_S[,j,ind]<-(stats::rlnorm(1, meanlog =mu[j], sdlog =sigma[j]))
        }

      }

    } else {return(biomass)}
  }



  return(B_S)
}
