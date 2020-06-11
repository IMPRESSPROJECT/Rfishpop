#' @title Projecting our Exploited Population on based of desired captures or efforts
#'
#'
#' @description This function allows us to extend our simulated Population through the years on based of the desired captures for such years (strategy="catch") or on based of the desired effort f (component of fishing mortality F = f * SEL) for such years (strategy="effort") .
#'
#'
#' @param Pop.Mod A list containing the components returned by Population.Modeling function (main function).
#' @param new.years The years for which our Population must be projected.
#' @param strategy The strategy for the next years, this can be "catch" if we want to fix the catches for the next years (argument my.catch must be used) or "effort" if we want to fix the value of effort f, component of fishing mortality F = f * SEL (argument my.effort must be used).
#' @param my.catch Weight of the captures for each of the years for which the projections are carried out.
#' @param my.effort A vector containing for the projecting years the value for which the effort of the last year must be multiplied to obtain the desired effort f (component of fishing mortality F = f * SEL) for the new years.
#' @param tol The level of tolerence in the optimization process. If tol=NULL the default value of 0.01 is used.
#' @param limit.f Maximum value of effort f (component of fishing mortality F = f * SEL) considered in the optimization process. If limit.f=NULL the default value of 4 is used.
#' @return The object Pop.Mod updated containing the information for the new years for which the population has been projected.
#' @details This package is an open and developing project, for this reason this function must still be optimized for good execution times.
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' }
#' @examples
#'
#'# First we introduce the basic parameters to define the population.
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
#' # We extract the F and SSB.
#' Sum.Pop.Mod(Pop.Mod,c("F","SSB"))
#' #Pop.Mod=Population.Modeling.Projections(Pop.Mod,new.years=2021,my.catch=1500,
#' #limit.f=NULL,tol=NULL,strategy="catch")
#' #Sum.Pop.Mod(Pop.Mod,c("C"))
#' Pop.Mod=Population.Modeling.Projections(Pop.Mod,new.years=2021,my.effort=0.6,
#' limit.f=NULL,tol=NULL,strategy="effort")
#' Sum.Pop.Mod(Pop.Mod,c("C"))
#' @export




Population.Modeling.Projections=function(Pop.Mod,new.years,my.catch,tol,limit.f,strategy,my.effort){

  if(strategy=="catch"){
  len.new.years=length(new.years)
  len.my.catch=length(my.catch)

  if(len.new.years!=len.my.catch){stop("The lengths of new.years and my.catch vectors must be the same")}
  for(kk in 1:len.new.years){

  Pop.Mod1=Population.Modeling.Projections.Internal(Pop.Mod,new.years=new.years[kk],my.catch=my.catch[kk],tol=tol,limit.f=limit.f)
  Pop.Mod=Pop.Mod1
  }}


  if(strategy=="effort"){
    len.new.years=length(new.years)
    len.my.catch=length(my.effort)

    if(len.new.years!=len.my.catch){stop("The lengths of new.years and my.effort vectors must be the same")}
    for(kk in 1:len.new.years){

      Pop.Mod1=Population.Modeling.Projections.Internal2(Pop.Mod,new.years=new.years[kk],my.effort=my.effort[kk])
      Pop.Mod=Pop.Mod1
    }}

  return(Pop.Mod1)

}


