#' @title Steepness of the Stock Recruitment Relationship
#'
#' @description Returns the steepness of the stock recruitment relationship if type="steepness" or the values of the parameters of the stock recruitment relationship for a given value of the steepness if type="parameters".
#'
#'
#' @param Pop.Mod A list containing the components returned by Population.Modeling function (main function).
#' @param Fish.years The number of recent years to estimate the mean of SEL (selectivity).
#' @param Bio.years The number of recent years to estimate the mean of M, Mat, WC, and W (natural mortality, maturity, stock weight and capture weight).
#' @param type The  desired approach which has two possibilities, type="steepness" or  type="parameters".
#' @param h The desired value of the steepness when type="parameters". In other case (type="steepness") this parameter is equal to NULL.
#' @param Method The procedure to obtain the age vector of weight (stock and captures), natural mortality, selectivity and maturity. By default is "mean" which means that the mean of the last "Bio.years" is used. The alternative option is "own", the user can introduce these elements.
#' @param par If Method="own" it is a list containing the matrices whose columns report for each iteration the age vector of weight (stock and captures), natural mortality, selectivity and maturity. In other case is equal to NULL.
#' @details The function returns the steepness of a stock recruitment relationship. Remember that the steepness is commonly defined as the fraction of recruitment from an unfished population obtained when the spawning stock biomass is 20 percentage of its unfished level. The other point of view and possibility of this function is to compute the parameters of a stock recruitment model for which the corresponding steepness is equal to a desired value.
#'
#' @return \itemize{
#' \item{h:}{if type="steepness" is an array whose third dimension is the number of iterations. For each iteration the value of the steepness is reported.}
#' \item{parameters:}{if type="parameters" is an array whose third dimension is the number of iterations. For each iteration the value of the parameters of the stock-recruitment relationship .}}
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
#' f=rep(0.5,number_years)
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
#'
#'
#'
#'
#' # Suppose that we are interesting in obtaining the values of the parameters of
#' #the stock-recruitment model for a given steepness (h=0.5). The following line
#' #allows us to do that.
#' parameters_value<-steepness(Pop.Mod,3,3,type="parameters",h=0.5,Method="mean",par=NULL)
#' # Suppose that we are interesting in obtaining the value of the steepness for
#' #the stock-recruitment model which is used. The following line
#' #allows us to do that.
#' steepness_value<-steepness(Pop.Mod,3,3,type="steepness",Method="mean",par=NULL)
#' @export
steepness<-function(Pop.Mod,Fish.years,Bio.years,type,h,Method,par){
  if(is.null(Method)){Method="mean"}

  SR<-Pop.Mod$Info$SR
  N<-Pop.Mod$Matrices$N;ages<-as.numeric(rownames(N))
  niter<-dim(N)[3]

  if (type=="steepness"){
    RE<-BYR.eq(Pop.Mod,0,Fish.years,Bio.years,FALSE,Method,par)
    A0<-array(0,dim=c(1,1,niter))
    R<-array(0,dim=c(1,1,niter))
    h<-array(0,dim=c(1,1,niter))
  for(ind2 in 1:niter){
  A0[1,1,ind2]<-RE$DPM[,7,ind2]
  SR<-Pop.Mod$Info$SR
  type<-SR$type
  if (type=="BH"){
    a_BH<-SR$par[1]
    b_BH<-SR$par[2]
    R[1,1,ind2]<-(A0[1,1,ind2]*a_BH)/(b_BH+A0[1,1,ind2])
    h[1,1,ind2]<-((0.2*A0[1,1,ind2]*a_BH)/(b_BH+0.2*A0[1,1,ind2]))*(1/R[1,1,ind2])
  }

  if (type=="RK"){
    a_RK<-SR$par[1]
    b_RK<-SR$par[2]
    R[1,1,ind2]<-(a_RK*A0[1,1,ind2]*exp(-b_RK*A0[1,1,ind2]))
    h[1,1,ind2]<-(a_RK*0.2*A0[1,1,ind2]*exp(-b_RK*0.2*A0[1,1,ind2]))*(1/R[1,1,ind2])
  }
  }
    return(h=h)}

  if (type=="parameters"){
    A0<-array(0,dim=c(1,1,niter))
    R0<-array(0,dim=c(1,1,niter))
    a_BH<-array(0,dim=c(1,1,niter))
    b_BH<-array(0,dim=c(1,1,niter))
    a_RK<-array(0,dim=c(1,1,niter))
    b_RK<-array(0,dim=c(1,1,niter))
    RE<-BYR.eq(Pop.Mod,0,Fish.years,Bio.years,FALSE,Method,par)
    for(ind2 in 1:niter){


    R0[1,1,ind2]<-RE$DPM[,5,ind2]
    A0[1,1,ind2]<-RE$DPM[,7,ind2]
    SR<-Pop.Mod$Info$SR
    type<-SR$type
    if (type=="BH"){
      a_BH[1,1,ind2]<-(4*h*R0[1,1,ind2])/(5*h-1)
      b_BH[1,1,ind2]<-(A0[1,1,ind2]*(1-h))/(5*h-1)

    }

    if (type=="RK"){
      b_RK[1,1,ind2]<-log(5*h)/(0.8*A0[1,1,ind2])
      a_RK[1,1,ind2]<-exp(b_RK[1,1,ind2]*A0[1,1,ind2])/(A0[1,1,ind2]/R0[1,1,ind2])

    }

    }
    if (type=="BH"){return(parameters=list(a_BH=a_BH,b_BH=b_BH))}
    if (type=="RK"){return(parameters=list(a_RK=a_RK,b_RK=b_RK))}
  }
}



