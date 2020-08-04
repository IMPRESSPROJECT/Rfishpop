#' @title Length Based Spawning Potential Ratio (LB-SPR)
#'
#' @description The function provides required information for computing Length Based Spawning Potential Ratio (LB-SPR): Length distributions of catches.
#'
#' @param Pop.Mod A list containing the components returned by Population.Modeling function (main function).
#' @param CV The coefficient of variation associated to the log-normal distribution used in Distribution.length function (see Details of such function).
#' @param RF.value The number of values generated for each age (given a year and an iteration) from the log-normal distribution used in Distribution.length function (see Details of such function). By default RF.value=1000.
#' @details The function reports the length distributions of catches for each year and iteration in our Pop.Mod object.
#'
#'
#' @return length: the length distributions of catches for each year and iteration.
#'
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' }
#' @examples
#' # The first step is to simulate the population as follows.
#' ctrPop<-list(years=seq(1980,2020,1),niter=1,N0=15000,ages=0:15,minFage=2,
#' maxFage=5,ts=0,tc=0.5,seed=NULL)
#' number_ages<-length(ctrPop$ages);number_years<-length(ctrPop$years)
#' Mvec=c(1,0.6,0.5,0.4,0.35,0.35,0.3,rep(0.3,9))
#' M<-matrix(rep(Mvec,number_years),ncol = number_years)
#' colnames(M)<-ctrPop$years
#' rownames(M)<-ctrPop$ages
#' ctrBio<-list(M=M,CV_M=0, L_inf=20, t0=0, k=0.3, CV_L=0, CV_LC=0, a=6*10^(-6), b=3,
#'              a50_Mat=4, ad_Mat=-0.2,CV_Mat=0)
#'# Logistic selectivity
#' ctrSEL<-list(type="Logistic", par=list(a50_Sel=2.3, ad_Sel=-0.2),CV_SEL=0)
#' f=matrix(rep(0.5,number_years),ncol=number_years,nrow=1,byrow=TRUE)
#' ctrFish<-list(f=f,ctrSEL=ctrSEL)
#' a_BH=15000; b_BH=50; CV_REC_BH=0
#' SR<-list(type="BH",par=c(a_BH,b_BH,CV_REC_BH))
#' Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)
#'
#' # Then, we use the function to obtain the catches length distribution.
#' # UNCOMMENT THE FOLLOWING LINES
#' # resul=Data.to.LB.SPR(Pop.Mod,CV=0.2)
#'
#' # Furthermore, than the data provided by Data.to.LB.SPR
#' # function the LB-SPR method also needs some life history parameters.
#' # For example:
#'
#' L_inf=Pop.Mod$Info$ctrBio$L_inf # von Bertalanffy asymptotic length (Linf)
#' k=Pop.Mod$Info$ctrBio$k
#' t0=Pop.Mod$Info$ctrBio$t0
#' x50=Pop.Mod$Info$ctrBio$a50_Mat
#' L50=Length_VB(L_inf, k, x50, t0) # Length at 50% maturity (L50)
#' xd=Pop.Mod$Info$ctrBio$ad_Mat
#' x95=(-xd)*log(19)+x50
#' L95=Length_VB(L_inf, k, x95, t0) # Length at 95% maturity (L95)
#' # Since M is a constant through the year but age dependent,
#' # the natural mortality divided by von Bertalanffy K coefficient
#' # is computed as follows, using the mean of the age vector.
#'   M.vec=Pop.Mod$Matrices$M[,1,1]
#'   MK <-mean(M.vec)/k
#' # All the information collected above must be used through the
#' # following lines of code to apply the LB-SPR model.
#'
#' # CREATE THE OBJECT LB_pars
#' #library(LBSPR)
#' #MyPars <- new("LB_pars")
#' #MyPars@Species <- "MySpecies"
#' #MyPars@Linf <- L_inf
#' #MyPars@L50 <-L50 #Length at 50% maturity (L50)
#' #MyPars@L95 <-L95 #Length at 95% maturity (L95)
#' #MyPars@MK <-MK  #the natural mortality divided by von Bertalanffy k coefficient
#' # Add such information to the length distribution.
#' #freq=resul[[1]]
#' #write.csv(freq, file="len.csv")
#' #Len <- new("LB_lengths", LB_pars=MyPars, file=paste0("len.csv"),
#' #dataType="freq", header=TRUE)
#' # Now, the distribution length has the required format to be
#' # introduced in the LBSPRfit().
#' # library(LBSPR)
#' # myFit<- LBSPRfit(MyPars,Len)
#' @export

Data.to.LB.SPR=function(Pop.Mod,CV,RF.value=1000){
  Len.list=list()


  niter=dim(Pop.Mod$Matrices$N)[3]
  number_years=dim(Pop.Mod$Matrices$N)[2]

  if(RF.value!=1000){
    L.D<-Distribution.length(Pop.Mod,CV=CV,Type="LengthC", RF.value=RF.value)}
  L.D<-Distribution.length(Pop.Mod,CV=CV,Type="LengthC")

  for (i in 1:niter){
    freq=L.D[,,i]
    freq=as.data.frame(freq)
    years=as.numeric(colnames(Pop.Mod$Matrices$N))
    colnames(freq)=c(years)
    rownames(freq)=as.numeric(rownames(L.D))+0.5
    Len.list[[i]]=freq

  }
  return(length=Len.list)
}
