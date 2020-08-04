#' @title Data for Length Based Indicators (LBI)
#'
#' @description The function provides required information for computing Length Based Indicators: Length distributions of catches and the corresponding average weight per length.
#'
#' @param Pop.Mod A list containing the components returned by Population.Modeling function (main function).
#' @param CV The coefficient of variation associated to the log-normal distribution used in Distribution.length function (see Details of such function).
#' @param RF.value The number of values generated for each age (given a year and an iteration) from the log-normal distribution used in Distribution.length function (see Details of such function). By default RF.value=1000.
#' @details The function reports the length distributions of catches for each year and iteration in our Pop.Mod object. Furthermore, the corresponding average weight per length (for each year and iteration) is also provided. The catches length distribution is computed using the Distribution.length function of this package.
#'
#'
#' @return A list containing the following components.\itemize{
#' \item{length:} the length distributions of catches for each year and iteration.
#' \item{weight:} the average weight per length for each year and iteration.
#'}
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' }
#' @examples
#'
#'# The first step is to simulate the population.
#' ctrPop<-list(years=seq(1980,2020,1),niter=1,N0=15000,ages=0:15,minFage=2,
#' maxFage=5,tc=0.5,seed=NULL)
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
#' # Then the function is used to obtain the length distributions and average
#' # weight per length.
#' # UNCOMMENT THE FOLLOWING LINES
#' #resul=Data.to.LBI(Pop.Mod,CV=0.2)
#' #freq=resul$length[[1]]
#' #wal=resul$weight[[1]]
#'
#' # Furthermore, than the data provided by Data.to.LBI
#' # function the LBI method also needs some life history parameters.
#' # For example:
#'
#' L_inf=Pop.Mod$Info$ctrBio$L_inf # von Bertalanffy asymptotic length (Linf)
#' k=Pop.Mod$Info$ctrBio$k
#' t0=Pop.Mod$Info$ctrBio$t0
#' x50=Pop.Mod$Info$ctrBio$a50_Mat
#' L50=Length_VB(L_inf, k, x50, t0) # Length at 50% maturity (L50)
#' # Since M is a constant through the year but age dependent,
#' # the natural mortality divided by von Bertalanffy k coefficient
#' # is computed as follows, using the mean of the age vector.
#' M.vec=Pop.Mod$Matrices$M[,1,1]
#' MK <-mean(M.vec)/k
#'
#' # Finally, running the following line of code after load the
#' # required code LBI are computed.
#' #LBI=lb_ind(data=freq,binwidth=3,linf=L_inf,lmat=L50,mk_ratio=MK,weight=wal)
#' @export

Data.to.LBI=function(Pop.Mod,CV, RF.value=1000){
  data.list=list()
  weight.list=list()
  a=Pop.Mod$Info$ctrBio$a; b=Pop.Mod$Info$ctrBio$b

  niter=dim(Pop.Mod$Matrices$N)[3]
  number_years=dim(Pop.Mod$Matrices$N)[2]

  if(RF.value!=1000){
  L.D<-Distribution.length(Pop.Mod,CV=CV,Type="LengthC", RF.value=RF.value)}
  L.D<-Distribution.length(Pop.Mod,CV=CV,Type="LengthC")

  for (i in 1:niter){

    freq=(cbind(as.numeric(rownames(L.D))+0.5,L.D[,,i]))
    freq=as.data.frame(freq)
    years=as.numeric(colnames(Pop.Mod$Matrices$N))
    colnames(freq)=c("MeanLength",years)
    wal=freq
    wal[,2:(number_years+1)]<-rep(a*freq[,1]^b,number_years)
    data.list[[i]]=freq
    weight.list[[i]]=wal}

  return(list(length=data.list,weight=weight.list))
}
