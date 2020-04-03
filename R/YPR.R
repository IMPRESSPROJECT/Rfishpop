#' @title Yield-per-Recruit
#'
#' @description Return yield-per-Recruit for each iteration.
#'
#'@param Pop.Mod A list containing the components returned by Population.Modeling function (main function).
#' @param f.grid A sequence of fishing efforts.
#' @param Fish.years The number of recent years to estimate the mean of SEL (selectivity).
#' @param Bio.years The number of recent years to estimate the mean of M and WC (natural mortality and capture weight).
#' @param Method The procedure to obtain the age vector of weight (captures), selectivity and natural mortality. By default is "mean" which means that the mean of the last "Bio.years" is used. The alternative option is "own", the user can introduce these elements.
#' @param par If Method="own" it is a list containing the matrices whose columns report for each iteration the age vector of weight (captures), natural mortality, and selectivity. In other case is equal to NULL.
#' @param plot A vector of two elements. The first one is a logical parameter. By default is equal to TRUE, which means that a biomass per recruit graph is done. The second element refers to which iteration must be plotted.
#' @details The function return the yield-per-recruit.
#'
#' @return An array whose third dimension corresponds to the iterations. For each iteration the array contains a matrix reporting the yield-per-recruit for a range of overall fishing mortalities.
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
#' f.grid<-seq(0.00,1.5,by=0.01)
#' YPR(Pop.Mod,f.grid,3,3,plot=c(TRUE,1), Method="mean",par=NULL)
#' # If par is not NULL must be something like (assuming that WC, M,
#' # and SEL are defined previously).
#' # par=list(); par$WC<-WC; par$SEL<-SEL; par$M<-M
#' @export





YPR<-function(Pop.Mod,f.grid,Fish.years,Bio.years,plot,Method,par){
  if(is.null(Method)){Method="mean"}
  if(is.null(plot)){plot=TRUE}
  min_age<-Pop.Mod$Info$minFage
  max_age<-Pop.Mod$Info$maxFage
  N<-Pop.Mod$Matrices$N;ages<-as.numeric(rownames(N))

  number_ages<-nrow(N)
  number_years<-ncol(N)
  niter<-dim(N)[3]

  if(Method=="mean"){
  W.ct<-(Sum.Pop.Mod(Pop.Mod,c("WC"))$WC)
  Mt<-Pop.Mod$Matrices$M
  st<-(Sum.Pop.Mod(Pop.Mod,c("SEL"))$SEL)

  W.c<-matrix(0,ncol=niter,nrow=number_ages);M<-W.c;s<-W.c
  for (ind2 in 1:niter){
  W.c[,ind2]<-apply(W.ct[,(number_years-Bio.years+1):number_years,ind2],1,mean)
  M[,ind2]<-apply(Mt[,(number_years-Bio.years+1):number_years,ind2],1,mean)
  s[,ind2]<-apply(st[,(number_years-Fish.years+1):number_years,ind2],1,mean)
  }}

  if(Method=="own"){
    W.c<-par$WC
    M<-par$M
    s<-par$SEL
  }

  l<-length(f.grid);F.vector<-1:l;ypr.vector<-1:l


  R<-array(0, dim=c(l, 2, niter))
  for (ind2 in 1:niter){
  a<-1:l
  for (j in 1:l){
    F<-f.grid[j]*s #Now is a vector but later will be a matrix and a mean will be necessary
    F<-matrix(F,ncol=niter)
    Z<-F+M
    N<-1
    ypr<-0


  for (i in 1:(number_ages-1)){

    ypr<-ypr+N*F[i,ind2]*W.c[i,ind2]*((1-exp(-Z[i,ind2]))/Z[i,ind2])
    N<-N*exp(-Z[i,ind2])
  }

  ypr<-ypr+(N*F[number_ages,ind2]*W.c[number_ages,ind2]/Z[number_ages,ind2])

  ypr.vector[j]<-ypr
  F.vector[j]<-mean(F[(min_age-ages[1]+1):(max_age-ages[1]+1),ind2])
  }
  ind<-which(ypr.vector<0);ypr.vector[ind]<-0
  results<-cbind(F.vector,ypr.vector)
  colnames(R)<-c("F","YPR")

  R[,,ind2]<-results}
  if(plot[1]==1){
    graphics::plot(R[,1,plot[2]],R[,2,plot[2]], main="Yield per recruit",xlab="F (fishing mortality)",ylab="YPR",type="l")
  }
  return(R)}
