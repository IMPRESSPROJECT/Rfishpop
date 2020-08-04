#' @title Yield-per-Recruit
#'
#' @description Return yield-per-Recruit for each iteration.
#'
#'@param Pop.Mod A list containing the components returned by Population.Modeling function (main function).
#' @param f.grid A sequence of fishing efforts.
#' @param Fish.years The number of recent years to estimate the mean of SEL (selectivity, see information about such elements in Sum.Pop.Mod function).
#' @param Bio.years The number of recent years to estimate the mean of M and WC (natural mortality and catch weight, see information about such elements in Population.Modeling function and Sum.Pop.Mod function).
#' @param Method The procedure to obtain the age vector of weight (catches), selectivity and natural mortality. By default is "mean" which means that the mean of the last "Bio.years" is used. The alternative option is "own", the user can introduce these elements.
#' @param par If Method="own" it is a list containing the matrices whose columns report for each iteration the age vector of weight (catches), natural mortality, and selectivity. In other case is equal to NULL.
#' @param plot A vector of two elements. The first one is a logical parameter. By default is equal to TRUE, which means that a biomass per recruit graph is done. The second element refers to which iteration must be plotted.
#' @details The function return the yield-per-recruit.
#'
#' @return An array whose third dimension corresponds to the iterations. For each iteration the array contains a matrix reporting the yield-per-recruit for a range of overall fishing mortalities.
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' \item{Maria Grazia Pennino}
#' }
#' @examples
#' ctrPop<-list(years=seq(1980,2020,by=1),niter=2,N0=15000,ages=0:15,minFage=2,
#'              maxFage=5,tc=0.5,seed=NULL)
#' number_ages<-length(ctrPop$ages);number_years<-length(ctrPop$years)
#'Mvec=c(1,0.6,0.5,0.4,0.35,0.35,0.3,rep(0.3,9))
#'M<-matrix(rep(Mvec,number_years),ncol = number_years)
#'colnames(M)<-ctrPop$years
#'rownames(M)<-ctrPop$ages
#'ctrBio<-list(M=M,CV_M=0.2, L_inf=20, t0=-0.25, k=0.3, CV_L=0, CV_LC=0, a=6*10^(-6), b=3,
#'            a50_Mat=1, ad_Mat=-0.5,CV_Mat=0)
#'ctrSEL<-list(type="cte", par=list(cte=0.5),CV_SEL=0)
#'f=matrix(rep(0.5,number_years),ncol=number_years,nrow=2,byrow=TRUE)
#'ctrFish<-list(f=f,ctrSEL=ctrSEL)
#'a_BH=15000; b_BH=50; CV_REC_BH=0
#'SR<-list(type="BH",par=c(a_BH,b_BH,CV_REC_BH))
#'Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)
#' f.grid<-seq(0.00,1.5,by=0.01)
#' YPR(Pop.Mod,f.grid,3,3,plot=c(TRUE,1), Method="mean",par=NULL)
#'
#' # The following commented lines refers to par argument.
#' # If par is not NULL must be something like (assuming that WC, M,
#' # SEL are defined previously).
#' # par=list(); par$WC<-WC; par$SEL<-SEL; par$M<-M
#' # YPR(Pop.Mod,f.grid,plot=c(TRUE,1),Method="own",par=par)
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
