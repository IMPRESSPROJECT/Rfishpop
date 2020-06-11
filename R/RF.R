#' @title Reference Fishery Mortalities
#'
#' @description Returns the reference fishery mortality which produces maximum YPR (FM_type="F_max"), the reference fishery mortality at which the slope of the YPR curve is reduced to 0.1 of that estimated at the origin (FM_type="F_0.1"), the reference fishery mortality at which the MSY is attained (FM_type="F_msy") and the reference fishery mortality which will drive the stock to extinction (FM_type="F_Crash"). Furthermore for each of these fishery mortalities the corresponding effort of fishing, YPR, BYR, B, Y and R in equilibrium is also returned.
#'
#' @param Pop.Mod A list containing the components returned by Population.Modeling function (main function).
#' @param Fish.years The number of recent years to estimate the mean of SEL (selectivity).
#' @param Bio.years The number of recent years to estimate the mean of M, Mat, WC, and W (natural mortality, maturity, stock weight and capture weight).
#' @param Method The procedure to obtain the age vector of weight (stock and captures), natural mortality, selectivity and maturity. By default is "mean" which means that the mean of the last "Bio.years" is used. The alternative option is "own", the user can introduce these elements.
#' @param par If Method="own" it is a list containing the age vector of weight (stock and captures), natural mortality, selectivity and maturity (for the first iteration). In other case is equal to NULL.
#' @param FM_type which of the four reference fishery mortalities must be computed. The possibilities have been described above: FM_type="F_max", FM_type="F_0.1", FM_type="F_msy" and FM_type="F_Crash".
#' @param iters A vector containing the iteration for which the reference fishery mortalities must be computed.
#' @details The function returns the reference fishery mortality which produces maximum YPR (FM_type="F_max"), the reference fishery mortality at which the slope of the YPR curve is reduced to 0.1 of that estimated at the origin (FM_type="F_0.1"), the reference fishery mortality at which the MSY is attained (FM_type="F_msy") and  the reference fishery mortality which will drive the stock to extinction (FM_type="F_Crash"). Furthermore for each of these fishery mortalities the corresponding effort of fishing, YPR, BYR, B, Y and R in equilibrium is also returned. If the fishing effort is equal to 10 can be that optimize process had not found the correct value in the default sequence (except for FM_type="F_Crash", it this case the sequence finish at 60).
#'
#' @return One of the three following elements depending the above selection:
#' \item{F_max:}{the value of F that produces maximum YPR with the corresponding effort of fishing, YPR, BYR, B, Y and R in equilibrium.}
#' \item{F_0.1:}{the value of F at which the slope of the YPR curve is reduced to 0.1 of that estimated at the origin with the corresponding effort of fishing, YPR, BYR, B, Y and R in equilibrium.}
#' \item{F_msy:}{the value of F at which the MSY is attained with the corresponding effort of fishing, YPR, BYR, B, Y and R in equilibrium.}
#' \item{F_Crash:}{the value of F which will drive the stock to extinction with the corresponding effort of fishing, YPR, BYR, B, Y and R in equilibrium.}
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' }
#' @examples
#' ctrPop<-list(years=seq(1980,2020,by=1),niter=2,N0=10000,ages=0:15,minFage=4,
#' maxFage=7,ts=0,tc=0.5,tseed=NULL)
#' number_ages<-length(ctrPop$ages);number_years<-length(ctrPop$years)
#' M<-matrix(rep(0.4,number_ages*number_years),ncol = number_years)
#' colnames(M)<-ctrPop$years
#' rownames(M)<-ctrPop$ages
#' ctrBio<-list(M=M,CV_M=0.2, L_inf=124.5, t0=0, k=0.164, CV_L=0.2, CV_LC=0.2,
#' a=4.5*10^(-6),b=3.1049,a50_Mat=3, ad_Mat=-0.5,CV_Mat=0.2)
#' ctrSEL<-list(type="Logistic", par=list(a50_Sel=1.5, ad_Sel=-1),CV_SEL=0.2)
#' f=matrix(rep(0.5,number_years),ncol=number_years,nrow=2,byrow=TRUE)
#' ctrFish<-list(f=f,ctrSEL=ctrSEL)
#' a_BH=10000; b_BH=400; CV_REC_BH=0.2
#' SR<-list(type="BH",par=c(a_BH,b_BH,CV_REC_BH))
#' Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)
#' # We compute the reference fishery mortalities.
#' # RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type="F_max",iters=1:2)
#' # RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type="F_0.1",iters=1:2)
#' # RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type="F_msy",iters=1:2)
#' # RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type="F_Crash",iters=1:2)
#' # If par is not NULL must be something like (assuming that W, WC,
#' # M, Mat and SEL are defined previously).
#' # par=list(); par$W<-W; par$WC<-WC; par$SEL<-SEL; par$Mat<-Mat; par$M<-M

#'
#' @export
#'
RF<-function(Pop.Mod, Fish.years,Bio.years,Method,par,FM_type,iters){
  if(is.null(Method)){Method="mean"}
  n=length(iters)


  if(n<=20){

    rest<-20-n
    iters_list<-as.list(c(iters,rep(1,rest)))
    Resul<-array(0,dim=c(1,7,n))
    Pop.Mod2<-Pop.Mod

    Pop.Mod2$Matrices$N<-Pop.Mod$Matrices$N[,,c(iters,rep(1,rest))]
    Pop.Mod2$Matrices$F<-Pop.Mod$Matrices$F[,,c(iters,rep(1,rest))]
    Pop.Mod2$Matrices$M<-Pop.Mod$Matrices$M[,,c(iters,rep(1,rest))]
    Pop.Mod2$Matrices$W<-Pop.Mod$Matrices$W[,,c(iters,rep(1,rest))]
    Pop.Mod2$Matrices$Mat<-Pop.Mod$Matrices$Mat[,,c(iters,rep(1,rest))]
    Pop.Mod2$Matrices$C_N<-Pop.Mod$Matrices$C_N[,,c(iters,rep(1,rest))]
    Pop.Mod2$Matrices$C_W<-Pop.Mod$Matrices$C_W[,,c(iters,rep(1,rest))]

    Resul[,,iters]<-RF_U(Pop.Mod2, Fish.years = Fish.years,Bio.years = Bio.years,Method=Method,par=par,FM_type=FM_type,iters =iters_list)[,,iters]

  }



  if (n>20){
    blo<-20
    block<-n%/%blo
    rest<-n-blo*block
    seq1<-c(0,(1:(block))*20)+1
    seq1<-seq1[-length(seq1)]
    seq2<-seq(blo,n-rest,by=blo)
    ll<-length(seq1)
    Resul<-array(0,dim=c(1,7,n))
    Pop.Mod2<-Pop.Mod
    for(i in 1:(ll)){
      Pop.Mod2$Matrices$N<-Pop.Mod$Matrices$N[,,seq1[i]:seq2[i]]
      Pop.Mod2$Matrices$F<-Pop.Mod$Matrices$F[,,seq1[i]:seq2[i]]
      Pop.Mod2$Matrices$M<-Pop.Mod$Matrices$M[,,seq1[i]:seq2[i]]
      Pop.Mod2$Matrices$W<-Pop.Mod$Matrices$W[,,seq1[i]:seq2[i]]
      Pop.Mod2$Matrices$Mat<-Pop.Mod$Matrices$Mat[,,seq1[i]:seq2[i]]
      Pop.Mod2$Matrices$C_N<-Pop.Mod$Matrices$C_N[,,seq1[i]:seq2[i]]
      Pop.Mod2$Matrices$C_W<-Pop.Mod$Matrices$C_W[,,seq1[i]:seq2[i]]
      Resul[,,seq1[i]:seq2[i]]<-RF_U(Pop.Mod2, Fish.years = Fish.years,Bio.years = Bio.years,Method=Method,par=par,FM_type=FM_type,iters =1:blo)
    }

    if(rest>0){
      i=ll
      a<-c((seq2[i]+1):n,1:(blo-rest))
      Pop.Mod2$Matrices$N<-Pop.Mod$Matrices$N[,,a]
      Pop.Mod2$Matrices$F<-Pop.Mod$Matrices$F[,,a]
      Pop.Mod2$Matrices$M<-Pop.Mod$Matrices$M[,,a]
      Pop.Mod2$Matrices$W<-Pop.Mod$Matrices$W[,,a]
      Pop.Mod2$Matrices$Mat<-Pop.Mod$Matrices$Mat[,,a]
      Pop.Mod2$Matrices$C_N<-Pop.Mod$Matrices$C_N[,,a]
      Pop.Mod2$Matrices$C_W<-Pop.Mod$Matrices$C_W[,,a]
      Resul[,,(seq2[i]+1):n]<-RF_U(Pop.Mod2, 3,3,Method="mean",par=NULL,FM_type=FM_type,iters =1:blo)[,,1:rest]
    }}
  colnames(Resul)=c("f","F","YPR","BPR","R","Y","B")
  return(Resul)

}
