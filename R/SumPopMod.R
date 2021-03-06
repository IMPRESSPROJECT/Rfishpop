#' @title Information of the Exploited Population (Structured by Age) simulated using Population.Modeling.
#'
#' @description This function allows us to extract additional information obtained in the simulation process of Population.Modeling (main function). The specified information that can be extracted is explained above.
#'
#'
#' @param Pop.Mod A list containing the components returned by Population.Modeling function (main function).
#' @param Elements A vector specifing which of the following elements must be reported by the function.\itemize{
#' \item{"Z": }{Third dimensional array containing the instantaneous mortality for each age, year and iteration.}
#' \item{"LS":}{Third dimensional arraycontaining the (stock) length for each age, year and iteration (at 1st of January).}
#' \item{"LC":}{Third dimensional array containing the length of the catches for each age, year and iteration (at tc).}
#' \item{"WS":}{Third dimensional array containing the population weight for each age, year and iteration.}
#' \item{"WSSB":}{Third dimensional array containing the weight of the mature population for each age, year and iteration.}
#' \item{"C":}{Weight of the catches for each year and iteration.}
#' \item{"SEL":}{Selectivity by age, for each iteration.}
#' \item{"BIO":}{Total biomass for each year and iteration.}
#' \item{"SSB":}{Maturity biomass for each year (spawning stock) and iteration. First age (0) is not contributing to the SSB, their contribution has been assumed negligible.}
#' \item{"REC":}{Population numbers at first age.}
#' \item{"F":}{Mean fishing mortality (only takes the values between the minFage and maxFage.}
#' \item{"WC":}{Third dimensional array containing the weight corresponding to the catch length for each age (at tc), year and iteration.}}
#' @return A list containing the the objects specified before using argument "Elements".
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
#' maxFage=7,tc=0.5,seed=NULL)
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
#' CV_REC_C=0.2
#' # If the spawning stock recruiment relationship is constant:
#' SR<-list(type="cte",par=c(CV_REC_C))
#' # If the spawning stock recruitment relationship is Beverton-Holt Recruitment Model:
#' SR<-list(type="BH",par=c(a_BH,b_BH,CV_REC_BH))
#' # If the spawning stock recruitment relationship is Ricker Recruitment Model:
#' SR<-list(type="RK",par=c(a_RK,b_RK,CV_REC_RK))
#'
#' # The following lines allow us to use the described function.
#' Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)
#' # We extract the F and SSB.
#' Sum.Pop.Mod(Pop.Mod,c("F","SSB"))
#' @export

Sum.Pop.Mod<-function(Pop.Mod,Elements){
  ctrFish=Pop.Mod$Info$ctrFish
  seed=Pop.Mod$Info$seed
  set.seed(seed)



  niter<-dim(Pop.Mod$Matrices$N)[3]

  #Matrices
  N<-Pop.Mod$Matrices$N
  F<-Pop.Mod$Matrices$F
  M<-Pop.Mod$Matrices$M
  W<-Pop.Mod$Matrices$W
  Mat<-Pop.Mod$Matrices$Mat
  C_W<-Pop.Mod$Matrices$C_W
  C_N<-Pop.Mod$Matrices$C_N
  f=Pop.Mod$Info$ctrFish$f
  # Bio
  a<-Pop.Mod$Info$ctrBio$a
  b<-Pop.Mod$Info$ctrBio$b
  min_age<-Pop.Mod$Info$minFage
  max_age<-Pop.Mod$Info$maxFage
  ages<-as.numeric(rownames(N))
  years<-as.numeric(colnames(N))
  number_years<-ncol(N)
  number_ages<-nrow(N)

  Z<-M+F

  # L

  L<-N

  for(i in 1:number_years){
    L[,i,]<-(W[,i,]/a)^(1/b)
  }


  # Lc


  W_c<-N

  for(i in 1:number_ages){
    for(j in 1:number_years){
      W_c[i,j,]<-C_W[i,j,]/C_N[i,j,]
    }
  }


  L_c<-N

  for(i in 1:number_years){
    L_c[,i,]<-(W_c[,i,]/a)^(1/b)
  }


  WMA<-N;WM<-N
  for(i in 1:number_ages){
    for(j in 1:number_years){
      WM[i,j,]<-N[i,j,]*W[i,j,]
      WMA[i,j,]<-WM[i,j,]*Mat[i,j,]
    }}

  column.names <- years
  row.names <- ""
  matrix.names <- 1:niter

  SSB<-array(rep(0,number_years), dim=c(1, number_years, niter),dimnames=list(row.names,column.names,matrix.names))

  for(j in 1:number_years){
    if(niter>1){ SSB[,j,]<-colSums(WMA[,j,])} else {SSB[,j,1]<-sum(WMA[,j,1])}
 }


  # SELECTIVITY: We need to generate stochastic values of a50_Sel (CV_SEL)

  s<-N
  for(kk in 1:niter){
    for(i in 1:number_ages){
      for (j in 1:number_years){
        s[i,j,kk]<-F[i,j,kk]/f[kk,j]
      }
    }
  }

  year_C_W<-array(rep(0,number_years), dim=c(1, number_years, niter),dimnames=list("",column.names,matrix.names))
  biomass<-array(rep(0,number_years), dim=c(1, number_years, niter),dimnames=list("",column.names,matrix.names))
  Recruiment<-array(rep(0,number_years), dim=c(1, number_years, niter),dimnames=list("",column.names,matrix.names))
  F_mean<-array(rep(0,number_years), dim=c(1, number_years, niter),dimnames=list("",column.names,matrix.names))

  for (ind in 1:niter){
  year_C_W[,,ind]<-colSums(C_W[,,ind])
  biomass[,,ind]<-colSums(WM[,,ind])
  Recruiment[,,ind]<-N[1,,ind]
  F_mean[,,ind]<-apply(F[(min_age-ages[1]+1):(max_age-ages[1]+1),,ind], 2, mean)}





  l_el<-length(Elements)
  el<-list()
  for (i in 1:l_el){
    if(Elements[i]=="Z"){el[[i]]<-Z}
    if(Elements[i]=="LS"){el[[i]]<-L}
    if(Elements[i]=="LC"){el[[i]]<-L_c}
    if(Elements[i]=="WS"){el[[i]]<-WM}

    if(Elements[i]=="C"){el[[i]]<-year_C_W}
    if(Elements[i]=="WSSB"){el[[i]]<-WMA}
    if(Elements[i]=="SEL"){el[[i]]<-s}

    if(Elements[i]=="BIO"){el[[i]]<-biomass}
    if(Elements[i]=="SSB"){el[[i]]<-SSB}
    if(Elements[i]=="REC"){el[[i]]<-Recruiment}
    if(Elements[i]=="F"){el[[i]]<-F_mean}
    if(Elements[i]=="WC"){el[[i]]<-W_c}
  }
  names(el)<-Elements
  return(el)
}
