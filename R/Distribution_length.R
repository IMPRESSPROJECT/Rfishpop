#' @title Stock Length and Catches Length Distribution for each year
#'
#' @description Return the stock length or catches length distribution for each year and iteration.
#'
#'
#' @param Pop.Mod A list containing the components returned by Population.Modeling function (main function).
#' @param CV The coefficient of variation associated to the log-normal distribution (see Details).
#' @param Type An indicator of which distribution length must be computed, length stock distribution (Type="LengthS") whereas length catches distribution (Type="LengthC").
#' @param RF.value The number of values generated for each age (given a year and an iteration) from the log-normal distribution (see details). By default RF.value=1000.
#' @return An array whose third dimension is the number of iterations, and the second one is the different years. Hence each column contains the distribution length (stock or catches) for each year.
#' @details The function returns the stochastic length distribution of the stock (Type="LengthS") or length catches distribution (Type="LengthC") for each year and iteration.
#' In the case of the stock length distribution it is computed generating for each age, year and iteration RF.value random values from a log-normal distribution centered in the corresponding stock length and whose variability comes from the given CV. For the catches length distribution the mean of the log-normal distribution is given by the corresponding catch length. For the stock length distribution the distribution obtain for each age (given a year and an iteration) is scaled using the corresponding stock number (N matrix), whereas in the catch distribution this role is for the catch matrix C.
#' @author
#'  \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' }
#' @examples
#'
#' # First we introduce the basic parameters to define the population.
#' # Note that N0 is equal to 10000 individuals, and hence below we are
#' # consistent with this unit when we introduce the biological and
#' # stock-recruitment parameters.
#' ctrPop<-list(years=seq(1980,2020,by=1),niter=2,N0=10000,ages=0:15,minFage=4,
#' maxFage=7,tc=0.5,tseed=NULL)
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
#'
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
#'
#'
#'
#'
#'
#' # We compute the catches length distribution:
#' # UNCOMMENT THE FOLLOWING LINES
#' #L.D<-Distribution.length(Pop.Mod,CV=0.2,Type="LengthC")
#' # We compute the stock length distribution:
#' # L.D<-Distribution.length(Pop.Mod,CV=0.2,Type="LengthS")
#' @export



Distribution.length<-function(Pop.Mod,CV,Type,RF.value=1000){
  message("It is working, please give us a while :)")


  a=list()
  N<-Pop.Mod$Matrices$N;niter<-dim(N)[3];number_years<-dim(N)[2]
  Nc=N
  C_Nc=Pop.Mod$Matrices$C_N
  LSc<-Sum.Pop.Mod(Pop.Mod,c("LS"))$LS
  LCc=Sum.Pop.Mod(Pop.Mod,c("LC"))$LC

  for (i in 1:niter){
    seed=Pop.Mod$Info$seed
    N<-Nc[,,i]
    C_N<-C_Nc[,,i]
    LS<-LSc[,,i]
    LC<-LCc[,,i]
    L_inf<-Pop.Mod$Info$ctrBio$L_inf

    a[[i]]=list(seed=seed,N=N,C_N=C_N,LS=LS,LC=LC,L_inf=L_inf)
  }

  if(Type=="LengthS"){
    for (i in 1:niter){
      seed=Pop.Mod$Info$seed
      N<-Nc[,,i]

      LS<-LSc[,,i]

      L_inf<-Pop.Mod$Info$ctrBio$L_inf

      a[[i]]=list(seed=seed,N=N,LS=LS,L_inf=L_inf)
    }

  }
  if(Type=="LengthC"){
    for (i in 1:niter){
      seed=Pop.Mod$Info$seed

      C_N<-C_Nc[,,i]

      LC<-LCc[,,i]
      L_inf<-Pop.Mod$Info$ctrBio$L_inf

      a[[i]]=list(seed=seed,N=C_N,LS=LC,L_inf=L_inf)
    }

  }


  numCores <- parallel::detectCores()

  cl <- parallel::makeCluster(numCores)

  parallel::clusterExport(cl,c("Distribution.length.aux"),
                          envir=environment())

  ret<-parallel::parLapply(cl, a,Distribution.length.aux, CV,RF.value)
  parallel::stopCluster(cl)


  matrix.names <- 1:niter
  max_length<-round(L_inf+(2*CV*L_inf))
  lengths<-0:max_length
  column.names <- colnames(N)
  res<-array(0, dim=c(max_length+1, number_years, niter),dimnames=list(lengths,column.names,matrix.names))

  for(i in 1:niter){
    res[,,i]=ret[[i]][,,1]
  }


return(res)



}




