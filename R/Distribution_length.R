#' @title Length distribution
#'
#' @description Return the stock length or captures length distribution for each year and iteration.
#'
#'
#' @param Pop.Mod A list containing the components returned by Population.Modeling function (main function).
#' @param CV The coefficient of variation associated to the log-normal distribution (see Details).
#' @param Type An indicator of which distribution length must be computed, length stock distribution (Type="LengthS") whereas length captures distribution (Type="LengthC").
#' @param scale A rescale parameter, the matrix N (Type="LengthS") or C (Type="LengthC") is divided by this value to avoid large times of computation. See details.
#' @return An array whose third dimension is the number of iterations, and the second one is the different years. Hence each column contains the distribution length (stock or captures) for each year.
#' @details The function returns the stochastic length distribution of the stock (Type="LengthS") or length captures distribution (Type="LengthC") for each year and iteration.
#' In the case of the stock length distribution it is computed generating for each age, year and iteration random values from a log-normal distribution centered in the corresponding stock length and whose variability comes from the given CV. The number of values generated in each case is determined by the population size N if scale=NULL, in other case the number of values is determined by N/scale. In this way the distribution can be approximated correctly without the need of generation of large sequences of values. For the captures length distribution the mean of the log-normal distribution is given by the corresponding capture length, and the number of random values is given by the corresponding number of captures, C if scale=NULL or by C/scale in other case.
#' @author
#' \itemize{
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
#' maxFage=7,ts=0,tc=0.5,tseed=NULL)
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
#' f=rep(0.5,number_years)
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
#' # We compute the captures length distribution:
#' L.D<-Distribution.length(Pop.Mod,CV=0.2,Type="LengthC",scale=NULL)
#' # We compute the stock length distribution:
#'  L.D<-Distribution.length(Pop.Mod,CV=0.2,Type="LengthS",scale=NULL)
#' @export




Distribution.length<-function(Pop.Mod,CV,Type,scale){
  if(CV==0){ stop("The value of CV must be strictly positive")}

  if(is.numeric(Pop.Mod$Info$seed)){set.seed(Pop.Mod$Info$seed)}
  N<-Pop.Mod$Matrices$N
  C<-Pop.Mod$Matrices$C_N
  if(is.numeric(scale)){N<-N/scale
  C<-C/scale}





  number_years<-dim(N)[2];number_ages<-dim(N)[1]
  Lm<-Sum.Pop.Mod(Pop.Mod,c("LS"))$LS
  L.cm<-Sum.Pop.Mod(Pop.Mod,c("LC"))$LC
  niter<-dim(N)[3]

  if(niter==1){
    if(Type=="LengthS"){
      Ls<-list()
      for(j in 1:number_years){
        Ls[[j]]<-0
      }
      for (j in 1:number_years) {
        L<-Lm[,j,1]

        for(i in 1:number_ages){

          v<-(CV*mean(L[i]))^2
          m<-mean(L[i])
          if (m>0){
            mu<-log(m^2/sqrt(m^2+v))}
          sigma<-sqrt(log(v/m^2+1))


          if(N[i,j,1]>1 &m>0&sigma>0){
            Ls[[j]]<-c(Ls[[j]],(stats::rlnorm((N[i,j,1]), meanlog =mu, sdlog =sigma)))}
        }}

      L_inf<-Pop.Mod$Info$ctrBio$L_inf
      max_length<-round(L_inf+(2*CV*L_inf))
      lengths<-0:max_length
      column.names <- colnames(N)
      RS<-array(0, dim=c(max_length+1, number_years, 1),dimnames=list(lengths,column.names,1))

      for (i in 1:number_years){
        aux<-c(Ls[[i]][-1],max_length)
        ind_aux<-which(aux>max_length);aux[ind_aux]<-max_length
        # trick +1 to count zero
        d<-stats::setNames(tabulate(floor(aux+1)), 0:max_length)
        # take_out artificial
        d[max_length+1]<-d[max_length+1]-1
        RS[,i,1]<-d
      }

      return(RS)}

    if(Type=="LengthC"){

      Lsc<-list()
      for(j in 1:number_years){
        Lsc[[j]]<-0
      }
      for (j in 1:number_years) {
        L.c<-L.cm[,j,1]
        for(i in 1:number_ages){

          v.c<-(CV*mean(L.c[i]))^2
          m.c<-mean(L.c[i])
          if (m.c>0){
            mu.c<-log(m.c^2/sqrt(m.c^2+v.c))}
          sigma.c<-sqrt(log(v.c/m.c^2+1))
          if(C[i,j,1]>1&m.c>0&sigma.c>0){
            Lsc[[j]]<-c(Lsc[[j]],(stats::rlnorm((C[i,j,1]), meanlog =mu.c, sdlog =sigma.c)))}
        }}

      L_inf<-Pop.Mod$Info$ctrBio$L_inf
      max_length<-round(L_inf+(2*CV*L_inf))
      lengths<-0:max_length
      column.names <- colnames(N)

      RC<-array(0, dim=c(max_length+1, number_years, 1),dimnames=list(lengths,column.names,1))

      for (i in 1:number_years){
        aux<-c(Lsc[[i]][-1],max_length)
        ind_aux<-which(aux>max_length);aux[ind_aux]<-max_length
        # trick +1 to count zero
        d<-stats::setNames(tabulate(floor(aux+1)), 0:max_length)
        # take_out artificial
        d[max_length+1]<-d[max_length+1]-1
        RC[,i,1]<-d
      }

      return(RC)}}




  if (niter>1){
  ### We need to compute the maximun posible length
  L_inf<-Pop.Mod$Info$ctrBio$L_inf
  max_length<-round(L_inf+(2*CV*L_inf))

  ### Vector of lengths
  lengths<-0:max_length
  column.names <- colnames(N)
  matrix.names <- 1:niter
  RLD<-array(0, dim=c(max_length+1, number_years, niter),dimnames=list(lengths,column.names,matrix.names))


  if(Type=="LengthS"){RS=f_RLD(Lm,N,number_ages,number_years,niter,column.names,max_length,CV,Pop.Mod)
  return(LD=RS)
  }

  if(Type=="LengthC"){RC=f_RLD(L.cm,C,number_ages,number_years,niter,column.names,max_length,CV,Pop.Mod)
  return(LD=RC)
  }}


    }


