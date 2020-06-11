#' @title Index of biomass (for each year) and abundance (for each year and age)
#'
#' @description Returns the indices of abundance for each year, age and iteration, and indices of biomass for each year and iteration.
#'
#'
#' @param Pop.Mod A list containing the components returned by Population.Modeling function (main function).
#' @param type indices type which can be "biomass" or "abundance".
#' @param par list of the parameters required of computed the selected index. \itemize{
#' \item  type="biomass": \itemize{
#' \item q_B which is the vector of annual catchability coefficients.
#' \item gamma which is the density dependent parameter.
#' \item CV which is the biomass coefficient of variation. Default value 0.}
#' \item  type="abundance": \itemize{
#' \item q_A which is the matrix of annual and age specific catchability coefficients.
#' \item gamma which is the density dependent parameter.
#' \item CV which is the biomass coefficient of variation. Default value 0.
#' }}
#'
#' @details The function returns the index of abundance for each year, age and iteration, and the index of biomass for each year and iteration. The biomass index for year t is \deqn{IB_t=q_B_t*BIO_t^{gamma}} where q_B_t is the catchability coefficient for year t and BIO_t is the biomass for year t when CV=0. If CV is different than 0 the biomass index for year t is \deqn{IB_t=q_B_t*BIO_t^{gamma}*epsilon_t} where q_B_t is the catchability coefficient, BIO_t is the biomass for year t, and epsilon_t is the residual generated from a log normal distribution center in zero and whose variability determined for the corresponding CV. The abundance index for year t and age i is \deqn{IA_it=q_A_it*N_it^{gamma}} where q_A_it is the catchability coefficient and N_it is the abundance for year t and age i when CV=0. If CV is different than 0 the abundance index for year t and age i is \deqn{IA_it=q_A_it*N_it^{gamma}*epsilon_it} where q_A_it is the catchability coefficient, N_it is the abundance for year t and age i, and epsilon_it is the residual generated from a log normal distribution center in zero and variability determined for the corresponding CV.
#'
#' @return An array containing the indices of abundance for each year, age, and iteration if "type=abundance" or the indices of biomass for each year and iteration if "type=biomass".
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
#'
#' # Now,we can compute the index of abundance or biomass.
#' # For biomass index:
#' q_B<-rep(0.01,41);gamma<-1;CV<-0.2; par<-list(q_B,gamma,CV)
#'  #I<-Sampling_Survey(Pop.Mod,type="biomass",par=par)
#'
#' # For abundance index:
#' q_A<-matrix(0.2,ncol=41,nrow=16);gamma<-1;CV<-0.2; par<-list(q_A,gamma,CV)
#'  #I<-Sampling_Survey(Pop.Mod,type="abundance",par=par)

#' @export






Sampling_Survey<-function(Pop.Mod,type,par){

  N<-Pop.Mod$Matrices$N
  number_ages<-nrow(N)
  number_years<-ncol(N)
  biomass<-Sum.Pop.Mod(Pop.Mod,c("BIO"))$BIO
  niter<-dim(Pop.Mod$Matrices$N)[3]
  con_list<-as.list(1:niter)

  for (i in 1:niter){
    con_list[[i]]<-list(N[,,i],biomass[,,i])
  }

  ####################################
  R_bio<-function(con_list,par){
    N<-con_list[[1]]
    biomass<-con_list[[2]]
    number_ages<-nrow(N)
    number_years<-ncol(N)
    q_B<-par[[1]]
    gamma<-par[[2]]
    CV<-par[[3]]
    if(is.null(CV)){CV=0}
    years<-colnames(N)

    Index_Biomass<-(q_B)*(biomass^gamma)
    names(Index_Biomass)<-years
    if(CV>0){
      v<-(CV*mean(biomass))^2
      m<-mean(biomass)

      mu<-0
      sigma<-sqrt(log(v/m^2+1))

      residuals<-stats::rlnorm(number_years, meanlog =mu, sdlog =sigma)
      Index_Biomass<-((q_B)*(biomass^gamma)*residuals)
      names(Index_Biomass)<-years
    }
    return(Index=Index_Biomass)}



  if (type=="biomass"){
    numCores <- parallel::detectCores()

    cl <- parallel::makeCluster(numCores)



    ret<-parallel::parLapply(cl,con_list,R_bio,par=par)

    parallel::stopCluster(cl)
    ret=array(as.numeric(unlist(ret)), dim=c(1, number_years, niter))
    colnames(ret)<-colnames(N)
  }
  ###################################

  R_abun<-function(con_list,par){
    N<-con_list[[1]]
    biomass<-con_list[[2]]
    number_ages<-nrow(N)
    number_years<-ncol(N)
    q_A<-par[[1]]
    gamma<-par[[2]]
    CV<-par[[3]]

    Index_Abundance<-N
    mean_N<-apply(N,1,mean)
    if(CV==0){
      for (i in 1:number_ages){
        for (j in 1:number_years){
          Index_Abundance[i,j]<-(q_A[i,j])*(N[i,j]^(gamma))
        }
      }
    }

    if(CV>0){
      for (i in 1:number_ages){
        v<-(CV*mean(N[i,]))^2
        m<-mean(N[i,])

        mu<-0
        sigma<-sqrt(log(v/m^2+1))

        residuals<-stats::rlnorm(number_years, meanlog =mu, sdlog =sigma)
        for (j in 1:number_years){
          Index_Abundance[i,j]<-(residuals[j]*(q_A[i,j])*(N[i,j]^(gamma)))
        }
      }
    }
    return(Index=Index_Abundance)
  }
  if (type=="abundance"){
    numCores <- parallel::detectCores()

    cl <- parallel::makeCluster(numCores)



    ret<-parallel::parLapply(cl,con_list,R_abun,par=par)

    parallel::stopCluster(cl)
    ret=array(as.numeric(unlist(ret)), dim=c(number_ages, number_years, niter))
    colnames(ret)<-colnames(N)
    rownames(ret)<-rownames(N)
  }


  return(ret)
}

