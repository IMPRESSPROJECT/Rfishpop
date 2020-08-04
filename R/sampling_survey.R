#' @title Indices of abundance (for each year and age) and biomass (for each year), and sampling stock length (for each year)
#'
#' @description Returns the indices of abundance for each year, age and iteration, the indices of biomass for each year and iteration, and a sample of the stock length (for each year and iteration) derived from the corresponding length distribution.
#'
#'
#' @param Pop.Mod A list containing the components returned by Population.Modeling function (main function).
#' @param type An indicator of which index/sample must be computed: index of biomass (type="biomass"), index of abundance (type="abundance"), or stock length sample (type="LengthS").
#' @param q_A The matrix of annual and age specific catchability coefficients.
#' @param gamma The density dependent parameter.
#' @param CV_A The coefficient of variation associated to the abundance index. Default value 0.
#' @param par A vector containing the required parameters if type="LengthS". \itemize{
#' \item{CV The coefficient of variation associated to the log-normal distribution used to compute the stock length distribution from which the length sample is obtained (see Details).}
#' \item{RF.value The number of values generated for each age (given a year and an iteration) from the log-normal distribution used to obtain the stock length distribution from which the length sample is obtained  (see Details). By default RF.value=1000.}
#' \item{sample.size The size of the desired sample.}}
#' @param tsampling Time of the year at which the sampling is carried out. This parameter takes a value between 0 and 1, since the year is considered as a unit. By default tsampling=0.
#'
#'
#' @details The function returns the index of abundance for each year, age and iteration, the index of biomass for each year and iteration or a stock length sample generating random values from the computed length distribution.
#'
#' ABUNDANCE INDEX:
#'
#' The abundance index for year t and age i is \deqn{IA_it=q_A_it*N_it^{gamma}} where q_A_it is the catchability coefficient, gamma is the density dependent parameter and N_it is the abundance for year t and age i when CV_A=0. If CV_A is different than 0 the abundance index for year t and age i is generated from a log normal distribution center in \deqn{q_A_it*N_it^{gamma}} and variability determined for the corresponding CV_A.
#' Note that the matrix N (stock numbers) in this function is computed at the time instant tsampling, i.e., it is not the matrix N computed in the main function Population.Modeling function instead we update such matrix to time instant tsampling.
#'
#'
#' BIOMASS INDEX:
#'
#' Firstly, we compute the W (third dimensional array containing the weight corresponding to the stock length for each age, year and iteration) at time tsampling. Then, we define WS (third dimensional array containing the weight for each age, year and iteration) multiplying the abundance index IA by the weight W. Finally, the index of biomass for each year is the sum by columns of matrix WS.
#'
#' SAMPLE LENGTH:
#'
#' A stock length sample is obtained generating random values from the computed stock length distribution.
#'
#'  How the stock distribution length is computed?
#'  The stock distribution length from which the sampling procedure is carried out is computed as follows.
#'  It is computed generating for each age, year and iteration RF.value random values from a log-normal distribution centered in the corresponding stock length and whose variability comes from the given CV. The distribution obtained for each age (given a year and an iteration) is scaled using the corresponding abundance index matrix (IA).
#'
#'  Remember that the matrix LS (stock length) is computed at the time instant tsampling.
#'
#' @return An array containing the indices of abundance for each year, age, and iteration. If type="biomass" an array containing the indices of biomass for each year and iteration is also reported. If type="LengthS" an array containing in each column (year) the corresponding length sample for each iteration (third dimension) is also provided.
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' \item{Maria Grazia Pennino}
#' }
#' @examples
#' library(Rfishpop)
#' ctrPop<-list(years=seq(1980,2020,by=1),niter=2,N0=10000,ages=0:15,minFage=4,
#' maxFage=7,ts=0,tc=0.5,tseed=NULL)
#' number_ages<-length(ctrPop$ages);number_years<-length(ctrPop$years)
#' M<-matrix(rep(0.4,number_ages*number_years),ncol = number_years)
#' colnames(M)<-ctrPop$years
#' rownames(M)<-ctrPop$ages
#' ctrBio<-list(M=M,CV_M=0, L_inf=124.5, t0=-0.25, k=0.164, CV_L=0, CV_LC=0, a=4.5*10^(-6), b=3.1049,
#' a50_Mat=3, ad_Mat=-0.5,CV_Mat=0)
#' ctrSEL<-list(type="Logistic", par=list(a50_Sel=1.5, ad_Sel=-1),CV_SEL=0)
#' f=matrix(rep(0.5,number_years),ncol=number_years,nrow=2,byrow=TRUE)
#' ctrFish<-list(f=f,ctrSEL=ctrSEL)
#' a_RK=10; b_RK=0.0002; CV_REC_RK=0.2
#' SR<-list(type="RK",par=c(a_RK,b_RK,CV_REC_RK))
#' Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)
#' # Uncomments the following lines.
#' # Abundance
#' #q_A<-matrix(1,ncol=41,nrow=16);gamma<-1;CV_A<-0
#' #IA=Sampling_Survey(Pop.Mod=Pop.Mod,type="abundance",q_A=q_A,
#' #gamma=gamma,CV_A=CV_A,tsampling=0)
#' # Biomass
#' #IB=Sampling_Survey(Pop.Mod=Pop.Mod,type="biomass",q_A=q_A,
#' #gamma=gamma,CV_A=CV_A,tsampling=0)
#' #IB$biomass
#' #IB$abundance
#' # Length
#' #par=c(0.2,1000,3000)
#' #IL=Sampling_Survey(Pop.Mod=Pop.Mod,type="LengthS",q_A=q_A,
#' #gamma=gamma,CV_A=CV_A,par=par,tsampling=0)
#' #IL$length
#' #IL$abundance
#'
#' @export






Sampling_Survey<-function(Pop.Mod,type,q_A,gamma,CV_A=0,par,tsampling=0){

  if (Pop.Mod$Info$ts==tsampling){
      N<-Pop.Mod$Matrices$N
      W<-Pop.Mod$Matrices$W
      LS=Sum.Pop.Mod(Pop.Mod,c("LS"))$LS
      } else {
      New.matrices=Population.Modeling.tsampling(Pop.Mod,tsampling)
      N=New.matrices$N
      LS=New.matrices$LS
      number_years<-ncol(N)
      a<-Pop.Mod$Info$ctrBio$a
      b<-Pop.Mod$Info$ctrBio$b
      W<-N

      for(i in 1:number_years){
        W[,i,]<-Weight(LS[,i,],a,b)
      }

    }

  number_ages<-nrow(N)
  number_years<-ncol(N)
  niter<-dim(N)[3]
  con_list<-as.list(1:niter)

  for (i in 1:niter){
    con_list[[i]]<-list(N[,,i])
  }


  R_abun<-function(con_list,q_A,gamma,CV_A){
    N<-con_list[[1]]
    number_ages<-nrow(N)
    number_years<-ncol(N)
    CV<-CV_A

    Index_Abundance<-N

    if(CV==0){
      for (i in 1:number_ages){
        for (j in 1:number_years){
          Index_Abundance[i,j]<-(q_A[i,j])*(N[i,j]^(gamma))
        }
      }
    }

    if(CV>0){
      for (i in 1:number_ages){
        for (j in 1:number_years){
        v<-(CV*(q_A[i,j])*(N[i,j]^(gamma)))^2
        m<-(q_A[i,j])*(N[i,j]^(gamma))

        mu<-log(m^2/sqrt(m^2+v))
        sigma<-sqrt(log(v/m^2+1))

          Index_Abundance[i,j]<-stats::rlnorm(1, meanlog =mu, sdlog =sigma)
        }
      }
    }
    return(Index=Index_Abundance)
  }


    numCores <- parallel::detectCores()

    cl <- parallel::makeCluster(numCores)



    ret<-parallel::parLapply(cl,con_list,R_abun,q_A=q_A,gamma=gamma,CV_A=CV_A)

    parallel::stopCluster(cl)
    ret=array(as.numeric(unlist(ret)), dim=c(number_ages, number_years, niter))
    colnames(ret)<-colnames(N)
    rownames(ret)<-rownames(N)

    # SAVE ABUNDANCE INDEX
    N=ret

    if (type=="abundance"){
      return(N)
    }

    if (type=="biomass"){
     WM<-N
      for(i in 1:number_ages){
        for(j in 1:number_years){
          WM[i,j,]<-N[i,j,]*W[i,j,]

        }}

     years<-as.numeric(colnames(N))
     column.names <- years
     row.names <- ""
     matrix.names <- 1:niter
     biomass<-array(rep(0,number_years), dim=c(1, number_years, niter),dimnames=list("",column.names,matrix.names))

     for (ind in 1:niter){
       biomass[,,ind]<-colSums(WM[,,ind])
       }
      return(list(biomass=biomass,abundance=N))
    }

      if(type=="LengthS"){
      our.sample<-Sampling_length(Pop.Mod,CV=par[1],
                                  sample.size=par[3],RF.value=par[2],N,LS)
      return(list(length=our.sample,abundance=N))
      }


}

