#' @title Modeling an Exploited Population (Structured by Age)
#'
#' @description Provides a flexible and generic operating model (OM) which simulates the real dynamics of the fishery system. The OM is formed by biological, fishery and control components. The stock is described as age structured population along the time.
#'
#'@param ctrPop A list containg the following elements:  \itemize{
#'   \item years=vector containing the years for which the dynamics of the fishery system is simulated.
#'   \item ages=vector containing the different ages presented in the stock. The oldest age defines a plus group. The first age must be 0, in the current version the OM is not adapted to treat with populations started in other ages.
#'   \item niter=number of iterations of the simulation process (see Details for more information about the stochastic and deterministic performance).
#'   \item N0=a number corresponding to the population size at first age and year, or a vector containing the population size at first year for all ages. It is worth to mention that the user needs to know in which units is the number introduced here. For example, it can refers to individuals or thousands of individuals.
#'   \item minFage=minimum age for which the corresponding fishing mortality is considered to compute the mean fishing mortality.
#'   \item maxFage=maximum age for which the corresponding fishing mortality is considered to compute the mean fishing mortality.
#'   \item tc=time of the year at which the catches are simulated. This parameter takes a value between 0 and 1, since the year is considered as a unit. By default  tc=0.5, the catches occurs at mid of the year (this is an approximation since the catches can occur during a period of the year not in a particular moment).
#'   \item seed=a numeric value to introduce into set.seed() function for having a reproducible result. By default is NULL, which means that the results are not  reproducible, that is, each run returns different results.}
#'@param ctrBio A list containg the following biological parameters of the population: \itemize{
#'   \item M=matrix containing the rates of instantaneous natural mortality for each year and age.
#'   \item CV_M=coefficient of variation associated to the natural mortality. In each stochastic iteration the rates of instantaneous natural mortality come from a log-normal distribution centered on the corresponding value of M and variability determined by CV_M.
#'   \item L_inf=asymptotic average maximum body size (parameter of Von Bertalanffy Growth Model, see Details). It is important to take care of the units of the introduced number, which must be in accordance with others parameters. The unit can be, for example, cm.
#'   \item t0=hypothetical age at which the species has zero length (parameter of Von Bertalanffy Growth Model, see Details).
#'   \item k=growth rate coefficient that determines how quickly the maximum is attained (parameter of Von Bertalanffy Growth Model, see Details)
#'   \item CV_L=coefficient of variation associated to the stock length.  In each stochastic iteration the length comes from a normal distribution centered on the corresponding length obtained from Von Bertalanffy Growth Model (for each age at ts=0 the time of the year at which the stock is simulated which has been fixed in 0 by default) and the corresponding variability is determined by CV_L.
#'   \item CV_LC=coefficient of variation associated to the captures length. In each stochastic iteration the length comes from a normal distribution centered on the corresponding length obtained from Von Bertalanffy Growth Model (for each age at tc) and the corresponding variability is determined by CV_LC.
#'   \item a=allometric growth parameter. It is important to introduce this parameter in a unit according to the units used to measure the length and weight.
#'   \item b=scaling constant of allometric growth. It is important to introduce this parameter in a unit according to the units used to measure the length and weight.
#'   \item a50_Mat=x-value of the sigmoid's midpoint of the logistic function used to generate the maturity matrix.
#'   \item ad_Mat=minus the inverse of the logistic growth rate (steepness of the curve) used to generate the maturity matrix. Zero or positive are NOT valid.
#'   \item CV_Mat=coefficient of variation associated to the a50_Mat parameter. In each stochastic iteration the maturity matrix values comes from a logistic function whose a50_Mat parameter is generated from a log-normal distribution centered on the given value of a50_Mat and variability detrmined by CV_Mat.}
#'@param ctrFish A list containg the following fishing parameters of the population: \itemize{
#' \item f=is the annual component of fishing mortality F = f * SEL. Can be different for each of the iterations, hence we introduce a matrix whose rows contain the annual vector for each iteration.
#' \item ctrSEL= list of two objects, "type" which specifies the selectivity function considered and "par" which contains the corresponding parameters. Below the different selectivity functions are described.\itemize{
#' \item type="cte", a constant selectivity function is used, which means that there is no dependence on the age of the the prey and the probability to be captured. In this case the element "par" contains:\itemize{
#'\item cte= the constant probability value (cte is a value in [0,1]).
#' }
#' \item type="Andersen", the Andersen selectivity function that is dependent on the ratio of a parameter to the prey age is considered, and the element "par" contains the values of its parameters (see Details):\itemize{
#' \item p1=define the value of \eqn{x} at which the transition between the two gaussian functions happens; that is, \eqn{x=p5/exp(p1)}. Negative values are NOT valid.
#' \item p3=ascending slope parameter. Negative values or zero are NOT valid.
#' \item p4=descending slope parameter. Negative values or zero are NOT valid.
#' \item p5=define the value of \eqn{x} at which the transition between the two gaussian functions happens; that is, \eqn{x=p5/exp(p1)}. Negative values or zero are NOT valid.
#' }
#' \item type="Gamma", the Gamma selectivity function which only depends on the prey age is considered, the element "par" contains the values of (see Details)\itemize{
#'  \item gamma=the size of the mesh used by the fleet (according to the units of the remaining parameters, for example, cm). Negative values or zero are NOT valid.
#'  \item beta=the scale parameter. Negative values or zero are NOT valid.
#'  \item alpha=the shape parameter. Negative values are NOT valid.
#'  }
#' \item type="Logistic", the Logistic selectivity function based on ages is considered (see Details), the element "par" contains the parameters a50_Sel and ad_Sel. \itemize{
#' \item a50_Sel=age at which the captured proportion is equal to 50 percent the x-value of the sigmoid's midpoint of the logistic function.
#' \item ad_Sel=minus the inverse of the logistic growth rate (steepness of the curve). Zero or positive are NOT valid.
#'}}
#' \item CV_SEL=coefficient of variation associated to the selectivity. \itemize{
#' \item for type="cte", the values of the selectivity matrix for each stochastic iteration comes from a uniform distribution with mean equal to "cte" and whose variability is determined by CV_SEL.
#' \item for type="Logistic", the value of the parameter a50_SEL in each stochastic iteration comes from a normal distribution centered on the given value of such parameter, and whose variability is determined by CV_SEL.
#' \item for type="Andersen", the value of the parameter p1 in each stochastic iteration comes from a normal distribution centered on the given value of such parameter, and whose variability is determined by CV_SEL.
#' \item for type="Gamma", the value of the parameter alpha in each stochastic iteration comes from a normal distribution centered on the given value of such parameter, and whose variability is determined by CV_SEL.
#' }}
#'@param SR A list containing the following parameters of spawning stock recruitment relationship: \itemize{
#' \item type=stock recruitment model options. Three choices are available: "cte" which means that a constant recruitment is used, "BH" which corresponds to Beverton-Holt Recruitment Model and "RK" which refers to Ricker Recruitment Model. The default option is "cte".
#' \item par= parameters of the stock recruitment model specified in the argument "type". In order to introduce the following parameters correctly we need to know which are the units of biomass and  N0 to be consistent.  \itemize{
#'   \item type="cte" then recruitment is equal to the population size at first age and year introduced previously in N0 parameter then par contains only CV_REC_C which is the associated coefficient of variation. For the stochastic iterations the recruitment is generated from a  log-normal distribution centered on N0 and whose variability is determined by CV_REC_C.
#'   \item type="BH" then par is equal to a vector containing the following elements: \itemize{
#'     \item a_BH=maximum number of recruitments produced (parameter of Beverton-Holt Recruitment Model, see Details).
#'     \item b_BH=spawning stock needed to produce recruitment equal to half maximum (parameter of Beverton-Holt Recruitment Model, see Details).
#'     \item CV_REC_BH=coefficient of variation associated to the Beverton-Holt Recruitment Model. In each stochastic iteration the deterministic equation of the model is multipled by log-normal residuals centered on 0 and whose variability is determined by CV_REC_BH.}
#'   \item type="RK" then par is equal to a vector containing the following elements: \itemize{
#'      \item a_RK=recruits-per-spawner at low stock levels (parameter of Ricker Recruitment Model, see Details).
#'      \item b_RK=relates to the rate of decreasing of recruits-per-spawner as SSB increases (parameter of Ricker Recruitment Model, see Details).
#'      \item CV_REC_RK=coefficient of variation associated to the Ricker Recruitment Model. In each  stochastic iteration the deterministic equation of the model is multipled by log-normal residuals centered on 0 and whose variability is determined by CV_REC_RK.}
#'      }}
#' @details STOCHASTIC AND DETERMINISTIC PERFORMANCE
#'
#'
#' The first iteration (niter=1) contains the results corresponding to the deterministic case, to wit, all the coefficients of variation (CV's) associated to the different biological and fishery components are zero. The next iterations contain the stochastic results.
#' Hence, if niter=1 (CV's=0), we only obtain the deterministic performance.
#' If niter=1 and  the CV's are different than 0, the function returns ERROR because the first iteration is the deterministic one and hence the CV's can not be used.
#' If niter>1, the CV's are 0 and the matrix of the annual component of fishing mortality has constant rows, the function returns ERROR because it does has sense to repeat a deterministic process (a process without variability).
#'
#' VON BERTALANFFY GROWTH MODEL
#'
#' Von Bertalanffy Growth Model equation is
#' \deqn{L(x)=L_inf*(1-exp(-k*(x-t0)))} where L_inf is the asymptotic average maximum body size, t0 is hypothetical age at which the species has zero length, and k is growth rate coefficient and x is the vector of ages where the function must be computed.
#'
#' LOGISTIC FUNCTION
#'
#' The logistic function equation is \deqn{L(x)<-1/(1+exp((x-x50)/xd))} where x50 is the x-value of the sigmoid's midpoint of the logistic function, and xd is minus the inverse of the logistic growth rate (steepness of the curve).
#'
#' SELECTIVITY FUNCTIONS
#'
#' Andersen selectivity function is \deqn{SA(x)<-p0+p2exp(-(ln(p5/x)-p1)^2/p4) if ln(p5/x)<=p1} and \deqn{SA(x)<-p0+p2exp(-(ln(p5/x)-p1)^2/p3) if ln(p5/x)>p1.} We fixed \eqn{p0=0}, which is the beginning size for the plateau, and \eqn{p2=1}, which is the maximum value attainable by the function. Remember that \eqn{p3} and \eqn{p4} are the ascending and descending slope parameters, respectively, whereas \eqn{p1} and \eqn{p5} define the value of \eqn{x} at which the transition between the two gaussian functions happens; that is, \eqn{x=p5/exp(p1)}.
#'
#' Gamma selectivity function is \deqn{SG(x)<-((x/((alpha-1) beta gamma))^(alpha-1))exp(alpha-1-(1/beta gamma)),} where gamma is the size of the mesh, alpha is the shape parameter and beta is the scale parameter.
#'
#' RECRUITMENT MODELS
#'
#' The formulation of the Beverton-Holt Recruitment Model equation is: \deqn{R=a_BH*SSB/(b_BH+SSB)} where SSB is the maturity biomass (spawning stock), a_BH is the maximum number of recruitments produced and b_BH is the spawning stock needed to produce recruitment equal to half maximum.
#'
#' The formulation of the Ricker equation is: \deqn{R=a_RK*SSB*exp(-b_RK*SSB)} where SSB is the maturity biomass (spawning stock), a_RK is the recruits-per-spawner at low stock levels and b_RK relates to the rate of decreasing of recruits-per-spawner as SSB increases.
#'
#' STOCK NUMBERS
#'
#' If N0 is a number then its is updated using the recruitment model evaluated at the SSB of the first year (such SSB assume negligible the contribution of age 0). More precisely, the maturity at first age (age 0 years) is assumed 0. Whereas if N0 is a vector such procedure of updated is not carried out.
#' Remember that: N0 is a number corresponding to the population size at first age and year, or a vector containing the population size at first year for all ages.
#'
#' Note that the time of the year at which the stock is simulated has been fixed in 0 by default, which means that the stock numbers correspond to January. No possibility of changing such default value is provided.
#' @return A list containing the following components.
#' \item{Matrices:}{Arrays of matrices containing the population values.}\itemize{
#' \item{N: }{Third dimensional array containing the population size for each age, year and iteration.}
#' \item{F:}{Third dimensional array containing the instantaneous fishing mortality for each age, year and iteration.}
#' \item{M:}{Third dimensional array containing the instantaneous natural mortality for each age, year and iteration.}
#' \item{W:}{Third dimensional array containing the weight corresponding to the stock length for each age, year and iteration (at ts the time of the year at which the stock is simulated which has been fixed in 0 by default).}
#' \item{Mat:}{Third dimensional array containing the proportion of mature at each age, year and iteration.}
#' \item{C_N:}{Third dimensional arraycontaining the number of captures for each age, year and iteration.}
#' \item{C_W:}{Third dimensional array containing the weight of captures for each age, year and iteration (at tc).}}
#' \item{Info:}{The information used to create the population.}\itemize{
#' \item{ctrFish:}{Argument of the function containing fishing information.}
#' \item{ctrBio:}{Argument of the function containing biological information.}
#' \item{SR:}{Argument of the function containing stock-recruitment relationship.}
#' \item{minFage:}{minimum age for which the corresponding fishing mortality is considered to compute the mean fishing mortality.}
#' \item{maxFage:}{maximum age for which the corresponding fishing mortality is considered to compute the mean fishing mortality.}
#' \item{ts:}{time of the year at which the stock is simulated. It has been fixed in 0 by default. No possibility of changing such default value is provided.}
#'  \item{tc:}{time of the year at which the catches are simulated. This parameter takes a value between 0 and 1, since the year is considered as a unit. By default is 0.5 half of the year.}
#' \item{seed:}{a numeric value to introduce into set.seed() function for having a reproducible result. By default is NULL, which means that the results are not  reproducible, that is, each run returns different results.}}
#'
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' \item{Maria Grazia Pennino}
#' }
#'
#' @examples
#' # First we introduce the basic parameters to define the population.
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
#' # To access to the matrices:
#' Pop.Mod$Matrices
#' # To access to the information
#' Pop.Mod$Info
#' @export

Population.Modeling<-function(ctrPop,ctrBio,ctrFish,SR){
  CV_REC_C=NULL; CV_REC_BH=NULL; CV_REC_RK=NULL
  ts<-0 # We have fixed it to avoid complications
  tc<-ctrPop$tc

  if(is.null(tc)){tc=0.5}

  seed=ctrPop$seed
  if(is.numeric(seed)){
  set.seed(seed)}

  Sel_type=ctrFish$ctrSEL$type

  type<-SR$type
  if(type=="BH"){CV_REC_BH<-SR$par[3]}
  if(type=="RK"){CV_REC_RK<-SR$par[3]}

  if(is.null(type)){type="cte"}
  if(type=="cte"){CV_REC_C<-SR$par[1]}

  CV_LC<-ctrBio$CV_LC
  CV_SEL<-ctrFish$ctrSEL$CV_SEL
  CV_L<-ctrBio$CV_L
  CV_M<-ctrBio$CV_M
  CV_Mat<-ctrBio$CV_Mat

  niter<-ctrPop$niter



  ### Establishing default values and checkings
  if(is.null(niter)){niter=1}
  if(is.null(CV_SEL)){CV_SEL=0}
  if(is.null(CV_M)){CV_M=0}
  if(is.null(CV_L)){CV_L=0}

  if(is.null(CV_Mat)){CV_Mat=0}
  if(is.null(CV_LC)){CV_LC=0}
  if(is.null(CV_REC_BH)){CV_REC_BH=0}
  if(is.null(CV_REC_RK)){CV_REC_RK=0}
  if(is.null(CV_REC_C)){CV_REC_C=0}

  check.sum<-sum(c(CV_SEL,CV_M,CV_L,CV_Mat,CV_LC,CV_REC_BH,CV_REC_RK,CV_REC_C))
  if(check.sum>0 & niter==1){ stop("niter is equal to 1 whereas the coefficients of variation are not zero. Note that the first iteration is the deterministic one, hence you need al least one iteration more to use the coefficients of variation")}

  f<-ctrFish$f
  nf=dim(f)[1]
  ll=length(apply(f,2,unique))
  if(check.sum==0 & niter>1 & ll==number_years){ stop("The coefficients of variation are zero whereas the number of iterations (niter) is greater than 1, and the matrix of the annual component of fishing mortality has constant rows")}

  if(niter!=nf){stop("The annual component of fishing mortality must be a matrix whose rows contain the annual vector for each iteration")}


  ### Taking out the values of the list
  years<-ctrPop$years
  ages<-ctrPop$ages
  if(ages[1]!=0){stop("The first age must be 0 as we explained in the help page")}

  N0<-ctrPop$N0
  a50_Sel<-ctrFish$a50_Sel
  ad_Sel<-ctrFish$ad_Sel
  f<-ctrFish$f
  M<-ctrBio$M
  L_inf<-ctrBio$L_inf
  t0<-ctrBio$t0
  k<-ctrBio$k
  a<-ctrBio$a
  b<-ctrBio$b
  a50_Mat<-ctrBio$a50_Mat
  ad_Mat<-ctrBio$ad_Mat
  min_age<-ctrPop$minFage
  max_age<-ctrPop$maxFage
  type<-SR$type
  if (type=="BH"){
    a_BH<-SR$par[1]
    b_BH<-SR$par[2]
  }

  if (type=="RK"){
    a_RK<-SR$par[1]
    b_RK<-SR$par[2]
  }
  ### New names

  number_years<-length(years)
  number_ages<-length(ages)

  lN0<-length(N0)
  if(lN0!=1 & lN0!=number_ages){
    stop("ERROR N0 must be a number or a vector of dimension number of ages")
  }



  column.names <- years
  row.names <- ages
  matrix.names <- 1:niter

  N<-array(rep(0,number_ages*number_years), dim=c(number_ages, number_years, niter),dimnames=list(row.names,column.names,
                                                                                                  matrix.names))
  # SELECTIVITY: We need to generate stochastic values of a50_Sel (CV_SEL)

  if(Sel_type=="Logistic"){

  a50_Sel<-ctrFish$ctrSEL$par$a50_Sel
  ad_Sel<-ctrFish$ctrSEL$par$ad_Sel

  ### Deterministic
  sd<-matrix(rep(Logistic(x=ages,x50=a50_Sel,xd=ad_Sel),number_years),ncol=number_years)
  colnames(sd)<-years
  rownames(sd)<-ages

  ### Stochastic (Log normal distribution)
  s<-N

  if(CV_SEL>0){

    for(j in 1:number_years){
      s[,j,-1]=stochastic_logistic_SEL_1(a50_Sel,ad_Sel,CV_SEL,niter,s,ages)
    }

  s[,,1]<-sd
  }}

  if(Sel_type=="cte"){

    ### Deterministic
    cte<-ctrFish$ctrSEL$par$cte
    sd<-matrix(rep(cte,number_years*number_ages),ncol=number_years)
    colnames(sd)<-years
    rownames(sd)<-ages

    ### Stochastic (Uniform)
    s<-N

    if(CV_SEL>0){
      for(j in 1:number_years){
        s[,j,-1]=stochastic_cte_SEL_1(cte,CV_SEL,niter,s,number_ages)
      }


      s[,,1]<-sd
    }}

  if(Sel_type=="Andersen"){
    ### Deterministic
    p1<-ctrFish$ctrSEL$par$p1;p3<-ctrFish$ctrSEL$par$p3;p4<-ctrFish$ctrSEL$par$p4;p5<-ctrFish$ctrSEL$par$p5
    sd<-matrix(rep(andersen(x=ages,p1=p1,p3=p3,p4=p4,p5=p5),number_years),ncol=number_years)
    colnames(sd)<-years
    rownames(sd)<-ages

    ### Stochastic (Uniform)
    s<-N

    if(CV_SEL>0){

      for(j in 1:number_years){
        s[,j,-1]=stochastic_andersen_SEL_1(p1=p1,p3=p3,p4=p4,p5=p5,CV_SEL,niter,s,ages)

      }
      s[,,1]<-sd
    }
  }

  if(Sel_type=="Gamma"){

    alpha<-ctrFish$ctrSEL$par$alpha
    gamma<-ctrFish$ctrSEL$par$gamma
    beta<-ctrFish$ctrSEL$par$beta


    ### Deterministic
    sd<-matrix(rep(gamma_SEL(x=ages,alpha=alpha,gamma=gamma,beta=beta),number_years),ncol=number_years)
    colnames(sd)<-years
    rownames(sd)<-ages

    ### Stochastic (Log normal distribution)
    s<-N

    if(CV_SEL>0){
      for(j in 1:number_years){
        s[,j,-1]=stochastic_gamma_SEL_1(alpha,beta,gamma,CV_SEL,niter,s,ages)

      }
      s[,,1]<-sd
    }}

  if(niter==1 & CV_SEL==0){s[,,1]<-sd}
  if(niter>1 & CV_SEL==0) {s[,,1:niter]<-sd}

  ### FISHING MORTALITY



  F<-N
  for(kk in 1:niter){
  for(i in 1:number_ages){
    for (j in 1:number_years){
      F[i,j,kk]<-s[i,j,kk]*f[kk,j]
    }
  }
  }
  ### NATURAL MORTALITY
  Md<-M


  ### Stochastic (Log normal distribution; CV_M)
  M<-N

  if(CV_M>0){
  for (i in 1:number_ages){
    for(j in 1:number_years){
    m<-Md[i,j]
    v<-(CV_M*m)^2
    mu<-log(m^2/sqrt(m^2+v))
    sigma<-sqrt(log(v/m^2+1))
    M[i,j,]<-stats::rlnorm(niter,mu,sigma)
    }}
    M[,,1]<- Md
  }
  ### Deterministic
  if(niter==1 & CV_M==0){M[,,1]<- Md}
  if(niter>1 & CV_M==0) {M[,,1:niter]<- Md}


  ### MORTALITY

  Z<-M+F

  if (length(N0)==1){
    N[1,1,]<-N0
    for (i in 2:(number_ages-1)){
      N[i,1,]<-(exp(-Z[i-1,1,]))*N[i-1,1,]
    }
    i=number_ages
    N[i,1,]<-((exp(-Z[i-1,1,]))*N[i-1,1,])/(1-exp(-Z[i,1,]))
    }

  if (length(N0)>1){
    N[,1,]<-N0
  }

  ### LENGTH
  L<-N

  Ld<-matrix(rep(Length_VB(L_inf,k,ages+ts,t0),number_years), ncol=number_years,nrow=number_ages)


  ### Stochastic (Normal CV_L)


  if(CV_L>0){
  for (i in 1:number_ages){
    for(j in 1:number_years){
      m<-Ld[i,j]
      v<-(CV_L*m)
      L[i,j,]<-stats::rnorm(niter,m,v)

    }}
    L[,,1]<-Ld
  }
  ### Deterministic
  if(niter==1 & CV_L==0){L[,,1]<-Ld}
  if(niter>1 & CV_L==0) {L[,,1:niter]<-Ld}


  ### WEIGHTS

  W<-N

  for(i in 1:number_years){
    W[,i,]<-Weight(L[,i,],a,b)
  }



  WM<-N

  for(i in 1:number_ages){
    WM[i,1,]<-N[i,1,]*W[i,1,]
  }



  ### MATURITY (Log normal distribution CV_Mat)
  Mat<-N
  ### Deterministic
  Matd<-matrix(rep(Logistic(x=ages,x50=a50_Mat,xd=ad_Mat),number_years),ncol = number_years)

  if(CV_Mat>0){
  m<-a50_Mat
  v<-(CV_Mat*m)^2
  mu<-log(m^2/sqrt(m^2+v))
  sigma<-sqrt(log(v/m^2+1))
  a50_Mat_niter<-matrix(stats::rlnorm((niter-1)*number_years,mu,sigma),ncol=number_years)



  for(j in 1:number_years){
  for (ind in 2:niter){
   Mat[,j,ind]<-Logistic(x=ages,x50=a50_Mat_niter[ind-1,j],xd=ad_Mat)
  }}
  Mat[,,1]<-Matd
  }


  if(niter==1 & CV_Mat==0){Mat[,,1]<-Matd}
  if(niter>1 & CV_Mat==0) {Mat[,,1:niter]<-Matd}



  ### We assume for now that the maturity at age 0 is 0. Then we make the required
  ### changes to have effectively 0 proportion of maturity at age 0.

  Mat[1,,]<-rep(0, number_years)


  WMA<-N

  for(i in 1:number_ages){
    WMA[i,1,]<-WM[i,1,]*Mat[i,1,]
  }


  row.names <- ""

  ### SSB

  SSB<-array(rep(0,number_years), dim=c(1, number_years, niter),dimnames=list(row.names,column.names,matrix.names))

  if (length(N0)==1){



  if(niter>1){ SSB[,1,]<-colSums(WMA[-1,1,])} else {SSB[,1,1]<-sum(WMA[-1,1,1])}
  if(type=="cte"){if(CV_REC_C>0){N[1,1,-1]=Log.normal(rep(N[1,1,1],niter-1),CV_REC_C)}
    }
  if(type=="BH"){
    if(CV_REC_BH>0){N[1,1,-1]<-Log.normal(RBH(SSB[,1,-1],a_BH,b_BH),CV_REC_BH)#Beverton-Holt Recruitment Model
    N[1,1,1]<-RBH(SSB[,1,1],a_BH,b_BH)}
    }
  if(type=="RK"){
    if(CV_REC_RK>0){N[1,1,-1]<-Log.normal(RRK(SSB[,1,-1],a_RK,b_RK),CV_REC_RK)# Ricker Recruitment Model
    N[1,1,1]<-RRK(SSB[,1,1],a_RK,b_RK)}
  }
  #if(niter==1 & CV_REC_C==0){if(type=="cte"){N[1,1,1]<-N[1,1,1]}}
  #if(niter>1 & CV_REC_C==0) {if(type=="cte"){N[1,1,]<-rep(N[1,1,1],niter)}}


  if(niter==1 & CV_REC_BH==0){if(type=="BH"){N[1,1,1]<-RBH(SSB[,1,1],a_BH,b_BH)}}
  if(niter>1 & CV_REC_BH==0) {if(type=="BH"){N[1,1,]<-RBH(SSB[,1,],a_BH,b_BH)}}

  if(niter==1 & CV_REC_RK==0){if(type=="RK"){N[1,1,1]<-RRK(SSB[,1,1],a_RK,b_RK)}}
  if(niter>1 & CV_REC_RK==0) {if(type=="RK"){N[1,1,]<-RRK(SSB[,1,],a_RK,b_RK)}}

}


  for(j in 2:number_years){
    for(i in 2:number_ages){

      if(i==number_ages){N[i,j,]<-N[i-1,j-1,]*exp(-Z[i-1,j-1,])+N[i,j-1,]*exp(-Z[i,j-1,])}else{N[i,j,]<-N[i-1,j-1,]*exp(-Z[i-1,j-1,])}


      WM[i,j,]<-N[i,j,]*W[i,j,]
      WMA[i,j,]<-WM[i,j,]*Mat[i,j,]
      if(niter>1){SSB[,j,]<-colSums(WMA[,j,])} else {SSB[,j,1]<-sum(WMA[,j,])}
      if(type=="cte"){
        if(CV_REC_C>0){N[1,j,-1]<-Log.normal(rep(N[1,1,1],niter-1),CV_REC_C)
        N[1,j,1]=N[1,1,1]
        }
        if(niter==1 & CV_REC_C==0){if(type=="cte"){N[1,j,1]<-N[1,1,1]}}
        if(niter>1 & CV_REC_C==0) {if(type=="cte"){N[1,j,1:niter]<-rep(N[1,1,1],niter)}}
      }
      if(type=="BH"){
        if(CV_REC_BH>0){N[1,j,-1]<-Log.normal(RBH(SSB[,j,-1],a_BH,b_BH),CV_REC_BH)
        N[1,j,1]<-RBH(SSB[,j,1],a_BH,b_BH)
        }
        if(niter==1 & CV_REC_BH==0){if(type=="BH"){N[1,j,1]<-RBH(SSB[,j,1],a_BH,b_BH)}}
        if(niter>1 & CV_REC_BH==0) {if(type=="BH"){N[1,j,1:niter]<-RBH(SSB[,j,1:niter],a_BH,b_BH)}}
      }
      if(type=="RK"){
        if(CV_REC_RK>0){ N[1,j,-1]<-Log.normal(RRK(SSB[,j,-1],a_RK,b_RK),CV_REC_RK)
        N[1,j,1]<-RRK(SSB[,j,1],a_RK,b_RK)
        }
        if(niter==1 & CV_REC_RK==0){if(type=="RK"){N[1,j,1]<-RRK(SSB[,j,1],a_RK,b_RK)}}
        if(niter>1 & CV_REC_RK==0) {if(type=="RK"){N[1,j,1:niter]<-RRK(SSB[,j,1:niter],a_RK,b_RK)}}
      }}}


  for(j in 1:number_years){
    i=1

  WM[i,j,]<-N[i,j,]*W[i,j,]
  }

### CAPTURES

  C_N<-N


  for(i in 1:number_ages){
    for(j in 1:number_years){
      C_N[i,j,]<-(F[i,j,]/Z[i,j,])*(N[i,j,]*(1-exp(-Z[i,j,])))
    }
  }


### LENGTH CAPTURES
  L_c<-N

### Stochastic (Normal CV_LC)

  L_cd<-matrix(rep(Length_VB(L_inf,k,ages+tc,t0),number_years), ncol=number_years,nrow=number_ages)

  if(CV_LC>0){
  for (i in 1:number_ages){
    for(j in 1:number_years){
      m<-L_cd[i,j]
      v<-(CV_LC*m)
      L_c[i,j,]<-stats::rnorm(niter,m,v)
    }}
    L_c[,,1]<-L_cd
  }

  if(niter==1 & CV_LC==0){L_c[,,1]<-L_cd}
  if(niter>1 & CV_LC==0) {L_c[,,1:niter]<-L_cd}



  W_c<-N

  for(i in 1:number_years){
    W_c[,i,]<-Weight(L_c[,i,],a,b)
  }


  C_W<-N

  for(i in 1:number_ages){
    for(j in 1:number_years){
      C_W[i,j,]<-C_N[i,j,]*W_c[i,j,]
    }
  }





  results<-list()
  results[[1]]<-list(N=N,M=M,F=F,W=W,Mat=Mat,C_N=C_N,C_W=C_W)
  results[[2]]<-list(ctrBio=ctrBio[-1],ctrFish=ctrFish,SR=SR,minFage=min_age,maxFage=max_age,seed=seed,ts=ts,tc=tc)
  names(results)<-c("Matrices","Info")
  return(results)}
