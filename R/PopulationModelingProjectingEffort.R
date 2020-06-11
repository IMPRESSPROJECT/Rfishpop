
Population.Modeling.Projecting.effort<-function(Pop.Mod,new.years,f.new,my.catch,iter){

  ts<-Pop.Mod$Info$ts
  tc<-Pop.Mod$Info$tc

  if(is.null(ts)){ts=0}
  if(is.null(tc)){tc=0.5}



  seed=Pop.Mod$Info$seed
  set.seed(Pop.Mod$Info$seed)

  Sel_type=Pop.Mod$Info$ctrFish$ctrSEL$type

  ### Taking out the values of the list
  type<-Pop.Mod$Info$SR$type

  CV_REC_BH<-Pop.Mod$Info$SR$par[3]



  CV_REC_RK<-Pop.Mod$Info$SR$par[3]


  if(is.null(type)){type="cte"}

  CV_LC<-Pop.Mod$Info$ctrBio$CV_LC
  CV_SEL<-Pop.Mod$Info$ctrFish$ctrSEL$CV_SEL
  CV_L<-Pop.Mod$Info$ctrBio$CV_L
  CV_M<-Pop.Mod$Info$ctrBio$CV_M
  CV_Mat<-Pop.Mod$Info$ctrBio$CV_Mat

  N=Pop.Mod$Matrices$N
  niter<-dim(N)[3]



  ### Establishing default values and checkings
  if(is.null(niter)){niter=1}
  if(is.null(CV_SEL)){CV_SEL=0}
  if(is.null(CV_M)){CV_M=0}
  if(is.null(CV_L)){CV_L=0}

  if(is.null(CV_Mat)){CV_Mat=0}
  if(is.null(CV_LC)){CV_LC=0}
  if(is.null(CV_REC_BH)){CV_REC_BH=0}
  if(is.null(CV_REC_RK)){CV_REC_RK=0}

  check.sum<-sum(c(CV_SEL,CV_M,CV_L,CV_Mat,CV_LC,CV_REC_BH,CV_REC_RK))


  ### Taking out the values of the list

  years<-as.numeric(colnames(N))
  ages<-as.numeric(rownames(N))
  #N0<-ctrPop$N0
  a50_Sel<-Pop.Mod$Info$ctrFish$ctrSEL$par$a50_Sel
  ad_Sel<-Pop.Mod$Info$ctrFish$ctrSEL$par$ad_Sel
  f<-Pop.Mod$Info$ctrFish$f
  M<-Pop.Mod$Matrices$M
  #M=ctrBio$M
  L_inf<-Pop.Mod$Info$ctrBio$L_inf
  t0<-Pop.Mod$Info$ctrBio$t0
  k<-Pop.Mod$Info$ctrBio$k
  a<-Pop.Mod$Info$ctrBio$a
  b<-Pop.Mod$Info$ctrBio$b
  a50_Mat<-Pop.Mod$Info$ctrBio$a50_Mat
  ad_Mat<-Pop.Mod$Info$ctrBio$ad_Mat
  min_age<-Pop.Mod$Info$minFage
  max_age<-Pop.Mod$Info$maxFage
  type<-Pop.Mod$Info$SR$type
  if (type=="BH"){
    a_BH<-Pop.Mod$Info$SR$par[1]
    b_BH<-Pop.Mod$Info$SR$par[2]
  }

  if (type=="RK"){
    a_RK<-Pop.Mod$Info$SR$par[1]
    b_RK<-Pop.Mod$Info$SR$par[2]
  }
  ### New names

  number_years<-length(years)
  number_ages<-length(ages)

  #lN0<-length(N0)
  #if(lN0!=1 & lN0!=number_ages){
  #  stop("ERROR N0 must be a number or a vector of dimension number of ages")
  #}

  #new.years=c(2021)
  number.new.years=length(new.years)

  column.names <- c(years,new.years)
  row.names <- ages
  matrix.names <- 1:niter

  N<-array(rep(0,number_ages*(number_years+number.new.years)), dim=c(number_ages, number_years+number.new.years, niter),dimnames=list(row.names,column.names,
                                                                                                                                      matrix.names))
  N[,1:number_years,]=Pop.Mod$Matrices$N


  # SELECTIVITY: We need to generate stochastic values of a50_Sel (CV_SEL)

  if(Sel_type=="Logistic"){

    a50_Sel<-Pop.Mod$Info$ctrFish$ctrSEL$par$a50_Sel
    ad_Sel<-Pop.Mod$Info$ctrFish$ctrSEL$par$ad_Sel

    ### Deterministic
    sd=Sum.Pop.Mod(Pop.Mod,c("SEL"))$SEL[,,1]
    m_aux=matrix(rep(Logistic(x=ages,x50=a50_Sel,xd=ad_Sel),number.new.years), ncol=number.new.years,nrow=number_ages)
    sd<-cbind(sd,m_aux)


    ### Stochastic (Log normal distribution)
    s<-N

    if(CV_SEL>0){
      s<-stochastic_logistic_SEL_1(a50_Sel,ad_Sel,CV_SEL,niter,s,ages,number_years+number.new.years,seed)

      s[,,1]<-sd
      aux_SEL=Sum.Pop.Mod(Pop.Mod,c("SEL"))$SEL
      s[,1:number_years,]=aux_SEL
    }}

  if(Sel_type=="cte"){

    ### Deterministic
    cte<-Pop.Mod$Info$ctrFish$ctrSEL$par$cte
    sd<-matrix(rep(cte,number_ages*(number_years+number.new.years)),ncol=(number_years+number.new.years))
    colnames(sd)<-c(years,new.years)
    rownames(sd)<-ages

    ### Stochastic (Uniform)
    s<-N

    if(CV_SEL>0){
      s<-stochastic_cte_SEL_1(cte,CV_SEL,niter,s,number_years+number.new.years,number_ages,seed)

      s[,,1]<-sd
      aux_SEL=Sum.Pop.Mod(Pop.Mod,c("SEL"))$SEL
      s[,1:number_years,]=aux_SEL
    }}

  if(Sel_type=="Andersen"){
    ### Deterministic
    p1<-Pop.Mod$Info$ctrFish$ctrSEL$par$p1;p3<-Pop.Mod$Info$ctrFish$ctrSEL$par$p3
    p4<-Pop.Mod$Info$ctrFish$ctrSEL$par$p4;p5<-Pop.Mod$Info$ctrFish$ctrSEL$par$p5

    sd=Sum.Pop.Mod(Pop.Mod,c("SEL"))$SEL[,,1]
    m_aux=matrix(rep(andersen(x=ages,p1=p1,p3=p3,p4=p4,p5=p5),number.new.years), ncol=number.new.years,nrow=number_ages)

    sd<-cbind(sd,m_aux)


    ### Stochastic (Uniform)
    s<-N

    if(CV_SEL>0){
      s<-stochastic_andersen_SEL_1(p1=p1,p3=p3,p4=p4,p5=p5,CV_SEL,niter,s,ages,number_years+number.new.years,seed)

      s[,,1]<-sd
      aux_SEL=Sum.Pop.Mod(Pop.Mod,c("SEL"))$SEL
      s[,1:number_years,]=aux_SEL
    }
  }

  if(Sel_type=="Gamma"){

    alpha<-Pop.Mod$Info$ctrFish$ctrSEL$par$alpha
    gamma<-Pop.Mod$Info$ctrFish$ctrSEL$par$gamma
    beta<-Pop.Mod$Info$ctrFish$ctrSEL$par$beta


    ### Deterministic
    sd=Sum.Pop.Mod(Pop.Mod,c("SEL"))$SEL[,,1]
    m_aux=matrix(rep(gamma_SEL(x=ages,alpha=alpha,gamma=gamma,beta=beta),number.new.years), ncol=number.new.years,nrow=number_ages)

    sd<-cbind(sd,m_aux)


    ### Stochastic (Log normal distribution)
    s<-N

    if(CV_SEL>0){
      s<-stochastic_gamma_SEL_1(alpha,beta,gamma,CV_SEL,niter,s,ages,number_years+number.new.years,seed)

      s[,,1]<-sd
      aux_SEL=Sum.Pop.Mod(Pop.Mod,c("SEL"))$SEL
      s[,1:number_years,]=aux_SEL
    }}

  if(niter==1 & CV_SEL==0){s[,,1]<-sd}
  if(niter>1 & CV_SEL==0) {s[,,1:niter]<-sd}

  ### FISHING MORTALITY



  F<-N
  ######################################################CUIDADO

  fm=matrix(0,ncol=number_years+number.new.years, nrow=niter)
  f.new=rep(f.new,niter)
  for(i in 1:niter){

    fm[i,]=c(f[i,],f.new[i])
  }
  f=fm

  F[,1:number_years,]=Pop.Mod$Matrices$F
  for(kk in 1:niter){
    for(i in 1:number_ages){
      for (j in (number_years+1):(number_years+number.new.years)){
        F[i,j,kk]<-s[i,j,kk]*f[kk,j]
      }
    }
  }
  ### NATURAL MORTALITY
  m_aux=matrix(rep(M[,number_years,1],number.new.years), ncol=number.new.years,nrow=number_ages)

  Md<-cbind(M[,,1],m_aux)


  ### Stochastic (Log normal distribution; CV_M)
  M<-N

  M[,1:number_years,]=Pop.Mod$Matrices$M
  if(CV_M>0){
    for (i in 1:number_ages){
      for(j in (number_years+1):(number_years+number.new.years)){
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

  # Completing first column of N (not necessary)

  ### LENGTH
  L<-N

  m_aux=matrix(rep(Length_VB(L_inf,k,ages+ts,t0),number.new.years), ncol=number.new.years,nrow=number_ages)


  Ld<-cbind(Sum.Pop.Mod(Pop.Mod,c("LS"))$LS[,,1],m_aux)


  ### Stochastic (Normal CV_L)
  aux_L=Sum.Pop.Mod(Pop.Mod,c("LS"))$LS
  L[,1:number_years,]=aux_L
  if(CV_L>0){
    for (i in 1:number_ages){
      for(j in (number_years+1):(number_years+number.new.years)){
        m<-Ld[i,j]
        v<-(CV_L*m)
        L[i,j,]<-rnorm.seed(niter,m,v,seed)
      }}
    L[,,1]<-Ld
  }
  ### Deterministic
  if(niter==1 & CV_L==0){L[,,1]<-Ld}
  if(niter>1 & CV_L==0) {L[,,1:niter]<-Ld}


  ### WEIGHTS

  W<-N
  W[,1:number_years,]=Pop.Mod$Matrices$W

  for(i in (number_years+1):(number_years+number.new.years)){
    W[,i,]<-Weight(L[,i,],a,b)
  }



  # WM first column



  ### MATURITY (Log normal distribution CV_Mat)
  Mat<-N
  ### Deterministic


  m_aux=matrix(rep(Logistic(x=ages,x50=a50_Mat,xd=ad_Mat),number.new.years), ncol=number.new.years,nrow=number_ages)
  Matd<-cbind(Pop.Mod$Matrices$Mat[,,1],m_aux)

  Mat[,1:number_years,]=Pop.Mod$Matrices$Mat
  if(CV_Mat>0){
    for (ind in 2:niter){
      Mat[,(number_years+1):(number_years+number.new.years),ind]<-matrix(rep(Mat[,1,ind],number.new.years), ncol=number.new.years,nrow=number_ages)
    }
    Mat[,,1]<-Matd
  }


  if(niter==1 & CV_Mat==0){Mat[,,1]<-Matd}
  if(niter>1 & CV_Mat==0) {Mat[,,1:niter]<-Matd}



  # Defining WMA first column not necessary.

  ### SSB
  row.names=""
  SSB<-array(rep(0,(number_years+number.new.years)), dim=c(1, number_years+number.new.years, niter),dimnames=list(row.names,column.names,matrix.names))

  SSB[,1:number_years,]=Sum.Pop.Mod(Pop.Mod,c("SSB"))$SSB

  WM=N;WMA=N

  WM[,1:number_years,]=Sum.Pop.Mod(Pop.Mod,c("WS"))$WS

  WMA[,1:number_years,]=Sum.Pop.Mod(Pop.Mod,c("WSSB"))$WSSB



  for(j in (number_years+1):(number_years+number.new.years)){
    for(i in 2:number_ages){
      N[i,j,]<-N[i-1,j-1,]*exp(-Z[i-1,j-1,])

      WM[i,j,]<-N[i,j,]*W[i,j,]
      WMA[i,j,]<-WM[i,j,]*Mat[i,j,]
      if(niter>1){SSB[,j,]<-colSums(WMA[,j,])} else {SSB[,j,1]<-sum(WMA[,j,])}
      if(type=="cte"){N[1,j,]<-N[1,1,]}
      if(type=="BH"){
        if(CV_REC_BH>0){ N[1,j,-1]<-Log.normal(RBH(SSB[,j,-1],a_BH,b_BH),CV_REC_BH,seed)
        N[1,j,1]<-RBH(SSB[,j,1],a_BH,b_BH)
        }
        if(niter==1 & CV_REC_BH==0){if(type=="BH"){N[1,j,1]<-RBH(SSB[,j,1],a_BH,b_BH)}}
        if(niter>1 & CV_REC_BH==0) {if(type=="BH"){N[1,j,1:niter]<-RBH(SSB[,j,1:niter],a_BH,b_BH)}}
      }
      if(type=="RK"){
        if(CV_REC_RK>0){ N[1,j,-1]<-Log.normal(RRK(SSB[,j,-1],a_RK,b_RK),CV_REC_RK,seed)
        N[1,j,1]<-RRK(SSB[,j,1],a_RK,b_RK)
        }
        if(niter==1 & CV_REC_RK==0){if(type=="RK"){N[1,j,1]<-RRK(SSB[,j,1],a_RK,b_RK)}}
        if(niter>1 & CV_REC_RK==0) {if(type=="RK"){N[1,j,1:niter]<-RRK(SSB[,j,1:niter],a_RK,b_RK)}}
      }}}


  ### CAPTURES

  C_N<-N
  C_N[,1:number_years,]=Pop.Mod$Matrices$C_N


  for(i in 1:number_ages){
    for(j in (number_years+1):(number_years+number.new.years)){
      C_N[i,j,]<-(F[i,j,]/Z[i,j,])*(N[i,j,]*(1-exp(-Z[i,j,])))
    }
  }


  ### LENGTH CAPTURES
  L_c<-N

  m_aux=matrix(rep(Length_VB(L_inf,k,ages+tc,t0),number.new.years), ncol=number.new.years,nrow=number_ages)


  L_cd<-cbind(Sum.Pop.Mod(Pop.Mod,c("LC"))$LC[,,1],m_aux)

  L_c[,1:number_years,]=Sum.Pop.Mod(Pop.Mod,c("LC"))$LC
  ### Stochastic (Normal CV_LC)

  if(CV_LC>0){
    for (i in 1:number_ages){
      for(j in (number_years+1):(number_years+number.new.years)){
        m<-L_cd[i,j]
        v<-(CV_LC*m)
        L_c[i,j,]<-rnorm.seed(niter,m,v,seed)
      }}
    L_c[,,1]<-L_cd
  }

  if(niter==1 & CV_LC==0){L_c[,,1]<-L_cd}
  if(niter>1 & CV_LC==0) {L_c[,,1:niter]<-L_cd}



  W_c<-N
  W_c[,1:number_years,]=Sum.Pop.Mod(Pop.Mod,c("WC"))$WC


  for(i in (number_years+1):(number_years+number.new.years)){
    W_c[,i,]<-Weight(L_c[,i,],a,b)
  }


  C_W<-N

  C_W[,1:number_years,]=Pop.Mod$Matrices$C_W

  for(i in 1:number_ages){
    for(j in (number_years+1):(number_years+number.new.years)){
      C_W[i,j,]<-C_N[i,j,]*W_c[i,j,]
    }
  }





  year_C_W<-array(rep(0,number_years+number.new.years), dim=c(1, number_years+number.new.years, niter),dimnames=list("",column.names,matrix.names))

  for (ind in 1:niter){
    year_C_W[,,ind]<-colSums(C_W[,,ind])
  }


return((year_C_W[,number_years+number.new.years,iter]-my.catch))}
