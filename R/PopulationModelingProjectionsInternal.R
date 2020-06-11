Population.Modeling.Projections.Internal=function(Pop.Mod,new.years,my.catch,limit.f,tol){
#
  ## Function required in due to parallelize



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

  Sum.Pop.Mod<-function(Pop.Mod,Elements){
    ts<-Pop.Mod$Info$ts
    tc<-Pop.Mod$Info$tc

    seed=Pop.Mod$Info$seed
    set.seed(Pop.Mod$Info$seed)

    Sel_type=Pop.Mod$Info$ctrFish$ctrSEL$type
    CV_SEL<-Pop.Mod$Info$ctrFish$ctrSEL$CV_SEL
    ### Stochastic Parameters

    CV_LC<-Pop.Mod$Info$ctrBio$CV_LC
    CV_L<-Pop.Mod$Info$ctrBio$CV_L


    niter<-dim(Pop.Mod$Matrices$N)[3]

    #Matrices
    N<-Pop.Mod$Matrices$N
    F<-Pop.Mod$Matrices$F
    M<-Pop.Mod$Matrices$M
    W<-Pop.Mod$Matrices$W
    Mat<-Pop.Mod$Matrices$Mat
    C_W<-Pop.Mod$Matrices$C_W

    # Bio
    L_inf<-Pop.Mod$Info$ctrBio$L_inf
    t0<-Pop.Mod$Info$ctrBio$t0
    k<-Pop.Mod$Info$ctrBio$k
    a<-Pop.Mod$Info$ctrBio$a
    b<-Pop.Mod$Info$ctrBio$b
    min_age<-Pop.Mod$Info$minFage
    max_age<-Pop.Mod$Info$maxFage
    ages<-as.numeric(rownames(N))
    years<-as.numeric(colnames(N))
    number_years<-ncol(N)
    number_ages<-nrow(N)

    Z<-M+F

    ### LENGTH
    L<-N

    Ld<-matrix(rep(Length_VB(L_inf,k,ages+ts,t0),number_years), ncol=number_years,nrow=number_ages)


    ### Stochastic (Normal CV_L)


    if(CV_L>0){
      for (i in 1:number_ages){
        for(j in 1:number_years){
          m<-Ld[i,j]
          v<-(CV_L*m)
          L[i,j,]<-rnorm.seed(niter,m,v,seed)
        }}
      L[,,1]<-Ld
    }
    ### Deterministic
    if(niter==1 & CV_L==0){L[,,1]<-Ld}
    if(niter>1 & CV_L==0) {L[,,1:niter]<-Ld}



    ### LENGTH CAPTURES
    L_c<-N

    ### Stochastic (Normal CV_LC)

    L_cd<-matrix(rep(Length_VB(L_inf,k,ages+tc,t0),number_years), ncol=number_years,nrow=number_ages)

    if(CV_LC>0){
      for (i in 1:number_ages){
        for(j in 1:number_years){
          m<-L_cd[i,j]
          v<-(CV_LC*m)
          L_c[i,j,]<-rnorm.seed(niter,m,v,seed)
        }}
      L_c[,,1]<-L_cd
    }

    if(niter==1 & CV_LC==0){L_c[,,1]<-L_cd}
    if(niter>1 & CV_LC==0) {L_c[,,1:niter]<-L_cd}




    W_c<-N

    for(i in 1:number_years){
      W_c[,i,]<-Weight(L_c[,i,],a,b)
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

    if(Sel_type=="Logistic"){

      a50_Sel<-Pop.Mod$Info$ctrFish$ctrSEL$par$a50_Sel
      ad_Sel<-Pop.Mod$Info$ctrFish$ctrSEL$par$ad_Sel

      ### Deterministic
      sd<-matrix(rep(Logistic(x=ages,x50=a50_Sel,xd=ad_Sel),number_years),ncol=number_years)
      colnames(sd)<-years
      rownames(sd)<-ages

      ### Stochastic (Log normal distribution)
      s<-N

      if(CV_SEL>0){
        s<-stochastic_logistic_SEL_1(a50_Sel,ad_Sel,CV_SEL,niter,s,ages,number_years,seed=seed)

        s[,,1]<-sd
      }}

    if(Sel_type=="cte"){

      ### Deterministic
      cte<-Pop.Mod$Info$ctrFish$ctrSEL$par$cte
      sd<-matrix(rep(cte,number_years*number_ages),ncol=number_years)
      colnames(sd)<-years
      rownames(sd)<-ages

      ### Stochastic (Uniform)
      s<-N

      if(CV_SEL>0){
        s<-stochastic_cte_SEL_1(cte,CV_SEL,niter,s,number_years,number_ages,seed=seed)

        s[,,1]<-sd
      }}

    if(Sel_type=="Andersen"){
      ### Deterministic
      p1<-Pop.Mod$Info$ctrFish$ctrSEL$par$p1;p3<-Pop.Mod$Info$ctrFish$ctrSEL$par$p3;p4<-Pop.Mod$Info$ctrFish$ctrSEL$par$p4;p5<-Pop.Mod$Info$ctrFish$ctrSEL$par$p5
      sd<-matrix(rep(andersen(x=ages,p1=p1,p3=p3,p4=p4,p5=p5),number_years),ncol=number_years)
      colnames(sd)<-years
      rownames(sd)<-ages

      ### Stochastic (Uniform)
      s<-N

      if(CV_SEL>0){
        s<-stochastic_andersen_SEL_1(p1=p1,p3=p3,p4=p4,p5=p5,CV_SEL,niter,s,ages,number_years,seed=seed)

        s[,,1]<-sd
      }
    }

    if(Sel_type=="Gamma"){

      alpha<-Pop.Mod$Info$ctrFish$ctrSEL$par$alpha
      gamma<-Pop.Mod$Info$ctrFish$ctrSEL$par$gamma
      beta<-Pop.Mod$Info$ctrFish$ctrSEL$par$beta


      ### Deterministic
      sd<-matrix(rep(gamma_SEL(x=ages,alpha=alpha,gamma=gamma,beta=beta),number_years),ncol=number_years)
      colnames(sd)<-years
      rownames(sd)<-ages

      ### Stochastic (Log normal distribution)
      s<-N

      if(CV_SEL>0){
        s<-stochastic_gamma_SEL_1(alpha,beta,gamma,CV_SEL,niter,s,ages,number_years,seed=seed)

        s[,,1]<-sd
      }}


    if(niter==1 & CV_SEL==0){s[,,1]<-sd}
    if(niter>1 & CV_SEL==0) {s[,,1:niter]<-sd}


    year_C_W<-array(rep(0,number_years), dim=c(1, number_years, niter),dimnames=list("",column.names,matrix.names))
    biomass<-array(rep(0,number_years), dim=c(1, number_years, niter),dimnames=list("",column.names,matrix.names))
    Recruiment<-array(rep(0,number_years), dim=c(1, number_years, niter),dimnames=list("",column.names,matrix.names))
    F_mean<-array(rep(0,number_years), dim=c(1, number_years, niter),dimnames=list("",column.names,matrix.names))

    for (ind in 1:niter){
      year_C_W[,,ind]<-colSums(C_W[,,ind])
      biomass[,,ind]<-colSums(WM[,,ind])
      Recruiment[,,ind]<-N[1,,ind]
      F_mean[,,ind]<-apply(F[(min_age+1):(max_age+1),,ind], 2, mean)}





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


  andersen<-function(x,p1,p3,p4,p5){
    if(x[1]==0){x[1]<-0.0001}
    p0<-0
    p2<-1

    number_ages<-length(x);sel_andersen<-1:number_ages

    for (i in 1:number_ages){
      if(log(p5/x[i])<=p1){
        sel_andersen[i]<-p0+p2*exp(-((log(p5/x[i])-p1)^2/p4))
      }


      if(log(p5/x[i])>p1){
        sel_andersen[i]<-p0+p2*exp(-((log(p5/x[i])-p1)^2/p3))
      }

    }
    return(SA=sel_andersen)
  }
  gamma_SEL<-function(x,alpha,beta,gamma){
    if(x[1]==0){x[1]<-0.0001}
    sel_gamma<-((x/((alpha-1)*beta*gamma))^(alpha-1))*exp(alpha-1-(x/(beta*gamma)))
    return(SG=sel_gamma)}
  Logistic<-function(x,x50,xd){
    if(x[1]==0){x[1]<-0.0001}
    Lo<-1/(1+exp((x-x50)/xd))
    return(Lo)
  }
  stochastic_andersen_SEL_1<-function(p1,p3,p4,p5,CV_SEL,niter,s,ages,number_years,seed){
    if(is.numeric(seed)){set.seed(seed)}
    m<-p1
    sigma<-(CV_SEL*m)


    p1_niter<-stats::rnorm(niter,m,sigma)
    if(sum(p1_niter<0)>0){ stop("A random value of p1 is negative, and this has not sense.")}

    for (ind in 2:niter){
      s[,,ind]<-matrix(rep(andersen(x=ages,p1=p1_niter[ind],p3=p3,p4=p4,p5=p5),number_years),ncol = number_years)
    }


    ### Log normal option
    #m<-p1
    #v<-(CV_SEL*m)^2
    #mu<-log(m^2/sqrt(m^2+v))
    #sigma<-sqrt(log(v/m^2+1))

    #p1_niter<-stats::rlnorm(niter,mu,sigma)


    #for (ind in 2:niter){
    #  s[,,ind]<-matrix(rep(andersen(x=ages,p1=p1_niter[ind],p3=p3,p4=p4,p5=p5),number_years),ncol = number_years)
    #}

    return(s)
  }


  stochastic_cte_SEL_1<-function(alpha,CV_SEL,niter,s,number_years,number_ages,seed){
    if(is.numeric(seed)){set.seed(seed)}
    ### We use a Uniform(a,b)
    b<-max(alpha+0.5*sqrt(12)*CV_SEL*alpha,0)
    a<-min(alpha-0.5*sqrt(12)*CV_SEL*alpha,1)




    for (ind in 2:niter){
      values_runif<-stats::runif(number_ages*number_years,min = a,max = b)
      if(sum(values_runif<0 & values_runif>1)>0){ stop("A random value of cte is outside of [0,1]")}

      s[,,ind]<-matrix(values_runif,ncol = number_years)
    }

    return(s)
  }

  stochastic_gamma_SEL_1<-function(alpha,beta,gamma,CV_SEL,niter,s,ages,number_years,seed){
    if(is.numeric(seed)){set.seed(seed)}
    m<-alpha
    sigma<-(CV_SEL*m)

    alpha_Sel_niter<-stats::rnorm(niter,m,sigma)


    for (ind in 2:niter){
      s[,,ind]<-matrix(rep(gamma_SEL(x=ages,alpha=alpha_Sel_niter[ind],beta=beta,gamma=gamma),number_years),ncol = number_years)
    }



    ### log normal
    #m<-alpha
    #v<-(CV_SEL*m)^2
    # mu<-log(m^2/sqrt(m^2+v))
    # sigma<-sqrt(log(v/m^2+1))
    #
    # alpha_Sel_niter<-stats::rlnorm(niter,mu,sigma)
    #
    #
    # for (ind in 2:niter){
    #   s[,,ind]<-matrix(rep(gamma(x=ages,alpha=alpha_Sel_niter[ind],beta=beta,gamma=gamma),number_years),ncol = number_years)
    # }

    return(s)
  }
  stochastic_logistic_SEL_1<-function(a50_Sel,ad_Sel,CV_SEL,niter,s,ages,number_years,seed){
    if(is.numeric(seed)){set.seed(seed)}
    m<-a50_Sel
    sigma<-(CV_SEL*m)


    a50_Sel_niter<-stats::rnorm(niter,m,sigma)

    for (ind in 2:niter){
      s[,,ind]<-matrix(rep(Logistic(x=ages,x50=a50_Sel_niter[ind],xd=ad_Sel),number_years),ncol = number_years)
    }

    ### log normal
    # m<-a50_Sel
    # v<-(CV_SEL*m)^2
    # mu<-log(m^2/sqrt(m^2+v))
    # sigma<-sqrt(log(v/m^2+1))
    #
    # a50_Sel_niter<-stats::rlnorm(niter,mu,sigma)
    #
    #
    # for (ind in 2:niter){
    #   s[,,ind]<-matrix(rep(Logistic(x=ages,x50=a50_Sel_niter[ind],xd=ad_Sel),number_years),ncol = number_years)
    # }

    return(s)
  }
  Length_VB<-function(L_inf,k,ages,t0){
    L<-L_inf*(1-exp(-k*(ages-t0)))
    return(L)
  }

  rnorm.seed=function(n, mean = 0, sd = 1,seed){
    set.seed(seed)
    return(stats::rnorm(n,mean,sd))}
  Weight<-function(L,a,b){
    number_ages<-length(L)
    av<-rep(a,number_ages)
    bv<-rep(b,number_ages)
    W<-av*L^bv
    return(W)
  }
  Log.normal<-function(mean,cv,seed){
    if(is.numeric(seed)){set.seed(seed)}
    m<-mean;n_values<-length(m)
    v<-(cv*m)^2
    mu<-log(m^2/sqrt(m^2+v))

    sigma<-sqrt(log(v/m^2+1))


    values.log.normal<-1:n_values
    for(i in 1:n_values){
      if(is.numeric(mu[i])){}else{stop("See Log-normal function, problem mu")}
      if(sigma[i]<0){stop("See Log-normal function, problem sigma")}
      values.log.normal[i]<-stats::rlnorm(1,mu[i],sigma[i])
    }

    return(values.log.normal)
  }

  RBH<-function(SSB,a,b){
    niter<-length(SSB)
    av<-rep(a,niter)
    bv<-rep(b,niter)
    R<-(av*SSB)/(bv+SSB)
    return(R)
  }
  RRK<-function(SSB,a,b){
    niter<-length(SSB)
    av<-rep(a,niter)
    bv<-rep(b,niter)
    R<-(a*SSB)*exp(-b*SSB)
    return(R)
  }


  # We need to obtain the effort associated to my.catch

  fun=function(f.new,Pop.Mod,new.years,my.catch,iter){return(Population.Modeling.Projecting.effort(Pop.Mod=Pop.Mod,new.years=new.years,f.new=f.new,my.catch=my.catch,iter=iter))}
  N=Pop.Mod$Matrices$N
  if(is.null(tol)){tol=0.01}
  if(is.null(limit.f)){limit.f=4}
  niter<-dim(N)[3]


  fun2=function(iter,Pop.Mod,new.years,my.catch,tol,limit.f){
    xx<- stats::uniroot(fun, c(0,limit.f),Pop.Mod=Pop.Mod,new.years=new.years,my.catch=my.catch, tol=tol,iter=iter)
    effort=xx$root
    return(effort)}

  iters_list<-as.list(1:niter)


  numCores <- parallel::detectCores()

  cl <- parallel::makeCluster(numCores)

  parallel::clusterExport(cl,c("fun2","fun","Population.Modeling.Projecting.effort","Sum.Pop.Mod",
                               "andersen","gamma_SEL","Logistic",
                              "stochastic_andersen_SEL_1","stochastic_cte_SEL_1",
                              "stochastic_gamma_SEL_1","stochastic_logistic_SEL_1",
                              "Length_VB","rnorm.seed","Weight","Log.normal","RRK","RBH"),
                          envir=environment())

  ret<-parallel::parLapply(cl, iters_list,fun2,Pop.Mod=Pop.Mod,new.years=new.years,my.catch=my.catch,tol=tol,limit.f=limit.f)
  parallel::stopCluster(cl)



  v.effort=unlist(ret)


  # We project the population

  Pop.Mod1=Population.Modeling.Projecting(Pop.Mod,new.years,f.new=v.effort)


  return(Pop.Mod1)

}
