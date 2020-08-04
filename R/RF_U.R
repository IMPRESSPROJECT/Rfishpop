
RF_U<-function(Pop.Mod, Fish.years,Bio.years,Method,par,FM_type,iters,N0,prop){
  if(is.null(Method)){Method="mean"}

# DUE TO PARALLELIZATION READ HERE SUM.POP.MOD


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


### USEFUL NEW FUNCTIONS



RF_aux<-function(Pop.Mod, Fish.years,Bio.years,Method,par,FM_type,ind2,N0,prop){

  if (FM_type=="F_BPR"){
    niter<-1
    xx<-array(0,dim=c(1,1,niter))
    FBPRp_resul<-array(0,dim = c(1,7,niter))
    colnames(FBPRp_resul)<-c("f","F","YPR","BPR","R","Y","B")

    # Computes max SPR
    BPR0=BPR_RF(Pop.Mod,0,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par, ind2=ind2)[,2]


    BPR_p <- function (f, Pop.Mod, Fish.years,Bio.years,plot,Method=Method,par=par,ind2=ind2,BPR0=BPR0,prop=prop) return(abs(BPR_RF(Pop.Mod,f,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par, ind2=ind2)[,2]-prop*BPR0))

    xx<- stats::optimize(BPR_p, c(0,10), Pop.Mod=Pop.Mod, Fish.years=Fish.years,Bio.years=Bio.years,plot=FALSE,Method=Method,par=par,ind2=ind2,BPR0=BPR0,prop=prop, tol=0.01)


    f.BPRp<-xx$minimum

    RE<-BYR.eq_RF(Pop.Mod,f.BPRp,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2=ind2,N0=N0)



    DPM<-RE$DPM; FBPRp_resul[,,1]<-DPM

    return(F_BRPp=FBPRp_resul)
  }



  if (FM_type=="F_Crash"){
  type=Pop.Mod$Info$SR$type
  a=Pop.Mod$Info$SR$par[1]
  b=Pop.Mod$Info$SR$par[2]

      if(type=="BH"){
        slopeSR=RBH(SSB=0.0000001, a, b)
      }

      if(type=="RK"){
        slopeSR=RRK(SSB=0.0000001, a, b)
      }

      if(type=="cte"){
        results=matrix(rep(NA,7), ncol=7)
        colnames(results)<-c("f","F","YPR","BPR","R","Y","B")
        Fc_resul=results
      }

  if(type=="RK" | type=="BH"){

      crash <- function (f, Pop.Mod,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2,slopeSR=slopeSR,N0=N0){
        RE<-BYR.eq_RF(Pop.Mod,f,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2=ind2,N0=N0)
        return(abs(RE$DPM[,4]-(0.0000001/slopeSR)))
      }

      #xx <- stats::uniroot(crash,c(0.0001,60), Pop.Mod=Pop.Mod, Fish.years=Fish.years,Bio.years=Bio.years,plot=FALSE,Method=Method,par=par,ind2=ind2, tol=0.01)
      #f.c<-xx$root

      xx <- stats::optimize(crash,c(0.0001,60), Pop.Mod=Pop.Mod, Fish.years=Fish.years,Bio.years=Bio.years,plot=FALSE,Method=Method,par=par,ind2=ind2, slopeSR=slopeSR, N0=N0,tol=0.01)
      f.c<-xx$minimum

      RE<-BYR.eq_RF(Pop.Mod,f.c,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2,N0=N0)

      DPM<-RE$DPM; Fc_resul<-DPM}

      return(F_Crash=Fc_resul)

    }

    if (FM_type=="F_max"){
      N<-Pop.Mod$Matrices$N;ages<-as.numeric(rownames(N))

      niter<-1
      # Find F that produces max YPR: Fmax
      xx<-array(0,dim=c(1,1,niter))
      Fmax_resul<-array(0,dim = c(1,7,niter))
      colnames(Fmax_resul)<-c("f","F","YPR","BPR","R","Y","B")

      YPR_op <- function (f, Pop.Mod, Fish.years,Bio.years,plot,Method=Method,par=par,ind2=ind2) YPR_RF(Pop.Mod,f,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par, ind2=ind2)[,2]

      xx<- stats::optimize(YPR_op, c(0,10), Pop.Mod=Pop.Mod, Fish.years=Fish.years,Bio.years=Bio.years,plot=FALSE,Method=Method,par=par,ind2=ind2, maximum=TRUE, tol=0.01)


      f.max<-xx$maximum

      RE<-BYR.eq_RF(Pop.Mod,f.max,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2=ind2,N0=N0)


      # All the values for the Fmax

      DPM<-RE$DPM; Fmax_resul[,,1]<-DPM

      return(F_max=Fmax_resul)
    }

    if (FM_type=="F_0.1"){

      # F 0.1 function
      slope<-0.1

      xx <-YPR_RF(Pop.Mod,0.00001,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2=ind2)
      trgSlope <- xx[,2]/xx[,1]*slope


      yprSlope <- function (f, Pop.Mod, Fish.years,Bio.years,plot=FALSE,trgSlope,Method=Method,par=par,ind2){
        xxx <- YPR_RF(Pop.Mod, c(f, f+0.0001),Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2)
        xSlope <- (xxx[2,2]-xxx[1,2])/(xxx[2,1]-xxx[1,1])
        return(trgSlope-xSlope)
      }

      xx <- stats::uniroot(yprSlope, c(0.0001,10), Pop.Mod=Pop.Mod, Fish.years=Fish.years,Bio.years=Bio.years,plot=FALSE,trgSlope,Method=Method,par=par,ind2=ind2, tol=0.01)
      f.01<-xx$root

      RE<-BYR.eq_RF(Pop.Mod,f.01,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2,N0=N0)


      # All the values for the F 0.1

      DPM<-RE$DPM; F01_resul<-DPM
      return(F_0.1=F01_resul)
    }

    if (FM_type=="F_msy"){
      ### Fmsy

      MSY <- function (f, Pop.Mod,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2,N0=N0){
        RR<-BYR.eq_RF(Pop.Mod,f,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2=ind2,N0=N0)
        RR$DPM[,6]
      }

      a<-MSY(seq(0.01,5,by=0.01), Pop.Mod,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2,N0=N0)
      ind<-which(a==0);if(length(ind)>0){b<-min(ind)};
      c<-seq(0.01,10,by=0.01)
      if(length(ind)==0){b=length(c)}
      # Return F and corresponding equilibrium values
      xx <- stats::optimize(MSY, c(0,c[b]), Pop.Mod=Pop.Mod,Fish.years=Fish.years,Bio.years=Bio.years,plot=FALSE,Method=Method,par=par,ind2=ind2,N0=N0,maximum=TRUE, tol=0.01)

      f.msy<-xx$maximum
      RE<-BYR.eq_RF(Pop.Mod,f.msy,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2,N0=N0)


      # All the values for the Fmsy

      DPM<-RE$DPM; Fmsy_resul<-DPM
      return(F_msy=Fmsy_resul)

    }
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



### REQUIRED FUNCTIONS FOR THE ABOVE CODE

BYR.eq_RF<-function(Pop.Mod,f.grid,Fish.years,Bio.years,plot,Method,par,ind2,N0){

    SR<-Pop.Mod$Info$SR
    min_age<-Pop.Mod$Info$minFage
    max_age<-Pop.Mod$Info$maxFage

    ypr<-YPR_RF(Pop.Mod,f.grid,Fish.years,Bio.years,plot,Method,par,ind2=ind2)
    bpr<-BPR_RF(Pop.Mod,f.grid,Fish.years,Bio.years,plot,Method,par,ind2=ind2)
    ### First steps
    F<-bpr[,1];l<-length(F)
    ind<-which(bpr[,2]<0);bpr[ind,2]<-0
    ind<-which(ypr[,2]<0);ypr[ind,2]<-0
    R<-1:l;B<-1:l;Y<-1:l


    type<-SR$type
    if (type=="BH"){
      a_BH<-SR$par[1]
      b_BH<-SR$par[2]
      for (i in 1:l){
        if(bpr[i,2]>0){
          R[i]<-(-b_BH+a_BH*bpr[i,2])/bpr[i,2]} else {R[i]<-0}
        B[i]<-bpr[i,2]*R[i]
        Y[i]<-ypr[i,2]*R[i]
      }
    }

    if (type=="RK"){
      a_RK<-SR$par[1]
      b_RK<-SR$par[2]
      for (i in 1:l){
        if(bpr[i,2]>0){
        R[i]<-log(1/(a_RK*bpr[i,2]))/(-b_RK*bpr[i,2])} else {R[i]<-0}
        B[i]<-bpr[i,2]*R[i]
        Y[i]<-ypr[i,2]*R[i]
      }
    }

    if (type=="cte"){

      for (i in 1:l){
        R[i]<-N0
        B[i]<-bpr[i,2]*R[i]
        Y[i]<-ypr[i,2]*R[i]
      }
    }
    ind<-which(R<0);R[ind]<-0
    ind<-which(B<0);B[ind]<-0
    ind<-which(Y<0);Y[ind]<-0
    results<-cbind(f=f.grid,F,YPR=ypr[,2],BPR=bpr[,2],R.eq=R,Y.eq=Y,B.eq=B)


    return(list(DPM=results))
  }








  BPR_RF<-function(Pop.Mod,f.grid,Fish.years,Bio.years,plot,Method,par,ind2){
    min_age<-Pop.Mod$Info$minFage
    max_age<-Pop.Mod$Info$maxFage


    N<-Pop.Mod$Matrices$N[,,ind2];ages<-as.numeric(rownames(N))
    number_years<-ncol(N)
    number_ages<-nrow(N)

    if(Method=="mean"){
      W<-Pop.Mod$Matrices$W[,,ind2]
      W<-apply(W[,(number_years-Bio.years+1):number_years],1,mean)
      M<-Pop.Mod$Matrices$M[,,ind2];M<-apply(M[,(number_years-Bio.years+1):number_years],1,mean)
      Mat<-Pop.Mod$Matrices$Mat[,,ind2];Mat<-apply(Mat[,(number_years-Bio.years+1):number_years],1,mean)
      s<-(Sum.Pop.Mod(Pop.Mod,c("SEL"))$SEL)[,,ind2]; s<-apply(s[,(number_years-Bio.years+1):number_years],1,mean)
    }

    if(Method=="own"){
      W<-par$W[,ind2]
      M<-par$M[,ind2]
      Mat<-par$Mat[,ind2]
      s<-par$SEL[,ind2]
    }


    l<-length(f.grid);F.vector<-1:l;bpr.vector<-1:l
    a<-1:l
    for (j in 1:l){
      F<-f.grid[j]*s #Now is a vector but later will be a matrix and a mean will be necessary
      Z<-F+M

      N<-1
      bpr<-0
      for (i in 1:(number_ages-1)){

        bpr<-bpr+N*W[i]*Mat[i]
        N<-N*exp(-Z[i])
      }

      bpr<-bpr+(N*W[number_ages]*Mat[number_ages]/(1-exp(-Z[number_ages])))

      bpr.vector[j]<-bpr
      F.vector[j]<-mean(F[(min_age-ages[1]+1):(max_age-ages[1]+1)])
    }
    ind<-which(bpr.vector<0);bpr.vector[ind]<-0
    results<-cbind(F.vector,bpr.vector)
    colnames(results)<-c("F","BPR")
    rownames(results)<-NULL
    if(plot==TRUE){
      graphics::plot(results[,1],results[,2], main="Biomass per recruit",xlab="F (fishing mortality)",ylab="BPR",type="l")
    }
    return(results)}









  YPR_RF<-function(Pop.Mod,f.grid,Fish.years,Bio.years,plot,Method,par,ind2){
    min_age<-Pop.Mod$Info$minFage
    max_age<-Pop.Mod$Info$maxFage
    N<-Pop.Mod$Matrices$N[,,ind2];ages<-as.numeric(rownames(N))

    number_ages<-nrow(N)
    number_years<-ncol(N)

    if(Method=="mean"){
      W.c<-(Sum.Pop.Mod(Pop.Mod,c("WC"))$WC)[,,ind2];W.c<-apply(W.c[,(number_years-Bio.years+1):number_years],1,mean)
      M<-Pop.Mod$Matrices$M[,,ind2];M<-apply(M[,(number_years-Bio.years+1):number_years],1,mean)
      s<-(Sum.Pop.Mod(Pop.Mod,c("SEL"))$SEL)[,,ind2];s<-apply(s[,(number_years-Bio.years+1):number_years],1,mean)
    }

    if(Method=="own"){
      W.c<-par$WC[,ind2]
      M<-par$M[,ind2]
      s<-par$SEL[,ind2]
    }

    l<-length(f.grid);F.vector<-1:l;ypr.vector<-1:l
    a<-1:l
    for (j in 1:l){
      F<-f.grid[j]*s #Now is a vector but later will be a matrix and a mean will be necessary
      Z<-F+M
      N<-1
      ypr<-0


      for (i in 1:(number_ages-1)){

        ypr<-ypr+N*F[i]*W.c[i]*((1-exp(-Z[i]))/Z[i])
        N<-N*exp(-Z[i])
      }

      ypr<-ypr+(N*F[number_ages]*W.c[number_ages]/Z[number_ages])

      ypr.vector[j]<-ypr
      F.vector[j]<-mean(F[(min_age-ages[1]+1):(max_age-ages[1]+1)])
    }
    ind<-which(ypr.vector<0);ypr.vector[ind]<-0
    results<-cbind(F.vector,ypr.vector)
    colnames(results)<-c("F","YPR")
    rownames(results)<-NULL
    if(plot==TRUE){
      graphics::plot(results[,1],results[,2], main="Yield per recruit",xlab="F (fishing mortality)",ylab="YPR",type="l")
    }
    return(results)}



 ### Now, we can use all the functions

  myfxn <- function(ind2,Fish.years, Pop.Mod,Bio.years, Method, par, FM_type,N0,prop){
    return(RF_aux(Pop.Mod, Fish.years,Bio.years,Method,par,FM_type,ind2,N0,prop))}


  n=length(iters)


  iters_list<-as.list(iters)


  numCores <- parallel::detectCores()

  cl <- parallel::makeCluster(numCores)

  parallel::clusterExport(cl,c("myfxn","RF_aux","BYR.eq_RF","YPR_RF","BPR_RF","Sum.Pop.Mod","RBH","RRK"),
                          envir=environment())

  ret<-parallel::parLapply(cl, iters_list,myfxn, Pop.Mod=Pop.Mod, Bio.years=Bio.years,Method=Method,par=par,FM_type=FM_type,Fish.years=Fish.years,N0=N0,prop=prop)
  parallel::stopCluster(cl)



  return(array(unlist(ret),dim = c(1,7,n)))
}
