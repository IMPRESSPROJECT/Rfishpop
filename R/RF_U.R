#' @title Reference Fishery Mortalities (alternative)
#'
#' @description Return the reference fishery mortality which produces maximum YPR (FM_type="F_max"), the reference fishery mortality at which the slope of the YPR curve is reduced to 0.1 of that estimated at the origin (FM_type="F_0.1") and the reference fishery mortality at which the MSY is attained (FM_type="F_msy"). Furthermore for each of these fishery mortalities the corresponding effort of fishing, YPR, BYR, B, Y and R in equilibrium is also returned. THIS FUNCTION IS EQUIVALENT TO RF FUNCTION, BUT THE IMPLEMENTATION IS A LITTLE DIFFERENT THEN DEPENDING ON THE NUMBERS OF ITERATIONS FOR WHICH THE REFERENCE FISHERY MORTALITIES MUST BE COMPUTED THIS RF_U CAN BE FASTER THAN THE RF FUNCTION IF YOU USE A CODE SIMILAR TO THE ONE IN THE EXAMPLE.
#'
#' @param Pop.Mod A list containing the components returned by Population.Modeling function (main function).
#' @param Fish.years The number of recent years to estimate the mean of SEL (selectivity).
#' @param Bio.years The number of recent years to estimate the mean of M, Mat, WC, and W (natural mortality, maturity, stock weight and capture weight).
#' @param Method The procedure to obtain the age vector of weight (stock and captures), natural mortality, selectivity and maturity. By default is "mean" which means that the mean of the last "Bio.years" is used. The alternative option is "own", the user can introduce these elements.
#' @param par If Method="own" it is a list containing the age vector of weight (stock and captures), natural mortality, selectivity and maturity (for the first iteration). In other case is equal to NULL.
#' @param FM_type which of the three reference fishery mortalities must be computed. The possibilities have been described above: FM_type="F_max", FM_type="F_0.1" and FM_type="F_msy".
#' @param iters A vector containing the iteration for which the reference fishery mortalities must be computed.
#' @details The function returns the reference fishery mortality which produces maximum YPR (FM_type="F_max"), the reference fishery mortality at which the slope of the YPR curve is reduced to 0.1 of that estimated at the origin (FM_type="F_0.1") and the reference fishery mortality at which the MSY is attained (FM_type="F_msy"). Furthermore for each of these fishery mortalities the corresponding effort of fishing, YPR, BYR, B, Y and R in equilibrium is also returned. If the fishing effort is equal to 10 can be that optimize process had not found the correct value in the default sequence.
#'
#' @return One of the three following elements depending the above selection:
#' \item{F_max:}{the value of F that produces maximum YPR with the corresponding effort of fishing, YPR, BYR, B, Y and R in equilibrium.}
#' \item{F_0.1:}{the value of F at which the slope of the YPR curve is reduced to 0.1 of that estimated at the origin with the corresponding effort of fishing, YPR, BYR, B, Y and R in equilibrium.}
#' \item{F_msy:}{the value of F at which the MSY is attained with the corresponding effort of fishing, YPR, BYR, B, Y and R in equilibrium.}
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' }
#' @examples
#' # If you desire to compute the reference fishery mortalities for a large number of iterations,
#' # example 1000, the following code using the RF_U function is faster than the one using the RF
#' # function directly. The basic idea of the code is to divide the task of computing the reference
#' # fishery mortalities for 1000 iterations in blocks of 20. This is the same that the RF function do
#' # internally but it is slower than do that outside of the package.
#' #library(abind)
#' #start_time <- Sys.time()
#' #n<-1000
#' #block<-n%/%20
#' #seq1<-c(0,(1:(block-1))*20)+1
#' #seq2<-seq(20,n,by=20)
#' #ll<-length(seq1)
#' #Resul<-array(0,dim=c(1,7,n))
#' #Pop.Mod2<-Pop.Mod
#' #for(i in 1:ll){
#'  # Pop.Mod2$Matrices$N<-Pop.Mod$Matrices$N[,,seq1[i]:seq2[i]]
#'  # Pop.Mod2$Matrices$F<-Pop.Mod$Matrices$F[,,seq1[i]:seq2[i]]
#'  # Pop.Mod2$Matrices$M<-Pop.Mod$Matrices$M[,,seq1[i]:seq2[i]]
#'  # Pop.Mod2$Matrices$W<-Pop.Mod$Matrices$W[,,seq1[i]:seq2[i]]
#'  # Pop.Mod2$Matrices$Mat<-Pop.Mod$Matrices$Mat[,,seq1[i]:seq2[i]]
#'  # Pop.Mod2$Matrices$C_N<-Pop.Mod$Matrices$C_N[,,seq1[i]:seq2[i]]
#'  # Pop.Mod2$Matrices$C_W<-Pop.Mod$Matrices$C_W[,,seq1[i]:seq2[i]]
#'  # Resul[,,seq1[i]:seq2[i]]<-RF_U(Pop.Mod2, 3,3,Method="mean",par=NULL,FM_type="F_max",iters =1:20)
#'  # }
#'  # end_time <- Sys.time()
#'
#'
#' @export



RF_U<-function(Pop.Mod, Fish.years,Bio.years,Method,par,FM_type,iters){
  if(is.null(Method)){Method="mean"}

  ### SUM.POP.MOD function

  Sum.Pop.Mod<-function(Pop.Mod,Elements){
    ts<-Pop.Mod$Info$ts
    tc<-Pop.Mod$Info$tc
    if(is.numeric(Pop.Mod$Info$seed)){set.seed(Pop.Mod$Info$seed)}
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
          L[i,j,]<-stats::rnorm(niter,m,v)
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
        s<-stochastic_logistic_SEL_1(a50_Sel,ad_Sel,CV_SEL,niter,s,ages,number_years,seed=Pop.Mod$Info$seed)

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
        s<-stochastic_cte_SEL_1(cte,CV_SEL,niter,s,number_years,number_ages,seed=Pop.Mod$Info$seed)

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
        s<-stochastic_andersen_SEL_1(p1=p1,p3=p3,p4=p4,p5=p5,CV_SEL,niter,s,ages,number_years,seed=Pop.Mod$Info$seed)

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
        s<-stochastic_gamma_SEL_1(alpha,beta,gamma,CV_SEL,niter,s,ages,number_years,seed=Pop.Mod$Info$seed)

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


















  ######################################################################## useful functions



  RF_aux<-function(Pop.Mod, Fish.years,Bio.years,Method,par,FM_type,ind2){

    if (FM_type=="F_max"){
      N<-Pop.Mod$Matrices$N;ages<-as.numeric(rownames(N))

      niter<-1
      # Find F that produces max YPR: Fmax
      xx<-array(0,dim=c(1,1,niter))
      Fmax_resul<-array(0,dim = c(1,7,niter))
      colnames(Fmax_resul)<-c("f","F","YPR","BPR","R","Y","B")

      YPR_op <- function (f, Pop.Mod, Fish.years,Bio.years,plot,Method=Method,par=par,ind2=ind2) YPR_RF_taking_iteration(Pop.Mod,f,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par, ind2=ind2)[,2]

      xx<- stats::optimize(YPR_op, c(0,10), Pop.Mod=Pop.Mod, Fish.years=Fish.years,Bio.years=Bio.years,plot=FALSE,Method=Method,par=par,ind2=ind2, maximum=TRUE, tol=0.01)


      f.max<-xx$maximum

      RE<-BYR.eq_RF(Pop.Mod,f.max,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2=ind2)


      # All the values for the Fmax

      DPM<-RE$DPM; Fmax_resul[,,1]<-DPM

      return(F_max=Fmax_resul)
    }

    if (FM_type=="F_0.1"){

      # F 0.1 function
      slope<-0.1

      xx <-YPR_RF_taking_iteration(Pop.Mod,0.00001,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2=ind2)
      trgSlope <- xx[,2]/xx[,1]*slope


      yprSlope <- function (f, Pop.Mod, Fish.years,Bio.years,plot=FALSE,trgSlope,Method=Method,par=par,ind2){
        xxx <- YPR_RF_taking_iteration(Pop.Mod, c(f, f+0.0001),Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2)
        xSlope <- (xxx[2,2]-xxx[1,2])/(xxx[2,1]-xxx[1,1])
        return(trgSlope-xSlope)
      }

      xx <- stats::uniroot(yprSlope, c(0.0001,10), Pop.Mod=Pop.Mod, Fish.years=Fish.years,Bio.years=Bio.years,plot=FALSE,trgSlope,Method=Method,par=par,ind2=ind2, tol=0.01)
      f.01<-xx$root

      RE<-BYR.eq_RF(Pop.Mod,f.01,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2)


      # All the values for the F 0.1

      DPM<-RE$DPM; F01_resul<-DPM
      return(F_0.1=F01_resul)
    }

    if (FM_type=="F_msy"){
      ### Fmsy

      MSY <- function (f, Pop.Mod,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2){
        RR<-BYR.eq_RF(Pop.Mod,f,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2)
        RR$DPM[,6]
      }

      a<-MSY(seq(0.01,5,by=0.01), Pop.Mod,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2)
      ind<-which(a==0);if(length(ind)>0){b<-min(ind)};
      c<-seq(0.01,10,by=0.01)
      if(length(ind)==0){b=length(c)}
      # Return F and corresponding equilibrium values
      xx <- stats::optimize(MSY, c(0,c[b]), Pop.Mod=Pop.Mod,Fish.years=Fish.years,Bio.years=Bio.years,plot=FALSE,Method=Method,par=par,ind2=ind2,maximum=TRUE, tol=0.01)

      f.msy<-xx$maximum
      RE<-BYR.eq_RF(Pop.Mod,f.msy,Fish.years,Bio.years,plot=FALSE,Method=Method,par=par,ind2)


      # All the values for the Fmsy

      DPM<-RE$DPM; Fmsy_resul<-DPM
      return(F_msy=Fmsy_resul)

    }
  }














  BYR.eq_RF<-function(Pop.Mod,f.grid,Fish.years,Bio.years,plot,Method,par,ind2){
    ### We define useful functions



    YPR_use<-function(Pop.Mod,f.grid,Fish.years,Bio.years,min_age,max_age,plot,Method,par,ind2=ind2){
      min_age<-Pop.Mod$Info$minFage
      max_age<-Pop.Mod$Info$maxFage
      N<-Pop.Mod$Matrices$N[,,ind2];ages<-as.numeric(rownames(N))

      number_ages<-nrow(N)
      number_years<-ncol(N)

      if(Method=="mean"){
        W.c<-(Sum.Pop.Mod(Pop.Mod,c("WC"))$WC)[,,ind2];W.c<-apply(W.c[,(number_years-Bio.years+1):number_years],1,mean)
        M<-Pop.Mod$Matrices$M[,,ind2];M<-apply(M[,(number_years-Bio.years+1):number_years],1,mean)
        s<-(Sum.Pop.Mod(Pop.Mod,c("SEL"))$SEL)[,,ind2]; s<-apply(s[,(number_years-Bio.years+1):number_years],1,mean)}
      if(Method=="own"){
        W.c<-par$WC
        M<-par$M
        s<-par$SEL
      }

      l<-length(f.grid);F.vector<-1:l;ypr.vector<-1:l;Nlist<-matrix(0,ncol=l,nrow=number_ages)
      a<-1:l
      for (j in 1:l){
        F<-f.grid[j]*s #Now is a vector but later will be a matrix and a mean will be necessary
        Z<-F+M
        N<-1
        ypr<-0

        N.vec<-1:number_ages;N.vec[1]<-1
        for (i in 1:(number_ages-1)){

          ypr<-ypr+N*F[i]*W.c[i]*((1-exp(-Z[i]))/Z[i])
          N<-N*exp(-Z[i]);N.vec[i+1]<-N
        }
        Nlist[,j]<-N.vec
        ypr<-ypr+(N*F[number_ages]*W.c[number_ages]/Z[number_ages])

        ypr.vector[j]<-ypr
        F.vector[j]<-mean(F[(min_age-ages[1]+1):(max_age-ages[1]+1)])
      }
      colnames(Nlist)<-F.vector
      rownames(Nlist)<-ages
      return(Nlist)}
    #####
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
    ind<-which(R<0);R[ind]<-0
    ind<-which(B<0);B[ind]<-0
    ind<-which(Y<0);Y[ind]<-0
    results<-cbind(f=f.grid,F,YPR=ypr[,2],BPR=bpr[,2],R.eq=R,Y.eq=Y,B.eq=B)

    Nlist<-YPR_use(Pop.Mod,f.grid,Fish.years,Bio.years,min_age,max_age,plot,Method=Method,par=par,ind2=ind2)

    for (i in 1:l){
      Nlist[,i]<-Nlist[,i]*R[i]
    }
    return(list(DPM=results,N=Nlist))
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
      W<-par$W
      M<-par$M
      Mat<-par$Mat
      s<-par$SEL
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
      W.c<-par$WC
      M<-par$M
      s<-par$SEL
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


  YPR_RF_taking_iteration<-function(Pop.Mod,f.grid,Fish.years,Bio.years,plot,Method,par,ind2){
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
      W.c<-par$WC
      M<-par$M
      s<-par$SEL
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
    return(results)}











  #################################################################################################

  myfxn <- function(ind2,Fish.years, Pop.Mod,Bio.years, Method, par, FM_type){
    return(RF_aux(Pop.Mod, Fish.years,Bio.years,Method,par,FM_type,ind2))}


  n=length(iters)


  iters_list<-as.list(iters)


  numCores <- parallel::detectCores()

  cl <- parallel::makeCluster(numCores)

  parallel::clusterExport(cl,c("Logistic","myfxn","YPR_RF_taking_iteration","RF_aux","BYR.eq_RF","YPR_RF","BPR_RF","Sum.Pop.Mod","Length_VB","Weight"),
                          envir=environment())

  ret<-parallel::parLapply(cl, iters_list,myfxn, Pop.Mod=Pop.Mod, Bio.years=Bio.years,Method="mean",par=NULL,FM_type=FM_type,Fish.years=Fish.years)
  parallel::stopCluster(cl)



  return(array(unlist(ret),dim = c(1,7,n)))
}
