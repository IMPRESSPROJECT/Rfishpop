


Pop.Mod.proj.internal=function(Pop.Mod,new.years,my.catch,tol,limit.f,M.new.year){
ts<-Pop.Mod$Info$ts
tc<-Pop.Mod$Info$tc

seed=Pop.Mod$Info$seed
if(is.numeric(seed)){
  set.seed(seed)}

Sel_type=Pop.Mod$Info$ctrFish$ctrSEL$type

type<-Pop.Mod$Info$SR$type

type<-Pop.Mod$Info$SR$type
if(type=="BH"){CV_REC_BH<-Pop.Mod$Info$SR$par[3]}
if(type=="RK"){CV_REC_RK<-Pop.Mod$Info$SR$par[3]}
if(type=="cte"){CV_REC_C<-Pop.Mod$Info$SR$par[1]}


CV_LC<-Pop.Mod$Info$ctrBio$CV_LC
CV_SEL<-Pop.Mod$Info$ctrFish$ctrSEL$CV_SEL
CV_L<-Pop.Mod$Info$ctrBio$CV_L
CV_M<-Pop.Mod$Info$ctrBio$CV_M
CV_Mat<-Pop.Mod$Info$ctrBio$CV_Mat

niter<-dim(Pop.Mod$Matrices$N)[3]


years<-as.numeric(colnames(Pop.Mod$Matrices$N))
ages<-as.numeric(rownames(Pop.Mod$Matrices$N))
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
if (type=="BH"){
  a_BH<-Pop.Mod$Info$SR$par[1]
  b_BH<-Pop.Mod$Info$SR$par[2]
}

if (type=="RK"){
  a_RK<-Pop.Mod$Info$SR$par[1]
  b_RK<-Pop.Mod$Info$SR$par[2]
}



number_years<-length(years)
number_ages<-length(ages)

column.names <- c(years,new.years)
row.names <- ages
matrix.names <- 1:niter

N<-array(rep(0,number_ages*(number_years+1)), dim=c(number_ages, number_years+1, niter),dimnames=list(row.names,column.names,
                                                                                                                                    matrix.names))
N[,1:number_years,]=Pop.Mod$Matrices$N


# SELECTIVITY

if(Sel_type=="Logistic"){

  a50_Sel<-Pop.Mod$Info$ctrFish$ctrSEL$par$a50_Sel
  ad_Sel<-Pop.Mod$Info$ctrFish$ctrSEL$par$ad_Sel

  ### Deterministic
  aux_SEL=Sum.Pop.Mod(Pop.Mod,c("SEL"))$SEL
  sd=aux_SEL[,,1]
  m_aux=matrix((Logistic(x=ages,x50=a50_Sel,xd=ad_Sel)), ncol=1,nrow=number_ages)
  sd<-cbind(sd,m_aux)

  s<-N

  if(CV_SEL>0){
    j=number_years+1
    s[,j,-1]=stochastic_logistic_SEL_1(a50_Sel,ad_Sel,CV_SEL,niter,s,ages)
    s[,1:number_years,]=aux_SEL
    s[,,1]<-sd

  }}

if(Sel_type=="cte"){

  cte<-Pop.Mod$Info$ctrFish$ctrSEL$par$cte

  ### Deterministic
  aux_SEL=Sum.Pop.Mod(Pop.Mod,c("SEL"))$SEL
  sd=aux_SEL[,,1]
  m_aux=matrix(rep(cte,number_ages), ncol=1,nrow=number_ages)
  sd<-cbind(sd,m_aux)

  s<-N

  if(CV_SEL>0){
    j=number_years+1
    s[,j,-1]=stochastic_cte_SEL_1(cte,CV_SEL,niter,s,number_ages)
    s[,1:number_years,]=aux_SEL
    s[,,1]<-sd

  }}

if(Sel_type=="Andersen"){

  p1<-Pop.Mod$Info$ctrFish$ctrSEL$par$p1;p3<-Pop.Mod$Info$ctrFish$ctrSEL$par$p3
  p4<-Pop.Mod$Info$ctrFish$ctrSEL$par$p4;p5<-Pop.Mod$Info$ctrFish$ctrSEL$par$p5

  ### Deterministic
  aux_SEL=Sum.Pop.Mod(Pop.Mod,c("SEL"))$SEL
  sd=aux_SEL[,,1]
  m_aux=matrix(andersen(x=ages,p1=p1,p3=p3,p4=p4,p5=p5), ncol=1,nrow=number_ages)
  sd<-cbind(sd,m_aux)

  s<-N

  if(CV_SEL>0){
    j=number_years+1
    s[,j,-1]=stochastic_andersen_SEL_1(p1=p1,p3=p3,p4=p4,p5=p5,CV_SEL,niter,s,ages)
    s[,1:number_years,]=aux_SEL
    s[,,1]<-sd

  }}

if(Sel_type=="Gamma"){

  alpha<-Pop.Mod$Info$ctrFish$ctrSEL$par$alpha
  gamma<-Pop.Mod$Info$ctrFish$ctrSEL$par$gamma
  beta<-Pop.Mod$Info$ctrFish$ctrSEL$par$beta

  ### Deterministic
  aux_SEL=Sum.Pop.Mod(Pop.Mod,c("SEL"))$SEL
  sd=aux_SEL[,,1]
  m_aux=matrix(gamma_SEL(x=ages,alpha=alpha,gamma=gamma,beta=beta), ncol=1,nrow=number_ages)
  sd<-cbind(sd,m_aux)

  s<-N

  if(CV_SEL>0){
    j=number_years+1
    s[,j,-1]=stochastic_gamma_SEL_1(alpha,beta,gamma,CV_SEL,niter,s,ages)
    s[,1:number_years,]=aux_SEL
    s[,,1]<-sd

  }}


if(niter==1 & CV_SEL==0){s[,,1]<-sd}
if(niter>1 & CV_SEL==0) {s[,,1:niter]<-sd}




### NATURAL MORTALITY
m_aux=matrix(M.new.year, ncol=1,nrow=number_ages)

Md<-cbind(M[,,1],m_aux)


### Stochastic (Log normal distribution; CV_M)
M<-N

M[,1:number_years,]=Pop.Mod$Matrices$M
if(CV_M>0){
  for (i in 1:number_ages){
    j=(number_years+1)
      m<-Md[i,j]
      v<-(CV_M*m)^2
      mu<-log(m^2/sqrt(m^2+v))
      sigma<-sqrt(log(v/m^2+1))
      M[i,j,]<-stats::rlnorm(niter,mu,sigma)
    }
  M[,,1]<- Md
}
### Deterministic
if(niter==1 & CV_M==0){M[,,1]<- Md}
if(niter>1 & CV_M==0) {M[,,1:niter]<- Md}


### MORTALITY


Z.aux=Sum.Pop.Mod(Pop.Mod,c("Z"))$Z
Z=N
Z[,1:number_years,]=Z.aux


### LENGTH
L<-N

m_aux=matrix((Length_VB(L_inf,k,ages+ts,t0)), ncol=1,nrow=number_ages)
aux_L=Sum.Pop.Mod(Pop.Mod,c("LS"))$LS

Ld<-cbind(aux_L[,,1],m_aux)


### Stochastic (Normal CV_L)

L[,1:number_years,]=aux_L
if(CV_L>0){
  for (i in 1:number_ages){
    j=(number_years+1)
      m<-Ld[i,j]
      v<-(CV_L*m)
      L[i,j,]<-stats::rnorm(niter,m,v)
    }
  L[,,1]<-Ld
}
### Deterministic
if(niter==1 & CV_L==0){L[,,1]<-Ld}
if(niter>1 & CV_L==0) {L[,,1:niter]<-Ld}


### WEIGHTS

W<-N
W[,1:number_years,]=Pop.Mod$Matrices$W

i=number_years+1
W[,i,]<-Weight(L[,i,],a,b)


### MATURITY (Log normal distribution CV_Mat)
Mat<-N
### Deterministic

m_aux=matrix(Logistic(x=ages,x50=a50_Mat,xd=ad_Mat), ncol=1,nrow=number_ages)
Matd<-cbind(Pop.Mod$Matrices$Mat[,,1],m_aux)


Mat[,1:number_years,]=Pop.Mod$Matrices$Mat

if(CV_Mat>0){
  m<-a50_Mat
  v<-(CV_Mat*m)^2
  mu<-log(m^2/sqrt(m^2+v))
  sigma<-sqrt(log(v/m^2+1))
  a50_Mat_niter<-matrix(stats::rlnorm((niter-1)*1,mu,sigma),ncol=1)



  j=number_years+1
    for (ind in 2:niter){
      Mat[,j,ind]<-Logistic(x=ages,x50=a50_Mat_niter[ind-1,1],xd=ad_Mat)
    }
  Mat[,,1]<-Matd
}


if(niter==1 & CV_Mat==0){Mat[,,1]<-Matd}
if(niter>1 & CV_Mat==0) {Mat[,,1:niter]<-Matd}

### We assume for now that the maturity at age 0 is 0. Then we make the required
### changes to have effectively 0 proportion of maturity at age 0.

Mat[1,,]<-rep(0, number_years+1)




### SSB
row.names=""
SSB<-array(rep(0,(number_years+1)), dim=c(1, number_years+1, niter),dimnames=list(row.names,column.names,matrix.names))

SSB[,1:number_years,]=Sum.Pop.Mod(Pop.Mod,c("SSB"))$SSB

WM=N;WMA=N

WM[,1:number_years,]=Sum.Pop.Mod(Pop.Mod,c("WS"))$WS

WMA[,1:number_years,]=Sum.Pop.Mod(Pop.Mod,c("WSSB"))$WSSB


j=number_years+1
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
      if(CV_REC_BH>0){ N[1,j,-1]<-Log.normal(RBH(SSB[,j,-1],a_BH,b_BH),CV_REC_BH)
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
    }}

j=number_years+1
i=1
WM[i,j,]<-N[i,j,]*W[i,j,]



### LENGTH CAPTURES
L_c<-N

m_aux=matrix((Length_VB(L_inf,k,ages+tc,t0)), ncol=1,nrow=number_ages)

LC.aux=Sum.Pop.Mod(Pop.Mod,c("LC"))$LC

L_cd<-cbind(LC.aux[,,1],m_aux)

L_c[,1:number_years,]=LC.aux
### Stochastic (Normal CV_LC)

if(CV_LC>0){
  for (i in 1:number_ages){
      j=(number_years+1)
      m<-L_cd[i,j]
      v<-(CV_LC*m)
      L_c[i,j,]<-stats::rnorm(niter,m,v)
    }
  L_c[,,1]<-L_cd
}

if(niter==1 & CV_LC==0){L_c[,,1]<-L_cd}
if(niter>1 & CV_LC==0) {L_c[,,1:niter]<-L_cd}



W_c<-N
W_c[,1:number_years,]=Sum.Pop.Mod(Pop.Mod,c("WC"))$WC


i=(number_years+1)
W_c[,i,]<-Weight(L_c[,i,],a,b)



#### Loading Pop.Mod.Pro.Effort
Pop.Mod.Pro.Effort=function(f.new,iter,effort,niter,number_years,number_ages,s,Pop.Mod,my.catch,W_c,M,N){
  f=effort
  F.aux=Pop.Mod$Matrices$F
  j=(number_years+1)
  AUXIL=M
  AUXIL[,j,]=0

  F<-AUXIL
  fm=matrix(0,ncol=number_years+1, nrow=niter)
  i=iter

  fm[i,]=c(f[i,],f.new)

  f=fm

  F[,1:number_years,]=F.aux
  kk=iter
  for(i in 1:number_ages){
    j=(number_years+1)
    F[i,j,kk]<-s[i,j,kk]*f[kk,j]
  }



  Z=M+F

  C_N<-AUXIL
  C_N[,1:number_years,]=Pop.Mod$Matrices$C_N


  for(i in 1:number_ages){
    j=number_years+1
    C_N[i,j,]<-(F[i,j,]/Z[i,j,])*(N[i,j,]*(1-exp(-Z[i,j,])))
  }


  C_W<-AUXIL

  C_W[,1:number_years,]=Pop.Mod$Matrices$C_W

  for(i in 1:number_ages){
    j=(number_years+1)
    C_W[i,j,]<-C_N[i,j,]*W_c[i,j,]
  }




  year_C_W<-colSums(C_W[,,iter])



  return(abs(year_C_W[number_years+1]-my.catch))}

### Looking effort

effort=f

fun2=function(iter,effort,niter,number_years,number_ages,s,Pop.Mod,my.catch,W_c,M,N,tol, limit.f){
  #xx<- stats::uniroot(fun, c(0,limit.f),Pop.Mod=Pop.Mod,new.years=new.years,my.catch=my.catch, tol=tol,iter=iter,SEL=SEL,LS=LS,SSB=SSB,WS=WS,WSSB=WSSB,LC=LC,WC=WC)
  #effort=xx$root
  xx<- stats::optimize(Pop.Mod.Pro.Effort, c(0,limit.f),iter=iter,effort=effort,niter=niter,number_years=number_years,number_ages=number_ages,s=s,Pop.Mod=Pop.Mod,my.catch=my.catch,W_c=W_c,M=M,N=N,tol=tol)
  effort=xx$minimum
  return(effort)}
if(is.null(tol)){tol=0.01}
if(is.null(limit.f)){limit.f=4}
iters_list<-as.list(1:niter)


numCores <- parallel::detectCores()

cl <- parallel::makeCluster(numCores)

parallel::clusterExport(cl,c("fun2","Pop.Mod.Pro.Effort"),
                        envir=environment())

ret<-parallel::parLapply(cl, iters_list,fun2,effort=effort,niter=niter,number_years=number_years,number_ages=number_ages,s=s,Pop.Mod=Pop.Mod,my.catch=my.catch,W_c=W_c,M=M,N=N,tol=tol, limit.f=limit.f)
parallel::stopCluster(cl)



v.effort=unlist(ret)


a=Pop.Mod.Proj(f.new=v.effort,effort,niter,F,N,number_years,number_ages,s,Pop.Mod,Z,W_c,M)

Pop.Mod$Info$ctrFish$f=a$f
results<-list()
results[[1]]<-list(N=N,M=M,F=a$F,W=W,Mat=Mat,C_N=a$C_N,C_W=a$C_W)
results[[2]]<-list(ctrBio=Pop.Mod$Info$ctrBio,ctrFish=Pop.Mod$Info$ctrFish,SR=Pop.Mod$Info$SR,minFage=min_age,maxFage=max_age,seed=seed,ts=ts,tc=tc)
names(results)<-c("Matrices","Info")
return(results)}
