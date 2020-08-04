

Population.Modeling.tsampling<-function(Pop.Mod,tsampling){

  seed=Pop.Mod$Info$seed
  if(is.numeric(seed)){
  set.seed(seed)}

  CV_L<-Pop.Mod$Info$ctrBio$CV_L

  old.N=Pop.Mod$Matrices$N
  niter<-dim(old.N)[3]
  years<-as.numeric(colnames(old.N))
  ages<-as.numeric(rownames(old.N))



  L_inf<-Pop.Mod$Info$ctrBio$L_inf
  t0<-Pop.Mod$Info$ctrBio$t0
  k<-Pop.Mod$Info$ctrBio$k

  number_years<-length(years)
  number_ages<-length(ages)


  ### FISHING MORTALITY

  F=Pop.Mod$Matrices$F

  ### NATURAL MORTALITY
  M=Pop.Mod$Matrices$M


  ### MORTALITY

  Z<-M+F

  N=old.N

  ### LENGTH
  L<-N

  Ld<-matrix(rep(Length_VB(L_inf,k,ages+tsampling,t0),number_years), ncol=number_years,nrow=number_ages)


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


  for(i in 1:number_ages){
    for(j in 1:number_years){
      N[i,j,]<-(old.N[i,j,]*(exp(-Z[i,j,]*tsampling)))
    }
  }


  list=list(N=N,LS=L)
  return(list)}
