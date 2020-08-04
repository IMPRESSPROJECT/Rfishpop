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
