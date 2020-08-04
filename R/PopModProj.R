Pop.Mod.Proj=function(f.new,effort,niter,F,N,number_years,number_ages,s,Pop.Mod,Z,W_c,M){
  f=effort
  F.aux=Pop.Mod$Matrices$F
  F<-N
  fm=matrix(0,ncol=number_years+1, nrow=niter)
  for(i in 1:niter){
    
    fm[i,]=c(f[i,],f.new[i])
  }
  f=fm
  
  F[,1:number_years,]=F.aux
  for(kk in 1:niter){
    for(i in 1:number_ages){
      j=(number_years+1)
      F[i,j,kk]<-s[i,j,kk]*f[kk,j]
    }
  }
  
  
  Z=M+F
  
  C_N<-N
  C_N[,1:number_years,]=Pop.Mod$Matrices$C_N
  
  
  for(i in 1:number_ages){
    j=number_years+1
    C_N[i,j,]<-(F[i,j,]/Z[i,j,])*(N[i,j,]*(1-exp(-Z[i,j,])))
  }
  
  
  C_W<-N
  
  C_W[,1:number_years,]=Pop.Mod$Matrices$C_W
  
  for(i in 1:number_ages){
    j=(number_years+1)
    C_W[i,j,]<-C_N[i,j,]*W_c[i,j,]
  }
  
  
 
  
  return(list(Z=Z,F=F,C_W=C_W,C_N=C_N,f=f))}

