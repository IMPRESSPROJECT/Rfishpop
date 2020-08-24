## -----------------------------------------------------------------------------
library(Rfishpop)
ctrPop<-list(years=seq(1980,2020,by=1),niter=1,N0=15000,ages=0:15,minFage=2,
maxFage=5,tc=0.5,seed=NULL)
number_ages<-length(ctrPop$ages);number_years<-length(ctrPop$years)

Mvec=c(1,0.6,0.5,0.4,0.35,0.35,0.3,rep(0.3,9))
M<-matrix(rep(Mvec,number_years),ncol = number_years)
colnames(M)<-ctrPop$years
rownames(M)<-ctrPop$ages

ctrBio<-list(M=M,CV_M=0, L_inf=20, t0=-0.25, k=0.3, CV_L=0, CV_LC=0, a=6*10^(-6), b=3,
           a50_Mat=1, ad_Mat=-0.5,CV_Mat=0)

ctrSEL<-list(type="cte", par=list(cte=0.5),CV_SEL=0)

f=matrix(rep(0.5,number_years),ncol=number_years,nrow=1,byrow=TRUE)

ctrFish<-list(f=f,ctrSEL=ctrSEL)

a_BH=15000; b_BH=50; CV_REC_BH=0

SR<-list(type="BH",par=c(a_BH,b_BH,CV_REC_BH))

Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)

## -----------------------------------------------------------------------------
RE<-BYR.eq(Pop.Mod,0,3,3,c(FALSE,1),Method="mean",par=NULL)
N_eq<-RE$N
N_eq

## -----------------------------------------------------------------------------
rf=RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type="F_msy",iters=1,plot=FALSE)
rf
fmsy=rf$F_msy[,1,1];fmsy

## -----------------------------------------------------------------------------
ctrPop<-list(years=seq(1980,2020,by=1),niter=2,N0=N_eq,ages=0:15,minFage=2,
maxFage=5,tc=0.5,tseed=NULL)

f=matrix(c(rep(fmsy*2,25),rep(fmsy,number_years-25)),ncol=number_years,nrow=2,byrow=TRUE)

ctrFish<-list(f=f,ctrSEL=ctrSEL)

ctrBio<-list(M=M,CV_M=0.2, L_inf=20, t0=-0.25, k=0.3, CV_L=0, CV_LC=0, a=6*10^(-6), b=3,
           a50_Mat=1, ad_Mat=-0.5,CV_Mat=0)

Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)

## -----------------------------------------------------------------------------
q_A<-matrix(1,ncol=41,nrow=16);gamma<-1;CV_A<-0
IA=Sampling_Survey(Pop.Mod=Pop.Mod,type="abundance",q_A=q_A,gamma=gamma,CV_A=CV_A,tsampling=0)

## -----------------------------------------------------------------------------
cbind(Pop.Mod$Matrices$N[,41,1],IA[,41,1])

## -----------------------------------------------------------------------------
q_A<-matrix(1,ncol=41,nrow=16);gamma<-1;CV_A<-0.2

IA=list()

for (i in 1:100){
IA[[i]]=Sampling_Survey(Pop.Mod=Pop.Mod,type="abundance",q_A=q_A,gamma=gamma,CV_A=CV_A,tsampling=0)
}


plot(1:16,Pop.Mod$Matrices$N[,41,1],type="l",xlab = "Ages",ylab="Abundance index")

for (i in 1:100){
lines(1:16,IA[[i]][,41,1],col="red")
}

lines(1:16,Pop.Mod$Matrices$N[,41,1])


## -----------------------------------------------------------------------------
a=matrix(0,ncol=100,nrow=16)

for (i in 1:100){
  a[,i]=IA[[i]][,41,1]
}

amean=rowMeans(a)

plot(1:16,Pop.Mod$Matrices$N[,41,1],type="l",xlab = "Ages",ylab="Abundance index")

lines(1:16,amean,col="red")

## -----------------------------------------------------------------------------
plot(1:16,Pop.Mod$Matrices$N[,41,2],type="l",xlab = "Ages",ylab="Abundance index")

for (i in 1:100){
  lines(1:16,IA[[i]][,41,2],col="red")
}

lines(1:16,Pop.Mod$Matrices$N[,41,2])

## -----------------------------------------------------------------------------
q_A<-matrix(1,ncol=41,nrow=16);gamma<-1;CV_A<-0
IA=Sampling_Survey(Pop.Mod=Pop.Mod,type="abundance",q_A=q_A,gamma=gamma,CV_A=CV_A,tsampling=0.4)
#IA

## -----------------------------------------------------------------------------
q_A<-matrix(1,ncol=41,nrow=16);gamma<-1;CV_A<-0
IA=Sampling_Survey(Pop.Mod=Pop.Mod,type="biomass",q_A=q_A,gamma=gamma,CV_A=CV_A,tsampling=0)

## -----------------------------------------------------------------------------
cbind(IA$biomass[,,1],Sum.Pop.Mod(Pop.Mod,c("BIO"))$BIO[,,1])
bio0=IA$biomass

## -----------------------------------------------------------------------------
q_A<-matrix(1,ncol=41,nrow=16);gamma<-1;CV_A<-0.2

IA=list()

for (i in 1:100){
  IA[[i]]=Sampling_Survey(Pop.Mod=Pop.Mod,type="biomass",q_A=q_A,gamma=gamma,CV_A=CV_A,tsampling=0)
}

plot(1:41,bio0[,,1],type="l",xlab = "Years",ylab="Biomass index")

for (i in 1:100){
  lines(1:41,IA[[i]]$biomass[,,1],col="red")
}

lines(1:41,bio0[,,1])

## -----------------------------------------------------------------------------
a=matrix(0,ncol=100,nrow=41)

for (i in 1:100){
  a[,i]=IA[[i]]$biomass[,,1]
}

amean=rowMeans(a)

plot(1:41,bio0[,,1],type="l",xlab = "Years",ylab="Biomass index")

lines(1:41,amean,col="red")

## -----------------------------------------------------------------------------
plot(1:41,bio0[,,2],type="l",xlab = "Years",ylab="Biomass index")

for (i in 1:100){
  lines(1:41,IA[[i]]$biomass[,,2],col="red")
}

lines(1:41,bio0[,,2])

## -----------------------------------------------------------------------------
q_A<-matrix(1,ncol=41,nrow=16);gamma<-1;CV_A<-0
IA=Sampling_Survey(Pop.Mod=Pop.Mod,type="biomass",q_A=q_A,gamma=gamma,CV_A=CV_A,tsampling=0.3)

IA$biomass

## -----------------------------------------------------------------------------
#IA$abundance

## -----------------------------------------------------------------------------
q_A<-matrix(1,ncol=41,nrow=16);gamma<-1;CV_A<-0
par=c(0.2,1000,3000)

## -----------------------------------------------------------------------------
our.sample1<-Sampling_Survey(Pop.Mod=Pop.Mod,type="LengthS",q_A=q_A,gamma=gamma,CV_A=CV_A,par=par,tsampling=0)

## -----------------------------------------------------------------------------
our.sample2<-Sampling_Survey(Pop.Mod=Pop.Mod,type="LengthS",q_A=q_A,gamma=gamma,CV_A=CV_A,par=par,tsampling=0.4)

## -----------------------------------------------------------------------------
plot(density(our.sample1$length[,41,1]), xlab="Length", main="Densities of length samples")
lines(density(our.sample2$length[,41,1]),col="red")

## -----------------------------------------------------------------------------
plot(density(our.sample1$length[,41,2]), xlab="Length", main="Densities of length samples")
lines(density(our.sample2$length[,41,2]),col="red")

## -----------------------------------------------------------------------------
#our.sample2$length

## -----------------------------------------------------------------------------
#our.sample2$abundance

## -----------------------------------------------------------------------------
IC=Sampling_Catch(Pop.Mod=Pop.Mod,type="catch numbers",CV_CN=0)

## -----------------------------------------------------------------------------
cbind(Pop.Mod$Matrices$C_N[,41,1],IC[,41,1])

## -----------------------------------------------------------------------------
IC=list()

for (i in 1:100){
IC[[i]]=Sampling_Catch(Pop.Mod=Pop.Mod,type="catch numbers",CV_CN=0.2)
}

plot(1:16,Pop.Mod$Matrices$C_N[,41,1],type="l",xlab = "Ages",ylab="Catch numbers")

for (i in 1:100){
lines(1:16,IC[[i]][,41,1],col="red")
}

lines(1:16,Pop.Mod$Matrices$C_N[,41,1])

## -----------------------------------------------------------------------------
a=matrix(0,ncol=100,nrow=16)

for (i in 1:100){
  a[,i]=IC[[i]][,41,1]
}

amean=rowMeans(a)

plot(1:16,Pop.Mod$Matrices$C_N[,41,1],type="l",xlab = "Ages",ylab="Catch numbers")

lines(1:16,amean,col="red")

## -----------------------------------------------------------------------------
plot(1:16,Pop.Mod$Matrices$C_N[,41,2],type="l",xlab = "Ages",ylab="Catch numbers")

for (i in 1:100){
  lines(1:16,IC[[i]][,41,2],col="red")
}

lines(1:16,Pop.Mod$Matrices$C_N[,41,2])

## -----------------------------------------------------------------------------
IC=Sampling_Catch(Pop.Mod=Pop.Mod,type="catch weight",CV_CN=0)

## -----------------------------------------------------------------------------
cbind(IC$weight[,,1],Sum.Pop.Mod(Pop.Mod,c("C"))$C[,,1])

## -----------------------------------------------------------------------------
IC=list()

for (i in 1:100){
  IC[[i]]=Sampling_Catch(Pop.Mod=Pop.Mod,type="catch weight",CV_CN=0.2)
}

plot(1:41,Sum.Pop.Mod(Pop.Mod,c("C"))$C[,,1],type="l",xlab = "years", ylab="Catch weight")

for (i in 1:100){
  lines(1:41,IC[[i]]$weight[,,1],col="red")
}

lines(1:41,Sum.Pop.Mod(Pop.Mod,c("C"))$C[,,1])

## -----------------------------------------------------------------------------
a=matrix(0,ncol=100,nrow=41)

for (i in 1:100){
  a[,i]=IC[[i]]$weight[,,1]
}

amean=rowMeans(a)

plot(1:41,Sum.Pop.Mod(Pop.Mod,c("C"))$C[,,1],type="l",xlab = "years", ylab="Catch weight")

lines(1:41,amean,col="red")

## -----------------------------------------------------------------------------
plot(1:41,Sum.Pop.Mod(Pop.Mod,c("C"))$C[,,2],type="l",xlab = "years", ylab="Catch weight")

for (i in 1:100){
  lines(1:41,IC[[i]]$weight[,,2],col="red")
}

lines(1:41,Sum.Pop.Mod(Pop.Mod,c("C"))$C[,,2])

## -----------------------------------------------------------------------------
IC=Sampling_Catch(Pop.Mod=Pop.Mod,type="catch weight",CV_CN=0)
#IC$weight
#IC$numbers

## -----------------------------------------------------------------------------
par=c(0.2,1000,3000)

## -----------------------------------------------------------------------------

IC=Sampling_Catch(Pop.Mod=Pop.Mod,type="LengthC",CV_CN=0,par=par)

q_A<-matrix(1,ncol=41,nrow=16);gamma<-1
IAS=Sampling_Survey(Pop.Mod=Pop.Mod,type="LengthS",q_A=q_A,gamma=gamma,CV_A=0,par=par,tsampling=0)

plot(density(IAS$length[,41,1]),xlab="Length", main="Densities of length samples")
lines(density(IC$length[,41,1]),col="red")

## -----------------------------------------------------------------------------
par=c(0.2,1000,3000)
IC=Sampling_Catch(Pop.Mod=Pop.Mod,type="LengthC",CV_CN=0.2,par=par)

IAS=Sampling_Survey(Pop.Mod=Pop.Mod,type="LengthS",q_A=q_A,gamma=gamma,CV_A=0.2,par=par,tsampling=0)

plot(density(IAS$length[,41,1]),xlab="Length", main="Densities of length samples")
lines(density(IC$length[,41,1]),col="red")

## -----------------------------------------------------------------------------
#IC$length

## -----------------------------------------------------------------------------
#IC$numbers

