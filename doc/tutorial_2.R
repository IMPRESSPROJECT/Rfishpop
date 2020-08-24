## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width=6, fig.height=6)
knitr::opts_chunk$set(cache=TRUE)

## -----------------------------------------------------------------------------
library(Rfishpop)
ctrPop<-list(years=seq(1980,2020,by=1),niter=2,N0=15000,ages=0:15,minFage=2,
maxFage=5,tc=0.5,seed=NULL)
number_ages<-length(ctrPop$ages);number_years<-length(ctrPop$years)

Mvec=c(1,0.6,0.5,0.4,0.35,0.35,0.3,rep(0.3,9))
M<-matrix(rep(Mvec,number_years),ncol = number_years)
colnames(M)<-ctrPop$years
rownames(M)<-ctrPop$ages

ctrBio<-list(M=M,CV_M=0.2, L_inf=20, t0=-0.25, k=0.3, CV_L=0, CV_LC=0, a=6*10^(-6), b=3,
           a50_Mat=1, ad_Mat=-0.5,CV_Mat=0)

ctrSEL<-list(type="cte", par=list(cte=0.5),CV_SEL=0)

f=matrix(rep(0.5,number_years),ncol=number_years,nrow=2,byrow=TRUE)

ctrFish<-list(f=f,ctrSEL=ctrSEL)

a_BH=15000; b_BH=50; CV_REC_BH=0

SR<-list(type="BH",par=c(a_BH,b_BH,CV_REC_BH))

Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)

## -----------------------------------------------------------------------------
f.grid<-seq(0.00,0.5,by=0.01)
bpr<-BPR(Pop.Mod,f.grid,Bio.years=3,Fish.years=3,plot=c(TRUE,1),Method="mean",par=NULL)
head(bpr[,,1])

## -----------------------------------------------------------------------------

W=Pop.Mod$Matrices$W[,41,]

SEL=Sum.Pop.Mod(Pop.Mod,"SEL" )$SEL[,41,]

Mat=Pop.Mod$Matrices$Mat[,41,]

M=Pop.Mod$Matrices$M[,41,]

par=list(); par$W<-W; par$SEL<-SEL; par$Mat<-Mat; par$M<-M

bpr=BPR(Pop.Mod,f.grid,plot=c(TRUE,1),Method="own",par=par)
head(bpr[,,1])
head(bpr[,,2])

## -----------------------------------------------------------------------------
ypr<-YPR(Pop.Mod,f.grid,3,3,plot=c(TRUE,1), Method="mean",par=NULL)
head(ypr[,,1])

## -----------------------------------------------------------------------------
WC=Sum.Pop.Mod(Pop.Mod,"WC" )$WC[,41,]

SEL=Sum.Pop.Mod(Pop.Mod,"SEL" )$SEL[,41,]

M=Pop.Mod$Matrices$M[,41,]

par=list(); par$WC<-WC; par$SEL<-SEL; par$M<-M
ypr=YPR(Pop.Mod,f.grid,plot=c(TRUE,1),Method="own",par=par)
head(ypr[,,1])
head(ypr[,,2])

## -----------------------------------------------------------------------------
RE<-BYR.eq(Pop.Mod,f.grid,3,3,FALSE,Method="mean",par=NULL)
N_eq<-RE$N
DPM<-RE$DPM
### First iteration
head(DPM[,,1])
### Second iteration
head(DPM[,,2])


## -----------------------------------------------------------------------------
par=list(); par$W<-W; par$SEL<-SEL; par$Mat<-Mat; par$M<-M; par$WC=WC
RE=BYR.eq(Pop.Mod,f.grid,plot=FALSE,Method="own",par=par)
N_eq<-RE$N
DPM<-RE$DPM
### First iteration
head(DPM[,,1])
### Second iteration
head(DPM[,,2])


## -----------------------------------------------------------------------------
RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type="F_max",iters=1:2,plot=TRUE)

## -----------------------------------------------------------------------------
RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type="F_0.1",iters=1:2,plot=TRUE)

## -----------------------------------------------------------------------------
RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type="F_msy",iters=1:2,plot=TRUE)
RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type="F_msy",iters=2,plot=TRUE)

## -----------------------------------------------------------------------------
RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type="F_Crash",iters=1:2,plot=TRUE)

## -----------------------------------------------------------------------------
RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type="F_BPR",iters=1:2,plot=TRUE,prop=0.4)

## -----------------------------------------------------------------------------
RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type=c("F_Crash","F_msy"),iters=1:2,plot=TRUE)

RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type=c("F_Crash","F_msy","F_0.1","F_max","F_BPR"),iters=1:2,plot=TRUE)

RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type=c("F_Crash","F_BPR"),iters=1:2,plot=TRUE,prop = c(0.2,0.4))

## -----------------------------------------------------------------------------
par=list(); par$W<-W; par$SEL<-SEL; par$Mat<-Mat; par$M<-M; par$WC=WC
RF.Result=RF(Pop.Mod,Method="own",par=par,FM_type=c("F_Crash","F_BPR"),iters=1:2,plot=TRUE)
RF.Result

## -----------------------------------------------------------------------------
plotRF(Pop.Mod,RF.result=RF.Result, iter=2)

## -----------------------------------------------------------------------------
L.D<-Distribution.length(Pop.Mod,CV=0.2,Type="LengthS")

## -----------------------------------------------------------------------------
L.D[,,1][,1]
plot(L.D[,,1][,1], type="b", pch=19, col="red", xlab="", ylab="",main = "Distribution of stock length year 1980 iteration 1")
LS<-Sum.Pop.Mod(Pop.Mod,c("LS"))

## -----------------------------------------------------------------------------
L.D<-Distribution.length(Pop.Mod,CV=0.2,Type="LengthC")

## -----------------------------------------------------------------------------
L.D[,,1][,1]
plot(L.D[,,1][,1], type="b", pch=19, col="red", xlab="", ylab="",main = "Distribution of capture length year 1980 iteration 1")
LC<-Sum.Pop.Mod(Pop.Mod,c("LC"))

