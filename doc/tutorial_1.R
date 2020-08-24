## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width=6, fig.height=6)

## -----------------------------------------------------------------------------
library(Rfishpop)
ctrPop<-list(years=seq(1980,2020,by=1),ages=0:15,niter=2,N0=10000,minFage=4,
maxFage=7,tc=0.5,seed=NULL)

## -----------------------------------------------------------------------------
number_ages<-length(ctrPop$ages);number_years<-length(ctrPop$years)
M<-matrix(rep(0.4,number_ages*number_years),ncol = number_years)
colnames(M)<-ctrPop$years
rownames(M)<-ctrPop$ages

## -----------------------------------------------------------------------------

ctrBio<-list(M=M,CV_M=0.2, L_inf=124.5, t0=0, k=0.164, CV_L=0.2, CV_LC=0.2, a=4.5*10^(-6), b=3.1049,
           a50_Mat=3, ad_Mat=-0.5,CV_Mat=0.2)

## -----------------------------------------------------------------------------
f=matrix(rep(0.5,number_years),ncol=number_years,nrow=2,byrow=TRUE)

## -----------------------------------------------------------------------------
ctrSEL<-list(type="cte", par=list(cte=0.5),CV_SEL=0.2)

## -----------------------------------------------------------------------------
ctrSEL<-list(type="Andersen", par=list(p1=2,p3=0.2,p4=0.2,p5=40),CV_SEL=0.05)

## -----------------------------------------------------------------------------
ctrSEL<-list(type="Gamma", par=list(gamma=10,alpha=15, beta=0.03),CV_SEL=0.05)

## -----------------------------------------------------------------------------
ctrSEL<-list(type="Logistic", par=list(a50_Sel=1.5, ad_Sel=-1),CV_SEL=0.2)

## -----------------------------------------------------------------------------

ctrFish<-list(f=f,ctrSEL=ctrSEL)

## -----------------------------------------------------------------------------
CV_REC_C=0.2
SR<-list(type="cte",par=c(CV_REC_C))

## -----------------------------------------------------------------------------
a_BH=10000; b_BH=400; CV_REC_BH=0.2
SR<-list(type="BH",par=c(a_BH,b_BH,CV_REC_BH))

## -----------------------------------------------------------------------------
a_RK=10; b_RK=0.0002; CV_REC_RK=0.2
SR<-list(type="RK",par=c(a_RK,b_RK,CV_REC_RK))

## -----------------------------------------------------------------------------
Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)


## -----------------------------------------------------------------------------
N=Pop.Mod$Matrices$N

## -----------------------------------------------------------------------------
F=Pop.Mod$Matrices$F

## -----------------------------------------------------------------------------
M=Pop.Mod$Matrices$M

## -----------------------------------------------------------------------------
W=Pop.Mod$Matrices$W

## -----------------------------------------------------------------------------
Mat=Pop.Mod$Matrices$Mat

## -----------------------------------------------------------------------------
C_N=Pop.Mod$Matrices$C_N

## -----------------------------------------------------------------------------
C_W=Pop.Mod$Matrices$C_W

## -----------------------------------------------------------------------------
Pop.Mod$Info

## -----------------------------------------------------------------------------
N_D<-N [,,1]

## -----------------------------------------------------------------------------
N_D<-N [,,2]

## ----eval=FALSE---------------------------------------------------------------
#  ctrPop<-list(years=seq(1980,2020,by=1),niter=1,N0=10000,ages=0:15,minFage=4,
#  maxFage=7,tc=0.5,seed=NULL)
#  
#  Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)
#  

## ----eval=FALSE---------------------------------------------------------------
#  ctrPop<-list(years=seq(1980,2020,by=1),niter=2,N0=10000,ages=0:15,minFage=4,
#  maxFage=7,tc=0.5,seed=NULL)
#  ctrBio<-list(M=M,CV_M=0, L_inf=124.5, t0=0, k=0.164, CV_L=0, CV_LC=0, a=4.5*10^(-6), b=3.1049,
#               a50_Mat=3, ad_Mat=-0.5,CV_Mat=0)
#  
#  ctrSEL<-list(type="Logistic", par=list(a50_Sel=1.5, ad_Sel=-1),CV_SEL=0)
#  
#  f=matrix(rep(0.5,number_years),ncol=number_years,nrow=2,byrow=TRUE)
#  ctrFish<-list(f=f,ctrSEL=ctrSEL)
#  
#  CV_REC_C=0
#  SR<-list(type="cte",par=c(CV_REC_C))
#  Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)

## ----eval=FALSE---------------------------------------------------------------
#  ctrPop<-list(years=seq(1980,2020,by=1),niter=2,N0=10000,ages=0:15,minFage=4,
#  maxFage=7,tc=0.5,seed=NULL)
#  ctrBio<-list(M=M,CV_M=0, L_inf=124.5, t0=0, k=0.164, CV_L=0, CV_LC=0, a=4.5*10^(-6), b=3.1049,
#               a50_Mat=3, ad_Mat=-0.5,CV_Mat=0)
#  
#  ctrSEL<-list(type="Logistic", par=list(a50_Sel=1.5, ad_Sel=-1),CV_SEL=0.2)
#  
#  f=matrix(rep(0.5,number_years),ncol=number_years,nrow=1,byrow=TRUE)
#  ctrFish<-list(f=f,ctrSEL=ctrSEL)
#  
#  SR<-list(type="cte",par=c(CV_REC_C))
#  Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)

## -----------------------------------------------------------------------------
Z<-Sum.Pop.Mod(Pop.Mod,c("Z"))
LS<-Sum.Pop.Mod(Pop.Mod,c("LS"))

LC<-Sum.Pop.Mod(Pop.Mod,c("LC"))
WS<-Sum.Pop.Mod(Pop.Mod,c("WS"))

WSSB<-Sum.Pop.Mod(Pop.Mod,c("WSSB"))
C<-Sum.Pop.Mod(Pop.Mod,c("C"))

SEL<-Sum.Pop.Mod(Pop.Mod,c("SEL"))
BIO<-Sum.Pop.Mod(Pop.Mod,c("BIO"))


SSB<-Sum.Pop.Mod(Pop.Mod,c("SSB"))
REC<-Sum.Pop.Mod(Pop.Mod,c("REC"))

F<-Sum.Pop.Mod(Pop.Mod,c("F"))
WC<-Sum.Pop.Mod(Pop.Mod,c("WC"))


## -----------------------------------------------------------------------------
E<-selecting_units(Pop.Mod,c("C","BIO","SSB"))

## -----------------------------------------------------------------------------
p1=2;p3=0.2;p4=0.2;p5=40
ages<-0:15
SA<-andersen(x=ages,p1=p1,p3=p3,p4=p4,p5=p5)
plot(ages,SA,type="b", pch=19, col="red", main="Andersen Selectivity function")

## -----------------------------------------------------------------------------
gamma=10;alpha=15; beta=0.03
ages<-seq(0, 15, by=0.1)
SG<-gamma_SEL(x=ages,alpha=alpha,gamma=gamma,beta=beta)
plot(ages,SG,type="b", pch=19, col="red", main="Gamma Selectivity function")

## -----------------------------------------------------------------------------
a50_Sel=1.5; ad_Sel=-1
ages<-0:15
LO<-Logistic(x=ages,x50=a50_Sel,xd=ad_Sel)
plot(ages,LO,type="b", pch=19, col="red", main="Logistic Selectivity function")

## -----------------------------------------------------------------------------
a50_Mat=3; ad_Mat=-0.5
ages<-0:15
Mat<-Logistic(x=ages,x50=a50_Mat,xd=ad_Mat)
plot(ages,Mat,type="b", pch=19, col="red", main="Maturity function")

## -----------------------------------------------------------------------------
L_inf=124.5; t0=0; k=0.164;ts=0 # (ts=0 is fixed by default)
ages<-0:15
LS<-Length_VB(L_inf,k,ages+ts,t0)
plot(ages,LS,type="b", pch=19, col="red", main="Stock Length")

## -----------------------------------------------------------------------------
a=4.5*10^(-6); b=3.1049
WS<-Weight(LS,a,b)

## -----------------------------------------------------------------------------
plot(ages,WS,type="b", pch=19, col="red", main="Stock Weight")

## -----------------------------------------------------------------------------
plot(LS,WS,type="b", pch=19, col="red", main="Length-Weight")

## -----------------------------------------------------------------------------
years=seq(1980,2020,by=1)
a_BH=10000; b_BH=400
R<-RBH(SSB$SSB[,,1],a_BH,b_BH)
plot(years,R,type="b", pch=19, col="red", main="Beverton-Holt Recruitment")

## -----------------------------------------------------------------------------
a_RK=10; b_RK=0.0002
R<-RRK(SSB$SSB[,,1],a_RK,b_RK)
plot(years,R,type="b", pch=19, col="red", main="Ricker Model")

## -----------------------------------------------------------------------------
a_BH=10000; b_BH=400; CV_REC_BH=0.2
SR<-list(type="BH",par=c(a_BH,b_BH,CV_REC_BH))
Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)
steepness_value<-steepness(Pop.Mod,Fish.years=3,Bio.years=3,type="steepness",Method="mean",par=NULL)
steepness_value

## -----------------------------------------------------------------------------
a_RK=10; b_RK=0.0002; CV_REC_RK=0.2
SR<-list(type="RK",par=c(a_RK,b_RK,CV_REC_RK))
Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)
steepness_value<-steepness(Pop.Mod,Fish.years=3,Bio.years=3,type="steepness",Method="mean",par=NULL)
steepness_value

## -----------------------------------------------------------------------------
a_BH=10000; b_BH=400; CV_REC_BH=0.2
SR<-list(type="BH",par=c(a_BH,b_BH,CV_REC_BH))
Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)
parameters_value<-steepness(Pop.Mod,3,3,type="parameters",h=0.93,Method="mean",par=NULL)
parameters_value

## -----------------------------------------------------------------------------
a_RK=10; b_RK=0.0002; CV_REC_RK=0.2
SR<-list(type="RK",par=c(a_RK,b_RK,CV_REC_RK))
Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)
parameters_value<-steepness(Pop.Mod,3,3,type="parameters",h=2.34,Method="mean",par=NULL)
parameters_value

## ----eval=FALSE---------------------------------------------------------------
#  par=list(); par$W<-W; par$SEL<-SEL; par$Mat<-Mat; par$M<-M
#  steepness(Pop.Mod,3,3,type="parameters",h=2.34,Method="own",par=par)

