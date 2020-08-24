## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width=6, fig.height=6)

## -----------------------------------------------------------------------------
library(Rfishpop)
ctrPop<-list(years=seq(1980,2020,by=1),niter=2,N0=10000,ages=0:15,minFage=4,
maxFage=7,ts=0.2,tc=0.5,tseed=NULL)

number_ages<-length(ctrPop$ages);number_years<-length(ctrPop$years)
M<-matrix(rep(0.4,number_ages*number_years),ncol = number_years)
colnames(M)<-ctrPop$years
rownames(M)<-ctrPop$ages
ctrBio<-list(M=M,CV_M=0.2, L_inf=124.5, t0=0, k=0.164, CV_L=0.2, CV_LC=0.2, a=4.5*10^(-6), b=3.1049,
           a50_Mat=3, ad_Mat=-0.5,CV_Mat=0.2)


ctrSEL<-list(type="Logistic", par=list(a50_Sel=1.5, ad_Sel=-1),CV_SEL=0.2)

f=matrix(rep(0.5,number_years),ncol=number_years,nrow=2,byrow=TRUE)
ctrFish<-list(f=f,ctrSEL=ctrSEL)

a_BH=10000; b_BH=400; CV_REC_BH=0.2; a_RK=10; b_RK=0.0002; CV_REC_RK=0.2

SR<-list(type="RK",par=c(a_RK,b_RK,CV_REC_RK))

Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)

## -----------------------------------------------------------------------------
Pop.Mod1=Population.Modeling.Projections(Pop.Mod,new.years=2021,my.catch=1500,limit.f=NULL,tol=NULL,strategy="catch")

## -----------------------------------------------------------------------------
Sum.Pop.Mod(Pop.Mod1,c("C"))

## -----------------------------------------------------------------------------
EF=Pop.Mod1$Info$ctrFish$f
colnames(EF)<-c(ctrPop$years,2021)
rownames(EF)<-c("Iter 1","Iter 2")
EF

## -----------------------------------------------------------------------------
Pop.Mod2=Population.Modeling.Projections(Pop.Mod,new.years=2021,my.effort=0.6,limit.f=NULL,tol=NULL,strategy="effort")

## -----------------------------------------------------------------------------
Sum.Pop.Mod(Pop.Mod2,c("C"))

## -----------------------------------------------------------------------------
EF=Pop.Mod2$Info$ctrFish$f
colnames(EF)<-c(ctrPop$years,2021)
rownames(EF)<-c("Iter 1","Iter 2")
EF

## -----------------------------------------------------------------------------
Pop.Mod3=Population.Modeling.Projections(Pop.Mod,new.years=c(2021,2022),my.catch=c(1500,1600),limit.f=NULL,tol=NULL,strategy="catch")

## -----------------------------------------------------------------------------
Sum.Pop.Mod(Pop.Mod3,c("C"))

## -----------------------------------------------------------------------------
Pop.Mod3b=Population.Modeling.Projections(Pop.Mod,new.years=c(2021,2022),my.catch=c(1500,1600),limit.f=NULL,tol=0.001,strategy="catch")

## -----------------------------------------------------------------------------
Sum.Pop.Mod(Pop.Mod3b,c("C"))

## -----------------------------------------------------------------------------
M<-matrix(rep(0.7,number_ages*number_years),ncol = number_years)
colnames(M)<-ctrPop$years
rownames(M)<-ctrPop$ages

Pop.Mod4=Population.Modeling.Projections(Pop.Mod,new.years=2021,my.catch=1500,limit.f=NULL,tol=NULL,strategy="catch",M.new.years=M)

## -----------------------------------------------------------------------------
Sum.Pop.Mod(Pop.Mod4,c("C"))

## -----------------------------------------------------------------------------
EF=Pop.Mod4$Info$ctrFish$f
colnames(EF)<-c(ctrPop$years,2021)
rownames(EF)<-c("Iter 1","Iter 2")
EF

