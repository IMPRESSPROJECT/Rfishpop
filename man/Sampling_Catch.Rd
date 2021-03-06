% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling_catch.R
\name{Sampling_Catch}
\alias{Sampling_Catch}
\title{Catch numbers (for each year and age), catch weight (for each year), and sampling catch length (for each year)}
\usage{
Sampling_Catch(Pop.Mod, type, CV_CN = 0, par)
}
\arguments{
\item{Pop.Mod}{A list containing the components returned by Population.Modeling function (main function).}

\item{type}{An indicator of which element must be computed: catch weight (type="catch weight"), catch numbers (type="catch numbers"), or length catch sample (type="LengthC").}

\item{CV_CN}{The coefficient of variation associated to the catch numbers. Default value 0.}

\item{par}{A vector containing the required parameters if type="LengthC". \itemize{
\item{CV The coefficient of variation associated to the log-normal distribution used to compute the catch length distribution from which the length sample is obtained (see Details).}
\item{RF.value The number of values generated for each age (given a year and an iteration) from the log-normal distribution used to obtain the catch length distribution from which the length sample is obtained  (see Details). By default RF.value=1000.}
\item{sample.size The size of the desired sample.}}}
}
\value{
An array containing the catch numbers for each year, age, and iteration. If type="catch weight" an array containing the catch weight for each year and iteration is also reported. If type="LengthC" an array containing in each column (year) the corresponding length sample for each iteration (third dimension) is also provided.
}
\description{
Returns the catch numbers for each year, age and iteration, the catch weight for each year and iteration, and a sample of the catch length (for each year and iteration) derived from the corresponding length distribution.
}
\details{
The function returns the catch numbers for each year, age and iteration, the catch weight for each year and iteration or a catch length sample (for each year and iteration) generating random values from the computed length distribution.

CATCH NUMBERS:

When CV_CN=0, the catch numbers are equal to the matrix C_N reported by the Population.Modeling function (main function).
Whereas if CV_CN>0, then the catch number for year t and age i (\eqn{{C_{N,S}}_{it}}) is generated from a log normal distribution center in \eqn{{C_N}_{it}} and variability determined for the corresponding CV_CN.

CATCH WEIGHT:

Firstly, we compute the WC (third dimensional array containing the weight corresponding to the catch length for each age, year and iteration) at time tc (see Population.Modeling function). Then, we define C_W (third dimensional array containing the catch weight for each age, year and iteration) multiplying the catch numbers (\eqn{C_{N,S}}) by the weight WC. Finally, the catch weight for each year is the sum by columns of matrix C_W.

SAMPLE LENGTH:

A catch length sample is obtained generating random values from the computed catch length distribution.

 How the catch distribution length is computed?
 The catch distribution length from which the sampling procedure is carried out is computed as follows.
 It is computed generating for each age, year and iteration RF.value random values from a log-normal distribution centered in the corresponding catch length (LC matrix) and whose variability comes from the given CV. The distribution obtained for each age (given a year and an iteration) is scaled using the corresponding catch numbers (\eqn{C_{N,S}}).

 Remember that the matrix LC (catch length) is computed at the time instant tc (see Population.Modeling function).
}
\examples{
library(Rfishpop)
ctrPop<-list(years=seq(1980,2020,by=1),niter=2,N0=10000,ages=0:15,minFage=4,
maxFage=7,ts=0,tc=0.5,tseed=NULL)
number_ages<-length(ctrPop$ages);number_years<-length(ctrPop$years)
M<-matrix(rep(0.4,number_ages*number_years),ncol = number_years)
colnames(M)<-ctrPop$years
rownames(M)<-ctrPop$ages
ctrBio<-list(M=M,CV_M=0, L_inf=124.5, t0=-0.25, k=0.164, CV_L=0, CV_LC=0, a=4.5*10^(-6), b=3.1049,
a50_Mat=3, ad_Mat=-0.5,CV_Mat=0)
ctrSEL<-list(type="Logistic", par=list(a50_Sel=1.5, ad_Sel=-1),CV_SEL=0)
f=matrix(rep(0.5,number_years),ncol=number_years,nrow=2,byrow=TRUE)
ctrFish<-list(f=f,ctrSEL=ctrSEL)
a_RK=10; b_RK=0.0002; CV_REC_RK=0.2
SR<-list(type="RK",par=c(a_RK,b_RK,CV_REC_RK))
Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)
# Uncomments the following lines.
# Catch numbers
#CV_CN<-0
#IN=Sampling_Catch(Pop.Mod=Pop.Mod,type="catch numbers",CV_CN=CV_CN)
# Catch weight
#IW=Sampling_Survey(Pop.Mod=Pop.Mod,type="catch weight",CV_CN=CV_CN)
#IW$weight
#IW$numbers
# Length
#par=c(0.2,1000,3000)
#IL=Sampling_Survey(Pop.Mod=Pop.Mod,type="LengthC",par=par,CV_CN=CV_CN)
#IL$length
#IL$numbers
}
\author{
\itemize{
\item{Marta Cousido-Rocha}
\item{Santiago Cerviño López}
\item{Maria Grazia Pennino}
}
}
