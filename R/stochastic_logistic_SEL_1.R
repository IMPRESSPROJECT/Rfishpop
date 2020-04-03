stochastic_logistic_SEL_1<-function(a50_Sel,ad_Sel,CV_SEL,niter,s,ages,number_years,seed){
  if(is.numeric(seed)){set.seed(seed)}
  m<-a50_Sel
  sigma<-(CV_SEL*m)


  a50_Sel_niter<-stats::rnorm(niter,m,sigma)

  for (ind in 2:niter){
    s[,,ind]<-matrix(rep(Logistic(x=ages,x50=a50_Sel_niter[ind],xd=ad_Sel),number_years),ncol = number_years)
  }

### log normal
# m<-a50_Sel
# v<-(CV_SEL*m)^2
# mu<-log(m^2/sqrt(m^2+v))
# sigma<-sqrt(log(v/m^2+1))
#
# a50_Sel_niter<-stats::rlnorm(niter,mu,sigma)
#
#
# for (ind in 2:niter){
#   s[,,ind]<-matrix(rep(Logistic(x=ages,x50=a50_Sel_niter[ind],xd=ad_Sel),number_years),ncol = number_years)
# }

return(s)
}
