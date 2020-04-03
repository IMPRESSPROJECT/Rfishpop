stochastic_gamma_SEL_1<-function(alpha,beta,gamma,CV_SEL,niter,s,ages,number_years,seed){
  if(is.numeric(seed)){set.seed(seed)}
   m<-alpha
  sigma<-(CV_SEL*m)

  alpha_Sel_niter<-stats::rnorm(niter,m,sigma)


  for (ind in 2:niter){
    s[,,ind]<-matrix(rep(gamma_SEL(x=ages,alpha=alpha_Sel_niter[ind],beta=beta,gamma=gamma),number_years),ncol = number_years)
  }



  ### log normal
  #m<-alpha
  #v<-(CV_SEL*m)^2
  # mu<-log(m^2/sqrt(m^2+v))
  # sigma<-sqrt(log(v/m^2+1))
  #
  # alpha_Sel_niter<-stats::rlnorm(niter,mu,sigma)
  #
  #
  # for (ind in 2:niter){
  #   s[,,ind]<-matrix(rep(gamma(x=ages,alpha=alpha_Sel_niter[ind],beta=beta,gamma=gamma),number_years),ncol = number_years)
  # }

  return(s)
}
