stochastic_gamma_SEL_1<-function(alpha,beta,gamma,CV_SEL,niter,s,ages){

   m<-alpha
  sigma<-(CV_SEL*m)

  alpha_Sel_niter<-stats::rnorm(niter-1,m,sigma)

  res=0
  for(ind in 1:(niter-1)){
    res=c(res,gamma_SEL(x=ages,alpha=alpha_Sel_niter[ind],beta=beta,gamma=gamma))
  }
  res=res[-1]

  return(res)
}
