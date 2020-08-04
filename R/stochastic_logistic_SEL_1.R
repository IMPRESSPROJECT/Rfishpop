stochastic_logistic_SEL_1<-function(a50_Sel,ad_Sel,CV_SEL,niter,s,ages){
  m<-a50_Sel
  sigma<-(CV_SEL*m)


  a50_Sel_niter<-stats::rnorm(niter-1,m,sigma)

  res=0
  for(ind in 1:(niter-1)){
  res=c(res,(Logistic(x=ages,x50=a50_Sel_niter[ind],xd=ad_Sel)))
  }
res=res[-1]
return(res)
}
