stochastic_andersen_SEL_1<-function(p1,p3,p4,p5,CV_SEL,niter,s,ages){
  m<-p1
  sigma<-(CV_SEL*m)


  p1_niter<-stats::rnorm(niter-1,m,sigma)
  if(sum(p1_niter<0)>0){ stop("A random value of p1 is negative, and this has not sense.")}

  res=0
  for(ind in 1:(niter-1)){
    res=c(res,andersen(x=ages,p1=p1_niter[ind],p3=p3,p4=p4,p5=p5))
  }
  res=res[-1]



  return(res)
}


