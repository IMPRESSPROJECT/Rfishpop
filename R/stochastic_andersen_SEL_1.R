stochastic_andersen_SEL_1<-function(p1,p3,p4,p5,CV_SEL,niter,s,ages,number_years,seed){
  if(is.numeric(seed)){set.seed(seed)}
  m<-p1
  sigma<-(CV_SEL*m)


  p1_niter<-stats::rnorm(niter,m,sigma)
  if(sum(p1_niter<0)>0){ stop("A random value of p1 is negative, and this has not sense.")}

  for (ind in 2:niter){
    s[,,ind]<-matrix(rep(andersen(x=ages,p1=p1_niter[ind],p3=p3,p4=p4,p5=p5),number_years),ncol = number_years)
  }


  ### Log normal option
  #m<-p1
  #v<-(CV_SEL*m)^2
  #mu<-log(m^2/sqrt(m^2+v))
  #sigma<-sqrt(log(v/m^2+1))

  #p1_niter<-stats::rlnorm(niter,mu,sigma)


  #for (ind in 2:niter){
  #  s[,,ind]<-matrix(rep(andersen(x=ages,p1=p1_niter[ind],p3=p3,p4=p4,p5=p5),number_years),ncol = number_years)
  #}

  return(s)
}


