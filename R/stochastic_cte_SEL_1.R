stochastic_cte_SEL_1<-function(alpha,CV_SEL,niter,s,number_years,number_ages,seed){
  if(is.numeric(seed)){set.seed(seed)}
  ### We use a Uniform(a,b)
  b<-max(alpha+0.5*sqrt(12)*CV_SEL*alpha,0)
  a<-min(alpha-0.5*sqrt(12)*CV_SEL*alpha,1)




  for (ind in 2:niter){
    values_runif<-stats::runif(number_ages*number_years,min = a,max = b)
    if(sum(values_runif<0 & values_runif>1)>0){ stop("A random value of cte is outside of [0,1]")}

    s[,,ind]<-matrix(values_runif,ncol = number_years)
  }

  return(s)
}
