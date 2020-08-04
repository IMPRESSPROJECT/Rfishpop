stochastic_cte_SEL_1<-function(alpha,CV_SEL,niter,s,number_ages){
  ### We use a Uniform(a,b)
  b<-max(alpha+0.5*sqrt(12)*CV_SEL*alpha,0)
  a<-min(alpha-0.5*sqrt(12)*CV_SEL*alpha,1)


    values_runif<-stats::runif(number_ages*(niter-1),min = a,max = b)
    if(sum(values_runif<0 & values_runif>1)>0){ stop("A random value of cte is outside of [0,1]")}


  return(values_runif)
}
