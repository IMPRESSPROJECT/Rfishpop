rnorm.seed=function(n, mean = 0, sd = 1,seed){
  set.seed(seed)
  return(stats::rnorm(n,mean,sd))}
