
Log.normal<-function(mean,cv,seed){
  if(is.numeric(seed)){set.seed(seed)}
  m<-mean;n_values<-length(m)
  v<-(cv*m)^2
  mu<-log(m^2/sqrt(m^2+v))

  sigma<-sqrt(log(v/m^2+1))


  values.log.normal<-1:n_values
  for(i in 1:n_values){
    if(is.numeric(mu[i])){}else{stop("See Log-normal function, problem mu")}
    if(sigma[i]<0){stop("See Log-normal function, problem sigma")}
    values.log.normal[i]<-stats::rlnorm(1,mu[i],sigma[i])
  }

  return(values.log.normal)
}

