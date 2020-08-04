log.normal.dis<-function(x,y,max_length,CV,seed){
  if(missing(seed)) {seed=NULL}
  if(is.numeric(seed)){set.seed(seed)}
  v<-(CV*mean(x))^2
  m<-mean(x)
  if (m>0){
    mu<-log(m^2/sqrt(m^2+v))}
  sigma<-sqrt(log(v/m^2+1))


  if(y>1 &m>0&sigma>0){
    aux<-(stats::rlnorm(y, meanlog =mu, sdlog =sigma))
    ind_aux<-which(aux>max_length);aux[ind_aux]<-max_length

    # trick to have table until max_length
    aux<-c(aux,max_length)

    # trick +1 to count zero
    d<-stats::setNames(tabulate(floor(aux+1)), 0:max_length)
    # take_out artificial
    d[max_length+1]<-d[max_length+1]-1
  } else {
    x<-c(x,max_length)
    d<-stats::setNames(tabulate(floor(x+1)), 0:max_length)
    d[max_length+1]<-d[max_length+1]-1}

  return(d)
}


