random.generator<-function(v,n){
  values<-as.numeric(names(v))

  repited<-as.numeric(v)

  created.vector<-rep(values,repited)

  SAMPLE.vec<-sample(created.vector,n)


  return(SAMPLE.vec)

}
