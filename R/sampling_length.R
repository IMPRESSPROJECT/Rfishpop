Sampling_length<-function(Pop.Mod,CV,sample.size,RF.value,N,LS){

  L.D=Distribution.length.sampling(Pop.Mod=Pop.Mod,CV=CV,RF.value=RF.value,N,LS)
  column.names <- colnames(L.D)
  number_years<-dim(L.D)[2]
  niter<-dim(L.D)[3]
  Our.sample<-array(0, dim=c(sample.size, number_years, niter),dimnames=list(1:sample.size,column.names,1:niter))


  for (i in 1:niter){
    Our.sample[,,i]<-apply(L.D[,,i],2,f<-function(v,n=sample.size){random.generator(v,n)})
  }

  return(Our.sample)




}
