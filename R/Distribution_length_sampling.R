Distribution.length.sampling<-function(Pop.Mod,CV,RF.value=1000,N,LS){

  message("It is working, please give us a while :)")

  a=list()
  niter=dim(N)[3]
  number_years=dim(N)[2]
  Nc=N
  LSc=LS

  for (i in 1:niter){
    seed=Pop.Mod$Info$seed
    N<-Nc[,,i]
    LS<-LSc[,,i]
    L_inf<-Pop.Mod$Info$ctrBio$L_inf

    a[[i]]=list(seed=seed,N=N,LS=LS,L_inf=L_inf)
  }



  numCores <- parallel::detectCores()

  cl <- parallel::makeCluster(numCores)

  parallel::clusterExport(cl,c("Distribution.length.aux"),
                          envir=environment())

  ret<-parallel::parLapply(cl, a,Distribution.length.aux, CV,RF.value)
  parallel::stopCluster(cl)


  matrix.names <- 1:niter
  max_length<-round(L_inf+(2*CV*L_inf))
  lengths<-0:max_length
  column.names <- colnames(N)
  res<-array(0, dim=c(max_length+1, number_years, niter),dimnames=list(lengths,column.names,matrix.names))

  for(i in 1:niter){
    res[,,i]=ret[[i]][,,1]
  }


  return(res)



}
