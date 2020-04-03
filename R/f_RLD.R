f_RLD<-function(Lm,N,number_ages,number_years,niter,column.names,max_length,CV,Pop.Mod){
  if(is.numeric(Pop.Mod$Info$seed)){set.seed(Pop.Mod$Info$seed)}
  if(is.numeric(Pop.Mod$Info$seed)){List_d<-mapply(log.normal.dis,Lm,N,max_length=max_length,CV=CV,seed=Pop.Mod$Info$seed)}
  else{List_d<-mapply(log.normal.dis,Lm,N,max_length=max_length,CV=CV)}
  lengths<-0:max_length
  matrix.names <- 1:niter
  all.values<-as.numeric(unlist(List_d))
  Array_d<-array(all.values, dim=c((max_length+1), number_ages*number_years, niter))
  RLD<-array(0, dim=c(max_length+1, number_years, niter),dimnames=list(lengths,column.names,matrix.names))

  fin<-number_ages*(1:number_years)
  init<-fin-(number_ages-1)
  for (i in 1:number_years) {
    RLD[,i,]<-apply(Array_d[,init[i]:fin[i],],c(1,3),sum)
  }
  return(RLD)
}



