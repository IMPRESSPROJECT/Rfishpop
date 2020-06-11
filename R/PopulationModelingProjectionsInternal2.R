Population.Modeling.Projections.Internal2=function(Pop.Mod,new.years,my.effort){

  # We need to obtain the effort associated to my.catch
  N=Pop.Mod$Matrices$N
  number_years=dim(N)[2]
  niter<-dim(N)[3]
  v.effort=1:niter

  f=Pop.Mod$Info$ctrFish$f


  for(i in 1:niter){
  aux.f=f[i,number_years]
  v.effort[i]=aux.f*my.effort
  }

  # We project the population

  Pop.Mod1=Population.Modeling.Projecting(Pop.Mod,new.years,f.new=v.effort)


  return(Pop.Mod1)

}
