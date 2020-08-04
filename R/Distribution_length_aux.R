
Distribution.length.aux<-function(a,CV,RF.value){

  if(CV==0){ stop("The value of CV must be strictly positive")}

  if(is.numeric(a$seed)){set.seed(a$seed)}


  N<-a$N;number_years<-dim(N)[2];number_ages<-dim(N)[1];niter<-dim(N)[3]

  Lm<-a$LS

    zero.count=length(which(Lm==0))
    if(zero.count>0){warning("Individuals of zero length are not taken into account. You can avoid individuals of length 0 changing the parameter t0 (Von Bertalanffy) and/or ts parameter")}

    L_inf<-a$L_inf
    max_length<-round(L_inf+(2*CV*L_inf))
    lengths<-0:max_length
    column.names <- colnames(N)
    RS<-array(0, dim=c(max_length+1, number_years, 1),dimnames=list(lengths,column.names,1))

    for (j in 1:number_years) {
      L<-Lm[,j]

      MA=matrix(0,nrow=max_length+1,ncol=number_ages)
      for(i in 1:number_ages){

        v<-(CV*mean(L[i]))^2
        m<-mean(L[i])
        if (m>0){
          mu<-log(m^2/sqrt(m^2+v))}
        sigma<-sqrt(log(v/m^2+1))


        if(m>0&sigma>0){
          a=(stats::rlnorm(RF.value, meanlog =mu, sdlog =sigma))
          aux<-c(a,max_length)
          ind_aux<-which(aux>max_length);aux[ind_aux]<-max_length
          # trick +1 to count zero
          d<-stats::setNames(tabulate(floor(aux+1)), 0:max_length)
          # take_out artificial
          d[max_length+1]<-d[max_length+1]-1
          cor=RF.value/N[i,j]
          MA[,i]=d/cor
        }

      }
      RS[,j,1]<-apply(MA,1,sum)

    }



    return(RS)

}






