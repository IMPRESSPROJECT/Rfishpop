#' @title Reference Points plots
#'
#' @description Returns a plot of the results provided by RF function (Reference Points)
#' @param Pop.Mod A list containing the components returned by Population.Modeling function (main function).
#' @param RF.result A list containing the components returned by RF function (Reference Points).
#' @param iter The iteration for which the plot must be carried out.
#' @details The function reports the plots: equilibrium biomass v. fishing mortality, equilibrium yield v. fishing mortality, equilibrium recruitment v. biomass,  and equilibrium yield v. biomass. The  reference fishery mortalities are also plotted in the previous curves. Note that only the fishery mortalities required by the argument FM_type in RF function are plotted. If prop is a vector of length greater than 1 only the values corresponding to the first one will be represented in the plot for simplicity. On the other hand, note that F_msy and F_max coincide when the recruitment relationship is constant and hence only one of both appears in the plot. The reference fishery mortalities F_max and F_Crash can be overlapped in equilibrium recruitment v. biomass and equilibrium yield v. biomass plots then only once can be shown in the plot.
#'
#'
#' @return The plot described above in Details.
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' \item{Maria Grazia Pennino}
#' }
#' @examples
#'
#' #RF.result=RF(Pop.Mod, 3,3,Method="mean",par=NULL,
#' #FM_type=c("F_Crash","F_msy","F_0.1","F_max","F_BPR")
#' #,iters=1:2,plot=TRUE)
#' # plotRF(Pop.Mod,RF.result, iter=2)
#' @export







plotRF=function(Pop.Mod,RF.result, iter){

  iters=as.numeric(dimnames(RF.result[[1]])[[3]])
  inditer=which(iters==iter)

  list_final=RF.result
  FM.aux=names(list_final)
  vec=c("F_Crash","F_msy" ,"F_0.1"  ,"F_max")
  ind=rep(0,4)
  for(i in 1:4){
    a=which(FM.aux==vec[i])
    if(length(a)>0){
    ind[i]=a}
  }
  ind=ind[ind>0]
  end_ind=ind
  l.ind=length(ind);if(l.ind==length(FM.aux)){FM_type=FM.aux}
  if(l.ind!=length(FM.aux)){FM_type=c(FM.aux[end_ind],"F_BPR")}

  ind=which(FM_type=="F_BPR")
  if(length(ind)>0){aux=FM_type[-ind]}else{aux=FM_type}
  l.aux=length(aux)


  # Saving information
  RFmax=list_final$F_max; ind_max=is.null(RFmax)
  RF01=list_final$F_0.1; ind_01=is.null(RF01)
  RFmsy=list_final$F_msy; ind_msy=is.null(RFmsy)
  RFcrash=list_final$F_Crash; ind_crash=is.null(RFcrash)

  ll=length(list_final)
  names(list_final)=rep("",ll)
  if(length(ind)>0){
    RFBPRpro=list_final[[l.aux+1]]; ind_BPR=is.null(RFBPRpro)} else {ind_BPR=TRUE;RFBPRpro=NULL}


  Fmax=RFmax[,2,inditer]
  F01=RF01[,2,inditer]
  Fmsy=RFmsy[,2,inditer]
  Fcrash=RFcrash[,2,inditer]
  FBRPpro=RFBPRpro[,2,inditer]


  fmax=RFmax[,1,inditer];if(is.null(fmax)){fmax=NA}
  f01=RF01[,1,inditer];if(is.null(f01)){f01=NA}
  fmsy=RFmsy[,1,inditer];if(is.null(fmsy)){fmsy=NA}
  fcrash=RFcrash[,1,inditer];if(is.null(fcrash)){fcrash=NA}
  fBRPpro=RFBPRpro[,1,inditer];if(is.null(fBRPpro)){fBRPpro=NA}

  bmax=RFmax[,7,inditer]
  b01=RF01[,7,inditer]
  bmsy=RFmsy[,7,inditer]
  bcrash=RFcrash[,7,inditer]
  bBRPpro=RFBPRpro[,7,inditer]

  n <- 4
  graphics::par(oma = c(4,1,1,1), mfrow = c(2, 2), mar = c(2, 2, 1, 1))



  f.ref=max(fmax,f01,fmsy,fcrash,fBRPpro,na.rm =TRUE)+0.2
  f.grid<-seq(0.00,f.ref,by=0.01)

  if(RF.result$INFO$Method=="mean"){
    Method=RF.result$INFO$Method
    Fish.years=RF.result$INFO$Fish.years
    Bio.years=RF.result$INFO$Bio.years
    par=NULL
  RE<-BYR.eq(Pop.Mod,f.grid,Fish.years=Fish.years,Bio.years=Bio.years,plot=FALSE,Method=Method,par=par)
  }

  if(RF.result$INFO$Method=="own"){
    Method=RF.result$INFO$Method
    Fish.years=3
    Bio.years=3
    par=RF.result$INFO$par
    RE<-BYR.eq(Pop.Mod,f.grid,Fish.years=Fish.years,Bio.years=Bio.years,plot=FALSE,Method=Method,par=par)
  }
  F=RE$DPM[,2,iter]



  # Plot 3


  graphics::plot(RE$DPM[,7,iter],RE$DPM[,5,iter],type="l", main="Recruitment v. biomass",xlab="",ylab="")
  type=Pop.Mod$Info$SR$type
  a=Pop.Mod$Info$SR$par[1]
  b=Pop.Mod$Info$SR$par[2]

  if(type=="BH"){
    slopeSR=RBH(SSB=0.0000001, a, b)
  }

  if(type=="RK"){
    slopeSR=RRK(SSB=0.0000001, a, b)
  }


  if(ind_max==FALSE){
    slope <- 1/RFmax[,4,inditer]
    intercept <- 0

    graphics::abline(intercept, slope,col="red")}

  if(ind_01==FALSE){

    slope <- 1/RF01[,4,inditer]
    intercept <- 0

    graphics::abline(intercept, slope,col="blue")}

  if(ind_msy==FALSE){

    slope <- 1/RFmsy[,4,inditer]
    intercept <- 0

    graphics::abline(intercept, slope,col="green")}

  if(ind_crash==FALSE){
    if(type!="cte"){
      x=c(0.0000001,bcrash);y=c( slopeSR,RFcrash[,5,inditer])
      slope <- diff(y)/diff(x)
      intercept <- y[1]-slope*x[1]

      graphics::abline(intercept, slope,col="brown")}}

  if(ind_BPR==FALSE){

    slope <- 1/RFBPRpro[,4,inditer]
    intercept <-0

    graphics::abline(intercept, slope,col="yellow")}

  graphics::points(bmax,RFmax[,5,inditer],col="red",cex=1.3,pch=19)
  graphics::points(b01,RF01[,5,inditer],col="blue",cex=1.3,pch=19)
  graphics::points(bmsy,RFmsy[,5,inditer],col="green",cex=1.3,pch=19)
  graphics::points(bcrash,RFcrash[,5,inditer],col="brown",cex=1.3,pch=19)
  graphics::points(bBRPpro,RFBPRpro[,5,inditer],col="yellow",cex=1.3,pch=19)

  # Plot 2

  DPM<-RE$DPM[,7,iter]
  graphics::plot(F,DPM,type="l", main="Biomass v. F",xlab="",ylab="")

  graphics::abline(v=Fmax,col="red")
  graphics::abline(v=F01,col="blue")
  graphics::abline(v=Fmsy,col="green")
  graphics::abline(v=Fcrash,col="brown")
  graphics::abline(v=FBRPpro,col="yellow")

  graphics::points(Fmax,RFmax[,7,inditer],col="red",cex=1.3,pch=19)
  graphics::points(F01,RF01[,7,inditer],col="blue",cex=1.3,pch=19)
  graphics::points(Fmsy,RFmsy[,7,inditer],col="green",cex=1.3,pch=19)
  graphics::points(Fcrash,RFcrash[,7,inditer],col="brown",cex=1.3,pch=19)
  graphics::points(FBRPpro,RFBPRpro[,7,inditer],col="yellow",cex=1.3,pch=19)





  # Plot 4


  graphics::plot(RE$DPM[,7,iter],RE$DPM[,6,iter],type="l", main="Yield v. biomass",xlab="",ylab="")

  graphics::abline(v=bmax,col="red")
  graphics::abline(v=b01,col="blue")
  graphics::abline(v=bmsy,col="green")
  graphics::abline(v=bcrash,col="brown")
  graphics::abline(v=bBRPpro,col="yellow")

  graphics::points(bmax,RFmax[,6,inditer],col="red",cex=1.3,pch=19)
  graphics::points(b01,RF01[,6,inditer],col="blue",cex=1.3,pch=19)
  graphics::points(bmsy,RFmsy[,6,inditer],col="green",cex=1.3,pch=19)
  graphics::points(bcrash,RFcrash[,6,inditer],col="brown",cex=1.3,pch=19)
  graphics::points(bBRPpro,RFBPRpro[,6,inditer],col="yellow",cex=1.3,pch=19)


  # plot 1
  DPM=RE$DPM[,6,iter]
  graphics::plot(F,DPM,type="l", main="Yield v. F",xlab="",ylab="")

  graphics::abline(v=Fmax,col="red")
  graphics::abline(v=F01,col="blue")
  graphics::abline(v=Fmsy,col="green")
  graphics::abline(v=Fcrash,col="brown")
  graphics::abline(v=FBRPpro,col="yellow")

  graphics::points(Fmax,RFmax[,6,inditer],col="red",cex=1.3,pch=19)
  graphics::points(F01,RF01[,6,inditer],col="blue",cex=1.3,pch=19)
  graphics::points(Fmsy,RFmsy[,6,inditer],col="green",cex=1.3,pch=19)
  graphics::points(Fcrash,RFcrash[,6,inditer],col="brown",cex=1.3,pch=19)
  graphics::points(FBRPpro,RFBPRpro[,6,inditer],col="yellow",cex=1.3,pch=19)

  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  graphics::legend('bottom',legend=c("F_max", "F_0.1","F_msy","F_Crash",
                                     FM_type[-end_ind]),
                   col=c("red", "blue","green","brown","yellow"), lwd = 5, xpd = TRUE, horiz = TRUE, cex = 0.8, seg.len=1, bty = 'n')





}



