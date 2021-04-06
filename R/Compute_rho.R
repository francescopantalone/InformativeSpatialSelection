#' @description computes the matrix \eqn{\Sigma_{\mathbf{x}',\mathbf{x}}}
#' @details, it is based on \code{RFcovmatrix} from the package  \pkg{RandomFields}.
#' @param x: an \eqn{n\times 2} matrix
#' @param xprime: an \eqn{n'\times 2}  matrix
#' @examples
#' .scale=1
#' sim1<-as.data.frame(RandomFields::RFsimulate(model = RandomFields::RMgauss() + 
#' RandomFields::RMtrend(mean = 10),
#' n=3000,
#' x=.scale*(0:1),
#' y=.scale*(0:1)))
#' cov(t(sim1))
#' RandomFields::RFcovmatrix(
#' x=.scale*(0:1),
#' y=.scale*(0:1),
#' model = RandomFields::RMgauss() + 
#' RandomFields::RMtrend(mean = 10))
#' Sigmaxprimex.f(cbind(c(0,0,1,1),c(0,1,0,1)),.scale=10)
#' x=matrix(runif(10),5,2);
#' xprime=matrix(runif(14),7,2)
#' Sigmaxprimex.f(x,xprime,1)
#' Sigmaxprimex.fcheck(x,xprime,1)
#' ell=sample(201^2,2);
#' ell=c(1,10);
#' cov(t(Y3000[ell,]));
#' Sigmaxprimex.f(X[ell,],.scale=.1);
#' Sigmaxprimex.fcheck(X[ell,],.scale=.1)
Sigmaxprimex.fcheck<-
  function(x,xprime=NULL,
    .scale=.1){
    if(!is.null(xprime)){
    RandomFields::RFcovmatrix(
    x=rbind(x,xprime)/.scale,
    model = RandomFields::RMgauss())[1:nrow(x),nrow(x)+(1:nrow(xprime))]}else{
      RandomFields::RFcovmatrix(
        x=x/.scale,
        model = RandomFields::RMgauss())}}

Sigmaxprimex.f<-
  function(xprime,x=xprime,
           .scale=.1){
    D<-apply(x,1,function(xx){
      apply(xprime,1,function(yy){
       sum((xx-yy)^2)
      })})
    exp(-D/(.scale^2))}


#' @description Computes \eqn{\rho(\mathbf{y}\mid\mathbf{x}}.
#' @param X a \eqn{2\times n} matrix
#' @param Y a \eqn{n\times p} matrix
#' @param ell a vector of integers
#' @param y a vector of same length than \code{ell}
#' @examples
#' x1<-generatex_i(201)
#' X<-generatesquaregridU(x1=x1,x2=x1)
#' Y<-generateYonsquaregridU(x1=x1, nrep=3)
#' ell=1:2
#' y=c(10,10)
#' updateY(X,Y,ell=1:2,y=c(0,0))
#' updateY(X,Y,ell=1,y=0)
updateY<-function(X,
                  Y,
                  y,
                  ell,
                  .scale=10,
                  sigmaxx=Sigmaxprimex.f(X[ell,],.scale=.scale),
                  sigmaxprimex=Sigmaxprimex.f(X,X[ell,],.scale=.scale),
                  A=sigmaxprimex%*%solve(sigmaxx)){
    Y+A%*%(y-Y[ell,])
}



#' @description Computes rho
#' @param x1 a
#' @param x2 a 
#' @param x1 a
#' @param x1 a
#' @param x1 a
#' 
#' @examples
#' library(InformativeSpatialSelection)
#' x1<-generatex_i(201)
#' x2=x1
#' X<-generatesquaregridU(x1=x1,x2=x1)
#' y=1:2
#' nrep=300
#' .scale=.1
#' data(package="InformativeSpatialSelection")
#' data(epsilon3000)
#' Y<-Y3000
#' epsilon<-epsilon3000
#' #' #---------------
#' ell=c(1,7)
#' sigmaxx=Sigmaxprimex.f(X[ell,],.scale=.scale)
#' sigmaxprimex=Sigmaxprimex.f(X,X[ell,],.scale=.scale)
#' Yc<-updateY(X,Y,ell=1:2,y=y,sigmaxx=sigmaxx,sigmaxprimex=sigmaxprimex)
#' 
#' plotRR<-function(x,Rr){
#' library(ggplot2)
#' X<-data.frame(cbind(x,t(Rr)))
#' names(X)<-c("x","mean","sd")
#' X$ymin=X$mean-qnorm(.975)*X$sd
#' X$ymax=X$mean+qnorm(.975)*X$sd
#' ggplot(data=X,aes(x=x,y=mean))+
#' geom_point()+
#' geom_line()+
#' scale_x_continuous(trans='log10')+
#' geom_ribbon(aes(ymin=ymin,ymax=ymax,fill="red",alpha=.3))}
#' 
#' 
#' 
#' .Gamma=c(outer(c(1,2,5),10^{-2:3},"*"))
#' .Beta=c(outer(c(1,2,5),10^{-1:2},"*"))
#' 
#' Rrgamma<-sapply(.Gamma,function(.gamma){
#' rho(x1=x1,x2=x2,X=X,ell=ell,y=y,nrep=nrep,
#' Y=Y,epsilon=epsilon,Yc=Yc,.gamma=.gamma,.beta=1)
#' })
#' 
#' RrBeta<-sapply(.Beta,function(.beta){
#' rho(x1=x1,x2=x2,X=X,ell=ell,y=y,nrep=nrep,Y=Y,epsilon=epsilon,Yc=Yc,.gamma=.5,.beta=.beta)
#' })
#' 
#' plotRR(.Gamma,Rrgamma[1,])
#' plotRR(.Beta,RrBeta[1,])
#' 
rho<-function(x1=generatex_i(201),
              x2=generatex_i(201),
              X=generateYonsquaregridU(x1=x1,x2=x2),
              ell,
              y,
              nrep=300,
              .scale=10,
              Y=generateYonsquaregridU(x1=x1,
                                       x2=x2,
                                       n=nrep,
                                       .scale=.scale,
                                       mu=0,
                                       .sigma=1),
              epsilon=generateYonsquaregridU(x1=x1,
                                       x2=x2,
                                       n=nrep,
                                       .scale=.scale,
                                       mu=0,
                                       .sigma=1),
              Yc=updateY(X,Y,ell,y),
              .beta=.1,
              .gamma=.1,
              .summary=TRUE){
  Zc<-exp(.beta*Yc+.gamma*epsilon)
  Z <-exp(.beta*Y +.gamma*epsilon)
  rho=(apply(Zc[ell,,drop=FALSE],2,prod)/(apply(Zc,2,sum)^(length(ell))))/
      (apply( Z[ell,,drop=FALSE],2,prod)/(apply( Z,2,sum)^(length(ell))))
  if(.summary){
    c(mean=mean(rho),
      sd=sqrt(var(rho)/length(rho)),
      quantile(rho,probs = c(.001,.01,.05,.95,.99,.999)))}else{rho}
}


compatible<-function(dd){
  dd<-dd[sapply(dd$y,length)==sapply(dd$ell,length),]
  uniqueyell<-unique(dd[c("ell","y")])
  uniqueyell$id<-1:nrow(uniqueyell)
  dd<-merge(dd,uniqueyell)
}

#' @description Computes rho
#' @param x1 a
#' @param x2 a
#' @param X  a
#' @param paramd  a dataframe with 4 named variables: ell, y, .Beta,.Gamma
#' @param .scale a positive real number
#' @param ncores an integer
#' @param nreppercore an integer
#' @examples
#' system.time(xxx<-parallelrho(ncores=2,nreppercore= 2,nrep= 2))
#' system.time(parallelrho(ncores=8,nreppercore= 2,nrep= 2))#user system elapsed 0.197   0.299  46.693 
#' system.time(parallelrho(ncores=2,nreppercore= 8,nrep= 2))#user system elapsed 0.213   0.130  98.344 
#' system.time(parallelrho(ncores=2,nreppercore= 2,nrep= 8))#user system elapsed 0.165   0.167  98.832 
#' system.time(parallelrho(ncores=8,nreppercore= 2,nrep=32))#user system elapsed 2.782   2.207 760.096 
#' system.time(parallelrho(ncores=8,nreppercore= 8,nrep= 8))#user system elapsed 2.214   2.205 728.399 
#' system.time(parallelrho(ncores=8,nreppercore=32,nrep= 2))#user system elapsed 2.168   2.351 756.155 
#' system.time(parallelrho(ncores=8,nreppercore=20,nrep=320))#user system elapsed  
#' system.time(parallelrho(ncores=8,nreppercore=80,nrep=80))#user system elapsed  
#' system.time(parallelrho(ncores=8,nreppercore=320,nrep=20))#user system elapsed
parallelrho<-function(x1=generatex_i(201),
                      x2=generatex_i(201),
                      X=generatesquaregridU(x1=x1,x2=x2),
                      .scale=10,
                      paramd=compatible(do.call(expand.grid,list(
                                        ell=list(1:2,3,1),
                                        y=list(0:1,0,1),
                                        .Beta=c(outer(c(1,2,5),10^{-2:3},"*")),
                                        .Gamma=c(outer(c(1,2,5),10^{-2:3},"*"))))),
                      ncores=parallel::detectCores(),
                      nreppercore=10,
                      nrep=10){
  library(doParallel)
  library(plyr)
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  ZZ<-plyr::adply(1:ncores,
              1,
              function(i,.paramd,.x1,.x2,.X,..scale,.nreppercore,.nrep){
    library(InformativeSpatialSelection)
    set.seed(i)
    plyr::rdply(.nreppercore,
                (function(..x1=.x1,..x2=.x2,..X=.X,...scale=..scale,..paramd=.paramd,..nrep=.nrep){
                  Y<-generateYonsquaregridU(x1=..x1,
                                            x2=..x2,
                                            n=..nrep,
                                            .scale=...scale,
                                            mu=0,
                                            .sigma=1);
                  epsilon=generateYonsquaregridU(x1=..x1,
                                                 x2=..x2,
                                                 n=..nrep,
                                                 .scale=...scale,
                                                 mu=0,
                                                 .sigma=1);
                  plyr::ddply(..paramd,
                              ~id,
                              function(...paramd,...X,.Y,.epsilon){
                                .ell<-unlist(unique(...paramd$ell))
                                .y  <-unlist(unique(...paramd$y))
                                Yc=updateY(X = ...X,Y = .Y,ell=.ell,y=.y)
                                plyr::ddply(...paramd,
                                            ~.Beta+.Gamma,
                                            function(....paramd,..ell){
                                              Zc<-exp(....paramd$.Beta*Yc +....paramd$.Gamma*.epsilon)
                                              Z <-exp(....paramd$.Beta*.Y +....paramd$.Gamma*.epsilon)
                                              num=(apply(Zc[..ell,,drop=FALSE],2,prod)/(apply(Zc,2,sum)^(length(..ell))))
                                              den=(apply( Z[..ell,,drop=FALSE],2,prod)/(apply(Z ,2,sum)^(length(..ell))))
                                              data.frame(num=num,den=den)},..ell=.ell)
                                },...X=..X,.Y=Y,.epsilon=epsilon)
                  })())
    },
    .paramd=paramd,
    .x1=x1,
    .x2=x2,
    .X=X,
    ..scale=.scale,
    .nreppercore=nreppercore,
    .nrep=nrep,
    .parallel=TRUE)
  merge(ZZ[c("id",".Beta",".Gamma","num","den")],unique(paramd[c("id","y","ell")]),by="id",all=T)
  }


rhofrommc<-function(Numden){
ZZ<-plyr::ddply(Numden,
            ~ id+.Beta+.Gamma,
            function(d){
              nas<-sum(is.na(d$num)|is.na(d$den))
              d<-d[with(d,!is.na(d$num)&!is.na(d$den)),]
              n<-nrow(d)
              num=mean(d$num,na.rm=TRUE)
              den=mean(d$den,na.rm=TRUE)
              rho=num/den
              vv<-var(cbind(d$num,d$den))/n
              f2<-c(1/den,-num/den)
              vratio<-t(f2)%*%vv%*%f2
              data.frame(n=n,
                         rho=rho,
                         vrho=vratio)}
)
ZZ<-merge(ZZ,unique(Numden[c("y","ell","id",".Beta",".Gamma")]),by=c("id",".Beta",".Gamma"),all=TRUE)
}
