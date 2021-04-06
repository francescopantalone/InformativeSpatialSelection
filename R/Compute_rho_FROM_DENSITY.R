#' @description Computes rho
#' @param Pop a data.frame or matrix containing the coordinates of the population.
#' @examples
#' library(InformativeSpatialSelection)
#' x1<-generatex_i(201)
#' x2=x1
#' X<-generatesquaregridU(x1=x1,x2=x1)
#' samplesize=1000
#' nrep=3000
#' .scale=.1
#'.scale=10
#'Y=generateYonsquaregridU(x1=x1,x2=x2,n=nrep,.scale=.scale,mu=0,.sigma=1)
#'epsilon=generateYonsquaregridU(x1=x1,x2=x2,n=nrep,.scale=.scale,mu=0,.sigma=1)
#'.beta=1
#'.gamma=1
#'.summary=TRUE
#' xx<-rho.binomial.fromdensity(Y=Y,epsilon=epsilon)
#' plot(xx$xx(),xx$yy,type='l')
#' xx<-rho.binomial.fromdensity(Y=Y,epsilon=epsilon,.beta=200,.gamma=0)
#' plot(xx$xx(),xx$yy,type='l')
rho.binomial.fromdensity<-
  function(x1=generatex_i(201),
           x2=generatex_i(201),
           X=generateYonsquaregridU(x1=x1,x2=x2),
           ell,
           samplesize=1000,
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
           .beta=.1,
           .gamma=.1,
           .summary=TRUE){
    Z <-exp(.beta*Y +.gamma*epsilon)
    N<-nrow(Z)
    sample.y<-c(plyr::aaply(1:nrep,1,function(i){
      Y[sample(N,size=samplesize,prob=Z[,i],replace=T),i]}))
    xx<-seq(min(sample.y),max(sample.y),length.out=300);
    yy<-ks::kde(x = sample.y,eval.points = xx)$estimate/dnorm(xx)
    list(xx=function(){seq(min(sample.y),max(sample.y),length.out=300)},yy=yy)
  }
  
