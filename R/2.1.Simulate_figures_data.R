#' Simulate data for figure 1.
#' @description
#' @return a dataframe.
Simulate_figure1_data <- function(){
  xx <- 201 
  x1<- x2<- generatex_i(xx)#This generates an arithmetic sequence  of length xx starting at 0 ending at 1.
  set.seed(42);
  theta1=1;
  theta2=c(.01,.1,1);
  plyr::ddply(data.frame(.n=1:3),.variables = ~.n,function(d){
    y = generateYonsquaregridU(length.out = xx,x1 = x1, x2 = x1, 
                               mu = 0, 
                               .sigma = theta1, 
                               .scale = theta2[unique(d$.n)])
    bindXYZpiI(x1 = x1,x2=x2,y=y,theta1=theta1,theta2=theta2[unique(d$.n)])})
} 

#' Simulate Z and samples for figure 2 data.
#' @description
#' Code that generates table  from package.
#' Simulates U and different replications of Y and Z.
#' @return a dataframe.
Simulate_figure2_data <- function(figure1_data=NULL){ 
  set.seed(9)
  xx <- 201 
  x1 <- x2<-generatex_i(xx)#This generates an arithmetic sequence  of length xx starting at 0 ending at 1.
  if(is.null(figure1_data)){figure1_data=get(data(figure1_data,package="InformativeSpatialSelection"))}
  y=figure1_data$y[figure1_data$.n==2]
  b=.4;cc=.3

  xi0 = c(log(10),log(10)-b^2/2-cc^2/2,log(10)-b^2/2-cc^2/2)
  xi1  = c(0,0,b)
  xi2 = c(0,sqrt(b^2+cc^2),cc)
  thetaparam=as.list(unique(figure1_data[figure1_data$.n==2,c("theta1","theta2")]))
  xi0text=c("log(10)","log(10)-.4^2/2-.3^2","log(10)-.4^2/2-.3^2/2")
  xi1text= paste0(character(0),xi1)
  xi2text= paste0(character(0),xi2)
  
  z <- generateZconditionnallyongridUY2(length.out = xx, 
                                        y = y, 
                                        xi0 = xi0,
                                        xi1  =xi1, 
                                        xi2 = xi2,
                                        .scale.epsilon=.1)
  UZ<-bindXYZpiI(x1,x2,y=y,y0=NA,z)
  UZ.long<-reshape2::melt(UZ,id=c("x1","x2","y"))
  #z<-10*z/mean(z)
  #range(z);mean(z)
  #range(log(z)/log(10));range((log(z)-.2*y)/log(10))
  #plot(y,log(z),type='l')
  n=10
  nstar<-rpois(1,n)
  #UZ<-plyr::rdply(3,data.frame(figure1_data[figure1_data$.n==1,c("x1","x2")],"z"=c(z)))
  #UZ[UZ$.n==1,"z"]<-10
  samplebpp2<-UZ[c("x1","x2")][pps(n,z[,2]),]
  sampleppp2<-spatstat::rpoispp(lambda =spatstat::as.im.data.frame(UZ[c("x1","x2","z2")]),nsim=1)
  samplebpp3<-UZ[c("x1","x2")][pps(n,z[,3]),]
  sampleppp3<-spatstat::rpoispp(lambda =spatstat::as.im.data.frame(UZ[c("x1","x2","z3")]),nsim=1)
  samples=
    rbind(data.frame(variable="z1",x1=runif(n),x2=runif(n),type="bpp"),
          data.frame(variable="z1",x1=runif(nstar),x2=runif(nstar),type="Ppp"),
          data.frame(variable="z2",samplebpp2,type="bpp"),
          data.frame(variable="z2",x1=sampleppp2$x,x2=sampleppp2$y,type="Ppp"),
          data.frame(variable="z3",samplebpp3,type="bpp"),
          data.frame(variable="z3",x1=sampleppp3$x,x2=sampleppp3$y,type="Ppp"))
  

  x<-UZ.long[UZ.long$x1==0&UZ.long$x2==0,]
  x$y<-min(figure1_data$y1)
  y<-x
  y$y<-max(figure1_data$y1)
  
UZ.long<-rbind(x,y,figure2_data$UZ.long)  

UZ.long$variable<-factor(UZ.long$variable,levels=c("z1","z2","z3","y0"),ordered=T)
samples$variable<-factor(samples$variable,levels=c("z1","z2","z3","y0"),ordered=T)

   figure2_data=list(UZ=UZ.long,samples=samples,
                     parameters=list(xi0=xi0,
                                     xi1=xi1,
                                     xi2=xi2,
                                     theta=thetaparam))
}

#' Simulate all X,Y,Z,the data for the project
#' @description
#' Simulates U and different replications of Y and Z.
#' @return a dataframe.

Simulate_figure3_data <- function(figure1_data=NULL){

  grid_size <- 201
  if(is.null(figure1_data)){figure1_data=get(data(figure1_data,package="InformativeSpatialSelection"))}
  y=figure1_data$y[figure1_data$.n==1]
  
  # Z <- generateZconditionnallyongridUY(length.out = grid_size, y = Y, 
  #                                      theta0 = 20, 
  #                                      xi1 = c(0, 0, 0.8), 
  #                                      .sigma = 5)
  b=.5;cc=.3
  Z <- generateZconditionnallyongridUY2(length.out = grid_size, 
                                        y = y, 
                                        theta0 = c(log(10),log(10)-b^2/2-cc^2/2,log(10)-b^2/2-cc^2/2),
                                        xi1  = c(0,0,b), 
                                        xi2 = c(0,sqrt(b^2+cc^2),cc),
                                        .scale.epsilon=.1)
 #apply(Z,2,range)
    figure3_data <- bindXYZpiI(x1 = generatex_i(grid_size), x2 = generatex_i(grid_size), y = y, z = Z)
}

#' Generates data for figure 4 of the paper
# 
#' @description
#' Data simulation for figure 4.
#' @return a list of objects necessary for figure3.
#' @usage 
#' Simulate_fig4()
Simulate_figure4_data <- function(figure1_data=NULL){
  # ------ library ------ #
  library(RandomFields)
  library(sampling)
  library(gstat)
  library(ggplot2)
  library(gridExtra)
  # --------------------- #
  # -- generate data -- #
  xx <- 201 
  x1 <- x2<-generatex_i(xx)#This generates an arithmetic sequence  of length xx starting at 0 ending at 1.
  if(is.null(figure1_data)){figure1_data=get(data(figure1_data,package="InformativeSpatialSelection"))}
  y=figure1_data$y1[figure1_data$.n==2]
  b=.4
  b=.4;cc=.3
  z <- generateZconditionnallyongridUY2(length.out = xx, 
                                        y = y, 
                                        xi0 = c(log(10)-b^2/2-cc^2/2,log(10)-b^2/2-cc^2/2),
                                        xi1  = c(0,b), 
                                        xi2 = c(sqrt(b^2+cc^2),cc),
                                        .scale.epsilon=0.1)
  XYZ<-bindXYZpiI(x1,x2,y=y,y0=NA,z)
  # ------------------- #
  # -- samples -- #
  N <- nrow(XYZ)
  n <- 100 # sample size
  M <- 1000 # number of replications
  s1_M <- matrix(data = 0, nrow = M, ncol = n)
  s2_M <- matrix(data = 0, nrow = M, ncol = n)
  s3_M <- matrix(data = 0, nrow = M, ncol = n)
  set.seed(42)
  for(i in 1:M)
  {
    s1_M[i, ] <- srs(n, XYZ, replace = F)
  }
  pi2 <- inclusionprobabilities(a = XYZ$z1, n = n)
  set.seed(42)
  for(i in 1:M)
  {
    s2_M[i, ] <- sample(x = nrow(XYZ), size = n, prob = pi2)
    # s2_M[i, ] <- pps(n, XYZ$z.z2, XYZ)
  }
  pi3 <- inclusionprobabilities(a = XYZ$z2, n = n)
  set.seed(42)
  for(i in 1:M)
  {
    s3_M[i, ] <- sample(x = nrow(XYZ), size = n, prob = pi3)
    # s3_M[i, ] <- pps(n, XYZ$z.z3, XYZ)
  }
  # ------------- #
  XYZ_so <- XYZ
  coordinates(XYZ_so) <- ~ x1 + x2 # created a spatial objects for gstat package
  # -- variogram -- #
  vgm_pop <- variogram(y ~ 1, XYZ_so) # population variogram, since we use all the field
  vgm_s1 <- list()
  length(vgm_s1) <- M
  vgm_s2 <- list()
  length(vgm_s2) <- M
  vgm_s3 <- list()
  length(vgm_s3) <- M
  for(i in 1:M) # sample variograms
  {
    vgm_s1[[i]] <- variogram(y ~ 1, XYZ_so[s1_M[i,], ])
    vgm_s2[[i]] <- variogram(y ~ 1, XYZ_so[s2_M[i,], ])
    vgm_s3[[i]] <- variogram(y ~ 1, XYZ_so[s3_M[i,], ])
  }
  # --------------- #
  # -- fit variogram -- #
  fit_vgm_pop <- fit.variogram(vgm_pop, model = vgm(1, "Gau", 1))
  fit_vgm_s1 <- list()
  length(fit_vgm_s1) <- M
  fit_vgm_s2 <- list()
  length(fit_vgm_s2) <- M
  fit_vgm_s3 <- list()
  length(fit_vgm_s3) <- M
  for(i in 1:M)
  {
    fit_vgm_s1[[i]] <- fit.variogram(vgm_s1[[i]], vgm(1, "Gau", 1))
    fit_vgm_s2[[i]] <- fit.variogram(vgm_s2[[i]], vgm(1, "Gau", 1))
    fit_vgm_s3[[i]] <- fit.variogram(vgm_s3[[i]], vgm(1, "Gau", 1))
  }
  # ------------------- #
  list(N = N, n = n, M = M, XYZ = XYZ, s1_M = s1_M, s2_M = s2_M, s3_M = s3_M, vgm_pop = vgm_pop, fit_vgm_pop = fit_vgm_pop, 
       vgm_s1 = vgm_s1, fit_vgm_s1 = fit_vgm_s1, vgm_s2 = vgm_s2,
       fit_vgm_s2 = fit_vgm_s2, vgm_s3 = vgm_s3, fit_vgm_s3 = fit_vgm_s3)
}


#' Simulate data for figure 6.
#' @description
#' @return a dataframe.
Simulate_figure6_data <- function(figure1_data=NULL){
  if(is.null(figure1_data)){figure1_data=get(data(figure1_data,package="InformativeSpatialSelection"))}
  y<-quantile(figure1_data$y,c(.75,.95,.99))
  yexp<-c(outer(c(1,2,5),10^{-3:3},"*"))
  yexp<-c(-yexp,yexp)
  Numden<-parallelrho(
      paramd=rbind(compatible(do.call(expand.grid,list(
        experiment="xi2",
        ell=20301,
        y=y,
        .Beta=1,
      .Gamma=c(outer(c(1,2,5),10^{-3:3},"*"))))),
      compatible(do.call(expand.grid,list(
      experiment="xi1",
      ell=20301,
      y=y,
      .Beta=c(outer(c(1,2,5),10^{-3:3},"*")),
      .Gamma=1))),
    compatible(do.call(expand.grid,list(
      experiment="y",
      ell=20301,
      y=yexp,
      .Beta=1,
      .Gamma=1)))),
    ncores=parallel::detectCores(),
    nreppercore=80,
    nrep= 100)
  rhofrommc(Numden)
} 




