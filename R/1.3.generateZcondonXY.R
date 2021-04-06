#' Generates a Gaussian spatial process with a Guassian semivariogram
#' @param  length.out number of x coordinates for the grid (default 100)
#' @param x1 a numeric vector (default: generatex_i(length.out))
#' @param x2 a numeric vector
#' @param y a numeric vector
#' @param xis a list of named numeric vectors. Each named vector has names "xi0","xi1","xi2",".scale.epsilon".
#' @param Variogram
#' @return a vector of length length.out 
#' @author DB
#' @example 
#' generateZconditionnallyongridUY(y=generateYonsquaregridU())
generateZconditionnallyongridUY<-
  function(length.out = 100, 
           x1=generatex_i(length.out),
           x2=x1,
           y,
           xis=list(c(xi0=5,
                     xi1=.5,
                     xi2=1,
                     .scale.epsilon=.1)),
           Variogram="Gaussian"){
    # z<-sapply(xi1,function(b){b*y+generateYonsquaregridU(length.out, x1,x2,Variogram,.scale=.scale,mu=xi0,.sigma=.sigma)})
    z<-sapply(xis,function(xi){exp(xi["xi0"] + xi["xi1"]* y + xi["xi2"] * 
        generateYonsquaregridU(length.out, x1,x2,Variogram,.scale=xi[".scale.epsilon"],mu=xi0,.sigma=1))})
    colnames(z) <- paste0("z", 1:length(xis))
    z
  }


#' Generates a Gaussian spatial process with a Guassian semivariogram
#' @param length.out number of x coordinates for the grid
#' @param y numeric vector
#' @param xi0 numeric vector
#' @param xi1 numeric vector of same length than xi0
#' @param xi2 numeric vector of same length than xi1 
#' @param .scale.epsilon numeric scale of the random process epsilon
#' @return a dataframe with two variables named x1 and x2. 
#' @author DB FP
#' @example 
#' generateZconditionnallyongridUY2(y=generateYonsquaregridU(),xi0=1:2,xi1=0:1,xi2=1:2)
generateZconditionnallyongridUY2<-
  function(length.out = 100, 
           x1=generatex_i(length.out),
           x2=x1,
           y,
           xi0=1,
           xi1=1,
           xi2=1,
           .scale.epsilon=.1,
           Variogram="Gaussian"){
    # z<-sapply(xi1,function(b){b*y+generateYonsquaregridU(length.out, x1,x2,Variogram,.scale=.scale,mu=xi0,.sigma=.sigma)})
      epsilon = generateYonsquaregridU(length.out, x1, x2,Variogram,
                                       .scale=.scale.epsilon,
                                       mu=0,
                                       .sigma=1)
      z<-exp(cbind(1,y,epsilon)%*%rbind(xi0,xi1,xi2))
    colnames(z) <- paste0("z", 1:length(xi1))
    z
  }

#' Generates a dataframe with the grid (U), Y and Z1 and Z2.
#' @param  length.out number of x coordinates for the grid
#' @return a dataframe with two variables named x1 and x2. 
#' @author DB
#' @example 
#' bindXYZpiI(x1=generatex_i(100),x2=generatex_i(100),y=generateYonsquaregridU())
bindXYZpiI<-function(x1,x2,...){
  do.call(cbind,c(list(generatesquaregridU(x1=x1, x2=x1),list(...))))}


if(FALSE){UY$Z[UY$.n == 2 & pmax(UY$x1, UY$x2) < .6] <- 2 / 0.36 # bottom-left quadrant
UY$Z[UY$.n == 2 & pmax(UY$x1, UY$x2) >= .6 & pmin(UY$x1, UY$x2) < .6] <- 2 / 0.24 # top-left quadrant & bottom-right quadrant
UY$Z[UY$.n == 2 & pmin(UY$x1, UY$x2) >= .6] <- 2 / .16 # top-right quadrant

UY$Z[UY$.n == 3] <- 10 * UY$y[UY$.n == 3] / mean(UY$y[UY$.n == 3])}
