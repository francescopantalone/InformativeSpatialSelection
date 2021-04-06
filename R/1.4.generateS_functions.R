#1: processes that creates subsets  of a discretized population
#'@examples
#' YZ=data.frame(x1=c(0:1,0:1),x2=c(0,0:1,1))
#' YZ[srs(10,YZ),]  
srs <- function(n, XYZ, replace = T)
{
  s <- sample(nrow(as.matrix(XYZ)), n, replace)
  return(s)
}

pps <- function(n, Z)
{
  pi <- sampling::inclusionprobabilities(Z, n)
  s <- which(sampling::UPmaxentropy(pi) == 1)
  return(s)
}



binomialpp<-function(n=10,UXZ,position=c("x1","x2"), z="Z",nsim){
  sample(1:nrow(UXZ),size = n,replace=T,prob=UXZ[["z"]])}

#2: processes that create coordinates.

#'@examples
#' YZ=data.frame(x1=c(0:1,0:1),x2=c(0,0:1,1))
#' uniform(10,YZ)  
uniform<-function(n,YZ,position=c("x1","x2")){matrix(
  c(runif(n,min(YZ[[position[1]]]),max(YZ[[position[1]]])),
    runif(n,min(YZ[[position[2]]]),max(YZ[[position[2]]]))), n,2)}

#' poisson process sampling spatstat::rpoispp
#' @author DB
#' @param UXZ: a dataframe
#' @param position : the names of the position variables in \code{UYZ}
#' @param z: the names of the  
#' @examples 
#' data(UY)
#' UXZ<-UY[UY$.n==3,c("x1","x2","Z")]
#' poissonppcontinuous(UXZ)

library("sp")
library("spatstat")
library("gstat")
poissonpp<-function(UXZ,position=c("x1","x2"), z="Z",nsim=4){
  #sp::coordinates(UXZ)<-as.formula(paste0("~",paste(position,collapse="+")))
  spatstat::rpoispp(lambda =spatstat::as.im.data.frame(UXZ),nsim=nsim)}

