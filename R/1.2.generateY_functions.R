#' Generates a Gaussian spatial process \eqn{Y} with a Guassian semivariogram
#' @description 
#' This function is based on \code{RandomFields::RFsimulate}. It uses the option 
#' \code{model = RandomFields::RMgauss() + RandomFields::RMtrend(mean = 0)}
#' @param length.x1 number of x coordinates for the grid, default to 100, which implies a grid of 100 * 100.
#' @param length.x2  number of x coordinates for the grid, default to 100, which implies a grid of 100 * 100.
#' @param x1 a numerical vector.
#' @param x2 a numerical vector.
#' @param n number of replications (default (1))
#' @param Variogram (default \code{"Gaussian"}): type of variogram.
#' @param .scale a numerical value (Default 10)
#' @param mu a numerical value (Default 20)
#' @param .sigma a positive numerical value (Default 1)
#' @return a vector. 
#' @author DB
#' @examples 
#' x1<-generatex_i(20)
#' Y<-generateYonsquaregridU(x1=x1)
#' UY<-cbind(x1,x1,y=Y)
#' UY2<-cbind(x1,x1,y=generateYonsquaregridU(x1=x1,n=2))
#' library(ggplot2)
#' ggplot(data.frame(generatesquaregridU(x1=x1,x2=x1),y= Y),aes(x1,x2)) + 
#'   ggplot2::geom_tile(ggplot2::aes(fill = y),size=2) + 
#'   ggplot2::xlab("")+
#'   ggplot2::ylab("")+
#'   ggplot2::scale_x_continuous(expand=c(0,0))+
#'   ggplot2::scale_y_continuous(expand=c(0,0))+
#'   ggplot2::theme(
#' axis.text.x=ggplot2::element_blank(),
#'     axis.ticks.x=ggplot2::element_blank(),
#'     axis.text.y=ggplot2::element_blank(),
#'     axis.ticks.y=ggplot2::element_blank(),
#'     panel.border = ggplot2::element_blank(),
#'     panel.grid.major = ggplot2::element_blank(), 
#'     panel.grid.minor = ggplot2::element_blank(),
#'     panel.background = ggplot2::element_blank())+
#'   ggplot2::scale_fill_gradient(low = "white",     high = "black")  +
#'   ggplot2::labs(fill = "Signal")
generateYonsquaregridU<-function(length.out = 100, 
                                 x1=generatex_i(length.out), 
                                 x2=x1,
                                 Variogram="Gaussian",
                                 .scale=.1,
                                 mu=10,
                                 .sigma=1,
                                 nrep=1){
  if(Variogram=="Gaussian"){
y <- RandomFields::RFsimulate(model = RandomFields::RMgauss() + 
                                RandomFields::RMtrend(mean = 0),
                              n=nrep,
                    x=x1/.scale,
                    y=x2/.scale)

if(class(y)=="array"){
  y<-matrix(y,c(prod(dim(y)[1:2]),dim(y)[3]))
}else{
  y<-as.matrix(y@data)
  colnames(y) <- paste0("y",1:ncol(y))}
  mu+.sigma*y}}

#' Plot Y on the grid by means of a heat map.
#' @param  UY dataframe with x1 on the first column, x2 on the second column and y on the third column.
#' @return heat map. 
#' @author DB
#' @example 
#' plotY(generateYonsquaregridU())
plotY<-function(UY){
ggplot(UY, aes(x1,x2)) + 
  geom_tile(aes(fill = y)) + 
  scale_fill_gradient(low = "white", high = "black")}

