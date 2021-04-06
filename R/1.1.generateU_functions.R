#' Function that generates a sequence of points between 0 and 1.
#' @param length.out number of x coordinates for the grid
#' @return  a dataframe with two variables named x1 and x2. 
#' @author DB
#' @example 
#' generatex_i(5)
generatex_i<-function(length.out){seq(0, 1, length.out = length.out)}

#' Simple function that generates U as a grid on the [0,1]^2 square.
#' @param length.x1 number of x coordinates for the grid (optional and ignored if x1 is given)
#' @param length.x2 number of y coordinates for the grid (optional and ignored if x2 is given)
#' @param x1 x coordinates for the grid (optional if length.x1 is given)
#' @param x2 y coordinates for the grid (optional if length.x1 is given)
#' @return a dataframe with two variables named x1 and x2. 
#' @author DB
#' @example 
#' generatesquaregridU(5)
generatesquaregridU<-function(length.x1=100,
                              length.x2=length.x1,
                              x1=generatex_i(length.x1),
                              x2=generatex_i(length.x2)){
  expand.grid(x1=x1, x2=x2)}
