#' Gaussian Covariogram function.
#' @param h distance
#' @param .c a 2 dimension positive real vector,  .c[1] is the standard deviation, .c[2] the scale 
#' @author DB
#' @example 
#' C_Gaussian(1)
C_Gaussian <- function(h, .c = c(c0 = 1, scale = 1)){
  return(.c[1]^2 * (exp(-h^2/(.c[2]^2))))}


#' Exponential Covariogram function.
#' @param h distance
#' @param .c a 2 dimension positive real vector,  .c[1] is the standard deviation, .c[2] the scale 
#' @example 
#' C_exponential(1)
C_exponential<-function(h, .c = c(c0 = 1, scale = 1)){
  return(.c[1] ^ 2 * (exp(-h / .c[2])))}

#' Covariogram function.
#' @param h distance
#' @param type type of covariogram (spherical, exponential, Gaussian)
#' @example 
#' Covariogram.f(1)
Covariogram.f <- function(h, type = "exponential", ...){
  if(type == "spherical"){C_spherical(h, ...)}
  if(type == "exponential"){C_exponential(h, ...)}
  if(type == "Gaussian"){C_gaussian(h, ...)}}

#' Semivariogram function.
#' @param h distance
#' @param type type of covariogram (spherical, exponential, Gaussian)
#' @example 
#' Semivariogram.f(1)
Semivariogram.f <- function(h, type = "exponential", ...){
  Covariogram.f(0, type, ...) - Covariogram.f(h, type, ...)
}