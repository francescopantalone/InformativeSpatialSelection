#' Simulation UY
#'
#' The dataset contains three simulated gaussian random fields, with
#' \eqn{\mu=20}, \eqn{\sigma=5}, 
#' and scale = 10.
#' 
#' @format A data frame with 121203 rows and 6 columns:
#' \describe{
#'   \item{.n}{identification number of gaussian random field}
#'   \item{x1}{first coordinate}
#'   \item{x2}{second coordinate}
#'   \item{y}{response variable}
#'   \item{Z}{proxy variable}
#'   \item{I}{indicator variable}
#'          }
"UY"

#' Simulation XYZ
#'
#' The dataset contains a simulated random gaussian field, with
#' \eqn{\mu=20}, \eqn{\sigma=5}, scale = 10, and corresponding proxy variables.
#' 
#' @format A data frame with 2500 rows and 6 columns:
#' \describe{
#'   \item{x1}{first coordinate}
#'   \item{x2}{second coordinate}
#'   \item{y}{response variable}
#'   \item{z.z1}{proxy variable}
#'   \item{z.z2}{proxy variable}
#'   \item{z.z3}{proxy variable}
#'          }
"XYZ"

#' Simulation fig4
#'
#' These .rda contains all the data used for the figure 4 of the paper.
#' In particular, we have a simulated random gaussian field
#' (as in XYZ data, ergo \eqn{\mu=20}, \eqn{\sigma=5}, scale = 10), 
#' and \eqn{M=100} samples and corresponding variograms. 
#'
"fig4"

#' Rho simulation (simulation 1)
#' 
#' \eqn{y\sim N(0, 1)}, \eqn{\epsilon\sim N(0,1)}, scale = 0.1, and
#' number of replications M=10000.
#' 
#' \code{parallelrho_server(Y = Y[[i]],
#' epsilon = epsilon[[i]],
#' .scale = 0.1,
#' paramd = unique(compatible(rbind(do.call(expand.grid,list(
#'   experiment="gamma",
#'   ell=20301,
#'   y=yexp2,
#'   .Beta=1,
#'   .Gamma=c(outer(c(1,2,5),10^{-3:3},"*")))),
#'   do.call(expand.grid,list(
#'     experiment="beta",
#'     ell=20301,
#'     y=yexp2,
#'     .Beta=c(outer(c(1,2,5),10^{-3:3},"*")),
#'     .Gamma=1)),
#'   do.call(expand.grid,list(
#'     experiment="y",
#'     ell=20301,
#'     y=yexp1,
#'     .Beta=1,
#'     .Gamma=1))))),
#' ncores = 1)}
"Numden1"
