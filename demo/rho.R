# --------------------------------- #
#1. Set parameters values

alpha <- 10
mu <- 10
c0 <- sqrt(5)     #standard deviation
.scale <- 0.1

# --------------------------------- #
# --------------------------------- #
#2. Specific functions for our example
#'
C_h <- function(h, .c){C_exponential(h, .c)}



#'@param h a vector
#'@param ... additional  arguments for C_h
#'@example 
#'A_h(1,.sd=sqrt(5),.scale=.1)
#'A_h(c(1,2),.sd=sqrt(5),.scale=.1)
A_h <- function(h, C_h = C_exponential, .drop = TRUE, ...)
{
  Chs <- C_h(h, ...)
  C0 <- C_h(0, ...)
  B <- matrix(c(1, 1, 1, -1), 2, 2)
  plyr::aaply(Chs, 1, function(cc){
    B %*% diag(sqrt(c(C0 + cc, C0 - cc))) %*% B / 2 
  }, .drop = .drop)
}

#'
#' @param h a vector
#' @param ... additional  arguments for C_h
#' @example
#' h = c(1, 2); Ah = A_h(h); Bh = B_h(h, Ah = Ah); Bh
B_h <- function(h, Ah = NULL, C_h = C_exponential, .drop=TRUE, ...){    
  if(is.null(Ah)){Ah <- A_h(h, C_h, .drop = FALSE, ...)}
  m <- matrix(c(1, 0, 0, -1), 2, 2)
  if(length(dim(Ah)) == 2){Ah %*% m %*% Ah}else{
    plyr::aaply(Ah, 1, function(A){A %*% m %*% A}, .drop = .drop)}}
if(FALSE){
h <- 0.5
G_h <- rep(0, length(H))
for(h in 1:length(H))
{
  G <- 0
  for(i in 1:niter)
  {
    G <- G + t(y) %*% B_h(h, c0, .scale) %*% y * (((A_h(h, c0, .scale) %*% y)[1] * (A_h(h, c0, .scale) %*% y)[2]) / (mu ^ 2))
  }
  G_h[h] <- G / (2 * niter)
}
}

################################
#3.? Definition of K
F_Gauss <- function(x, .c2 = c2){(pnorm(x / .c2) - 1 / 2) * sqrt(2 * pi)}
K1 <- function(x, .c2 = c2){F_Gauss(x, .c2) + F_Gauss(1 - x, .c2)}
K_Gauss <- function(x_ell, .c0 = c0, .c2 = c2){.c0 * K1(x_ell[1], .c2) * K1(x_ell[2], .c2)}



################################
#4.1  Plot of K_Gauss
x <- seq(0, 1, length = 100)
y <- seq(0, 1, length = 100)
x1 <- rep(x, length(y))
x2 <- rep(y, each = length(x))
grid <- data.frame(x1, x2)
c0 <- 5
c2 <- 0.1

par(mfrow=c(1,4))
plot3D::scatter3D(grid$x1, grid$x2, 
                  plyr::daply(grid, ~ x1 + x2, function(K){K_Gauss(c(K$x1, K$x2), 1, .001)}))
plot3D::scatter3D(grid$x1, grid$x2, plyr::daply(grid, ~ x1 + x2, function(K){K_Gauss(c(K$x1, K$x2), 1, .1)}))
plot3D::scatter3D(grid$x1, grid$x2, plyr::daply(grid, ~ x1 + x2, function(K){K_Gauss(c(K$x1, K$x2), 1, 1)}))
plot3D::scatter3D(grid$x1, grid$x2, plyr::daply(grid, ~ x1 + x2, function(K){K_Gauss(c(K$x1, K$x2), 1, 1000)}))
par(mfrow = c(1, 1))

# grid_h <- matrix(0, 1, 2)
# h <- 0.5
# for(i in 1:nrow(grid))
# {
#   for(j in 1:nrow(grid))
#   {
#     # if((sqrt((grid[i, 1] - grid[j, 1]) ^ 2 + (grid[i, 2] - grid[j, 2]) ^ 2)) == h)
#     if(all.equal((sqrt((grid[i, 1] - grid[j, 1]) ^ 2 + (grid[i, 2] - grid[j, 2]) ^ 2)), h))
#     {
#       print("ok")
#       grid_h[i, ] <- rbind(cbind(grid[i, 1], grid[i, 2]))
#     }
#   }
# }

h <- 0.5
.scale <- 0.01
sigma_m <- matrix(c(C_h(0, c(c0, .scale)), -C_h(h, c(c0, .scale)), -C_h(h, c(c0, .scale)), C_h(0, c(c0, .scale))), nrow = 2, ncol = 2) / (C_h(0, c(c0, .scale)) ^ 2 - C_h(h, c(c0, .scale)) ^ 2)

# rho12 approximation number 1
#' @param y1 signal first unit
#' @param y2 signal second unit
#' @param beta model parameter
#' @param mu model parameter
#'@examples
#' y1=11;y2=9;alpha=10;beta=10; mu=10;
#' rho12.approx1(x1, x2, y1, y2, alpha, beta, mu)
rho12.approx1 <- function(y1, y2, alpha, beta, mu)
{
  num <- (alpha + beta * y1) * (alpha + beta * y2)
  den <- alpha ^ 2 + 2*alpha * beta * mu + (beta * mu) ^ 2
  num / den
}

#' # rho approximation number 2
#' @param x1 coordinates first unit
#' @param x2 coordinates second unit
#' @param y1 signal first unit
#' @param y2 signal second unit
#' @param beta model parameter
#' @param mu model parameter
#' @param .c vector of c0 (standard deviation) and .scale (scale), for covariogram computation
#' @examples
#' x1=rep(1/2,2)-h/(2*sqrt(2))*c(1,1)
#' x2=rep(1/2,2)+h/(2*sqrt(2))*c(1,1)
#' y1=11;y2=9;alpha=10;beta=10; mu=10;c0=sqrt(5);.scale=.1;
#' rho12(x1, x2, y1, y2, alpha, beta, mu, .c = c(c0,.scale))
# rho2 <- function(x1, x2, y1, y2, alpha, beta, mu, .c)
# {
#   h=sqrt(sum((x1-x2)^2))
#   inv.sigma <- matrix(c(C_h(0, .c), -C_h(h, .c), -C_h(h, .c), C_h(0, .c)), nrow = 2, ncol = 2) / (C_h(0, .c) ^ 2 - C_h(h, .c) ^ 2)
#   num <- (alpha + beta * y1) * (alpha + beta * y2)
#   den <- (alpha + beta *(mu + c(K_Gauss(c(1/2, 1/2), c0, c2), K_Gauss(c(1/2, 1/2), c0, c2)) %*% inv.sigma  %*% (c(y1, y2) - mu)) ^ 2)
#   num / den
# }

# rho12 approximation number 2
#' Specific to our case.
#' @param x1 coordinates first unit
#' @param x2 coordinates second unit
#' @param y1 signal first unit
#' @param y2 signal second unit
#' @param beta model parameter
#' @param mu model parameter
#' @param .c vector of c0 (standard deviation) and .scale (scale), for covariogram computation
#' @examples
#' x1=rep(1/2,2)-h/(2*sqrt(2))*c(1,1)
#' x2=rep(1/2,2)+h/(2*sqrt(2))*c(1,1)
#' y1=11;y2=9;alpha=10;beta=10; mu=10;c0=sqrt(5);.scale=.1;
#' rho12(x1, x2, y1, y2, alpha, beta, mu, .c = c(c0,.scale))
rho12.approx2 <- function(x1, x2, y1, y2, alpha, beta, mu,.c)
{
  h <- sqrt(sum((x1-x2)^2))
  inv.sigma <- matrix(c(C_h(0, .c), -C_h(h, .c), -C_h(h, .c), C_h(0, .c)), nrow = 2, ncol = 2) / (C_h(0, .c) ^ 2 - C_h(h, .c) ^ 2)
  inv.sigma <- solve(matrix(c(C_Gaussian(0, .c), 
                              C_Gaussian(h, .c), 
                              C_Gaussian(h, .c), 
                              C_Gaussian(0, .c)), nrow = 2, ncol = 2))
  num <- (alpha + beta * y1) * (alpha + beta * y2)
  den <- (alpha + beta *(mu + t(c(K_Gauss(x1, c0, c2), 
                                  K_Gauss(x2, c0, c2))) %*% inv.sigma  %*% (c(y1, y2) - mu))) ^ 2
  num / den
}

#den <- alpha + beta *(mu + t(c(K_Gauss(x1,c0,c2), K_Gauss(x2,c0,c2))) %*% inv.sigma  %*% (c(y2, y2) - mu))

# rho12 approximation number 2: denominator
#'@examples
#' @param x1 coordinates first unit
#' @param x2 coordinates second unit
#' @param y1 signal first unit
#' @param y2 signal second unit
#' @param beta model parameter
#' @param mu model parameter
#' @param .c vector of c0 (standard deviation) and .scale (scale), for covariogram computation
#' @examples
#' x1=rep(1/2,2)-h/(2*sqrt(2))*c(1,1)
#' x2=rep(1/2,2)+h/(2*sqrt(2))*c(1,1)
#' y1=11;y2=9;alpha=10;beta=10; mu=10;c0=sqrt(5);.scale=.1;
#' denom12.approx2(x1, x2, y1, y2, alpha, beta, mu, .c = c(c0,.scale))
denom12.approx2 <- function(x1, x2, y1, y2, alpha, beta, mu, .c)
{
  h <- sqrt(sum((x1 - x2) ^ 2))
  inv.sigma <- solve(matrix(c(C_Gaussian(0, .c), 
                              C_Gaussian(h, .c), 
                              C_Gaussian(h, .c), 
                              C_Gaussian(0, .c)), nrow = 2, ncol = 2))
  den <- (alpha + beta *(mu + t(c(K_Gauss(x1, c0, c2), 
                                  K_Gauss(x2, c0, c2))) 
                         %*% inv.sigma  
                         %*% (c(y1, y2) - mu))) ^ 2
  den
}

# -- #
h <- 0.5
y1 <- rnorm(1, 10)
y2 <- rnorm(1, 10)
beta <- 10
rho12_a1_p <- rho12.approx1(y1, y2, alpha, beta, mu)

# rho12_a1_p <- plyr::raply(1000,
#                       (function(){
#                         x1 <- runif(2)
#                         x2 <- runif(2)
#                         rho12.approx1(y1, y2, alpha, beta, mu)
#                       })())
# -- #
# -- #
K <- function(x_1, x_2, c0, c2)
{
  c0 * (1 - sqrt(2 * pi) * pnorm (x_1 / c2) + sqrt(2 * pi) * pnorm ((1 - x_1) / c2)) * (sqrt(2 * pi) * pnorm (x_2 / c2) + sqrt(2 * pi) * pnorm ((1 - x_2) / c2))
}

y1 <- rnorm(1, 10)
y2 <- rnorm(1, 10)
beta <- 10
c0 <- 5
c2 <- 0.1
rho12_a2_p <- plyr::raply(1000,
                      (function(){
                        x1 <- runif(2)
                        x2 <- runif(2)
                        h <- sqrt((x1[1] - x2[1])^2 + (x1[2] - x2[2])^2)
                        rho12.approx2(x1, x2, y1, y2, alpha, beta, mu, .c = c(c0, .scale))
                      })())

hist(rho12_a2_p, xlab = "rho approximation number 2", main = "")
# -- #
h <- .25
rho12_a2_p2 <- plyr::raply(1000,
                      (function(){
                        x1 <- runif(2)
                        x <- FALSE
                        while(x){theta <- runif(2 * pi)
                        x2 <- x1 + c(cos(theta), sin(theta)) * h
                        x <- all(x2 > 0) & all(x1 < 0)}
                        #                        h <- sqrt((x1[1] - x2[1])^2 + (x1[2] - x2[2])^2)
                        y <- matrix(c(5, Covariogram.f(h), Covariogram.f(h), 5), 2, 2) %*% rnorm(2, 10)
                        c(y, rho12.approx2(x1, x2, y[1], y[2], alpha, beta, mu, .c = c(c0, .scale)))
                      })())


hist(rho12_a2_p2, xlab = "rho", main = "")
plot3D::scatter3D(rho12_a2_p2[,1], rho12_a2_p2[,2], rho12_a2_p2[,3])
h <- .25
y <- qnorm(seq(.01, .99, length.out = 99),mean = mu)
x1 <- runif(2)
x <- FALSE
while(x){theta<-runif(2*pi)
x2 <- x1 + c(cos(theta), sin(theta)) * h
x <- all(x2 > 0) & all(x1 < 0)}
x1 <- rep(1/2, 2) - h / (2 * sqrt(2)) * c(1, 1)
x2 <- rep(1/2, 2) + h / (2 * sqrt(2)) * c(1, 1)

rho43_p <- plyr::aaply(y,1,
                       function(y1){
                         plyr::aaply(y,1,function(y2){
                           rho12(x1, x2, y1, y2, alpha, beta, mu, c(c0, .scale))           })})

rho43_d <- plyr::adply(y,1,
                       function(y1)
                       {
                         plyr::aaply(y,1,function(y2){
                           c(y1=y1,y2=y2,rho=rho12(x1, x2, y1, y2, alpha, beta, mu, c(c0, .scale))           )})})
plot3D::scatter3D(rho43_d$y1,rho43_d$y2,rho43_d$rho)

rho43_denom <- plyr::adply(y,1,
                           function(y1){
                             plyr::aaply(y,1,function(y2){
                               c(y1=y1,y2=y2,rho=denom12(x1, x2, y1, y2, alpha, beta, mu, .c = c(c0,.scale))           )})})
plot3D::scatter3D(rho43_d$y1,rho43_d$y2,rho43_denom$rho)



rho43_dapprox <- plyr::adply(y,1,
                             function(y1){
                               plyr::aaply(y,1,function(y2){
                                 c(y1=y1,y2=y2,rho=rho1( y1, y2, alpha, beta, mu,c0,.scale)           )})})
plot3D::scatter3D(rho43_d$y1,rho43_d$y2,rho43_dapprox$rho)
plot(rho43_d$y1,rho43_d$y2)
plot3D::scatter3D(rho43_d$y1,rho43_d$y2,rho43_d$rho)


heatmap(rho43_p,Rowv=NA,Colv=NA)
hist(rho1_p, xlab = "rho", main = "")
plot3D::scatter3D(rho1_p[, 1], rho1_p[, 2], rho1_p[, 3])

totoh <- plyr::raply(1000,
                     (function(){
                       x1 <- runif(2)
                       alpha <- runif(1, min = 0, max = 2 * pi)
                       x2 <- x1 + h * c(cos(alpha), sin(alpha))
                       h <- sqrt((x1[1] - x2[1])^2 + (x1[2] - x2[2])^2)
                       if((max(x2) > 1) | min(x2) < 0){NA}else{
                         rho2(x1, x2, y1, y2, alpha, beta, mu, .c = c(c0, .scale))}
                     })())

totoh <- totoh[!is.na(totoh)]
hist(totoh)


x2 <- c(.5, .5)
xx <- seq(0, 1, length.out = 100)
grid <- expand.grid(xx,xx)
tutu <- plyr::aaply(grid,1,
                    function(x1){
                      h <- sqrt((x1[1] - x2[1])^2 + (x1[2] - x2[2])^2)
                      rho2(x1, x2, y1, y2, alpha, beta, mu, h)
                    })
plot3D::scatter3D(grid[,1],grid[,2],tutu)


x <- plyr::rdply(1000,
                 (function(){
                   x <- runif(4)
                   h <- sqrt(sum((x[1:2] - x[3:4]) ^ 2))
                   rh <- rho2(x[1:2], x[3:4], y1, y2, alpha, beta, mu, h)
                   c(h = h, rh = rh)
                 })())

plot(x$h, x$rh)



# ex1

# h<-.5
#           
#           
# rho<-function(x1,x2,y1,y2)
#   {
#               x1[1]*x2[2]*(y1-10)*(y2-10)
#           }
#           
#           h=.5
#           y1<-rnorm(1,10)
#           y2<-rnorm(1,10)
#           
#           
#           toto<-plyr::raply(1000,
#                             (function(){
#                               x1<-runif(2)
#                               x2<-runif(2)
#                               rho(x1,x2,y1,y2)
#                             })())
#           
#           hist(toto)
#           
#           h<-.5
#           
#           totoh<-plyr::raply(1000,
#                              (function(){
#                                x1<-runif(2)
#                                alpha=runif(1,min = 0,max=2*pi)
#                                x2<-x1+h*c(cos(alpha),sin(alpha))
#                                if((max(x2)>1)|min(x2)<0){NA}else{
#                                  rho(x1,x2,y1,y2)}
#                              })())
#           
#           totoh<-totoh[!is.na(totoh)]
#           hist(totoh)
#           
#           x2<-c(.5,.5)
#           
#           xx<-seq(0,1,length.out = 100)
#           grid<-expand.grid(xx,xx)
#           tutu<-plyr::aaply(grid,1,
#                             function(x1){
#                               rho(x1,x2,y1,y2)
#                             })
#           plot3D::scatter3D(grid[,1],grid[,2],tutu)
#           
#           # ex2
#           
#           x=plyr::rdply(1000,
#                         (function(){
#                           x=runif(4)
#                           h=sqrt(sum((x[1:2]-x[3:4])^2))
#                           rh=rho(x[1:2],x[3:4],y1,y2)
#                           c(h=h,rh=rh)
#                         })())
#           
#           )
#                   plot(x$h,x$rh)
#                   
#                   
#                   
#                   set.seed(42)
#                   prc_y <- RFsimulate(model = RMgauss(c0 = c0, .scale = 0.1) + RMtrend(mean =mu), x = seq(0, 1, length = 100), y = seq(0, 1, length = 100))
#                   prc_z <- prc_y
#                   prc_z@data <- prc_y@data * alpha
#                   
#                   x <- seq(0, 1, length = 100)
#                   y <- seq(0, 1, length = 100)
#                   x1 <- rep(x, length(y))
#                   x2 <- rep(y, each = length(x))
#                   prc <- data.frame(x1, x2, y = prc_y@data)
#                   prc_grid <- expand.grid(seq(0, 1, length = 20), seq(0, 1, length = 20)) # grid
#                   colnames(prc_grid) <- c("x1", "x2")
#                   
#                   # d2 <- 
#                   # y <- mvrnorm(n = 100, mu = c(mu, mu), d2)
#                   
#                   niter <- 1000
#                   set.seed(42)
#                   y <- matrix(nrow = 2, ncol = niter, rnorm(2 * niter))
#                   H <- c(seq(0, 0.7, length.out = 30), seq(0.71, 1, length.out = 20))