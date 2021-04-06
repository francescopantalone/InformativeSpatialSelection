## ------------------------------ ##
# rho_{12}                         #
# nrep = 1000                      #
# mu = 0, .sigma = 1, .scale = 10  #
# y values seq(-100, 100, by = 5)  #
# varying y                        #
## ------------------------------ ##
# -- server -- #
setwd("/home/francesco/ISS")
# ------------ #
# -- library -- #
library(RandomFields)
library(sampling)
library(gstat)
library(doParallel)
library(InformativeSpatialSelection)
# ------------- #
num_cores <- 1
nrep <- 1000
# x1 <- generatex_i(201)
# x2 <- generatex_i(201)
# X <- generateYonsquaregridU(x1 = x1, x2 = x2)
x1 <- generatex_i(201)
x2 <- x1
X <- generatesquaregridU(x1 = x1, x2 = x1)
Y <- list()
registerDoParallel(cores = num_cores)
Y <- foreach(i = 1:num_cores) %dopar%
  {
    set.seed(i + 42)
    return(generateYonsquaregridU(x1 = x1,
                                  x2 = x2,
                                  n = nrep,
                                  .scale = 10,
                                  mu = 0,
                                  .sigma = 1))
  }
# save(Y, file = "Y_scale01.RData")
epsilon <- list()
registerDoParallel(cores = num_cores)
epsilon <- foreach(i = 1:num_cores) %dopar%
  {
    set.seed(i + 12)
    return(generateYonsquaregridU(x1 = x1,
                                  x2 = x2,
                                  n = nrep,
                                  .scale = 10,
                                  mu = 0,
                                  .sigma = 1))
  }
# save(epsilon, file = "epsilon_scale01.RData")
#
# system("echo \"Compute_rho_server.R, data generation done.\" | mail -s \"Process update\" \"pantalone.fra@gmail.com\"")
#


position_as_function_of_i=function(i){
  N=201
  c(N*(((N-1)/2)-i)+(((N+1)/2)+i),
    N*(((N+1)/2)+i)+(((N+1)/2)-i))}

parallelrho12_server <-function(x1 = generatex_i(201),
                                x2 = generatex_i(201),
                                X = generatesquaregridU(x1=x1,x2=x2),
                                Y, epsilon,
                                .scale=10,
                                paramd = compatible(do.call(expand.grid, list(
                                  ell = lapply(c(1,5,25,50),position_as_function_of_i),
                                  y = plyr::alply(as.matrix(expand.grid(c(-100,-10,1,0,1,10,100),c(-100,-10,1,0,1,10,100))),1,identity),
                                  .Beta = 1,
                                  .Gamma = 1))),
                                ncores = 1){
  library(doParallel)
  library(plyr)
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  ZZ<-plyr::ddply(paramd,
                  ~id,
                  function(.paramd, .X, .Y, .epsilon){
                    library(InformativeSpatialSelection)
                    .ell <- unlist(unique(.paramd$ell))
                    .y <- unlist(unique(.paramd$y))
                    Yc <- updateY(X = .X, Y = .Y, ell = .ell, y = .y)
                    plyr::ddply(.paramd,
                                ~.Beta + .Gamma,
                                function(..paramd, ..ell){
                                  Zc <- exp(..paramd$.Beta * Yc + ..paramd$.Gamma * .epsilon)
                                  Z <- exp(..paramd$.Beta * .Y + ..paramd$.Gamma * .epsilon)
                                  num <- (apply(Zc[..ell, ,drop = FALSE], 2, prod) / (apply(Zc, 2, sum) ^ (length(..ell))))
                                  den <- (apply( Z[..ell, ,drop = FALSE], 2, prod) / (apply(Z , 2, sum) ^ (length(..ell))))
                                  data.frame(num = num, den = den)},..ell = .ell)
                  }, .X = X, .Y = Y, .epsilon = epsilon, .parallel = TRUE)
  merge(ZZ[c("id", ".Beta", ".Gamma", "num", "den")], unique(paramd[c("id", "y", "ell")]), by = "id", all = T)
}

Numden <- list()
for(i in 1:num_cores)
{
  print(paste(i, "start", Sys.time()))
  Numden[[i]] <- parallelrho12_server(Y = Y[[i]],
                                      epsilon = epsilon[[i]],
                                      .scale = 10,
                                      paramd = compatible(do.call(expand.grid, list(
                                        ell = lapply(c(1,5,25,50),position_as_function_of_i),
                                        y = plyr::alply(as.matrix(expand.grid(seq(-100, 100, by = 5),
                                                                              seq(-100, 100, by = 5))),
                                                        1,identity),
                                        .Beta = 1,
                                        .Gamma = 1))),
                                      ncores = 1)
  save(Numden, file = "CR_9.RData")
  print(paste(i, "end", Sys.time()))
}
system("echo \"CR_9.R done.\" | mail -s \"Process done\" \"pantalone.fra@gmail.com\"")


# load("/Users/francesco/Google Drive/InformativeSpatialSelection/CR_9.RData")
# zzz_merge <- rhofrommc(Numden[[1]])
# zzz_merge <- na.omit(zzz_merge)
# zzz_merge$cimin <- sqrt(zzz_merge$vrho) * qnorm(0.025) + zzz_merge$rho
# zzz_merge$cimax <- sqrt(zzz_merge$vrho) * qnorm(0.975) + zzz_merge$rho
# zzz_merge2 <- cbind(zzz_merge, t(as.data.frame(zzz_merge$y)), t(as.data.frame(zzz_merge$ell)))
# colnames(zzz_merge2)[11:14] <- c("y1", "y2", "ell1", "ell2")
# # dist_fun <- function(i){
# #   return((2 * i * sqrt(2)) / 201)
# # }
# dist_fun <- function(X, ell1, ell2){
#   return(sqrt(sum((X[ell1,] - X[ell2,]) ^ 2)))
# }
# x1 <- generatex_i(201)
# x2 <- x1
# X <- generatesquaregridU(x1 = x1, x2 = x1)
# d <- rep(0, nrow(zzz_merge2))
# for(i in 1:nrow(zzz_merge2))
# {
#   d[i] <- dist_fun(X, zzz_merge2$ell1[i], zzz_merge2$ell2[i])
# }
# zzz_merge2 <- cbind(zzz_merge2, d)
# 
# pdf(file = "/Users/francesco/Google Drive/InformativeSpatialSelection/rho_12_2.pdf")
# par(mfrow = c(2, 2))
# for(i in 1:length(unique(zzz_merge2$d)))
# {
#   zzz <- subset(zzz_merge2, d == unique(zzz_merge2$d)[i])
#   scatterplot3d::scatterplot3d(zzz$y1, zzz$y2, zzz$rho, pch = 19,
#                                xlab = TeX("y_{1}"), ylab = TeX("y_{2}"), zlab = TeX("rho"),
#                                main = paste("dist = ", round(unique(zzz_merge2$d)[i], 3)))
# }
# dev.off()
# par(mfrow = c(1, 1))
# 
# pdf(file = "/Users/francesco/Google Drive/InformativeSpatialSelection/rho_12_heatmap.pdf")
# par(mfrow = c(2, 2))
# for(i in 1:length(unique(zzz_merge2$d)))
# {
#   zzz <- subset(zzz_merge2, d == unique(zzz_merge2$d)[i])
#   print(ggplot2::ggplot(zzz, ggplot2::aes(y1, y2)) + 
#     ggplot2::geom_tile(ggplot2::aes(fill = rho),size=2))
# }
# dev.off()
# par(mfrow = c(1, 1))
