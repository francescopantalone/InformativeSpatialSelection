## ------------------------------ ##
# rho_{1}                          #
# nrep = 10000                     #
# mu = 0, .sigma = 1, .scale = 0.1 #
# varying beta, gamma, y           #
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
num_cores <- 5
nrep <- 2000
# x1 <- generatex_i(201)
# x2 <- generatex_i(201)
# X <- generateYonsquaregridU(x1 = x1, x2 = x2)
x1 <- generatex_i(201)
x2 <- x1
X <- generatesquaregridU(x1 = x1, x2 = x1)
Y <- list()
registerDoParallel(cores = num_cores)
set.seed(42)
Y <- foreach(i = 1:num_cores) %dopar%
  {
    return(generateYonsquaregridU(x1 = x1,
                                  x2 = x2,
                                  n = nrep,
                                  .scale = 0.1,
                                  mu = 0,
                                  .sigma = 1))
  }
# save(Y, file = "Y_scale01.RData")
epsilon <- list()
registerDoParallel(cores = num_cores)
set.seed(12)
epsilon <- foreach(i = 1:num_cores) %dopar%
  {
    return(generateYonsquaregridU(x1 = x1,
                                  x2 = x2,
                                  n = nrep,
                                  .scale = 0.1,
                                  mu = 0,
                                  .sigma = 1))
  }
# save(epsilon, file = "epsilon_scale01.RData")
#
# system("echo \"Compute_rho_server.R, data generation done.\" | mail -s \"Process update\" \"pantalone.fra@gmail.com\"")
#

parallelrho_server <-function(x1 = generatex_i(201),
                              x2 = generatex_i(201),
                              X = generatesquaregridU(x1=x1,x2=x2),
                              Y, epsilon,
                              .scale=10,
                              paramd = compatible(do.call(expand.grid, list(
                                ell = list(1:2, 3, 1),
                                y = list(0:1, 0, 1),
                                .Beta = c(outer(c(1, 2, 5), 10 ^ {-2:3}, "*")),
                                .Gamma = c(outer(c(1, 2, 5), 10 ^ {-2:3}, "*"))))),
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
                    Yc <- updateY(X = .X, Y = .Y, ell = .ell, y = .y, .scale = .scale)
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
  Numden[[i]] <- parallelrho_server(Y = Y[[i]],
                                    epsilon = epsilon[[i]],
                                    .scale = 0.1,
                                    paramd = unique(compatible(rbind(do.call(expand.grid,list(
                                      experiment="gamma",
                                      ell=20201,
                                      y = c(0, outer(c(outer(c(-1,1), c(1,2,5),"*")), 10^(-1:2),"*")),
                                      .Beta=1,
                                      .Gamma = c(0, outer( c(1,2,5), 10^(-1:2),"*")))),
                                      do.call(expand.grid,list(
                                        experiment="beta",
                                        ell=20201,
                                        y = c(0, outer(c(outer(c(-1,1), c(1,2,5),"*")), 10^(-1:2),"*")),
                                        .Beta = c(0, outer(c(outer(c(-1,1), c(1,2,5),"*")), 10^(-1:2),"*")),
                                        .Gamma=1)),
                                      do.call(expand.grid,list(
                                        experiment="y",
                                        ell=20201,
                                        y = c(0, outer(c(outer(c(-1,1), c(1,2,5),"*")), 10^(-1:2),"*")),
                                        .Beta=1,
                                        .Gamma=1))))),
                                    ncores = 1)
  save(Numden, file = "CR_1.RData")
  print(paste(i, "end", Sys.time()))
}
system("echo \"CR_1.R done.\" | mail -s \"Process done\" \"pantalone.fra@gmail.com\"")



# load("/Users/francesco/Google Drive/InformativeSpatialSelection/CR_1.RData")
# rho_1 <- rhofrommc(rbind(Numden[[1]], Numden[[2]], Numden[[3]], Numden[[4]], Numden[[5]]))
# rho_1 <- na.omit(rho_1)
# rho_1$cimin <- sqrt(rho_1$vrho) * qnorm(0.025) + rho_1$rho
# rho_1$cimax <- sqrt(rho_1$vrho) * qnorm(0.975) + rho_1$rho
# rho_1 <- rho_1[is.finite(rowSums(rho_1)),]
# save(rho_1, file = "figure_rho1_data.rda")



# load("/Users/francesco/Google Drive/InformativeSpatialSelection/CR_1.RData")
# zzz_merge <- rhofrommc(rbind(Numden[[1]], Numden[[2]], Numden[[3]], Numden[[4]], Numden[[5]]))
# zzz_merge <- na.omit(zzz_merge)
# zzz_merge$cimin <- sqrt(zzz_merge$vrho) * qnorm(0.025) + zzz_merge$rho
# zzz_merge$cimax <- sqrt(zzz_merge$vrho) * qnorm(0.975) + zzz_merge$rho
# yplot <- c(-1, -0.5, -0.1, 0.1, 0.5, 1)
# dev.off()
# pdf(file = "/Users/francesco/Google Drive/InformativeSpatialSelection/rho_1_beta.pdf")
# zzz_merge1 <- zzz_merge[is.finite(rowSums(zzz_merge)),]
# ggplot(data = subset(zzz_merge1, .Gamma == 1 & y %in% yplot),
#        aes(x = .Beta, y=rho, ymin = cimin, ymax = cimax, fill = "red", alpha = .3)) +
#   geom_point() +
#   geom_line() +
#   scale_x_continuous(trans = 'log10') +
#   geom_ribbon() +
#   ggplot2::xlab(TeX("$\\xi_{1}$"))+
#   ggplot2::ylab(TeX("$\\rho_{1}$"))+
#   facet_wrap(~factor(y)) +
#   ggplot2::theme(legend.position = "none")
# dev.off()
# pdf(file = "/Users/francesco/Google Drive/InformativeSpatialSelection/rho_1_gamma.pdf")
# ggplot(data = subset(zzz_merge, .Beta == 1 & y %in% yplot ),
#        aes(x = .Gamma, y = rho, ymin = cimin, ymax = cimax, fill = "red", alpha = .3)) +
#   geom_point() +
#   geom_line() +
#   scale_x_continuous(trans = 'log10') +
#   geom_ribbon() +
#   ggplot2::xlab(TeX("$\\xi_{2}$"))+
#   ggplot2::ylab(TeX("$\\rho_{1}$"))+
#   facet_wrap(~factor(y)) +
#   ggplot2::theme(legend.position = "none")
# dev.off()
# pdf(file = "/Users/francesco/Google Drive/InformativeSpatialSelection/rho_1_y.pdf")
# ggplot(data = subset(zzz_merge, .Gamma == 1 & .Beta == 1 & y < 10),
#        aes(x = y, y = rho, ymin = cimin, ymax = cimax, fill = "red", alpha = .3)) +
#   geom_point() +
#   geom_line() +
#   scale_x_continuous(trans = 'log10') +
#   geom_ribbon()
# dev.off()






# load("/Users/francesco/Google Drive/InformativeSpatialSelection/CR_1.RData")
# zzz_merge <- rhofrommc(rbind(Numden[[1]], Numden[[2]], Numden[[3]], Numden[[4]], Numden[[5]]))
# zzz_merge <- na.omit(zzz_merge)
# zzz_merge$cimin <- sqrt(zzz_merge$vrho) * qnorm(0.025) + zzz_merge$rho
# zzz_merge$cimax <- sqrt(zzz_merge$vrho) * qnorm(0.975) + zzz_merge$rho
# yplot <- c(-1, -0.5, -0.1, 0.1, 0.5, 1)
# zzz_merge1 <- zzz_merge[is.finite(rowSums(zzz_merge)),]
# zzz_merge2 <- zzz_merge1
# zzz_merge2 <-cbind(zzz_merge, "y")
# colnames(zzz_merge2)[11] <- "scenario"
# zzz_merge2[which(zzz_merge2$y == -1 | zzz_merge2$y == 1), "scenario"] <- "-1,1"
# zzz_merge2[which(zzz_merge2$y == -0.5 | zzz_merge2$y == 0.5), "scenario"] <- "-0.5,0.5"
# zzz_merge2[which(zzz_merge2$y == -0.1 | zzz_merge2$y == 0.1), "scenario"] <- "-0.1,0.1"
# zzz_merge3 <- subset(zzz_merge2, y %in% yplot)
# zzz_merge3 <- cbind(zzz_merge3, "+")
# colnames(zzz_merge3)[12] <- "id"
# zzz_merge3[which(zzz_merge3$y %in% -c(0.1, 0.5, 1)), "id"] <- "-"
# 
# ggplot(data = subset(zzz_merge3, .Gamma == 1 & y %in% yplot),
#        aes(x = .Beta, y=rho, ymin = cimin, ymax = cimax, linetype = id, fill = "red", alpha = .3)) +
#   geom_point() +
#   geom_line() +
#   scale_x_continuous(trans = 'log10') +
#   geom_ribbon() +
#   ggplot2::xlab(TeX("$\\xi_{1}$"))+
#   ggplot2::ylab(TeX("$\\rho_{1}$"))+
#   facet_wrap(~factor(scenario)) +
#   ggplot2::theme(legend.position = "none")
# 
# 
# 
# ggplot(data = subset(zzz_merge3, .Beta == 1 & y %in% yplot ),
#        aes(x = .Gamma, y = rho, ymin = cimin, ymax = cimax, linetype = id, fill = "red", alpha = .3)) +
#   geom_point() +
#   geom_line() +
#   scale_x_continuous(trans = 'log10') +
#   geom_ribbon() +
#   ggplot2::xlab(TeX("$\\xi_{2}$"))+
#   ggplot2::ylab(TeX("$\\rho_{1}$"))+
#   facet_wrap(~factor(scenario)) +
#   ggplot2::theme(legend.position = "none")
