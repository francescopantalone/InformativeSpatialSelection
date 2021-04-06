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
                                  .scale = 1,
                                  mu = 10,
                                  .sigma = 2))
  }
save(Y, file = "Y_mu10scale1.RData")
epsilon <- list()
registerDoParallel(cores = num_cores)
set.seed(12)
epsilon <- foreach(i = 1:num_cores) %dopar%
  {
    return(generateYonsquaregridU(x1 = x1,
                                  x2 = x2,
                                  n = nrep,
                                  .scale = 1,
                                  mu = 10,
                                  .sigma = 2))
  }
save(epsilon, file = "epsilon_mu10scale1.RData")
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


# load("Y_scale1.RData")
# load("epsilon_scale1.RData")
# load("Numden.RData")
data(figure1_data)
y <- quantile(figure1_data$y, c(.75, .95, .99))
yexp1 <- c(outer(c(1, 2, 5), 10 ^ {-3:3}, "*"))
yexp1 <- c(-yexp1, yexp1)
yexp2 <- c( 0.1, 0.5, 1, 2, 10, 100)
yexp2 <- c(-yexp2, yexp2)
Numden <- list()
for(i in 1:num_cores)
{
  print(paste(i, "start", Sys.time()))
  Numden[[i]] <- parallelrho_server(Y = Y[[i]],
                                    epsilon = epsilon[[i]],
                                    .scale = 1,
                                    paramd = unique(compatible(rbind(do.call(expand.grid,list(
                                      experiment="gamma",
                                      ell=20301,
                                      y=yexp2,
                                      .Beta=1,
                                      .Gamma=c(outer(c(1,2,5),10^{-3:3},"*")))),
                                      do.call(expand.grid,list(
                                        experiment="beta",
                                        ell=20301,
                                        y=yexp2,
                                        .Beta=c(outer(c(1,2,5),10^{-3:3},"*")),
                                        .Gamma=1)),
                                      do.call(expand.grid,list(
                                        experiment="y",
                                        ell=20301,
                                        y=yexp1,
                                        .Beta=1,
                                        .Gamma=1))))),
                                    ncores = 1)
  save(Numden, file = "Numden_mu10scale1.RData")
  print(paste(i, "end", Sys.time()))
}
system("echo \"Compute_rho_server_mu10scale1.R done.\" | mail -s \"Process done\" \"pantalone.fra@gmail.com\"")


load("/Users/francesco/Google Drive/InformativeSpatialSelection/Numden_mu10scale1.RData")
zzz_merge <- rhofrommc(Numden[[1]])
# # Nwna <- nrow(zzz_merge)
zzz_merge <- na.omit(zzz_merge)
# # Nwona <- nrow(zzz_merge)
# # Nwona / Nwna
zzz_merge$cimin <- sqrt(zzz_merge$vrho) * qnorm(0.025) + zzz_merge$rho
zzz_merge$cimax <- sqrt(zzz_merge$vrho) * qnorm(0.975) + zzz_merge$rho
yexp1 <- c(outer(c(1, 2, 5), 10 ^ {-3:3}, "*"))
yexp1 <- c(-yexp1, yexp1)
yexp2 <- c(0.1, 0.5, 1, 2, 10, 100)
yexp2 <- c(-yexp2, yexp2)


ggplot(data = subset(zzz_merge, .Beta == 1 & .Gamma < 10 & is.element(y, yexp2)),
       aes(x = .Gamma, y = rho, ymin = cimin, ymax = cimax, fill = "red", alpha = .3)) + 
  geom_point() + geom_line() + geom_ribbon() +
  geom_hline(yintercept = 1, col = "red", linetype = "dashed") + facet_wrap(~factor(y))


ggplot(data = subset(zzz_merge, .Gamma == 1 & .Beta < 10 & is.element(y, yexp2)),
       aes(x = .Beta, y=rho, ymin = cimin, ymax = cimax, fill = "red", alpha = .3)) +
  geom_point() + geom_line() + geom_ribbon() + 
  geom_hline(yintercept = 1, col = "red", linetype = "dashed") + facet_wrap(~factor(y))


ggplot(data = subset(zzz_merge, .Gamma == 1 & .Beta == 1),
       aes(x = y, y = rho, ymin = cimin, ymax = cimax, fill = "red", alpha = .3)) +
  geom_point() + geom_line() + geom_ribbon()

