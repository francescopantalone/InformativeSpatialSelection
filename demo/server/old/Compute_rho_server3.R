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
num_cores <- 10
nrep <- 250
# # x1 <- generatex_i(201)
# # x2 <- generatex_i(201)
# # X <- generateYonsquaregridU(x1 = x1, x2 = x2)
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
                                  .scale = .1,
                                  mu = 0,
                                  .sigma = 1))
  }
epsilon <- list()
registerDoParallel(cores = num_cores)
set.seed(12)
epsilon <- foreach(i = 1:num_cores) %dopar%
  {
    return(generateYonsquaregridU(x1 = x1,
                                  x2 = x2,
                                  n = nrep,
                                  .scale = .1,
                                  mu = 0,
                                  .sigma = 1))
  }
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


# load("Y.RData")
# load("epsilon.RData")
# load("Numden.RData")
# data(figure1_data)
# y <- quantile(figure1_data$y, c(.75, .95, .99))
yexp1 <- c(outer(c(1, 2, 5), 10 ^ {-3:3}, "*"))
yexp1 <- c(-yexp1, yexp1)
yexp2 <- c(0.001, 0.1, 1, 2, 10, 100)
yexp2 <- c(-yexp2, yexp2)
Numden3 <- list()
for(i in 1:num_cores)
{
  print(paste(i, "start", Sys.time()))
  Numden3[[i]] <- parallelrho_server(Y = Y[[i]],
                                     epsilon = epsilon[[i]],
                                     paramd = rbind(compatible(do.call(expand.grid,list(
                                       experiment="gamma",
                                       ell=20301,
                                       y=yexp2,
                                       .Beta=1,
                                       .Gamma=c(0.01, 0.1, 0.5, 1, 10, 100)))),
                                       compatible(do.call(expand.grid,list(
                                         experiment="beta",
                                         ell=20301,
                                         y=yexp2,
                                         .Beta=c(0.01, 0.1, 0.5, 1, 10, 100),
                                         .Gamma=1))),
                                       compatible(do.call(expand.grid,list(
                                         experiment="y",
                                         ell=20301,
                                         y=yexp1,
                                         .Beta=1,
                                         .Gamma=1)))),
                                     ncores = 1)
  save(Numden3, file = "Numden3.RData")
  print(paste(i, "end", Sys.time()))
}
system("echo \"Compute_rho_server.R done.\" | mail -s \"Process done\" \"pantalone.fra@gmail.com\"")



# -- plot -- #

# zzz_merge <- rhofrommc(rbind(Numden[[1]], Numden[[2]], Numden[[3]], Numden[[4]], Numden[[5]],
#                              Numden[[6]], Numden[[7]], Numden[[8]], Numden[[9]], Numden[[10]]))
# zzz_merge <- na.omit(zzz_merge)
# zzz_merge$cimin <- sqrt(zzz_merge$vrho) * qnorm(0.025) + zzz_merge$rho
# zzz_merge$cimax <- sqrt(zzz_merge$vrho) * qnorm(0.975) + zzz_merge$rho

zzz_merge <- rhofrommc(rbind(Numden3[[1]], Numden3[[2]], Numden3[[3]], Numden3[[4]], Numden3[[5]],
                             Numden3[[6]], Numden3[[7]], Numden3[[8]], Numden3[[9]], Numden3[[10]]))
# Nwna <- nrow(zzz_merge)
zzz_merge <- na.omit(zzz_merge)
# Nwona <- nrow(zzz_merge)
# Nwona / Nwna
zzz_merge$cimin <- sqrt(zzz_merge$vrho) * qnorm(0.025) + zzz_merge$rho
zzz_merge$cimax <- sqrt(zzz_merge$vrho) * qnorm(0.975) + zzz_merge$rho

yexp1 <- c(outer(c(1, 2, 5), 10 ^ {-3:3}, "*"))
yexp1 <- c(-yexp1, yexp1)
yexp2 <- c(0.001, 0.1, 1, 2, 10, 100)
yexp2 <- c(-yexp2, yexp2)
#
figureXXgamma <- ggplot(data=subset(zzz_merge,.Beta==1 & is.element(y, yexp2)),
                        aes(x=.Gamma,
                            y=rho,
                            ymin=cimin,
                            ymax=cimax,fill="red",alpha=.3))+
  geom_point()+
  geom_line()+
  scale_x_continuous(trans='log10')+
  geom_ribbon()+facet_wrap(~factor(y))


a <- subset(zzz_merge,.Beta==1)
plot3D::scatter3D(x = a$.Gamma, y = a$y, z = a$rho)


#
#
# # Check if something is off
# # Consider bigger values of y
figureXXbeta <-
  ggplot(data=subset(zzz_merge,
                     .Gamma==1 & is.element(y, yexp2)),
         aes(x=.Beta,
             y=rho,
             ymin=cimin,ymax=cimax,fill="red",alpha=.3))+
  geom_point()+
  geom_line()+
  scale_x_continuous(trans='log10')+
  geom_ribbon()+facet_wrap(~factor(y))


# scatterplot3d::scatterplot3d(x = zzz_merge$.Beta, y = zzz_merge$y, z = zzz_merge$rho)
a <- subset(zzz_merge, .Gamma==1)
plot3D::scatter3D(x = a$.Beta, y = a$y, z = a$rho)


## Vary values of beta and gamma. We want to see how the relationship varies at different values of beta and gamma
figureXXy<-
  ggplot(data=subset(zzz_merge,
                     .Gamma==1&.Beta==1),
         aes(x=y,
             y=rho,
             ymin=cimin,ymax=cimax,fill="red",alpha=.3))+
  geom_point()+
  geom_line()+
  scale_x_continuous(trans='log10')+
  scale_y_continuous(trans='log10')+
  geom_ribbon()

#

