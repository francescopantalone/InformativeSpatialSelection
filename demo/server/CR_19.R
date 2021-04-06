## ------------------------------ ##
# rho_{12}                         #
# nrep = 10000                     #
# mu = 0, .sigma = 1, .scale = 0.1 #
# y values seq(-1, 1, by = 0.1)    #
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
Y <- foreach(i = 1:num_cores) %dopar%
  {
    set.seed(i + 42)
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
epsilon <- foreach(i = 1:num_cores) %dopar%
  {
    set.seed(i + 12)
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


# load("Y_scale1.RData")
# load("epsilon_scale1.RData")
# load("Numden.RData")
# load("Y_scale01.RData")
# load("epsilon_scale01.RData")

Numden <- list()
for(i in 1:num_cores)
{
  print(paste(i, "start", Sys.time()))
  Numden[[i]] <- parallelrho12_server(Y = Y[[i]],
                                      epsilon = epsilon[[i]],
                                      .scale = 0.1,
                                      paramd = compatible(do.call(expand.grid, list(
                                        ell = lapply(c(1,5,25,50),position_as_function_of_i),
                                        y = plyr::alply(as.matrix(expand.grid(seq(-1, 1, by = 0.1) ,
                                                                              seq(-1, 1, by = 0.1))),
                                                        1,identity),
                                        .Beta = 1,
                                        .Gamma = 1))),
                                      ncores = 1)
  save(Numden, file = "CR_19.RData")
  print(paste(i, "end", Sys.time()))
}
system("echo \"CR_19.R done.\" | mail -s \"Process done\" \"pantalone.fra@gmail.com\"")


# load("/Users/francesco/Google Drive/InformativeSpatialSelection/CR_19.RData")
# rho_12 <- rhofrommc(rbind(Numden[[1]], Numden[[2]]))
# rho_12 <- na.omit(rho_12)
# rho_12$cimin <- sqrt(rho_12$vrho) * qnorm(0.025) + rho_12$rho
# rho_12$cimax <- sqrt(rho_12$vrho) * qnorm(0.975) + rho_12$rho
# rho_12 <- cbind(rho_12, t(as.data.frame(rho_12$y)), t(as.data.frame(rho_12$ell)))
# colnames(rho_12)[11:14] <- c("y1", "y2", "ell1", "ell2")
# dist_fun <- function(X, ell1, ell2){
#   return(sqrt(sum((X[ell1,] - X[ell2,]) ^ 2)))
# }
# x1 <- generatex_i(201)
# x2 <- x1
# X <- generatesquaregridU(x1 = x1, x2 = x1)
# d <- rep(0, nrow(rho_12))
# for(i in 1:nrow(rho_12))
# {
#   d[i] <- dist_fun(X, rho_12$ell1[i], rho_12$ell2[i])
# }
# rho_12 <- cbind(rho_12, d)
# data("figure1_data")
# scenario1 <- cbind(figure1_data[figure1_data$.n==2,], 20001, 20602, "Scenario 1")
# colnames(scenario1)[7:9] <- c("ell1", "ell2", "scenario")
# scenario2 <- cbind(figure1_data[figure1_data$.n==2,], 19201, 21402, "Scenario 2")
# colnames(scenario2)[7:9] <- c("ell1", "ell2", "scenario")
# scenario3 <- cbind(figure1_data[figure1_data$.n==2,], 15201, 25402, "Scenario 3")
# colnames(scenario3)[7:9] <- c("ell1", "ell2", "scenario")
# scenario4 <- cbind(figure1_data[figure1_data$.n==2,], 10201, 30402, "Scenario 4")
# colnames(scenario4)[7:9] <- c("ell1", "ell2", "scenario")
# scenario <- rbind(scenario1, scenario2, scenario3, scenario4)
# save(rho_12, scenario, file = "figure_rho_12_data.rda")




# load("/Users/francesco/Google Drive/InformativeSpatialSelection/CR_19.RData")
# zzz_merge <- rhofrommc(rbind(Numden[[1]], Numden[[2]]))
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
# data("figure1_data")
# scenario1 <- cbind(figure1_data[figure1_data$.n==2,], 20001, 20602, "Scenario 1")
# colnames(scenario1)[7:9] <- c("ell1", "ell2", "scenario")
# scenario2 <- cbind(figure1_data[figure1_data$.n==2,], 19201, 21402, "Scenario 2")
# colnames(scenario2)[7:9] <- c("ell1", "ell2", "scenario")
# scenario3 <- cbind(figure1_data[figure1_data$.n==2,], 15201, 25402, "Scenario 3")
# colnames(scenario3)[7:9] <- c("ell1", "ell2", "scenario")
# scenario4 <- cbind(figure1_data[figure1_data$.n==2,], 10201, 30402, "Scenario 4")
# colnames(scenario4)[7:9] <- c("ell1", "ell2", "scenario")
# scenario <- rbind(scenario1, scenario2, scenario3, scenario4)
# # library(ggplot2)
# # library(latex2exp)
# grid_plot <- ggplot2::ggplot(scenario, ggplot2::aes(x1,x2)) + 
#   ggplot2::geom_tile(ggplot2::aes(fill = y1),size=2) + 
#   ggplot2::geom_point(ggplot2::aes(x = x1[ell1], y = x2[ell1]),
#                       fill = "white", color = "black", size = 1.5, stroke = .9, shape = 25) +
#   ggplot2::geom_point(ggplot2::aes(x = x1[ell2], y = x2[ell2]),
#                       fill = "white", color = "black", size = 1.5, stroke = .9, shape = 25) +
#   ggplot2::xlab("")+
#   ggplot2::ylab("")+
#   ggplot2::scale_x_continuous(expand=c(0,0))+
#   ggplot2::scale_y_continuous(expand=c(0,0))+
#   ggplot2::theme(
#     axis.ticks.x=ggplot2::element_blank(),
#     axis.ticks.y=ggplot2::element_blank(),
#     strip.text.x = element_text(angle = 0, hjust = 0),
#     panel.border = ggplot2::element_blank(),
#     panel.grid.major = ggplot2::element_blank(), 
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank())+
#   ggplot2::scale_fill_gradient(low = "white",     
#                                high = "black"#,guide=guide_colorbar(frame.colour=c("black"),frame.linewidth=2)
#   )  +
#   ggplot2::facet_wrap(~scenario, nrow = 1)+ 
#   ggplot2::labs(fill = "signal") + ggplot2::theme(legend.position="bottom")
# 
# contour_plot <- ggplot(data = zzz_merge2, aes(x = y1, y = y2, z = rho)) + geom_contour_filled() +
#   ggplot2::xlab(TeX("$y_{1}")) +
#   ggplot2::ylab(TeX("$y_{2}")) +
#   ggplot2::scale_x_continuous(expand = c(0, 0)) +
#   ggplot2::scale_y_continuous(expand = c(0, 0)) +
#   ggplot2::theme(
#     axis.ticks.x = ggplot2::element_blank(),
#     axis.ticks.y = ggplot2::element_blank(),
#     panel.border = ggplot2::element_blank(),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank()) +
#   facet_wrap(~round(d, 3), nrow = 1) +
#   ggplot2::labs(fill = TeX("$\\rho_{1,2}$")) + ggplot2::theme(legend.position="bottom")
# 
# dev.off()
# # pdf(file = "/Users/francesco/Google Drive/InformativeSpatialSelection/rho_12_contour.pdf")
# png(file = "/Users/francesco/Google Drive/InformativeSpatialSelection/rho_12_contour.png")
# grid.arrange(grid_plot, contour_plot, nrow = 2)
# dev.off()









# library(ggplot2)
# library(latex2exp)
# dev.off()
# pdf(file = "/Users/francesco/Google Drive/InformativeSpatialSelection/rho_12_contour.pdf")
# ggplot(data = zzz_merge2, aes(x = y1, y = y2, z = rho)) + geom_contour_filled() +
#   ggplot2::xlab(TeX("$y_{1}")) +
#   ggplot2::ylab(TeX("$y_{2}")) +
#   ggplot2::scale_x_continuous(expand = c(0, 0)) +
#   ggplot2::scale_y_continuous(expand = c(0, 0)) +
#   ggplot2::theme(
#     axis.ticks.x = ggplot2::element_blank(),
#     axis.ticks.y = ggplot2::element_blank(),
#     panel.border = ggplot2::element_blank(),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank()) +
#   facet_wrap(~round(d, 3)) +
#   ggplot2::labs(fill = TeX("$\\rho_{1,2}$"))
# dev.off()
# 
# pdf(file = "/Users/francesco/Google Drive/InformativeSpatialSelection/rho_12_heatmap.pdf")
# ggplot2::ggplot(zzz_merge2, ggplot2::aes(y1, y2)) +
#   ggplot2::geom_tile(ggplot2::aes(fill = rho),size=2) +
#   ggplot2::xlab(TeX("$y_{1}")) +
#   ggplot2::ylab(TeX("$y_{2}")) +
#   ggplot2::scale_x_continuous(expand = c(0, 0)) +
#   ggplot2::scale_y_continuous(expand = c(0, 0)) +
#   ggplot2::theme(
#     axis.ticks.x = ggplot2::element_blank(),
#     axis.ticks.y = ggplot2::element_blank(),
#     panel.border = ggplot2::element_blank(),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank()) +
#   facet_wrap(~round(d, 3)) +
#   ggplot2::labs(fill = TeX("$\\rho_{1,2}$"))
# dev.off()



