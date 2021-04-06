function(){
# ------ libraries ------ #
library(RandomFields)
library(sampling)
library(gstat)
library(ggplot2)
library(gridExtra)
# ----------------------- #
# -- generate data -- #
# grid_size <- 50
# set.seed(42)
# Y <- generateYonsquaregridU(length.out = grid_size, mu = 10, 
#                             .sigma = 1)
# # Z <- generateZconditionnallyongridUY(length.out = grid_size, y = Y, 
# #                                      alpha = 20, 
# #                                      beta = c(0, 0, 0.8), 
# #                                      .sigma = 5)
# Z <- generateZconditionnallyongridUY2(length.out = grid_size, y = Y, 
#                                       alpha = 1, .sigma.epsilon = 1,
#                                       beta = c(0, 0, 0.6))
# XYZ <- bindXYZpiI(x1 = generatex_i(grid_size), x2 = generatex_i(grid_size), y = Y, z = Z)
# plotY(XYZ)
# ------------------- #
# -- #
load("figure4_data.RData")
N <- nrow(XYZ)
f <- 0.1
n <- round(f * N)
s_srs <- srs(n = n, XYZ = XYZ, replace = FALSE)
s_pps <- pps(n = n, Z = XYZ$z.z3, XYZ = XYZ)
XYZ_so <- XYZ
coordinates(XYZ_so) <- ~ x1 + x2 # created a spatial objects for gstat package
krig_s1 <- list()
length(krig_s1) <- M
krig_s2 <- list()
length(krig_s1) <- M
krig_s3 <- list()
length(krig_s1) <- M
for(i in 1:M)
{
  krig_s1[[i]] <- krige(y ~ 1, XYZ_so[s1_M[i, ],], XYZ_so, model = fit_vgm_s1[[i]])
  krig_s2[[i]] <- krige(y ~ 1, XYZ_so[s2_M[i, ],], XYZ_so, model = fit_vgm_s2[[i]])
  krig_s3[[i]] <- krige(y ~ 1, XYZ_so[s3_M[i, ],], XYZ_so, model = fit_vgm_s3[[i]])
}
krig_s1_bias <- rep(0, M)
krig_s2_bias <- rep(0, M)
krig_s3_bias <- rep(0, M)
for(i in 1:M)
{
  krig_s1_bias[i] <- (sum(krig_s1[[i]]$var1.pred) - sum(XYZ_so$y)) / sum(XYZ_so$y)
  krig_s2_bias[i] <- (sum(krig_s2[[i]]$var1.pred) - sum(XYZ_so$y)) / sum(XYZ_so$y)
  krig_s3_bias[i] <- (sum(krig_s3[[i]]$var1.pred) - sum(XYZ_so$y)) / sum(XYZ_so$y)
}
mean(krig_s1_bias)
mean(krig_s2_bias)
mean(krig_s3_bias)
krig_s1_mse <- rep(0, M)
krig_s2_mse <- rep(0, M)
krig_s3_mse <- rep(0, M)
for(i in 1:M)
{
  krig_s1_mse[i] <- sum((krig_s1[[i]]$var1.pred - XYZ_so$y) ^ 2) / N
  krig_s2_mse[i] <- sum((krig_s2[[i]]$var1.pred - XYZ_so$y) ^ 2) / N
  krig_s3_mse[i] <- sum((krig_s3[[i]]$var1.pred - XYZ_so$y) ^ 2) / N
}
mean(krig_s1_mse)
mean(krig_s2_mse)
mean(krig_s3_mse)
# emp_var_srs <- variogram(y ~ 1, XYZ_so[s_srs, ])
# fit_var_srs <- fit.variogram(emp_var, vgm(1, "Gau", 1))
# krig_srs <- krige(y ~ 1, XYZ_so[s_srs,], XYZ_so, model = fit_var)
# krig_srs %>% as.data.frame %>% ggplot(aes(x = x1, y = x2)) +
#   geom_tile(aes(fill=var1.pred)) + coord_equal() + 
#   scale_fill_gradient(low = "yellow", high="red") + theme_bw()
# 
# emp_var_pps <- variogram(y ~ 1, XYZ_so[s_pps, ])
# fit_var_pps <- fit.variogram(emp_var, vgm(1, "Gau", 1))
# krig_pps <- krige(y ~ 1, XYZ_so[s_pps,], XYZ_so, model = fit_var)
# krig_pps %>% as.data.frame %>% ggplot(aes(x = x1, y = x2)) +
#   geom_tile(aes(fill=var1.pred)) + coord_equal() + 
#   scale_fill_gradient(low = "yellow", high="red") + theme_bw()
}