function(){
# "#U (x1,x2)
# # Y
# 
# Y: population of size N. (200*200) minimal size s.t 
# the sample looks good N*sampling rate large enough and only the high values are selected.
# 
# N2: minimal size such that the plot looks good. N2<N
# 
# 
# dataframe(x1,x2,y,z3,I1...I1000,ypred1,...,ypred1000)<-size N2
# 
# 
# 
# dataframe(x1,x2,y,ystar,z3,I1...I3,ypred1,...,ypred3) <- size N
# ystar~ rho f
# 
# 
# 
# #first graph show that the process we predict looks like a rho f process
# ggplot :
# subgraph 1: x1,x2,y
# subgraph 2: x1,x2,ypred1 <- naive 1
# subgraph 3: x1,x2,ypred2 <- naive 2
# subgraph 5: x1,x2, ystar <- independent on y
#  
# #second graph: show correction
#  ggplot 
#  subgraph 1: x1,x2,y
#  subgraph 5: x1,x2, ypred1
#  subgraph 6: x1,x2, ypred1 corrected <- do not have yet.
#  
# ggplot 
# "

# ------ libraries ------ #
library(InformativeSpatialSelection)
library(RandomFields)
library(sampling)
library(gstat)
library(ggplot2)
library(gridExtra)
library(grid)
# ----------------------- #





# -- work in progress -- #
xx <- 200
x1 <- x2 <- generatex_i(xx)#This generates an arithmetic sequence  of length xx starting at 0 ending at 1.
set.seed(42)
y <- generateYonsquaregridU(length.out = xx, x1 = x1, x2 = x1, mu = 0, .sigma = 1, .scale = 10)
b <- 2
cc <- .3
z <- generateZconditionnallyongridUY2(length.out = xx, 
                                      y = y, 
                                      alpha = c(log(100), log(100) - b ^ 2 / 2 - cc ^ 2 / 2,
                                                log(100) - b ^ 2 / 2 - cc ^ 2 / 2),
                                      beta  = c(0, 0, b), 
                                      gamma = c(0, sqrt(b ^ 2 + cc ^ 2), cc),
                                      .scale.epsilon = 10)
XYZ <- bindXYZpiI(x1 = generatex_i(xx), x2 = generatex_i(xx), y = y, z = z)
# ------------------- #
# -- just one shot -- #
N <- nrow(XYZ)
n <- 100
set.seed(12)
s1 <- srs(XYZ = XYZ, n = n, replace = FALSE)
s2_x1x2 <- spatstat::rpoispp(lambda = spatstat::as.im.data.frame(XYZ[c("x1", "x2", "z.z2")]), nsim = 1)
s3_x1x2 <- spatstat::rpoispp(lambda = spatstat::as.im.data.frame(XYZ[c("x1", "x2", "z.z3")]), nsim = 1)
s2_x1x2 <- data.frame(x1 = s2_x1x2$x, x2 = s2_x1x2$y) 
s3_x1x2 <- data.frame(x1 = s3_x1x2$x, x2 = s3_x1x2$y)

# s1 <- 0
# for(i in 1:nrow(s1_x1x2))
# {
#   append(s1, which(round(XYZ$x1, 2) == round(s1_x1x2$x1[i], 2) & round(XYZ$x2, 2) == round(s1_x1x2$x2[i], 2)))
#   # append() is used to face the possibility to have the same point sampled more than once
#   # round() is used, otherwise points don't match
# }
# s1 <- s1[-1]

s2 <- rep(0, nrow(s2_x1x2))
for(i in 1:nrow(s2_x1x2))
{
  s2[i] <- which(round(XYZ$x1, 2) == round(s2_x1x2$x1[i], 2) & round(XYZ$x2, 2) == round(s2_x1x2$x2[i], 2))[1]
  # round() is used, otherwise points don't match
}
s2 <- unique(s2)

s3 <- rep(0, nrow(s3_x1x2))
for(i in 1:nrow(s3_x1x2))
{
  s3[i] <- which(round(XYZ$x1, 2) == round(s3_x1x2$x1[i], 2) & round(XYZ$x2, 2) == round(s3_x1x2$x2[i], 2))[1]
  # round() is used, otherwise points don't match
}
s3 <- unique(s3)

XYZ_so <- XYZ
coordinates(XYZ_so) <- ~ x1 + x2 # created a spatial objects for gstat package
# vgm_pop <- variogram(y ~ 1, XYZ_so) # population variogram, since we use all the field
vgm_s1 <- variogram(y ~ 1, XYZ_so[s1, ])
vgm_s2 <- variogram(y ~ 1, XYZ_so[s2, ])
vgm_s3 <- variogram(y ~ 1, XYZ_so[s3, ])
# -- fit variogram -- #
fit_vgm_s1 <- fit.variogram(vgm_s1, vgm("Gau"))
fit_vgm_s2 <- fit.variogram(vgm_s2, vgm("Gau"))
fit_vgm_s3 <- fit.variogram(vgm_s3, vgm("Gau"))
krig_s1 <- krige(y ~ 1, XYZ_so[s1, ], XYZ_so, model = fit_vgm_s1)
krig_s2 <- krige(y ~ 1, XYZ_so[s2, ], XYZ_so, model = fit_vgm_s2)
krig_s3 <- krige(y ~ 1, XYZ_so[s3, ], XYZ_so, model = fit_vgm_s3)
I_s1 <- rep(0, nrow(XYZ))
I_s1[s1] <- 1
I_s2 <- rep(0, nrow(XYZ))
I_s2[s2] <- 1
I_s3 <- rep(0, nrow(XYZ))
I_s3[s3] <- 1
# I_s3 <- rep(0, nrow(XYZ))
# I_s3[s3] <- 1
samples <- data.frame(c("y", rep("Scenario 1", length(s1)), rep("Scenario 2", length(s2)), rep("Scenario 3", length(s3))),
                      c(0, XYZ$x1[s1], XYZ$x1[s2], XYZ$x1[s3]),
                      c(0, XYZ$x2[s1], XYZ$x2[s2], XYZ$x2[s3]),
                      c(0, krig_s1$var1.pred[s1], krig_s2$var1.pred[s2], krig_s3$var1.pred[s3]))
colnames(samples) <- c("Scenario", "x1", "x2", "y_pred")
samples$Scenario <- factor(x = samples$Scenario, levels = c("y", "Scenario 1", "Scenario 2", "Scenario 3"))
# pred_data <- data.frame(c(rep("Scenario 1", nrow(XYZ)), rep("Scenario 2", nrow(XYZ))),
#                         c(XYZ$x1, XYZ$x1),
#                         c(XYZ$x2, XYZ$x2),
#                         c(XYZ$y, XYZ$y),
#                         c(XYZ$z.z1, XYZ$z.z2),
#                         c(I_s1, I_s2),
#                         c(krig_s1$var1.pred, krig_s2$var1.pred),
#                         c(abs(krig_s1$var1.var), abs(krig_s2$var1.var)))
# colnames(pred_data) <- c("Scenario", "x1", "x2", "y", "z", "s", "y_pred","var")

pred_data <- data.frame(c(rep("y", nrow(XYZ)), rep("Scenario 1", nrow(XYZ)), rep("Scenario 2", nrow(XYZ)), rep("Scenario 3", nrow(XYZ))),
                        c(XYZ$x1, XYZ$x1, XYZ$x1, XYZ$x1),
                        c(XYZ$x2, XYZ$x2, XYZ$x2, XYZ$x2),
                        c(XYZ$y, XYZ$y, XYZ$y, XYZ$y),
                        c(rep(0, nrow(XYZ)), rep(0, nrow(XYZ)), XYZ$z.z1, XYZ$z.z2),
                        c(rep(1, nrow(XYZ)), I_s1, I_s2, I_s3),
                        c(XYZ$y, krig_s1$var1.pred, krig_s2$var1.pred, krig_s3$var1.pred),
                        c(rep(0, nrow(XYZ)), abs(krig_s1$var1.var), abs(krig_s2$var1.var), abs(krig_s3$var1.var)))
colnames(pred_data) <- c("Scenario", "x1", "x2", "y", "z", "s", "y_pred","var")
pred_data$Scenario <- factor(x = pred_data$Scenario, levels = c("y", "Scenario 1", "Scenario 2", "Scenario 3"))

pred_data$error <- pred_data$y - pred_data$y_pred
pred_data$pvalue <- 2 * (pnorm(abs(pred_data$error),
                               mean = 0, sd = sqrt(pred_data$var))) - 1
# pred_data$pvalue2 <- .98 * (pred_data$pvalue <= .98) + pred_data$pvalue * (pred_data$pvalue > .98)
pred_data$pvalue2 <- .95 * (pred_data$pvalue <= .95) + pred_data$pvalue * (pred_data$pvalue > .95)



plot_prediction <- ggplot(data = pred_data, aes(x = x1, y = x2, fill = y_pred)) + geom_tile() +
  ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::scale_x_continuous(expand=c(0,0)) +
  ggplot2::scale_y_continuous(expand=c(0,0)) + 
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()) +
  ggplot2::scale_fill_gradient(low = "white", high = "black") +
  geom_point(data = samples, shape = 6, color = "black", fill = "white", size = 2, stroke = 1.2) +
  facet_wrap(~Scenario, nrow = 1) +
  ggplot2::labs(fill = "signal") + 
  ggplot2::theme(legend.position="bottom",
                 legend.justification="right",
                 legend.margin=ggplot2::margin(0,0,0,0),
                 legend.box.margin=ggplot2::margin(-20,#top
                                                   18,#right
                                                   10,#bottom
                                                   10))#left
print(plot_prediction)
plot_error <- ggplot(pred_data, aes(x = x1, y = x2, fill = (y - y_pred))) + geom_tile() +
  ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::scale_x_continuous(expand=c(0,0)) +
  ggplot2::scale_y_continuous(expand=c(0,0)) + 
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()) +
  ggplot2::scale_fill_gradient(low = "white", high = "black") +
  geom_point(data = samples, shape = 6, color = "black", fill = "white", size = 2, stroke = 1.2) + 
  facet_wrap(~Scenario, nrow = 1) +
  ggplot2::labs(fill = "error") +
  ggplot2::theme(legend.position="bottom",
                 legend.justification="right",
                 legend.margin=ggplot2::margin(0,0,0,0),
                 legend.box.margin=ggplot2::margin(-20,#top
                                                   18,#right
                                                   10,#bottom
                                                   10))#left
print(plot_error)
# samplepos<-
# pred_data$disttosample<-apply(pred_data[c("x1","x2")],1,
#                               function(x){min(dist(x,))})


# pred_data$pvalue2<-pmax(pred_data$pvalue,.5)
plot_q <- ggplot(pred_data, aes(x = x1, y = x2, fill = pvalue2)) + geom_tile() +
  ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::scale_x_continuous(expand=c(0,0)) +
  ggplot2::scale_y_continuous(expand=c(0,0)) + 
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()) +
  ggplot2::scale_fill_gradient2(low = "white", mid="white",midpoint=.95,high = "black") +
  geom_point(data = samples, shape = 6, color = "black", fill = "white", size = 2, stroke = 1.2) + 
  facet_wrap(~Scenario, nrow = 1) +
  ggplot2::labs(fill = "standardised error") +
  ggplot2::theme(legend.position="bottom",
                 legend.justification="right",
                 legend.margin=ggplot2::margin(0,0,0,0),
                 legend.box.margin=ggplot2::margin(-20,#top
                                                   18,#right
                                                   10,#bottom
                                                   10))#left

print(plot_q)

plot_var <- ggplot(pred_data, aes(x = x1, y = x2, fill = var)) + geom_tile() +
  ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::scale_x_continuous(expand=c(0,0)) +
  ggplot2::scale_y_continuous(expand=c(0,0)) + 
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()) +
  ggplot2::scale_fill_gradient(low = "white", high = "black") +
  geom_point(data = samples, shape = 6, color = "black", fill = "white", size = 2, stroke = 1.2) + 
  facet_wrap(~Scenario, nrow = 1) +
  ggplot2::labs(fill = "estimated prediction variance") +
  ggplot2::theme(legend.position="bottom",
                 legend.justification="right",
                 legend.margin=ggplot2::margin(0,0,0,0),
                 legend.box.margin=ggplot2::margin(-20,#top
                                                   18,#right
                                                   10,#bottom
                                                   10))#left

print(plot_var)
# ---------------------- #

# -- old version -- #
# -- generate data -- #
grid_size <- 200
set.seed(42)
Y <- generateYonsquaregridU(length.out = grid_size, mu = 10,
                            .sigma = 1)
# Z <- generateZconditionnallyongridUY(length.out = grid_size, y = Y,
#                                      alpha = 20,
#                                      beta = c(0, 0, 0.8),
#                                      .sigma = 5)
Z <- generateZconditionnallyongridUY2(length.out = grid_size, y = Y,
                                      alpha = 1, .scale.epsilon = 1,
                                      beta = c(0, 0, 0.6))
XYZ <- bindXYZpiI(x1 = generatex_i(grid_size), x2 = generatex_i(grid_size), y = Y, z = Z)
# plotY(XYZ)
# ------------------- #
# -- just one shot -- #
N <- nrow(XYZ)
n <- 100 # sample size
set.seed(42)
s1 <- sample(x = nrow(XYZ), size = n)
pi2 <- inclusionprobabilities(a = XYZ$z.z2, n = n)
s2 <- sample(x = nrow(XYZ), size = n, prob = pi2)
pi3 <- inclusionprobabilities(a = XYZ$z.z3, n = n)
s3 <- sample(x = nrow(XYZ), size = n, prob = pi3)
XYZ_so <- XYZ
coordinates(XYZ_so) <- ~ x1 + x2 # created a spatial objects for gstat package
# vgm_pop <- variogram(y ~ 1, XYZ_so) # population variogram, since we use all the field
vgm_s1 <- variogram(y ~ 1, XYZ_so[s1, ])
vgm_s2 <- variogram(y ~ 1, XYZ_so[s2, ])
vgm_s3 <- variogram(y ~ 1, XYZ_so[s3, ])
# -- fit variogram -- #
# fit_vgm_s1 <- fit.variogram(vgm_s1, vgm(1, "Gau", 1, nugget = 0.1))
# fit_vgm_s2 <- fit.variogram(vgm_s2, vgm(1, "Gau", 1, nugget = 0.1))
# fit_vgm_s3 <- fit.variogram(vgm_s3, vgm(1, "Gau", 1, nugget = 0.1))
fit_vgm_s1 <- fit.variogram(vgm_s1, vgm("Gau"))
fit_vgm_s2 <- fit.variogram(vgm_s2, vgm("Gau"))
fit_vgm_s3 <- fit.variogram(vgm_s3, vgm("Gau"))
# fit_vgm_s1 <- fit.variogram(vgm_s1, vgm(1, "Gau", 1))
# fit_vgm_s2 <- fit.variogram(vgm_s2, vgm(1, "Gau", 1))
# fit_vgm_s3 <- fit.variogram(vgm_s3, vgm(1, "Gau", 1))
# fit_vgm_s1 <- fit.variogram(vgm_s1, vgm("Gau"))
# fit_vgm_s2 <- fit.variogram(vgm_s2, vgm("Gau"))
# fit_vgm_s3 <- fit.variogram(vgm_s3, vgm("Gau"))
krig_s1 <- krige(y ~ 1, XYZ_so[s1, ], XYZ_so, model = fit_vgm_s1)
krig_s2 <- krige(y ~ 1, XYZ_so[s2, ], XYZ_so, model = fit_vgm_s2)
krig_s3 <- krige(y ~ 1, XYZ_so[s3, ], XYZ_so, model = fit_vgm_s3)
I_s1 <- rep(0, nrow(XYZ))
I_s1[s1] <- 1
I_s2 <- rep(0, nrow(XYZ))
I_s2[s2] <- 1
I_s3 <- rep(0, nrow(XYZ))
I_s3[s3] <- 1
samples <- data.frame(c(rep("Scenario 1", n), rep("Scenario 2", n), rep("Scenario 3", n), rep("Scenario 4", n)),
                      c(rep(0, n), XYZ$x1[s1], XYZ$x1[s2], XYZ$x1[s3]),
                      c(rep(0, n), XYZ$x2[s1], XYZ$x2[s2], XYZ$x2[s3]),
                      c(XYZ$y[s1], krig_s1$var1.pred[s1], krig_s2$var1.pred[s2], krig_s3$var1.pred[s3]))
colnames(samples) <- c("Scenario", "x1", "x2", "y_pred")
pred_data <- data.frame(c(rep("Scenario 1", nrow(XYZ)), rep("Scenario 2", nrow(XYZ)), rep("Scenario 3", nrow(XYZ)), rep("Scenario 4", nrow(XYZ))),
                        c(XYZ$x1, XYZ$x1, XYZ$x1, XYZ$x1),
                        c(XYZ$x2, XYZ$x2, XYZ$x2, XYZ$x2),
                        c(XYZ$y, XYZ$y, XYZ$y, XYZ$y),
                        c(rep(0, nrow(XYZ)), XYZ$z.z1, XYZ$z.z2, XYZ$z.z3),
                        c(rep(1, nrow(XYZ)), I_s1, I_s2, I_s3),
                        c(XYZ$y, krig_s1$var1.pred, krig_s2$var1.pred, krig_s3$var1.pred),
                        c(rep(1, nrow(XYZ)), abs(krig_s1$var1.var), abs(krig_s2$var1.var),
                          abs(krig_s3$var1.var)))
colnames(pred_data) <- c("Scenario", "x1", "x2", "y", "z", "s", "y_pred","var")

pred_data$error <- pred_data$y - pred_data$y_pred
pred_data$pvalue <- 2 * (pnorm(abs(pred_data$error),
                               mean = 0, sd = sqrt(pred_data$var))) - 1
pred_data$pvalue2 <- .95 * (pred_data$pvalue <= .95) + pred_data$pvalue * (pred_data$pvalue > .95)



plot_prediction <- ggplot(data = pred_data, aes(x = x1, y = x2, fill = y_pred)) + geom_tile() +
  ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::scale_x_continuous(expand=c(0,0)) +
  ggplot2::scale_y_continuous(expand=c(0,0)) + 
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()) +
  ggplot2::scale_fill_gradient(low = "white", high = "black") +
  geom_point(data = samples, shape = 6, color = "black", fill = "white", size = 2, stroke = 1.2) +
  facet_wrap(~Scenario, nrow = 1) +
  ggplot2::labs(fill = "signal") + 
  ggplot2::theme(legend.position="bottom",
                 legend.justification="right",
                 legend.margin=ggplot2::margin(0,0,0,0),
                 legend.box.margin=ggplot2::margin(-20,#top
                                                   18,#right
                                                   10,#bottom
                                                   10))#left
print(plot_prediction)
plot_error <- ggplot(pred_data, aes(x = x1, y = x2, fill = (y - y_pred))) + geom_tile() +
  ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::scale_x_continuous(expand=c(0,0)) +
  ggplot2::scale_y_continuous(expand=c(0,0)) + 
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()) +
  ggplot2::scale_fill_gradient(low = "white", high = "black") +
  geom_point(data = samples, shape = 6, color = "black", fill = "white", size = 2, stroke = 1.2) + 
  facet_wrap(~Scenario, nrow = 1) +
  ggplot2::labs(fill = "error") +
  ggplot2::theme(legend.position="bottom",
                 legend.justification="right",
                 legend.margin=ggplot2::margin(0,0,0,0),
                 legend.box.margin=ggplot2::margin(-20,#top
                                                   18,#right
                                                   10,#bottom
                                                   10))#left
print(plot_error)
# samplepos<-
# pred_data$disttosample<-apply(pred_data[c("x1","x2")],1,
#                               function(x){min(dist(x,))})


# pred_data$pvalue2<-pmax(pred_data$pvalue,.5)
plot_q <- ggplot(pred_data, aes(x = x1, y = x2, fill = pvalue2)) + geom_tile() +
  ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::scale_x_continuous(expand=c(0,0)) +
  ggplot2::scale_y_continuous(expand=c(0,0)) + 
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()) +
  ggplot2::scale_fill_gradient2(low = "white", mid="white",midpoint=.95,high = "black") +
  geom_point(data = samples, shape = 6, color = "black", fill = "white", size = 2, stroke = 1.2) + 
  facet_wrap(~Scenario, nrow = 1) +
  ggplot2::labs(fill = "standardised error") +
  ggplot2::theme(legend.position="bottom",
                 legend.justification="right",
                 legend.margin=ggplot2::margin(0,0,0,0),
                 legend.box.margin=ggplot2::margin(-20,#top
                                                   18,#right
                                                   10,#bottom
                                                   10))#left

print(plot_q)

plot_var <- ggplot(pred_data, aes(x = x1, y = x2, fill = var)) + geom_tile() +
  ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::scale_x_continuous(expand=c(0,0)) +
  ggplot2::scale_y_continuous(expand=c(0,0)) + 
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()) +
  ggplot2::scale_fill_gradient(low = "white", high = "black") +
  geom_point(data = samples, shape = 6, color = "black", fill = "white", size = 2, stroke = 1.2) + 
  facet_wrap(~Scenario, nrow = 1) +
  ggplot2::labs(fill = "estimated prediction variance") +
  ggplot2::theme(legend.position="bottom",
                 legend.justification="right",
                 legend.margin=ggplot2::margin(0,0,0,0),
                 legend.box.margin=ggplot2::margin(-20,#top
                                                   18,#right
                                                   10,#bottom
                                                   10))#left

print(plot_var)
# ----------------- #



# # -- samples -- #
# # - just for now we take here the samples (we can play with the sample size n)
# # - and compute variograms, afterwards we generate them just once.
# N <- nrow(XYZ)
# n <- 100 # sample size
# M <- 100 # number of replications
# s1_M <- matrix(data = 0, nrow = M, ncol = n)
# s2_M <- matrix(data = 0, nrow = M, ncol = n)
# s3_M <- matrix(data = 0, nrow = M, ncol = n)
# set.seed(42)
# for(i in 1:M)
# {
#         s1_M[i, ] <- srs(n, XYZ, replace = F)
# }
# set.seed(42)
# for(i in 1:M)
# {
#         # Just for now we take pps samples by means of the function sample()
#         # and parameter prob, which is the wrong way to do it (sample does not
#         # respect perfectly that probabilities). We do this way for now to avoid
#         # long computing time.
#         # s2_M[i, ] <- pps(n, XYZ$z.z3, XYZ)
#         # s2_M[i, ] <- pps(n, XYZ$z.z2, XYZ)
#         pi <- inclusionprobabilities(a = XYZ$z.z2, n = n)
#         s2_M[i, ] <- sample(x = nrow(XYZ), size = n, prob = pi)
# }
# set.seed(42)
# for(i in 1:M)
# {
#         # Just for now we take pps samples by means of the function sample()
#         # and parameter prob, which is the wrong way to do it (sample does not
#         # respect perfectly that probabilities). We do this way for now to avoid
#         # long computing time.
#         # s3_M[i, ] <- pps(n, XYZ$z.z3, XYZ)
#         pi <- inclusionprobabilities(a = XYZ$z.z3, n = n)
#         s3_M[i, ] <- sample(x = nrow(XYZ), size = n, prob = pi)
# }
# # ------------- #
# XYZ_so <- XYZ
# coordinates(XYZ_so) <- ~ x1 + x2 # created a spatial objects for gstat package
# # -- variogram -- #
# # vgm_pop <- variogram(y ~ 1, XYZ_so) # population variogram, since we use all the field
# vgm_s1 <- list()
# length(vgm_s1) <- M
# vgm_s2 <- list()
# length(vgm_s2) <- M
# vgm_s3 <- list()
# length(vgm_s3) <- M
# for(i in 1:M) # sample variograms
# {
#         vgm_s1[[i]] <- variogram(y ~ 1, XYZ_so[s1_M[i,], ])
#         vgm_s2[[i]] <- variogram(y ~ 1, XYZ_so[s2_M[i,], ])
#         vgm_s3[[i]] <- variogram(y ~ 1, XYZ_so[s3_M[i,], ])
# }
# # --------------- #
# # -- fit variogram -- #
# fit_vgm_pop <- fit.variogram(vgm_pop, model = vgm(1, "Gau", 1))
# fit_vgm_s1 <- list()
# length(fit_vgm_s1) <- M
# fit_vgm_s2 <- list()
# length(fit_vgm_s2) <- M
# fit_vgm_s3 <- list()
# length(fit_vgm_s3) <- M
# for(i in 1:M)
# {
#         fit_vgm_s1[[i]] <- fit.variogram(vgm_s1[[i]], vgm(1, "Gau", 1))
#         fit_vgm_s2[[i]] <- fit.variogram(vgm_s2[[i]], vgm(1, "Gau", 1))
#         fit_vgm_s3[[i]] <- fit.variogram(vgm_s3[[i]], vgm(1, "Gau", 1))
# }
# # ------------- #
# # -- kriging -- #
# 
# # -- over the whole grid -- #
# krig_s1 <- list()
# length(krig_s1) <- M
# krig_s2 <- list()
# length(krig_s1) <- M
# krig_s3 <- list()
# length(krig_s1) <- M
# for(i in 1:M)
# {
#         krig_s1[[i]] <- krige(y ~ 1, XYZ_so[s1_M[i, ],], XYZ_so, model = fit_vgm_s1[[i]])
#         krig_s2[[i]] <- krige(y ~ 1, XYZ_so[s2_M[i, ],], XYZ_so, model = fit_vgm_s2[[i]])
#         krig_s3[[i]] <- krige(y ~ 1, XYZ_so[s3_M[i, ],], XYZ_so, model = fit_vgm_s3[[i]])
# }
# pred_whole_grid <- data.frame(XYZ$x1, XYZ$x2, XYZ$z.z3, t())
# krig_s1_bias <- rep(0, M)
# krig_s2_bias <- rep(0, M)
# krig_s3_bias <- rep(0, M)
# for(i in 1:M)
# {
#         krig_s1_bias[i] <- (sum(krig_s1[[i]]$var1.pred) - sum(XYZ_so$y)) / sum(XYZ_so$y)
#         krig_s2_bias[i] <- (sum(krig_s2[[i]]$var1.pred) - sum(XYZ_so$y)) / sum(XYZ_so$y)
#         krig_s3_bias[i] <- (sum(krig_s3[[i]]$var1.pred) - sum(XYZ_so$y)) / sum(XYZ_so$y)
# }
# mean(krig_s1_bias)
# mean(krig_s2_bias)
# mean(krig_s3_bias)
# krig_s1_mse <- rep(0, M)
# krig_s2_mse <- rep(0, M)
# krig_s3_mse <- rep(0, M)
# for(i in 1:M)
# {
#         krig_s1_mse[i] <- sum((krig_s1[[i]]$var1.pred - XYZ_so$y) ^ 2) / N
#         krig_s2_mse[i] <- sum((krig_s2[[i]]$var1.pred - XYZ_so$y) ^ 2) / N
#         krig_s3_mse[i] <- sum((krig_s3[[i]]$var1.pred - XYZ_so$y) ^ 2) / N
# }
# mean(krig_s1_mse)
# mean(krig_s2_mse)
# mean(krig_s3_mse)
# # ------------------------- #
# # ------------- #
}