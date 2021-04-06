# ------ library ------ #
library(RandomFields)
library(sampling)
library(gstat)
library(ggplot2)
# --------------------- #
# ------ generate fields ------ #
grid_size <- 30
set.seed(42)
Y <- generateYonsquaregridU(length.out = grid_size)
Z <- generateZconditionnallyongridUY(length.out = grid_size, y = Y)
XYZ <- bindXYZpiI(x1 = generatex_i(grid_size), x2 = generatex_i(grid_size), y = Y, z = Z)
# ----------------------------- #
# ------------ simulation ------------ #
# -- samples -- #
n <- 200 # sample size
M <- 100 # number of replications
s_srs <- matrix(data = 0, nrow = M, ncol = n)
s_pps1 <- matrix(data = 0, nrow = M, ncol = n)
s_pps2 <- matrix(data = 0, nrow = M, ncol = n)
set.seed(42)
for(i in 1:M)
{
  s_srs[i, ] <- which(srswor(n, nrow(XYZ)) == 1) 
}
pik_pps1 <- inclusionprobabilities(XYZ$z.z1, n)
set.seed(42)
for(i in 1:M)
{
  s_pps1[i, ] <- which(UPmaxentropy(pik_pps1) == 1) # UPmaxentropy() looks faster than UPtille()
}
pik_pps2 <- inclusionprobabilities(XYZ$z.z2, n)
set.seed(42)
for(i in 1:M)
{
  s_pps2[i, ] <- which(UPmaxentropy(pik_pps2) == 1)
}
# ------------- #
# ------------ kriging ------------ #
XYZ_so <- XYZ
coordinates(XYZ_so) <- ~ x1 + x2 # created a spatial objects for gstat package
# -- variogram -- #
vgm_pop <- variogram(y ~ 1, XYZ_so) # population variogram, since we use all the field
vgm_srs <- list()
length(vgm_srs) <- M
vgm_pps1 <- list()
length(vgm_pps1) <- M
vgm_pps2 <- list()
length(vgm_pps2) <- M
for(i in 1:M) # sample variograms
{
  vgm_srs[[i]] <- variogram(y ~ 1, XYZ_so[s_srs[i,], ])
  vgm_pps1[[i]] <- variogram(y ~ 1, XYZ_so[s_pps1[i,], ])
  vgm_pps2[[i]] <- variogram(y ~ 1, XYZ_so[s_pps2[i,], ])
}
# --------------- #
# -- fit variogram -- #
# fit_vgm_pop <- fit.variogram(vgm_pop, model = vgm(1, "Gau", 300, 1))
fit_vgm_pop <- fit.variogram(vgm_pop, model = vgm(1, "Gau", 1))
fit_vgm_srs <- list()
length(fit_vgm_srs) <- M
fit_vgm_pps1 <- list()
length(fit_vgm_pps1) <- M
fit_vgm_pps2 <- list()
length(fit_vgm_pps2) <- M
for(i in 1:M)
{
  # fit_vgm_srs[[i]] <- fit.variogram(vgm_srs[[i]], vgm(1, "Gau", 300, 1))
  # fit_vgm_pps1[[i]] <- fit.variogram(vgm_pps1[[i]], vgm(1, "Gau", 300, 1))
  # fit_vgm_pps2[[i]] <- fit.variogram(vgm_pps2[[i]], vgm(1, "Gau", 300, 1))
  fit_vgm_srs[[i]] <- fit.variogram(vgm_srs[[i]], vgm(1, "Gau", 1))
  fit_vgm_pps1[[i]] <- fit.variogram(vgm_pps1[[i]], vgm(1, "Gau", 1))
  fit_vgm_pps2[[i]] <- fit.variogram(vgm_pps2[[i]], vgm(1, "Gau", 1))
}
# ------------------- #
# -- prediction -- #
# XYZ_grid <- XYZ[, 1:2]
# coordinates(XYZ_grid) <- ~ x1 + x2
# pred_srs <- list()
# length(pred_srs) <- M
# pred_pps1 <- list()
# length(pred_pps1) <- M
# pred_pps2 <- list()
# length(pred_pps2) <- M
# for(i in 1:M)
# {
#   pred_srs[[i]] <- krige(formula = y ~ 1, locations = XYZ_so[s_srs[i,], ], newdata = XYZ_grid, model = fit_vgm_srs[[i]])
#   pred_pps1[[i]] <- krige(formula = y ~ 1, locations = XYZ_so[s_pps1[i,], ], newdata = XYZ_grid, model = fit_vgm_pps1[[i]])
#   pred_pps2[[i]] <- krige(formula = y ~ 1, locations = XYZ_so[s_pps2[i,], ], newdata = XYZ_grid, model = fit_vgm_pps2[[i]])
# }
# ---------------- #
# -- bias should be here -- #

# ------------ plot ------------ #
find_xy_plot <- function(vgm, fit_vgm)
{
  # this function is used to find the range of x and y, which is useful for have a nice plot
  # I guess there is a better way to do that, but for now it is handy.
  # vgm is a list, which contains all the sample variograms (M variograms)
  # fit_vgm is a list, which contains all the variogram fits (M variogram fits)
  M <- length(fit_vgm)
  xmin <- 100
  xmax <- 0
  ymin <- 100
  ymax <- 0
  for(i in 1:M)
  {
    if(min(variogramLine(fit_vgm[[i]], maxdist = max(vgm[[i]]$dist))$dist) < xmin)
    {
      xmin <- min(variogramLine(fit_vgm[[i]], maxdist = max(vgm[[i]]$dist))$dist)
    }
    if(max(variogramLine(fit_vgm[[i]], maxdist = max(vgm[[i]]$dist))$dist) > xmax)
    {
      xmax <- max(variogramLine(fit_vgm[[i]], maxdist = max(vgm[[i]]$dist))$dist)
    }
    if(min(variogramLine(fit_vgm[[i]], maxdist = max(vgm[[i]]$dist))$gamma) < ymin)
    {
      ymin <- min(variogramLine(fit_vgm[[i]], maxdist = max(vgm[[i]]$dist))$gamma)
    }
    if(max(variogramLine(fit_vgm[[i]], maxdist = max(vgm[[i]]$dist))$gamma) > ymax)
    {
      ymax <- max(variogramLine(fit_vgm[[i]], maxdist = max(vgm[[i]]$dist))$gamma)
    }
  }
  xy_plot <- c(xmin, xmax, ymin, ymax)
  names(xy_plot) <- c("xmin", "xmax", "ymin", "ymax")
  return(as.data.frame(t(xy_plot)))
}



par(mfrow = c(1, 3), cex = 0.5)
srs_mean <- 0
xy <- find_xy_plot(vgm_srs, fit_vgm_srs)
plot(variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist)), xlim = c(xy$xmin, xy$xmax), ylim = c(xy$ymin, xy$ymax), 
     type = "l", lwd = 3, main = "SRS")
for(i in 1:M)
{
  lines(variogramLine(fit_vgm_srs[[i]], maxdist = max(vgm_srs[[i]]$dist)), type = "l", lwd = 0.25)
  srs_mean <- srs_mean + variogramLine(fit_vgm_srs[[i]], maxdist = max(vgm_srs[[i]]$dist))
}
srs_mean <- srs_mean / M
lines(srs_mean, lty = 2, lwd = 3)
legend("bottomright", legend = c("pop", "repl", "repl mean"), lwd = c(3, 0.25, 3), lty = c(1, 1, 2), cex = 0.6)
pps1_mean <- 0
xy <- find_xy_plot(vgm_pps1, fit_vgm_pps1)
plot(variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist)), xlim = c(xy$xmin, xy$xmax), ylim = c(xy$ymin, xy$ymax), 
     type = "l", lwd = 3, main = "PPS % Z1")
for(i in 1:M)
{
  lines(variogramLine(fit_vgm_pps1[[i]], maxdist = max(vgm_pps1[[i]]$dist)), type = "l", lwd = 0.25)
  pps1_mean <- pps1_mean + variogramLine(fit_vgm_pps1[[i]], maxdist = max(vgm_pps1[[i]]$dist))
}
pps1_mean <- pps1_mean / M
lines(pps1_mean, lty = 2, lwd = 3)
legend("bottomright", legend = c("pop", "repl", "repl mean"), lwd = c(3, 0.25, 3), lty = c(1, 1, 2), cex = 0.6)
pps2_mean <- 0
xy <- find_xy_plot(vgm_pps2, fit_vgm_pps2)
plot(variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist)), xlim = c(xy$xmin, xy$xmax), ylim = c(xy$ymin, xy$ymax), 
     type = "l", lwd = 3, main = "PPS % Z2")
for(i in 1:M)
{
  lines(variogramLine(fit_vgm_pps2[[i]], maxdist = max(vgm_pps2[[i]]$dist)), type = "l", lwd = 0.25)
  pps2_mean <- pps2_mean + variogramLine(fit_vgm_pps2[[i]], maxdist = max(vgm_pps2[[i]]$dist))
}
pps2_mean <- pps2_mean / M
lines(pps2_mean, lty = 2, lwd = 3)
legend("bottomright", legend = c("pop", "repl", "repl mean"), lwd = c(3, 0.25, 3), lty = c(1, 1, 2), cex = 0.6)
par(mfrow = c(1, 1), cex = 1)
# ------------------------------ #




