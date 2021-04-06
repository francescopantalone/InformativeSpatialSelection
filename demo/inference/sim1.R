# ------------ library and data ------------ #
library(sampling)
library(samplingbook)
library(sp)
library(gstat)
library(ggplot2)
data(meuse, package = "sp")
# ------------------------------------------ #
# ------------ draw samples ------------ #
npop <- nrow(meuse)
n <- 60
M <- 100
# -- srs samples -- #
s_srs <- matrix(data = 0, M, n)
set.seed(42)
for(i in 1:M)
{
  s_srs[i, ] <- which(srswor(n, npop) == 1)
}
# ----------------- #
# -- pps samples -- #
s_pps <- matrix(data = 0, M, n)
x <- meuse$copper
set.seed(42)
for(i in 1:M)
{
  s_pps[i, ] <- pps.sampling(meuse$copper, n, method = "tille")$sample
}
# ----------------- #
# -------------------------------------- #
# ------------ kriging ------------ #
coordinates(meuse) <- ~ x + y
# -- variogram -- #
vgm_pop <- variogram(log(zinc) ~ 1, meuse)
vgm_srs <- list()
length(vgm_srs) <- M
vgm_pps <- list()
length(vgm_pps) <- M
for(i in 1:M)
{
  vgm_srs[[i]] <- variogram(log(zinc) ~ 1, meuse[s_srs[i,], ])
  vgm_pps[[i]] <- variogram(log(zinc) ~ 1, meuse[s_pps[i,], ])
}
# --------------- #
# -- fit variogram -- #
fit_vgm_pop <- fit.variogram(vgm_pop, vgm(1, "Gau", 300, 1))
fit_vgm_srs <- list()
length(fit_vgm_srs) <- M
fit_vgm_pps <- list()
length(fit_vgm_pps) <- M
for(i in 1:M)
{
  fit_vgm_srs[[i]] <- fit.variogram(vgm_srs[[i]], vgm(1, "Gau", 300, 1))
  fit_vgm_pps[[i]] <- fit.variogram(vgm_pps[[i]], vgm(1, "Gau", 300, 1))
}
# ------------------- #
# --------------------------------- #
# ------------ plot ------------ #
par(mfrow = c(1, 2), cex = 0.5)
srs_mean <- 0
plot(variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist)), ylim = c(0, 1), type = "l", col = "red", lwd = 1.42, main = "SRS")
for(i in 1:M)
{
  lines(variogramLine(fit_vgm_srs[[i]], maxdist = max(vgm_srs[[i]]$dist)), type = "l", lwd = 0.24)
  srs_mean <- srs_mean + variogramLine(fit_vgm_srs[[i]], maxdist = max(vgm_srs[[i]]$dist))
}
srs_mean <- srs_mean / M
lines(srs_mean, col = "blue", lwd = 1.42)
legend("topleft", legend = c("pop", "repl mean", "repl"), col = c("red", "blue", "black"), lty = 1)
pps_mean <- 0
plot(variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist)), ylim = c(0, 1), type = "l", col = "red", lwd = 1.42, main = "PPS % copper")
for(i in 1:M)
{
  lines(variogramLine(fit_vgm_pps[[i]], maxdist = max(vgm_pps[[i]]$dist)), type = "l", lwd = 0.24)
  pps_mean <- pps_mean + variogramLine(fit_vgm_pps[[i]], maxdist = max(vgm_pps[[i]]$dist))
}
pps_mean <- pps_mean / M
lines(pps_mean, col = "blue", lwd = 1.42)
legend("topleft", legend = c("pop", "repl mean", "repl"), col = c("red", "blue", "black"), lty = 1)
par(mfrow = c(1, 1), cex = 1)
# ------------------------------ #