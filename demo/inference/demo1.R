# ------------ library ------------ #
library(sampling)
library(samplingbook)
library(sp)
library(lattice)
library(gstat)
library(Spbsampling)
library(ggplot2)
data(meuse)
# --------------------------------- #
# ------------ draw samples ------------ #
npop <- nrow(meuse)
nsample <- 30
# -- srs sample -- #
set.seed(42)
srs_sample <- which(srswor(nsample, npop) == 1)
# ---------------- #
# -- pps sample -- #
set.seed(42)
pps_sample <- pps.sampling(meuse$copper, nsample, method = "tille")$sample
# ---------------- #
# -- pwd sample -- #
dis <- as.matrix(dist(cbind(meuse$x, meuse$y)))
std_dis <- stprod(dis ^ 10, rep(0, npop))
set.seed(42)
pwd_sample <- pwd(dis = std_dis, nsamp = nsample)
# ---------------- #
# -------------------------------------- #
# ------------ kriging ------------ #
# -- variogram -- #
coordinates(meuse) <- ~ x + y # converting to a spatial data frame
vgm_pop <- variogram(log(zinc) ~ 1, meuse)
vgm_srs <- variogram(log(zinc) ~ 1, meuse[srs_sample, ])
vgm_pps <- variogram(log(zinc) ~ 1, meuse[pps_sample, ])
vgm_pwd <- variogram(log(zinc) ~ 1, meuse[as.vector(pwd_sample), ])
# --------------- #
# -- fit variogram -- #
fit_vgm_pop <- fit.variogram(vgm_pop, vgm(1, "Gau", 300, 1))
fit_vgm_srs <- fit.variogram(vgm_srs, vgm(1, "Gau", 300, 1))
fit_vgm_pps <- fit.variogram(vgm_pps, vgm(1, "Gau", 300, 1))
fit_vgm_pwd <- fit.variogram(vgm_pwd, vgm(1, "Gau", 300, 1))
# ------------------- #
# -- predictions -- #
pred_pop <- variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist))
pred_srs <- variogramLine(fit_vgm_srs, maxdist = max(vgm_srs$dist))
pred_pps <- variogramLine(fit_vgm_pps, maxdist = max(vgm_pps$dist))
pred_pwd <- variogramLine(fit_vgm_pwd, maxdist = max(vgm_pwd$dist))
# ----------------- #
# --------------------------------- #
# ------------ plot ------------ #
vgm_pop1 <- vgm_pop
vgm_pop1$id <- "pop"
vgm_srs1 <- vgm_srs
vgm_srs1$id <- "srs"
vgm_pps1 <- vgm_pps
vgm_pps1$id <- "pps"
vgm_pwd1 <- vgm_pwd
vgm_pwd1$id <- "pwd"
vgLine <- rbind(cbind(variogramLine(fit_vgm_pop, maxdist = max(vgm_pop1$dist)), id ="pop"), 
  cbind(variogramLine(fit_vgm_srs, maxdist = max(vgm_srs1$dist)), id ="srs"),
  cbind(variogramLine(fit_vgm_pps, maxdist = max(vgm_pps1$dist)), id ="pps"),
  cbind(variogramLine(fit_vgm_pwd, maxdist = max(vgm_pwd1$dist)), id ="pwd")) 
ggplot(rbind(vgm_pop1, vgm_srs1, vgm_pps1, vgm_pwd1), aes(x = dist, y = gamma, colour = id)) + geom_line(data = vgLine) + 
  geom_point() 
# ------------------------------ #
