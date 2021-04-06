# ------------ library and data ------------ #
library(sampling)
library(samplingbook)
library(sp)
library(gstat)
library(geoR)
library(Spbsampling)
library(ggplot2)
data(meuse, package="sp")
# ------------------------------------------ #
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
# -- variogram -- #
coordinates(meuse) <- ~ x + y # converting to a spatial data frame
vgm_pop <- variogram(log(zinc) ~ 1, meuse)
vgm_srs <- variogram(log(zinc) ~ 1, meuse[srs_sample, ])
vgm_pps <- variogram(log(zinc) ~ 1, meuse[pps_sample, ])
vgm_pwd <- variogram(log(zinc) ~ 1, meuse[as.vector(pwd_sample), ])
# --------------- #
# -- fit variogram -- #
fit_vgm_pop <- fit.variogram.reml(log(zinc) ~ 1, meuse, model = vgm(1, "Gau", 300, 1))
fit_vgm_srs <- fit.variogram.reml(log(zinc) ~ 1, meuse[srs_sample, ], model = vgm(1, "Gau", 300, 1))
fit_vgm_pps <- fit.variogram.reml(log(zinc) ~ 1, meuse[pps_sample, ], model = vgm(1, "Gau", 300, 1))
fit_vgm_pwd <- fit.variogram.reml(log(zinc) ~ 1, meuse[as.vector(pwd_sample), ], model = vgm(1, "Gau", 300, 1))
# ------------------- #
# ------------ plot ------------ #
vgm_pop1 <- vgm_pop
vgm_pop1$id <- "pop"
vgm_srs1 <- vgm_srs
vgm_srs1$id <- "srs"
vgm_pps1 <- vgm_pps
vgm_pps1$id <- "pps"
vgm_pwd1 <- vgm_pwd
vgm_pwd1$id <- "pwd"
vgLine <- rbind( 
  cbind(variogramLine(fit_vgm_pop, maxdist = max(vgm_pop1$dist)), id = 
          "pop"), 
  cbind(variogramLine(fit_vgm_srs, maxdist = max(vgm_srs1$dist)), id = 
          "srs"),
  cbind(variogramLine(fit_vgm_pps, maxdist = max(vgm_pps1$dist)), id = 
          "pps"),
  cbind(variogramLine(fit_vgm_pwd, maxdist = max(vgm_pwd1$dist)), id = 
          "pwd")
) 
ggplot(rbind(vgm_pop1, vgm_srs1, vgm_pps1, vgm_pwd1), aes(x = dist, y = gamma, colour = id)) + geom_line(data = vgLine) + 
  geom_point()
# ------------------------------ #





