# ------------ library and data ------------ #
library(sampling)
library(samplingbook)
library(sp)
library(gstat)
library(Spbsampling)
library(ggplot2)
data(meuse, package = "sp")
# ------------------------------------------ #
# ------------ draw samples ------------ #
npop <- nrow(meuse)
nsample <- 60
niter <- 100
# -- srs samples -- #
srs_sample <- matrix(data = 0, niter, nsample)
set.seed(42)
for(i in 1:niter)
{
  srs_sample[i, ] <- which(srswor(nsample, npop) == 1)
}
# ----------------- #
# -- pps samples -- #
pps_sample <- matrix(data = 0, niter, nsample)
x <- meuse$copper
set.seed(42)
for(i in 1:niter)
{
  pps_sample[i, ] <- pps.sampling(meuse$copper, nsample, method = "tille")$sample
}
# ----------------- #
# -- PWD spread samples -- #
pwd_s_sample <- matrix(data = 0, niter, nsample)
dis <- as.matrix(dist(cbind(meuse$x, meuse$y)))
std_dis <- stprod(mat = dis, con = rep(0, npop))
set.seed(42)
pwd_s_sample <- pwd(dis = std_dis, n = nsample, nrepl = niter)
# ----------------- #
# -- PWD cluster samples -- #
pwd_c_sample <- matrix(data = 0, niter, nsample)
dis <- as.matrix(dist(cbind(meuse$x, meuse$y)))
std_dis <- stprod(mat = dis, con = rep(0, npop))
set.seed(42)
pwd_c_sample <- pwd(dis = std_dis, n = nsample, beta = -10, nrepl = niter)
# ----------------- #
# -------------------------------------- #
# ------------ kriging ------------ #
coordinates(meuse) <- ~ x + y
# -- variogram -- #
vgm_pop <- variogram(log(zinc) ~ 1, meuse)
vgm_srs <- list()
length(vgm_srs) <- niter
vgm_pps <- list()
length(vgm_pps) <- niter
vgm_s_pwd <- list()
length(vgm_s_pwd) <- niter
vgm_c_pwd <- list()
length(vgm_c_pwd) <- niter
for(i in 1:niter)
{
  vgm_srs[[i]] <- variogram(log(zinc) ~ 1, meuse[srs_sample[i,], ])
  vgm_pps[[i]] <- variogram(log(zinc) ~ 1, meuse[pps_sample[i,], ])
  vgm_s_pwd[[i]] <- variogram(log(zinc) ~ 1, meuse[pwd_s_sample[i,], ])
  vgm_c_pwd[[i]] <- variogram(log(zinc) ~ 1, meuse[pwd_c_sample[i,], ])
}
# --------------- #
# -- fit variogram -- #
fit_vgm_pop <- fit.variogram.reml(log(zinc) ~ 1, meuse, model = vgm(1, "Gau", 300, 1))
fit_vgm_srs <- list()
length(fit_vgm_srs) <- niter
fit_vgm_pps <- list()
length(fit_vgm_pps) <- niter
fit_vgm_s_pwd <- list()
length(fit_vgm_s_pwd) <- niter
fit_vgm_c_pwd <- list()
length(fit_vgm_c_pwd) <- niter
for(i in 1:niter)
{
  fit_vgm_srs[[i]] <- fit.variogram.reml(log(zinc) ~ 1, meuse[srs_sample[i,], ], model = vgm(1, "Gau", 300, 1))
  fit_vgm_pps[[i]] <- fit.variogram.reml(log(zinc) ~ 1, meuse[pps_sample[i,], ], model = vgm(1, "Gau", 300, 1))
  fit_vgm_s_pwd[[i]] <- fit.variogram.reml(log(zinc) ~ 1, meuse[pwd_s_sample[i,], ], model = vgm(1, "Gau", 300, 1))
  fit_vgm_c_pwd[[i]] <- fit.variogram.reml(log(zinc) ~ 1, meuse[pwd_c_sample[i,], ], model = vgm(1, "Gau", 300, 1))
}
# ------------------- #
# --------------------------------- #
# ------------ plot ------------ #
vgm_pop1 <- vgm_pop
vgm_pop1$id <- "pop"
vgm_srs1 <- list()
length(vgm_srs1) <- niter
vgm_srs1 <- vgm_srs
vgm_pps1 <- list()
length(vgm_pps1) <- niter
vgm_pps1 <- vgm_pps
vgm_s_pwd1 <- list()
length(vgm_s_pwd1) <- niter
vgm_s_pwd1 <- vgm_s_pwd
vgm_c_pwd1 <- list()
length(vgm_c_pwd1) <- niter
vgm_c_pwd1 <- vgm_s_pwd
vgLine <- list()
length(vgLine) <- niter

for(i in 1:niter)
{
  vgm_srs1[[i]]$id <- "srs"
  vgm_pps1[[i]]$id <- "pps"
  vgm_s_pwd1[[i]]$id <- "pwd_s"
  vgm_c_pwd1[[i]]$id <- "pwd_c"
  vgLine[[i]] <- rbind( 
    cbind(variogramLine(fit_vgm_pop, maxdist = max(vgm_pop1$dist)), id = 
            "pop"), 
    cbind(variogramLine(fit_vgm_srs[[i]], maxdist = max(vgm_srs1[[i]]$dist)), id = 
            "srs"),
    cbind(variogramLine(fit_vgm_pps[[i]], maxdist = max(vgm_pps1[[i]]$dist)), id = 
            "pps"),
    cbind(variogramLine(fit_vgm_s_pwd[[i]], maxdist = max(vgm_s_pwd1[[i]]$dist)), id = 
            "pwd_s"),
    cbind(variogramLine(fit_vgm_c_pwd[[i]], maxdist = max(vgm_c_pwd1[[i]]$dist)), id = 
            "pwd_c")
  ) 
}
vgm_pop11 <- list()
length(vgm_pop11) <- niter
for(i in 1:niter)
{
  vgm_pop11[[i]] <- vgm_pop1
}

par(mfrow = c(2,2), cex = 0.50)

srs_mean <- matrix(0, nrow(subset(vgLine[[1]], vgLine[[1]]$id == "srs")[,1:2]), ncol(subset(vgLine[[1]], vgLine[[1]]$id == "srs")[,1:2]))
srs_mean <- 0
plot(subset(vgLine[[1]], vgLine[[1]]$id == "pop")[,1:2], ylim = c(0, 1), type = "l", col = "red", lwd = 1.42, main = "SRS")
for(i in 1:niter)
{
  lines(subset(vgLine[[i]], vgLine[[i]]$id == "srs")[,1:2], type = "l", lwd = 0.24)
  srs_mean <- srs_mean + subset(vgLine[[i]], vgLine[[i]]$id == "srs")[,1:2]
}
srs_mean <- srs_mean / niter
lines(srs_mean, col = "blue", lwd = 1.42)
legend("topleft", legend = c("pop", "repl mean", "repl"), col = c("red", "blue", "black"), lty = 1)

pps_mean <- matrix(0, nrow(subset(vgLine[[1]], vgLine[[1]]$id == "pps")[,1:2]), ncol(subset(vgLine[[1]], vgLine[[1]]$id == "pps")[,1:2]))
pps_mean <- 0
plot(subset(vgLine[[1]], vgLine[[1]]$id == "pop")[,1:2], ylim = c(0, 1), type = "l", col = "red", lwd = 1.42, main = "PPS % copper")
for(i in 1:niter)
{
  lines(subset(vgLine[[i]], vgLine[[i]]$id == "pps")[,1:2], type = "l", lwd = 0.24)
  pps_mean <- pps_mean + subset(vgLine[[i]], vgLine[[i]]$id == "pps")[,1:2]
}
pps_mean <- pps_mean / niter
lines(pps_mean, col = "blue", lwd = 1.42)
legend("topleft", legend = c("pop", "repl mean", "repl"), col = c("red", "blue", "black"), lty = 1)

pwd_s_mean <- matrix(0, nrow(subset(vgLine[[1]], vgLine[[1]]$id == "pwd_s")[,1:2]), ncol(subset(vgLine[[1]], vgLine[[1]]$id == "pwd_s")[,1:2]))
pwd_s_mean <- 0
plot(subset(vgLine[[1]], vgLine[[1]]$id == "pop")[,1:2], ylim = c(0, 1), type = "l", col = "red", lwd = 1.42, main = "spread PWD")
for(i in 1:niter)
{
  lines(subset(vgLine[[i]], vgLine[[i]]$id == "pwd_s")[,1:2], type = "l", lwd = 0.24)
  pwd_s_mean <- pwd_s_mean + subset(vgLine[[i]], vgLine[[i]]$id == "pwd_s")[,1:2]
}
pwd_s_mean <- pwd_s_mean / niter
lines(pwd_s_mean, col = "blue", lwd = 1.42)
legend("topleft", legend = c("pop", "repl mean", "repl"), col = c("red", "blue", "black"), lty = 1)

pwd_c_mean <- matrix(0, nrow(subset(vgLine[[1]], vgLine[[1]]$id == "pwd_c")[,1:2]), ncol(subset(vgLine[[1]], vgLine[[1]]$id == "pwd_c")[,1:2]))
pwd_c_mean <- 0
plot(subset(vgLine[[1]], vgLine[[1]]$id == "pop")[,1:2], ylim = c(0, 1), type = "l", col = "red", lwd = 1.42, main = "cluster PWD")
for(i in 1:niter)
{
  lines(subset(vgLine[[i]], vgLine[[i]]$id == "pwd_c")[,1:2], type = "l", lwd = 0.24)
  pwd_c_mean <- pwd_c_mean + subset(vgLine[[i]], vgLine[[i]]$id == "pwd_c")[,1:2]
}
pwd_c_mean <- pwd_c_mean / niter
lines(pwd_c_mean, col = "blue", lwd = 1.42)
legend("topleft", legend = c("pop", "repl mean", "repl"), col = c("red", "blue", "black"), lty = 1)

par(mfrow = c(1,1), cex = 1)
# ------------------------------ #

