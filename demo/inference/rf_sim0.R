# ------------ library ------------ #
library(sampling)
library(samplingbook)
library(sp)
library(lattice)
library(gstat)
library(Spbsampling)
library(ggplot2)
library(RandomFields)
library(dplyr)
# --------------------------------- #
# ------------ simulation gaussian random field ------------ #
set.seed(42)
length.x=200
x1 <- generatex_i(length.x)
g_m_y <- RFsimulate(model = RMgauss() + RMtrend(mean = 10), x = 10*x1,y=10*x1)
grf <- cbind(generatesquaregridU(x1=x1, x2=x1), y = g_m_y@data[[1]])

ggplot(grf, aes(x1,x2)) + 
  geom_tile(aes(fill = y)) + 
  scale_fill_gradient(low = "white",     high = "black")


z <- grf$y * 2 # proxy variable
grf <- cbind(grf, z)

# ---------------------------------------------------------- #
# ------------ draw samples ------------ #
npop <- nrow(grf)
nsample <- 150
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
set.seed(42)
for(i in 1:niter)
{
  pps_sample[i, ] <- pps.sampling(grf$z, nsample, method = "tille")$sample
}
# ----------------- #
# -- PWD spread samples -- #
pwd_s_sample <- matrix(data = 0, niter, nsample)
dis <- as.matrix(dist(cbind(grf$x1, grf$x2)))
std_dis <- stprod(dis ^ 10, rep(0, npop))
set.seed(42)
pwd_s_sample <- pwd(dis = std_dis, nsamp = nsample, nrepl = niter)
# ----------------- #
# -- PWD cluster samples -- #
pwd_c_sample <- matrix(data = 0, niter, nsample)
dis <- as.matrix(dist(cbind(grf$x1, grf$x2)))
std_dis <- stprod(dis ^ -10, rep(0, npop))
set.seed(42)
pwd_c_sample <- pwd(dis = std_dis, nsamp = nsample, nrepl = niter)
# ----------------- #
# -------------------------------------- #
# ------------ kriging ------------ #
coordinates(grf) <- ~ x1 + x2
# -- variogram -- #
vgm_pop <- variogram(y ~ 1, grf)
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
  vgm_srs[[i]] <- variogram(y ~ 1, grf[srs_sample[i,], ])
  vgm_pps[[i]] <- variogram(y ~ 1, grf[pps_sample[i,], ])
  vgm_s_pwd[[i]] <- variogram(y ~ 1, grf[pwd_s_sample[i,], ])
  vgm_c_pwd[[i]] <- variogram(y ~ 1, grf[pwd_c_sample[i,], ])
}
# --------------- #
# -- fit variogram -- #
fit_vgm_pop <- fit.variogram.reml(y ~ 1, grf, model = vgm(1, "Gau", 300, 1)) # too computationally demanding when we use big grid
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
  fit_vgm_srs[[i]] <- fit.variogram.reml(y ~ 1, grf[srs_sample[i,], ], model = vgm(1, "Gau", 300, 1))
  fit_vgm_pps[[i]] <- fit.variogram.reml(y ~ 1, grf[pps_sample[i,], ], model = vgm(1, "Gau", 300, 1))
  fit_vgm_s_pwd[[i]] <- fit.variogram.reml(y ~ 1, grf[pwd_s_sample[i,], ], model = vgm(1, "Gau", 300, 1))
  fit_vgm_c_pwd[[i]] <- fit.variogram.reml(y ~ 1, grf[pwd_c_sample[i,], ], model = vgm(1, "Gau", 300, 1)) # with cluster sample, singular matrix obtained in the fit "singular V matrix in calc_VinvIminAw"
}
# ------------------- #
# -- prediction -- #
coordinates(grf_grid) <- ~ x1 + x2
pred_srs <- list()
length(pred_srs) <- niter
pred_pps <- list()
length(pred_pps) <- niter
pred_s_pwd <- list()
length(pred_s_pwd) <- niter
pred_c_pwd <- list()
length(pred_c_pwd) <- niter
for(i in 1:niter)
{
  pred_srs[[i]] <- krige(y ~ 1, grf, grf_grid, model = fit_vgm_srs[[i]])
  pred_pps[[i]] <- krige(y ~ 1, grf, grf_grid, model = fit_vgm_pps[[i]])
  pred_s_pwd[[i]] <- krige(y ~ 1, grf, grf_grid, model = fit_vgm_s_pwd[[i]])
  pred_c_pwd[[i]] <- krige(y ~ 1, grf, grf_grid, model = fit_vgm_c_pwd[[i]]) # with cluster sample, "In predict.gstat(g, newdata = newdata, block = block,  ... :Covariance matrix singular at location [0.473684,0.105263,0]: skipping..."
}
# ---------------- #
# -- bias -- #
mean_srs_pre <- rep(0, nrow(as.data.frame(grf_grid)))
for(i in 1:nrow(as.data.frame(grf_grid)))
{
  for(j in 1:niter)
  {
    mean_srs_pre[i] <- mean_srs_pre[i] + pred_srs[[j]]$var1.pred[i]  
  }
}
mean_srs_pre <- mean_srs_pre / niter

bias_srs_pre <- rep(0, nrow(as.data.frame(grf_grid)))
for(i in 1:nrow(as.data.frame(grf_grid)))
{
  bias_srs_pre[i] <- (mean_srs_pre[i] - grf$y[i]) / grf$y[i] # attention: grf here is a spatial object (and not a matrix)
}

mean_pps_pre <- rep(0, nrow(as.data.frame(grf_grid)))
for(i in 1:nrow(as.data.frame(grf_grid)))
{
  for(j in 1:niter)
  {
    mean_pps_pre[i] <- mean_pps_pre[i] + pred_pps[[j]]$var1.pred[i]  
  }
}
mean_pps_pre <- mean_pps_pre / niter

bias_pps_pre <- rep(0, nrow(as.data.frame(grf_grid)))
for(i in 1:nrow(as.data.frame(grf_grid)))
{
  bias_pps_pre[i] <- (mean_pps_pre[i] - grf$y[i]) / grf$y[i] # attention: grf here is a spatial object (and not a matrix)
}

mean_s_pwd_pre <- rep(0, nrow(as.data.frame(grf_grid)))
for(i in 1:nrow(as.data.frame(grf_grid)))
{
  for(j in 1:niter)
  {
    mean_s_pwd_pre[i] <- mean_s_pwd_pre[i] + pred_s_pwd[[j]]$var1.pred[i]  
  }
}
mean_s_pwd_pre <- mean_s_pwd_pre / niter

bias_s_pwd_pre <- rep(0, nrow(as.data.frame(grf_grid)))
for(i in 1:nrow(as.data.frame(grf_grid)))
{
  bias_s_pwd_pre[i] <- (mean_s_pwd_pre[i] - grf$y[i]) / grf$y[i] # attention: grf here is a spatial object (and not a matrix)
}


mean_c_pwd_pre <- rep(0, nrow(as.data.frame(grf_grid)))
for(i in 1:nrow(as.data.frame(grf_grid)))
{
  for(j in 1:niter)
  {
    mean_c_pwd_pre[i] <- mean_c_pwd_pre[i] + pred_s_pwd[[j]]$var1.pred[i]  
  }
}
mean_c_pwd_pre <- mean_c_pwd_pre / niter

bias_c_pwd_pre <- rep(0, nrow(as.data.frame(grf_grid)))
for(i in 1:nrow(as.data.frame(grf_grid)))
{
  bias_c_pwd_pre[i] <- (mean_c_pwd_pre[i] - grf$y[i]) / grf$y[i] # attention: grf here is a spatial object (and not a matrix)
}
# ---------- #
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
plot(subset(vgLine[[1]], vgLine[[1]]$id == "pop")[,1:2], type = "l", col = "red", lwd = 1.42, main = "SRS")
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
plot(subset(vgLine[[1]], vgLine[[1]]$id == "pop")[,1:2], type = "l", col = "red", lwd = 1.42, main = "PPS % z")
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
plot(subset(vgLine[[1]], vgLine[[1]]$id == "pop")[,1:2], type = "l", col = "red", lwd = 1.42, main = "spread PWD")
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
plot(subset(vgLine[[1]], vgLine[[1]]$id == "pop")[,1:2], type = "l", col = "red", lwd = 1.42, main = "cluster PWD")
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

