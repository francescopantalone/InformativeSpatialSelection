# ------------ library ------------ #
library(sampling)
library(samplingbook)
library(sp)
library(lattice)
library(gstat)
#library(Spbsampling)
library(ggplot2)
library(dplyr)
library(RandomFields)
library(pps)
# --------------------------------- #
# ------------ simulation gaussian random field ------------ #
set.seed(42)
prc_y <- RFsimulate(model = RMgauss(var = 5, scale = 0.1) + RMtrend(mean = 10), x = seq(0, 1, length = 100), y = seq(0, 1, length = 100))
prc_y2 <- RFsimulate(model = RMgauss(var = 5, scale = 0.1) + RMtrend(mean = 10), x = seq(0, 1, length = 100), y = seq(0, 1, length = 100))
alpha <- 0.5
prc_z <- prc_y2
prc_z@data <- alpha * prc_y@data + (1 - alpha) * prc_y2@data

x <- seq(0, 1, length = 100)
y <- seq(0, 1, length = 100)
x1 <- rep(x, length(y))
x2 <- rep(y, each = length(x))
prc <- data.frame(x1, x2, y = prc_y@data[[1]])
prc_grid <- expand.grid(seq(0, 1, length = 20), seq(0, 1, length = 20)) # grid
colnames(prc_grid) <- c("x1", "x2")
# prc %>% as.data.frame %>%
#   ggplot(aes(x1, x2)) + geom_point(aes(size = y), color = "blue", alpha = 3/4) +
#   ggtitle("Random field") + coord_equal() + theme_bw()
# ---------------------------------------------------------- #
# ------------ draw samples ------------ #
npop <- nrow(prc)
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
  toto<-runif(length(c(prc_z$variable1)))
  oo<-order(toto)
  pps_sample[i, ] <- oo[ppss(prc_z$variable1[oo], nsample)]
}
# ----------------- #
# -------------------------------------- #
# ------------ kriging ------------ #
set.seed(42)
# -- theoretical variogram -- #
prc_y_tv <- RFvariogram(RMgauss(var = 5, scale = 0.1) + RMtrend(mean = 10), distances = seq(0, 1, length.out = 100))
# --------------------------- #
coordinates(prc) <- ~ x1 + x2
# -- empirical variogram -- #
#vgm_pop <- variogram(y ~ 1, prc)
vgm_srs <- list()
length(vgm_srs) <- niter
vgm_pps <- list()
length(vgm_pps) <- niter
for(i in 1:niter)
{
  vgm_srs[[i]] <- variogram(y ~ 1, prc[srs_sample[i,], ])
  vgm_pps[[i]] <- variogram(y ~ 1, prc[pps_sample[i,], ])
}
# x<-seq(0,1,length.out=100)
# fit_vgm_srs[[1]] <- fit.variogram.reml(data=prc[srs_sample[1],],y~1,model =  vgm(model= "Gau",range =  300))
# y<-variogramLine(fit_vgm_srs[[1]],dist_vector=x)
# y<-cbind(y,prc_y_tv)
# ggplot(data=reshape2::melt(y,id.vars="dist"),aes(x=dist,y=value,color=variable))+
#   geom_line()+geom_point(data=vgm_srs[[1]],mapping=aes(x=dist,y=gamma,color="red"))
# ------------------------- #
# -- fit variogram -- #
#fit_vgm_pop <- fit.variogram.reml(y ~ 1, prc, model = vgm(1, "Gau", 300, 1)) # too computationally demanding when we use big grid
fit_vgm_srs <- list()
length(fit_vgm_srs) <- niter
fit_vgm_pps <- list()
length(fit_vgm_pps) <- niter
for(i in 1:niter)
{
  # fit_vgm_srs[[i]] <- fit.variogram(vgm_srs[[i]], vgm("Gau", 300), fit.sills = TRUE, fit.ranges = TRUE)
  # fit_vgm_pps[[i]] <- fit.variogram(vgm_pps[[i]], vgm("Gau", 300), fit.sills = TRUE, fit.ranges = TRUE)
  fit_vgm_srs[[i]] <- fit.variogram(vgm_srs[[i]], vgm("Gau")) # linear model has singular covariance matrix
  fit_vgm_pps[[i]] <- fit.variogram(vgm_pps[[i]], vgm("Gau"))
}
fit_vgm_srs1 <- list()
length(fit_vgm_srs1) <- niter
fit_vgm_pps1 <- list()
length(fit_vgm_pps1) <- niter
z_srs <- 0
z_pps <- 0
for(i in 1:niter)
{
  if(fit_vgm_srs[[i]]$range[2] >= 0)
  {
    z_srs <- z_srs + 1
    fit_vgm_srs1[[z_srs]] <- fit_vgm_srs[[i]]
  }
  if(fit_vgm_pps[[i]]$range[2] >= 0)
  {
    z_pps <- z_pps + 1 
    fit_vgm_pps1[[z_pps]] <- fit_vgm_pps[[i]]
  }
}
# ------------------- #
# -- prediction -- #
coordinates(prc_grid) <- ~ x1 + x2
pred_srs <- list()
length(pred_srs) <- niter
pred_pps <- list()
length(pred_pps) <- niter
for(i in 2:niter)
{
  pred_srs[[i]] <- krige(y ~ 1, prc, prc_grid, model = fit_vgm_srs1[[i]])
  pred_pps[[i]] <- krige(y ~ 1, prc, prc_grid, model = fit_vgm_pps1[[i]])
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
# vgm_pop1 <- vgm_pop
# vgm_pop1$id <- "pop"
vgm_srs1 <- list()
length(vgm_srs1) <- niter
vgm_srs1 <- vgm_srs
vgm_pps1 <- list()
length(vgm_pps1) <- niter
vgm_pps1 <- vgm_pps
vgLine <- list()
length(vgLine) <- niter
for(i in 1:min(z_srs, z_pps))
{
  vgm_srs1[[i]]$id <- "srs"
  vgm_pps1[[i]]$id <- "pps"
  vgLine[[i]] <- rbind( 
    cbind(variogramLine(fit_vgm_srs1[[i]], maxdist = max(vgm_srs1[[i]]$dist)), id = 
            "srs"),
    cbind(variogramLine(fit_vgm_pps1[[i]], maxdist = max(vgm_pps1[[i]]$dist)), id =
            "pps")
  ) 
}


par(mfrow = c(1, 2), cex = 0.50)

srs_mean <- matrix(0, nrow(subset(vgLine[[1]], vgLine[[1]]$id == "srs")[,1:2]), ncol(subset(vgLine[[1]], vgLine[[1]]$id == "srs")[,1:2]))
plot(x = seq(0, 1, length = length(prc_y_tv)), y = prc_y_tv, ylim = c(0, 7), type = "l", col = "red", lwd = 1.42, main = "SRS")
srs_mean <- 0
for(i in 1:min(z_srs, z_pps))
{
  lines(subset(vgLine[[i]], vgLine[[i]]$id == "srs")[, 1:2], type = "l", lwd = 0.24)
  srs_mean <- srs_mean + subset(vgLine[[i]], vgLine[[i]]$id == "srs")[,1:2]
}
srs_mean <- srs_mean / niter
lines(srs_mean, col = "blue", lwd = 1.42)
legend("bottomright", legend = c("theoretical", "repl mean", "repl"), col = c("red", "blue", "black"), lty = 1)

pps_mean <- matrix(0, nrow(subset(vgLine[[1]], vgLine[[1]]$id == "pps")[,1:2]), ncol(subset(vgLine[[1]], vgLine[[1]]$id == "pps")[,1:2]))
plot(x = seq(0, 1, length = length(prc_y_tv)), y = prc_y_tv, ylim = c(0, 7), type = "l", col = "red", lwd = 1.42, main = "PPS % z")
pps_mean <- 0
for(i in 1:min(z_srs, z_pps))
{
  lines(subset(vgLine[[i]], vgLine[[i]]$id == "pps")[,1:2], type = "l", lwd = 0.24)
  pps_mean <- pps_mean + subset(vgLine[[i]], vgLine[[i]]$id == "pps")[,1:2]
}
pps_mean <- pps_mean / niter
lines(pps_mean, col = "blue", lwd = 1.42, type = "l")
legend("bottomright", legend = c("pop", "repl mean", "repl"), col = c("red", "blue", "black"), lty = 1)

par(mfrow = c(1,1), cex = 1)
# ------------------------------ #
