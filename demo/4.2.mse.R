mse4.2<-function(){
  # ------ library ------ #
  library(RandomFields)
  library(sampling)
  library(gstat)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  # --------------------- #
  # -- generate data -- #
  set.seed(9)
  xx <- 201 
  x1 <- x2<-generatex_i(xx)#This generates an arithmetic sequence  of length xx starting at 0 ending at 1.
  # data(figure1_data)
  figure1_data <- Simulate_figure1_data()
  y=figure1_data$y1[figure1_data$.n==1]
  b=.5
  b=.5;cc=.3
  z <- generateZconditionnallyongridUY2(length.out = xx, 
                                        y = y, 
                                        alpha = c(log(10)-b^2/2-cc^2/2,log(10)-b^2/2-cc^2/2),
                                        beta  = c(0,b), 
                                        gamma = c(sqrt(b^2+cc^2),cc),
                                        .scale.epsilon=10)
  XYZ<-bindXYZpiI(x1,x2,y=y,y0=NA,z)
  # ------------------- #
  # -- samples -- #
  N <- nrow(XYZ)
  n <- 100
  M <- 1000 # number of replications
  s1_M <- matrix(data = 0, nrow = M, ncol = n)
  s2_M <- matrix(data = 0, nrow = M, ncol = n)
  s3_M <- matrix(data = 0, nrow = M, ncol = n)
  set.seed(42)
  for(i in 1:M)
  {
    s1_M[i, ] <- srs(n, XYZ, replace = F)
  }
  pi2 <- inclusionprobabilities(a = XYZ$z1, n = n)
  set.seed(42)
  for(i in 1:M)
  {
    s2_M[i, ] <- sample(x = nrow(XYZ), size = n, prob = pi2)
    # s2_M[i, ] <- pps(n, XYZ$z1)
  }
  pi3 <- inclusionprobabilities(a = XYZ$z2, n = n)
  set.seed(42)
  for(i in 1:M)
  {
    s3_M[i, ] <- sample(x = nrow(XYZ), size = n, prob = pi3)
    # s3_M[i, ] <- pps(n, XYZ$z2)
  }
  # ------------- #
  XYZ_so <- XYZ
  coordinates(XYZ_so) <- ~ x1 + x2 # created a spatial objects for gstat package
  # -- variogram -- #
  vgm_pop <- variogram(y ~ 1, XYZ_so) # population variogram, since we use all the field
  vgm_s1 <- list()
  length(vgm_s1) <- M
  vgm_s2 <- list()
  length(vgm_s2) <- M
  vgm_s3 <- list()
  length(vgm_s3) <- M
  for(i in 1:M) # sample variograms
  {
    vgm_s1[[i]] <- variogram(y ~ 1, XYZ_so[s1_M[i,], ])
    vgm_s2[[i]] <- variogram(y ~ 1, XYZ_so[s2_M[i,], ])
    vgm_s3[[i]] <- variogram(y ~ 1, XYZ_so[s3_M[i,], ])
  }
  # --------------- #
  # -- fit variogram -- #
  fit_vgm_pop <- fit.variogram(vgm_pop, model = vgm(1, "Gau", 1))
  fit_vgm_s1 <- list()
  length(fit_vgm_s1) <- M
  fit_vgm_s2 <- list()
  length(fit_vgm_s2) <- M
  fit_vgm_s3 <- list()
  length(fit_vgm_s3) <- M
  for(i in 1:M)
  {
    fit_vgm_s1[[i]] <- fit.variogram(vgm_s1[[i]], vgm(1, "Gau", 1))
    fit_vgm_s2[[i]] <- fit.variogram(vgm_s2[[i]], vgm(1, "Gau", 1))
    fit_vgm_s3[[i]] <- fit.variogram(vgm_s3[[i]], vgm(1, "Gau", 1))
  }
  # ------------------- #
  mse_s1 <- matrix(0, nrow = M, ncol = length(variogramLine(fit_vgm_s1[[i]], maxdist = max(vgm_s1[[i]]$dist))$gamma))
  for(i in 1:M)
  {
    mse_s1[i, ] <- (variogramLine(fit_vgm_s1[[i]], maxdist = max(vgm_s1[[i]]$dist))$gamma
                    - variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist))$gamma) ^ 2
  }
  # mean(colSums(mse_s1) / M)
  # plot(variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist))$gamma, colSums(mse_s1) / M, type = "l")
  mse_s2 <- matrix(0, nrow = M, ncol = length(variogramLine(fit_vgm_s2[[i]], maxdist = max(vgm_s2[[i]]$dist))$gamma))
  for(i in 1:M)
  {
    mse_s2[i, ] <- (variogramLine(fit_vgm_s2[[i]], maxdist = max(vgm_s2[[i]]$dist))$gamma
                    - variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist))$gamma) ^ 2
  }
  # mean(colSums(mse_s2) / M)
  # plot(variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist))$gamma, colSums(mse_s2) / M, type = "l")
  mse_s3 <- matrix(0, nrow = M, ncol = length(variogramLine(fit_vgm_s3[[i]], maxdist = max(vgm_s3[[i]]$dist))$gamma))
  for(i in 1:M)
  {
    mse_s3[i, ] <- (variogramLine(fit_vgm_s3[[i]], maxdist = max(vgm_s3[[i]]$dist))$gamma
                    - variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist))$gamma) ^ 2
  }
  # mean(colSums(mse_s3) / M)
  # plot(variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist))$gamma, colSums(mse_s3) / M, type = "l")
  mse_s1_plot <- as.data.frame(cbind(colSums(mse_s1) / M,
                                     variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist))$gamma))
  colnames(mse_s1_plot) <- c("mse", "distance")
  mse_s2_plot <- as.data.frame(cbind(colSums(mse_s2) / M,
                                     variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist))$gamma))
  colnames(mse_s2_plot) <- c("mse", "distance")
  mse_s3_plot <- as.data.frame(cbind(colSums(mse_s3) / M,
                                     variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist))$gamma))
  colnames(mse_s3_plot) <- c("mse", "distance")
  ggplot() + geom_line(data = mse_s1_plot, aes(x = distance, y = mse, linetype = "srs")) +
    geom_line(data = mse_s2_plot, aes(x = distance, y = mse, linetype = "bpp(z1, 100)")) +
    geom_line(data = mse_s3_plot, aes(x = distance, y = mse, linetype = "bpp(z2, 100)")) + 
    scale_linetype_manual(values = c(3, 2, 1)) + theme(legend.position = "bottom") + 
    labs(linetype = "Sampling scheme")
  
  
  
  rmse_s2_plot <- as.data.frame(cbind(mse_s2_plot[-c(1:10), 1] / mse_s1_plot[-c(1:10), 1], mse_s2_plot[-c(1:10), 2]))
  colnames(rmse_s2_plot) <- c("rmse", "distance")
  rmse_s2_plot$Sampling <- "bpp(z1, 100)"
  rmse_s3_plot <- as.data.frame(cbind(mse_s3_plot[-c(1:10), 1] / mse_s1_plot[-c(1:10), 1], mse_s3_plot[-c(1:10), 2]))
  colnames(rmse_s3_plot) <- c("rmse", "distance")
  rmse_s3_plot$Sampling <- "bpp(z2, 100)"
  
  # rmse_plot <- rbind(rmse_s2_plot, rmse_s3_plot)
  # ggplot() +
  #   geom_line(data = rmse_plot, aes(x = distance, y = rmse)) + facet_wrap(~Sampling, scales = "free") +
  #   scale_linetype_manual(values = c(2, 1)) + theme(legend.position = "bottom") + 
  #   labs(linetype = "Sampling scheme")
  
  ggplot() +
    geom_line(data = rmse_s2_plot, aes(x = distance, y = rmse, linetype = "bpp(z1, 100)")) +
    geom_line(data = rmse_s3_plot, aes(x = distance, y = rmse, linetype = "bpp(z2, 100)")) + 
    scale_linetype_manual(values = c(2, 1)) + theme(legend.position = "bottom") + 
    labs(linetype = "Sampling scheme")
  
  
  
  
  
  
  
  # mse_s1_mean <- apply(mse_s1, 2, mean)
  # mse_s1_sd <- apply(mse_s1, 2, sd)
  # s1_low <- mse_s1_mean - qnorm(.975) * mse_s1_sd
  # s1_up <- mse_s1_mean + qnorm(.975) * mse_s1_sd
  # mse_s2_mean <- apply(mse_s2, 2, mean)
  # mse_s2_sd <- apply(mse_s2, 2, sd)
  # s2_low <- mse_s2_mean - qnorm(.975) * mse_s2_sd
  # s2_up <- mse_s2_mean + qnorm(.975) * mse_s2_sd
  # mse_s3_mean <- apply(mse_s3, 2, mean)
  # mse_s3_sd <- apply(mse_s3, 2, sd)
  # s3_low <- mse_s3_mean - qnorm(.975) * mse_s3_sd
  # s3_up <- mse_s3_mean + qnorm(.975) * mse_s3_sd
  # 
  # mse_s1_data <- cbind(mse_s1_plot, s1_low, s1_up)
  # mse_s2_data <- cbind(mse_s2_plot, s1_low, s1_up)
  # mse_s3_data <- cbind(mse_s3_plot, s1_low, s1_up)
  # ggplot() + geom_line(data = mse_s1_data, aes(x = gamma, y = mse, linetype = "srs")) +
  #   geom_ribbon(data = mse_s1_data, aes(x = gamma, y = mse, ymin = s1_low, ymax = s1_up, fill = "band", alpha = 0.3)) + 
  #   geom_line(data = mse_s2_plot, aes(x = gamma, y = mse, linetype = "pps")) +
  #   geom_ribbon(data = mse_s2_data, aes(x = gamma, y = mse, ymin = s2_low, ymax = s2_up, fill = "band", alpha = 0.3)) + 
  #   geom_line(data = mse_s3_plot, aes(x = gamma, y = mse, linetype = "pps2")) + 
  #   geom_ribbon(data = mse_s3_data, aes(x = gamma, y = mse, ymin = s3_low, ymax = s3_up, fill = "band", alpha = 0.3)) + 
  #   scale_linetype_manual(values = c(3, 2, 1))
  # 
  # 
  # 
  # 
  # mse_s1_mean <- apply(mse_s1, 2, mean)
  # mse_s1_sd <- apply(mse_s1, 2, sd)
  # s1_low <- mse_s1_mean - qnorm(.975) * mse_s1_sd
  # s1_up <- mse_s1_mean + qnorm(.975) * mse_s1_sd
  
  #' X$ymin=X$mean-qnorm(.975)*X$sd
  #' X$ymax=X$mean+qnorm(.975)*X$sd
  
  # plot(variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist))$gamma, colSums(mse_s1) / M, type = "l", lty = 1)
  # lines(variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist))$gamma, colSums(mse_s2) / M, type = "l", lty = 2)
  # lines(variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist))$gamma, colSums(mse_s3) / M, type = "l", lty = 3)
  
  # mean(colSums(mse_s2) / M) / mean(colSums(mse_s1) / M)
  # mean(colSums(mse_s3) / M) / mean(colSums(mse_s1) / M)
}
