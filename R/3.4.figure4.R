#' Generates figure 4 of the paper
#' 
#' @param pdfexportfolder a character string indicating a path where to save the figure in pdf format. If NULL (default), no pdf export is created. 
#' @param pngconvert a boolean (FALSE by default) indicating whether we should convet the pdf file into png. If TRUE, a png will be created in the pdfexportfolder. If no pdf export folder is given, no png will be created
#' @return an object of classes "gg" and "ggplot": the first figure of the paper. 
#' @author FP
#' @seealso InformativeSpatialSelection::generateYonsquaregridU
#' @example 
#' figure4()
figure4 <- function(pdfexportfolder = NULL, pngconvert = FALSE){
  library(RandomFields)
  library(sampling)
  library(gstat)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  data(figure4_data)
  attach(figure4_data)
  # ------------ plot ------------ #
  # -- samples -- #
  # set.seed(42)
  # s1 <- srs(n, XYZ, replace = F)
  # pi2 <- inclusionprobabilities(a = XYZ$z1, n = n)
  # s2 <- sample(x = nrow(XYZ), size = n, prob = pi2)
  # pi3 <- inclusionprobabilities(a = XYZ$z2, n = n)
  # s3 <- sample(x = nrow(XYZ), size = n, prob = pi3)
  # XYZ_s1 <- XYZ[, c(1, 2, 3)]
  # XYZ_s1 <- as.data.frame(XYZ_s1)
  # XYZ_s1$I <- rep(0, nrow(XYZ_s1))
  # XYZ_s1[s1, 4] <- 1
  # XYZ_s1$Sampling <- rep("srs", nrow(XYZ_s1))
  # colnames(XYZ_s1) <- c("x1", "x2", "y", "I", "Sampling")
  # XYZ_s2 <- XYZ[, c(1, 2, 3)]
  # XYZ_s2 <- as.data.frame(XYZ_s2)
  # XYZ_s2$I <- rep(0, nrow(XYZ_s2))
  # XYZ_s2[s2, 4] <- 1
  # XYZ_s2$Sampling <- rep("bpp(z1, 100)", nrow(XYZ_s2))
  # colnames(XYZ_s2) <- c("x1", "x2", "y", "I", "Sampling")
  # XYZ_s3 <- XYZ[, c(1, 2, 3)]
  # XYZ_s3 <- as.data.frame(XYZ_s2)
  # XYZ_s3$I <- rep(0, nrow(XYZ_s2))
  # XYZ_s3[s3, 4] <- 1
  # XYZ_s3$Sampling <- rep("bpp(z2, 100)", nrow(XYZ_s2))
  # colnames(XYZ_s2) <- c("x1", "x2", "y", "I", "Sampling")
  # XYZ_map <- rbind(XYZ_s1, XYZ_s2, XYZ_s3)
  # XYZ_map$Sampling <- factor(XYZ_map$Sampling, levels = c("srs", "bpp(z1, 100)", "bpp(z2, 100)"))
  # samples_plot <- ggplot(data = XYZ_map, aes(x = x1, y = x2, fill = y)) + geom_tile() +
  #   ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::scale_x_continuous(expand=c(0,0)) +
  #   ggplot2::scale_y_continuous(expand=c(0,0)) + geom_point(data = XYZ_map[which(XYZ_map[, 4] == 1),], aes(x = x1, y = x2)) +
  #   facet_wrap(~Sampling, nrow = 1) +
  #   ggplot2::theme(
  #     axis.text.x = ggplot2::element_blank(),
  #     axis.ticks.x = ggplot2::element_blank(),
  #     axis.text.y = ggplot2::element_blank(),
  #     axis.ticks.y = ggplot2::element_blank(),
  #     panel.border = ggplot2::element_blank(),
  #     panel.grid.major = ggplot2::element_blank(), 
  #     panel.grid.minor = ggplot2::element_blank(),
  #     panel.background = ggplot2::element_blank()) +
  #   ggplot2::scale_fill_gradient(low = "white", high = "black")
  # ------------- #
  # -- densities -- #
  den_y_pop <- as.data.frame(cbind(density(XYZ$y)$x, density(XYZ$y)$y))
  colnames(den_y_pop) <- c("x", "y")
  den_y_s1_mean <- 0
  for(i in 1:M)
  {
    den_y_s1_mean <- den_y_s1_mean + cbind(density(XYZ[s1_M[i, ], "y"])$x, density(XYZ[s1_M[i, ], "y"])$y)
  }
  den_y_s1_mean <- as.data.frame(den_y_s1_mean / M)
  den_y_s1_mean$Sampling <- rep("bpp(1, 100)", nrow(den_y_s1_mean))
  colnames(den_y_s1_mean) <- c("x", "y", "Sampling")
  den_y_s2_mean <- 0
  for(i in 1:M)
  {
    den_y_s2_mean <- den_y_s2_mean + cbind(density(XYZ[s2_M[i, ], "y"])$x, density(XYZ[s2_M[i, ], "y"])$y)
  }
  den_y_s2_mean <- as.data.frame(den_y_s2_mean / M)
  den_y_s2_mean$Sampling <- rep("bpp(z1, 100)", nrow(den_y_s1_mean))
  colnames(den_y_s2_mean) <- c("x", "y", "Sampling")
  den_y_s3_mean <- 0
  for(i in 1:M)
  {
    den_y_s3_mean <- den_y_s3_mean + cbind(density(XYZ[s3_M[i, ], "y"])$x, density(XYZ[s3_M[i, ], "y"])$y)
  }
  den_y_s3_mean <- as.data.frame(den_y_s3_mean / M)
  den_y_s3_mean$Sampling <- rep("bpp(z2, 100)", nrow(den_y_s1_mean))
  colnames(den_y_s3_mean) <- c("x", "y", "Sampling")
  
  densities <- rbind(den_y_s1_mean, den_y_s2_mean, den_y_s3_mean)
  densities$Sampling <- factor(densities$Sampling, levels = c("bpp(1, 100)", "bpp(z1, 100)", "bpp(z2, 100)"))
  densities_plot <- ggplot(data = densities, aes(x = x, y = y)) + geom_line(aes(linetype = "mean")) +
    facet_wrap(~Sampling, nrow = 1) +
    geom_line(data = as.data.frame(den_y_pop), aes(x = x, y = y, linetype = "pop")) +
    scale_linetype_manual(values=c(2, 1)) + labs(x = "y", y = "density", linetype = "")
  # --------------- #
  # -- variograms -- #
  var_y_pop <- as.data.frame(variogramLine(fit_vgm_pop, maxdist = max(vgm_pop$dist)))
  colnames(var_y_pop) <- c("dist", "gamma")
  
  var_y_s1_mean <- 0
  for(i in 1:M)
  {
    var_y_s1_mean <- var_y_s1_mean + variogramLine(fit_vgm_s1[[i]], maxdist = max(vgm_s1[[i]]$dist))
  }
  var_y_s1_mean <- as.data.frame(var_y_s1_mean / M)
  var_y_s1_mean$Sampling <- rep("bpp(1, 100)", nrow(var_y_s1_mean))
  var_y_s2_mean <- 0
  for(i in 1:M)
  {
    var_y_s2_mean <- var_y_s2_mean + variogramLine(fit_vgm_s2[[i]], maxdist = max(vgm_s2[[i]]$dist))
  }
  var_y_s2_mean <- as.data.frame(var_y_s2_mean / M)
  var_y_s2_mean$Sampling <- rep("bpp(z1, 100)", nrow(var_y_s2_mean))
  var_y_s3_mean <- 0
  for(i in 1:M)
  {
    var_y_s3_mean <- var_y_s3_mean + variogramLine(fit_vgm_s3[[i]], maxdist = max(vgm_s3[[i]]$dist))
  }
  var_y_s3_mean <- as.data.frame(var_y_s3_mean / M)
  var_y_s3_mean$Sampling <- rep("bpp(z2, 100)", nrow(var_y_s3_mean))
  variograms <- rbind(var_y_s1_mean, var_y_s2_mean, var_y_s3_mean)
  variograms$Sampling <- factor(variograms$Sampling, levels = c("bpp(1, 100)", "bpp(z1, 100)", "bpp(z2, 100)"))
  
  variograms_plot <- ggplot(data = variograms, aes(x = dist, y = gamma)) + geom_line(aes(linetype = "mean")) +
    facet_wrap(~Sampling) +
    geom_line(data = var_y_pop, aes(x = dist, y = gamma, linetype = "pop")) + 
    scale_linetype_manual(values=c(2, 1)) + labs(x = "distance", y = "variogram", linetype = "")
  # ---------------- #
  # figure4 <- grid.arrange(samples_plot, densities_plot, variograms_plot, nrow = 3)
  figure4 <- grid.arrange(densities_plot, variograms_plot, nrow = 2)
  
  if(!is.null(pdfexportfolder)){
    figbasename="figure4.pdf"
    output=file.path(pdfexportfolder,figbasename)
    SweaveLst::graph2pdffile("print(figure4)",widthe = 6.5,height=12,output = output,
                             modify = function(x){
                               #x[5]<-"\\usepackage{amsfonts,mathrsfs}\n"
                               x<-gsub("signal","$\\\\\\\\mathbf{y}$",x)
                               x<-gsub("gamma","$\\\\\\\\gamma$",x)
                               x<-gsub(", scale=  1.10","",x)# to remove automatic rescaling the legend title would be better if done
                               rectangles<-grep("rectangle",x)
                               x[rectangles]<-gsub(paste0("\\path\\[","fill="),"\\fill\\[",x[rectangles])
                               x
                             })
    if(pngconvert){    
      temppngfile<-tempfile(pattern = ".png")
      pdftools::pdf_convert(output,"png",dpi = 1200,filenames = temppngfile)
      file.copy(from=temppngfile,to=pdfexportfolder,overwrite = TRUE)
    }
  }
  figure4
}