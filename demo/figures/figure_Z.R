#' Generates the figure Z of the paper
#' 
#' @return an object of classes "gg" and "ggplot": the first figure of the paper. 
#' @author FP
#' @example 
#' figureZ()
figureZ<-function(){
# -- library -- #
library(RandomFields)
library(sampling)
library(ggplot2)
library(gridExtra)
# ------------- #
# -- grid and processes Zs -- #
grid_size <- 30
set.seed(42)
Y <- generateYonsquaregridU(length.out = grid_size)
Z <- generateZconditionnallyongridUY(length.out = grid_size, y = Y)
XYZ <- bindXYZpiI(x1 = generatex_i(grid_size), x2 = generatex_i(grid_size), y = Y, z = Z)
n <- 200 # sample size
# --------------------------- #
z1_pps <- inclusionprobabilities(XYZ$z.z1, n)
set.seed(42)
s_z1_pps <- which(UPmaxentropy(z1_pps) == 1)
plot_z1 <- ggplot(data = XYZ, aes(x = x1, y = x2, fill = z.z1)) + geom_tile() +
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
  ggplot2::scale_fill_gradient(low = "white", high = "black")  +
  ggplot2::labs(fill = "Z1") + geom_point(data = XYZ[s_z1_pps, ], aes(x = x1, y = x2))
z2_pps <- inclusionprobabilities(XYZ$z.z2, n)
set.seed(42)
s_z2_pps <- which(UPmaxentropy(z2_pps) == 1)
plot_z2 <- ggplot(data = XYZ, aes(x = x1, y = x2, fill = z.z2)) + geom_tile() +
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
  ggplot2::scale_fill_gradient(low = "white", high = "black")  +
  ggplot2::labs(fill = "Z2") + geom_point(data = XYZ[s_z2_pps, ], aes(x = x1, y = x2))
grid.arrange(plot_z1, plot_z2, nrow = 1, ncol = 2)
}
