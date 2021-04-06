#' Generates the rho_{12} function of the paper
#' 
#' This function 
#'  - generates the rho_{12} function of the paper. 
#' @param light a boolean (FALSE by default) indicating whether we should produce a light (with less pixels) version of the file. 
#' @return an object of classes "gg" and "ggplot": the figure. 
#' @author FP DB
#' @seealso InformativeSpatialSelection::generateYonsquaregridU
figure_rho1_beta <- function(pdfexportfolder = NULL, pngconvert = FALSE){
  library(ggplot2)
  library(latex2exp)
  data(figure_rho12_data)
  attach(figure_rho12_data)
  yplot <- c(-1, -0.5, -0.1, 0.1, 0.5, 1)
  
  grid_plot <- ggplot2::ggplot(scenario, ggplot2::aes(x1,x2)) +
    ggplot2::geom_tile(ggplot2::aes(fill = y1),size=2) +
    ggplot2::geom_point(ggplot2::aes(x = x1[ell1], y = x2[ell1]),
                        fill = "white", color = "black", size = 1.5, stroke = .9, shape = 25) +
    ggplot2::geom_point(ggplot2::aes(x = x1[ell2], y = x2[ell2]),
                        fill = "white", color = "black", size = 1.5, stroke = .9, shape = 25) +
    ggplot2::xlab("")+
    ggplot2::ylab("")+
    ggplot2::scale_x_continuous(expand=c(0,0))+
    ggplot2::scale_y_continuous(expand=c(0,0))+
    ggplot2::theme(
      axis.ticks.x=ggplot2::element_blank(),
      axis.ticks.y=ggplot2::element_blank(),
      strip.text.x = element_text(angle = 0, hjust = 0),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank())+
    ggplot2::scale_fill_gradient(low = "white",
                                 high = "black"#,guide=guide_colorbar(frame.colour=c("black"),frame.linewidth=2)
    )  +
    ggplot2::facet_wrap(~scenario, nrow = 1)+
    ggplot2::labs(fill = "signal") + ggplot2::theme(legend.position="bottom")

  contour_plot <- ggplot(data = rho_12, aes(x = y1, y = y2, z = rho)) + geom_contour_filled() +
    ggplot2::xlab(TeX("$y_{1}")) +
    ggplot2::ylab(TeX("$y_{2}")) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()) +
    facet_wrap(~round(d, 3), nrow = 1) +
    ggplot2::labs(fill = TeX("$\\rho_{1,2}$")) + ggplot2::theme(legend.position="bottom")
  
  figureXX <- grid.arrange(grid_plot, contour_plot, nrow = 2)
  
  if(!is.null(pdfexportfolder)){    
    figbasename=paste0("figure_rho12.pdf")
    output=file.path(pdfexportfolder,figbasename)
    SweaveLst::graph2pdffile("print(figureXX +scale_shape_manual(values = c(21, 25),labels=c('PppPppPpp','bppbppbpp'),guide=guide_legend(order=1)))",
                             widthe = 8.5,
                             height = 2.88,
                             output = output,
                             modify = function(x){
                               #x[5]<-"\\usepackage{amsfonts,mathrsfs}\n"
                               x <- gsub("                      y0","$\\\\\\\\mathbf{y}$", x)
                               x <- gsub("PppPppPpp","$\\\\\\\\mathrm{Ppp}(\\\\\\\\mathbf{z})$", x)
                               x <- gsub("bppbppbpp","$\\\\\\\\mathrm{bpp}(\\\\\\\\mathbf{z},10)$", x)
                               x <- gsub("Sample","$~~~\\\\\\\\mathbf{d}$", x)
                               x <- gsub("z1","$\\\\\\\\mathbf{z}$", x)
                               x <- gsub("log(10)","\\\\log(10)$", x, fixed = T)
                               x <- gsub("log(z)=","$\\\\log(\\\\mathbf{z})=", x, fixed = T)
                               x <- gsub("log(z)","$\\\\log(\\\\mathbf{z})$", x, fixed = T)
                               x <- gsub("xi0",   "\\\\\\\\xi_0", x)
                               x <- gsub(", 25.67) {",   ", 32.05) {", x, fixed = T)
                               x <- gsub("xi1 y",  "\\\\\\\\xi_1\\\\\\\\mathbf{y}", x)
                               x <- gsub("xi2 e", "\\\\\\\\xi_2\\\\\\\\mathbf{e}$", x)
                               x <- gsub("2.d. y", "2.d. $\\\\mathbf{y}$", x, fixed = T)
                               x <- gsub(", scale=  1.10","",x)# to remove automatic rescaling the legend title would be better if done
                               rectangles <- grep("rectangle", x)
                               x[rectangles] <- gsub(paste0("\\path\\[","fill="),"\\fill\\[", x[rectangles])
                               x
                             })
    fs::file_show(file.path(getwd(),"latex/fig/figure_rho12.pdf"))
    if(pngconvert){
      temppngfile<-tempfile(pattern = ".png")
      pdftools::pdf_convert(output, "png", dpi = 1200,filenames = temppngfile)
      file.copy(from = temppngfile, to = file.path(pdfexportfolder, "figure_rho12.png"), overwrite = TRUE)
    }
  }
  figureXX
}

