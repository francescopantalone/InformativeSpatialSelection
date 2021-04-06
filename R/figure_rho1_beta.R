#' Generates the rho_{1} function of the paper, when \eqn{beta} is varying
#' 
#' This function 
#'  - generates the rho_{1} function of the paper, when \eqn{beta} is varying. 
#' @param light a boolean (FALSE by default) indicating whether we should produce a light (with less pixels) version of the file. 
#' @return an object of classes "gg" and "ggplot": the figure. 
#' @author FP DB
#' @seealso InformativeSpatialSelection::generateYonsquaregridU
figure_rho1_beta <- function(pdfexportfolder = NULL, pngconvert = FALSE){
  library(ggplot2)
  library(latex2exp)
  data(figure_rho1_data)
  yplot <- c(-1, -0.5, -0.1, 0.1, 0.5, 1)
  figureXX <- ggplot(data = subset(rho_1, .Gamma == 1 & y %in% yplot),
         aes(x = .Beta, y=rho, ymin = cimin, ymax = cimax, fill = "red", alpha = .3)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(trans = 'log10') +
    geom_ribbon() +
    ggplot2::xlab(TeX("$\\xi_{1}$"))+
    ggplot2::ylab(TeX("$\\rho_{1}$"))+
    facet_wrap(~factor(y)) +
    ggplot2::theme(legend.position = "none")

  if(!is.null(pdfexportfolder)){    
    figbasename=paste0("figure_rho1_beta.pdf")
    output=file.path(pdfexportfolder,figbasename)
    SweaveLst::graph2pdffile("print(figureXX +scale_shape_manual(values = c(21, 25),labels=c('PppPppPpp','bppbppbpp'),guide=guide_legend(order=1)))",
                             widthe = 8.5,
                             height = 2.88,
                             output = output,
                             modify = function(x){
                               #x[5]<-"\\usepackage{amsfonts,mathrsfs}\n"
                               x <- gsub("                      y0","$\\\\\\\\mathbf{y}$",x)
                               x <- gsub("PppPppPpp","$\\\\\\\\mathrm{Ppp}(\\\\\\\\mathbf{z})$",x)
                               x <- gsub("bppbppbpp","$\\\\\\\\mathrm{bpp}(\\\\\\\\mathbf{z},10)$",x)
                               x <- gsub("Sample","$~~~\\\\\\\\mathbf{d}$", x)
                               x <- gsub("z1","$\\\\\\\\mathbf{z}$", x)
                               x <- gsub("log(10)","\\\\log(10)$", x, fixed = T)
                               x <- gsub("log(z)=","$\\\\log(\\\\mathbf{z})=", x,fixed = T)
                               x <- gsub("log(z)","$\\\\log(\\\\mathbf{z})$", x, fixed = T)
                               x <- gsub("xi0",   "\\\\\\\\xi_0",x)
                               x <- gsub(", 25.67) {",   ", 32.05) {", x, fixed = T)
                               x <- gsub("xi1 y",  "\\\\\\\\xi_1\\\\\\\\mathbf{y}", x)
                               x <- gsub("xi2 e","\\\\\\\\xi_2\\\\\\\\mathbf{e}$", x)
                               x <- gsub("2.d. y","2.d. $\\\\mathbf{y}$", x, fixed = T)
                               x <- gsub(", scale=  1.10","",x)# to remove automatic rescaling the legend title would be better if done
                               rectangles<-grep("rectangle",x)
                               x[rectangles] <- gsub(paste0("\\path\\[", "fill="),"\\fill\\[", x[rectangles])
                               x
                             })
    fs::file_show(file.path(getwd(),"latex/fig/figure_rho1_beta.pdf"))
    if(pngconvert){
      temppngfile<-tempfile(pattern = ".png")
      pdftools::pdf_convert(output, "png", dpi = 1200, filenames = temppngfile)
      file.copy(from=temppngfile, to = file.path(pdfexportfolder, "figure_rho1_beta.png"), overwrite = TRUE)
    }
  }
  figureXX
}

