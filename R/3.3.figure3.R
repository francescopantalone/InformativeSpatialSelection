#' Generates figure 3 of the paper
#' 
#' @param tikzexportfolder a character string indicating a path where to save the figure in tikz format. If NULL (default), no tikz export is done. 
#' @param light a boolean (FALSE by default) indicating whether we should produce a light (with less pixels) version of the file. 
#' @param pngconvert a boolean (FALSE by default) indicating whether we should convet the pdf file into png. If yes, a png will be created  
#' @return an object of classes "gg" and "ggplot": the first figure of the paper. 
#' @author FP
#' @seealso InformativeSpatialSelection::generateYonsquaregridU
#' @examples
#' attach(figure3())
#' grid.arrange(p, ydensity,
#'              zdensity, blankPlot , 
#'              ncol = 2, nrow = 2, 
#'              widths = c(4,1.4), heights = c(4, 1.4))

figure3 <- function(pdfexportfolder = NULL, light = FALSE, pngconvert = FALSE){
  library(RandomFields)
  library(sampling)
  library(gstat)
  library(ggplot2)
  library(gridExtra)
  grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right"), which_legend = 1) {
    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[which_legend]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position="none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(position,
                       "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                              legend,
                                              ncol = 1,
                                              heights = grid::unit.c(unit(1, "npc") - lheight, lheight)),
                       "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                             legend,
                                             ncol = 2,
                                             widths = grid::unit.c(unit(1, "npc") - lwidth, lwidth)))
    
    grid::grid.newpage()
    grid::grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
  }
  # -- generate data -- #
  # grid_size <- 50
  # set.seed(42)
  # Y <- generateYonsquaregridU(length.out = grid_size, mu = 20, .sigma = 5)
  # Z <- generateZconditionnallyongridUY(length.out = grid_size, y = Y, alpha = 20, beta = c(0, 0, 0.8), .sigma = 5)
  # XYZ <- bindXYZpiI(x1 = generatex_i(grid_size), x2 = generatex_i(grid_size), y = Y, z = Z)
  # ------------------- #
  data(figure3_data)
  N <- nrow(figure3_data)
  n <- round(N * 0.1)
  set.seed(42)
  # s1 <- srs(n, XYZ, replace = F)
  #  s2 <- pps(n, XYZ$z.z2, XYZ)
  #  s3 <- pps(n, XYZ$z.z3, XYZ)
  XYZ2 <- reshape2::melt(figure3_data, id = c("x1","x2","y"))
  XYZ2 <- XYZ2[is.element(XYZ2$variable, c("z.z1", "z.z2", "z.z3")), ]
  XYZ2$variable <- droplevels(XYZ2$variable)
  levels(XYZ2$variable) <- c("$Z\\perp Y;\\beta=\\gamma=0;$", "$Z\\perp Y$; $\\beta=0$", "$Z \\not\\perp Y$; $\\beta\\neq 0$")
  p <- ggplot(data = XYZ2, aes(x = value, y = y)) + 
    #geom_density_2d_filled() + 
    #geom_line(alpha=.2,size=.2)+
    geom_point(alpha=.2,size=.2)+
    ggplot2::scale_x_continuous(trans="log10",
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme(legend.position = "none",
                   panel.background = element_rect(fill= "gray94",colour="gray94",size = 0.5, linetype = "solid"),
                   panel.border = ggplot2::element_blank(), 
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()#, 
                   #panel.background = ggplot2::element_blank()
                   ) +
    scale_fill_grey(start = 0.8, end = 0.2) + 
    facet_grid(~variable) +
    xlab("") + ylab("$\\mathbf{y}$")
  
  p_temp <- ggplot(data = XYZ2, aes(x = value, y = y)) + geom_density_2d_filled() + 
    ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme(panel.border = ggplot2::element_blank(), 
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(), 
                   panel.background = ggplot2::element_blank()) +
    scale_fill_grey(start = 0.8, end = 0.2) + facet_grid(~variable) +
    xlab("") + ylab("$\\mathbf{y}$")
  
  l <- lemon::g_legend(p_temp) 
  
  # Marginal density plot of x (top panel)
  ydensity <- ggplot(XYZ2, aes(x = y)) + geom_density(alpha = .5) + 
    coord_flip() + xlab("") + ylab("") +
    theme(axis.text = element_blank(), panel.background = element_blank(),
          axis.ticks = element_blank(), axis.line = element_blank())
  # Marginal density plot of y (right panel)
  zdensity <- ggplot(XYZ2, aes(x = value)) + geom_density(alpha = .5) + 
    facet_grid(~variable) + theme(strip.background = element_blank(), 
                                  strip.text.x = element_blank(),
                                  axis.text = element_blank(),
                                  panel.background = element_blank(),
                                  axis.ticks = element_blank(),
                                  axis.line = element_blank()) +
    scale_y_reverse() + xlab("") + ylab("")+ scale_x_continuous(trans = "log10")
  
  
  blankPlot <- ggplot() + geom_blank(aes(1, 1)) +
    theme(
      plot.background = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(), 
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )
  grid.arrange(p, ydensity,
               zdensity, blankPlot , 
               ncol = 2, nrow = 2, widths = c(4,1.4), heights = c(4, 1.4))
  
  
  if(!is.null(pdfexportfolder)){
    figbasename="figure3.pdf"
    output=file.path(pdfexportfolder,figbasename)
    SweaveLst::graph2pdffile("grid.arrange(p, ydensity,
               zdensity, blankPlot , 
               ncol = 2, nrow = 2, widths = c(4,1.4), heights = c(4, 1.4))",
                             widthe = 8.2,height=5,output = output,
                             modify = function(x){
                               #x[5]<-"\\usepackage{amsfonts,mathrsfs}\n"
                               x<-gsub("signal","$\\\\\\\\mathbf{y}$",x)
                               x<-gsub(", scale=  1.10","",x)# to remove automatic rescaling the legend title would be better if done
                               rectangles<-grep("rectangle",x)
                               x[rectangles]<-gsub(paste0("\\path\\[","fill="),"\\fill\\[",x[rectangles])
                               x
                             })
    #fs::file_show(output)
    if(pngconvert){    
      temppngfile<-tempfile(pattern = ".png")
      pdftools::pdf_convert(output,"png",dpi = 1200,filenames = temppngfile)
      file.copy(from=temppngfile,to=file.path(pdfexportfolder,"figure3.png"),overwrite = TRUE)
    }}
  list(p=p, ydensity=ydensity,zdensity=zdensity, blankPlot=blankPlot )
}

