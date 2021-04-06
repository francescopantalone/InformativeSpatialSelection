#' Generates the first figure of the paper
#' 
#' This function 
#'  - generates a regular grid of size 201 on the \eqn{[0,1]^2} square.
#'  - uses \code{InformativeSpatialSelection::generateYonsquaregridU} to generate the values of \eqn{Y}.
#'  - replicates the process 3 times 
#'  - generates a heatmap
#' @param pdfexportfolder a character string indicating a path where to save the figure in pdf format. If NULL (default), no tikz export is done. 
#' @param pngconvert a boolean indicating if a conversion to png. 
#' @param figure1_data a dataframe (default: \code{get(data(figure1_data,package="InformativeSpatialSelection"))})
#' @return an object of classes "gg" and "ggplot": the first figure of the paper. 
#' @author DB
#' @seealso InformativeSpatialSelection::generateYonsquaregridU
#' @examples
#' figure1()
#' figure1(pdfexportfolder=tempdir())
#' fs::file_show(file.path(tempdir(),"figure1standalone.pdf"))
figure1<-function(pdfexportfolder=NULL,pngconvert=FALSE,
                  figure1_data=NULL){
  if(is.null(figure1_data)){
    figure1_data=get(data(figure1_data,package="InformativeSpatialSelection"))}
  # xx<-201 
  # x1<-generatex_i(xx)#This generates an arithmetic sequence  of length xx starting at 0 ending at 1.
  # set.seed(42)
  # UY<-plyr::rdply(.n=3,function(){bindXYZpiI(x1,x1,y=generateYonsquaregridU(x1,x1))})
  figure1_data$title<-factor(paste0("1.",letters[figure1_data$.n]," $\\theta_1=",figure1_data$theta1,"$, ",
                                   "$\\theta_2=",figure1_data$theta2,"$"))
  figureXX<-ggplot2::ggplot(figure1_data, ggplot2::aes(x1,x2)) + 
    ggplot2::geom_tile(ggplot2::aes(fill = y1),size=2) + 
    ggplot2::xlab("")+
    ggplot2::ylab("")+
    ggplot2::scale_x_continuous(expand=c(0,0))+
    ggplot2::scale_y_continuous(expand=c(0,0))+
    ggplot2::theme(
      axis.text.x=ggplot2::element_blank(),
      axis.ticks.x=ggplot2::element_blank(),
      axis.text.y=ggplot2::element_blank(),
      axis.ticks.y=ggplot2::element_blank(),
      strip.text.x = element_text(angle = 0, hjust = 0),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(), 
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank())+
    ggplot2::scale_fill_gradient(low = "white",     
                                 high = "black"#,guide=guide_colorbar(frame.colour=c("black"),frame.linewidth=2)
                                 )  +
    ggplot2::facet_wrap(~title)+ 
    ggplot2::labs(fill = "signal")+ 
    ggplot2::theme(legend.position="bottom",
                   legend.justification="right",
                   legend.margin=ggplot2::margin(0,0,0,0),
                   legend.box.margin=ggplot2::margin(-20,#top
                                                     18,#right
                                                     10,#bottom
                                                     10))#left
  if(!is.null(pdfexportfolder)){
    figbasename=paste0("figure1.pdf")
    output=file.path(pdfexportfolder,figbasename)
    SweaveLst::graph2pdffile("print(figureXX)",widthe = 8.5,height=3.7,output = output,
                             modify = function(x){
                               #x[5]<-"\\usepackage{amsfonts,mathrsfs}\n"
                               x<-gsub("610.68,160.89","612,180",x)
                               x<-gsub("signal","$\\\\\\\\mathbf{y}$",x)
                               x<-gsub(", scale=  1.10","",x)# to remove automatic rescaling the legend title would be better if done
                               rectangles<-grep("rectangle",x)
                               x[rectangles]<-gsub(paste0("\\path\\[","fill="),"\\fill\\[",x[rectangles])
                               x
                             })
    if(pngconvert){    
      temppngfile<-tempfile(pattern = ".png")
      pdftools::pdf_convert(output,"png",dpi = 1200,filenames = temppngfile)
      file.copy(from=temppngfile,to=file.path(pdfexportfolder,"figure1.png"),overwrite = TRUE)
    }
  }
  figureXX
}
