#' Generates the second figure of the paper
#' 
#' This function 
#'  - generates a regular grid of size 201 on the \eqn{[0,1]^2} square.
#'  - uses \code{InformativeSpatialSelection::generateYonsquaregridU} to generate the values of \eqn{Y}. 
#' @param light a boolean (FALSE by default) indicating whether we should produce a light (with less pixels) version of the file. 
#' @return an object of classes "gg" and "ggplot": the first figure of the paper. 
#' @author DB
#' @seealso InformativeSpatialSelection::generateYonsquaregridU
#' @examples
#' figure2(light=TRUE)
#' figure2()
figure2<-function(pdfexportfolder=NULL,
                  pngconvert=FALSE,
                  figure2_data=NULL){
  library(sampling)
  library(ggplot2)
  if(is.null(figure2_data)){
    figure2_data=get(data(figure2_data,package="InformativeSpatialSelection"))}
  graphsubfigureslevels<-c("2.a. z=exp(xi0)=10",
                           "2.b. z=exp(xi0+ xi2 e)",
                           "2.c. z=exp(xi0+xi1 y+ xi2 e)",
                           "2.d. y")
  UZ$.t<-ordered(UZ$variable);
  levels(UZ$.t)<-  graphsubfigureslevels
  samples$.t<-ordered(samples$variable);levels(samples$.t)<-  graphsubfigureslevels
  levels(samples$type)<-paste0(levels(samples$type),"        ")
  #  UZ$.t<-factor(c("$\\\\\\\mathbf{z}=1$",
  #                  "",
  #                  "")Design variable ",1:3))[UZ$.n]
  
  #devtools::install_github("eliocamp/ggnewscale")
  library(ggnewscale)
  breakpoints<-c(1,5,10,20,50)
  figureXX<-ggplot(data=UZ, aes(x1,x2)) + 
    geom_tile(aes(fill = y),size=2) + 
    scale_fill_gradient(na.value="transparent",
                        low = "white",     
                        high = "black",
                        guide=guide_colorbar(order=3)) +
    labs(fill = "                      y0")+
    new_scale("fill") +
    geom_tile(aes(fill = value),size=2)  +
    labs(fill = "design variable")+ 
    xlab("")+
    ylab("")+
    #geom_point(data=UZ[figure2_data$I==1,],size=2,color="black",fill="white")+
    #geom_point(data=samples[samples$type=="bpp",],size=2)+
    geom_point(data=samples,mapping = aes(shape=type),fill="white",color="black",size=1.5,stroke=.9)+#shape=21+
    scale_shape_manual(values = c(21, 25),guide=guide_legend(order=1)) +
    labs(shape = "Sample")+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    theme(
      strip.text.x = element_text(angle = 0, hjust = 0),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      panel.border = element_blank(),
      panel.background=element_rect(fill="transparent",colour=NA),
      plot.background=element_rect(fill="transparent",colour=NA),
      legend.background = element_rect(fill = "transparent", colour = "transparent"),
      legend.key  = element_rect(fill = "transparent", colour = "transparent"),
      legend.spacing = unit(1, 'cm'),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()#,
      #panel.background = element_blank()
    )+
    scale_fill_gradient(na.value="transparent",
                        low = "white",     
                        high = "black",
                        breaks=breakpoints,
                        labels=breakpoints,
                        trans = 'log10',
                        guide=guide_colorbar(order=2))+#frame.colour=c("black"),frame.linewidth=2)
    
    labs(fill = "z1")+
    facet_wrap(~.t, nrow=1)+
    theme(legend.position="bottom",
          legend.justification="right",
          legend.margin=margin(0,0,0,0),
          legend.title=element_text(size=10),
          legend.box.margin=margin(-20,#top
                                   18,#right
                                   10,#bottom
                                   10))#left
  
  
  
  
  if(!is.null(pdfexportfolder)){    
    figbasename=paste0("figure2.pdf")
    output=file.path(pdfexportfolder,figbasename)
    SweaveLst::graph2pdffile("print(figureXX +scale_shape_manual(values = c(21, 25),labels=c('PppPppPpp','bppbppbpp'),guide=guide_legend(order=1)))",
                             widthe = 8.5,
                             height=2.88,
                             output = output,
                             modify = function(x){
                               #x[5]<-"\\usepackage{amsfonts,mathrsfs}\n"
                               x<-gsub("                      y0","$\\\\\\\\mathbf{y}$",x)
                               x<-gsub("PppPppPpp","$\\\\\\\\mathrm{Ppp}(\\\\\\\\mathbf{z})$",x)
                               x<-gsub("bppbppbpp","$\\\\\\\\mathrm{bpp}(\\\\\\\\mathbf{z},10)$",x)
                               x<-gsub("Sample","$~~~\\\\\\\\mathbf{d}$",x)
                               x<-gsub("z1","$\\\\\\\\mathbf{z}$",x)
                               #x<-gsub("log(10)","\\\\log(10)$",x,fixed=T)
                               #x<-gsub("log(z)=","$\\\\log(\\\\mathbf{z})=",x,fixed=T)
                               #x<-gsub("log(z)=","$\\\\log(\\\\mathbf{z})=",x,fixed=T)
                               #x<-gsub("log(z)","$\\\\log(\\\\mathbf{z})$",x,fixed=T)
                               x<-gsub("xi0",   "\\\\\\\\xi_0",x)
                               x<-gsub("=10",   "=10$",x)
                               x<-gsub("=exp",   "=\\\\\\\\exp",x)
                               x<-gsub(", 25.67) {",   ", 32.05) {",x,fixed=T)
                               x<-gsub("xi1 y",  "\\\\\\\\xi_1\\\\\\\\mathbf{y}",x)
                               x<-gsub("xi2 e)","\\\\\\\\xi_2\\\\\\\\mathbf{e})$",x,fixed=TRUE)
                               x<-gsub("2.d. y","2.d. $\\\\mathbf{y}$",x,fixed=T)
                               x<-gsub(", scale=  1.10","",x)# to remove automatic rescaling the legend title would be better if done
                               rectangles<-grep("rectangle",x)
                               x[rectangles]<-gsub(paste0("\\path\\[","fill="),"\\fill\\[",x[rectangles])
                               x
                             })
    fs::file_show(file.path(getwd(),"latex/fig/figure2.pdf"))
    if(pngconvert){
      temppngfile<-tempfile(pattern = ".png")
      pdftools::pdf_convert(output,"png",dpi = 1200,filenames = temppngfile)
      file.copy(from=temppngfile,to=file.path(pdfexportfolder,"figure2.png"),overwrite = TRUE)
    }
  }
  figureXX
}

