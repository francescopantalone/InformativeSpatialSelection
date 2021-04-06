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
#' figure6()
figure6<-function(pdfexportfolder=NULL,pngconvert=FALSE){
  data("figure6_data",package="InformativeSpatialSelection")
  figure6_data$cimin<-with(figure6_data,sqrt(vrho)*qnorm(0.025)+rho)
  figure6_data$cimax<-with(figure6_data,sqrt(vrho)*qnorm(0.975)+rho)
  
  
  library(ggplot2)
    figureXXgamma<-
      ggplot(data=subset(figure6_data,.Beta==1&!is.element(y,c(outer(c(1,2,5),10^{-3:3},"*")))),
                     aes(x=.Gamma,
                         y=rho,
                         ymin=cimin,
                         ymax=cimax,fill="red",alpha=.3))+
      geom_point()+
      geom_line()+
      scale_x_continuous(trans='log10')+
      geom_ribbon()+facet_wrap(~factor(y))
    
    figureXXbeta<-
      ggplot(data=subset(figure6_data,
                         .Gamma==1&!is.element(y,c(outer(c(1,2,5),10^{-3:3},"*")))),
             aes(x=.Beta,
                 y=rho,
                 ymin=cimin,ymax=cimax,fill="red",alpha=.3))+
      geom_point()+
      geom_line()+
      scale_x_continuous(trans='log10')+
      geom_ribbon()+facet_wrap(~factor(y))
    
    figureXXy<-
      ggplot(data=subset(figure6_data,
                         .Gamma==1&.Beta==1&is.element(y,c(outer(c(1,2,5),10^{-3:3},"*")))),
             aes(x=y,
                 y=rho,
                 ymin=cimin,ymax=cimax,fill="red",alpha=.3))+
      geom_point()+
      geom_line()+
      scale_x_continuous(trans='log10')+
      scale_y_continuous(trans='log10')+
      geom_ribbon()
    
  
    ggplot(data=figure6_data[figure6_data$.Gamma==1,],
           aes(x=.Beta,y=rho,color=factor(id),
               group=factor(id),
               ymin=cimin,ymax=cimax,fill="red",alpha=.3))+
      geom_point()+
      geom_line()+
      scale_x_continuous(trans='log10')+
      #scale_y_continuous(trans='log10')+
      geom_ribbon()+
      facet_wrap(~id)
    
    
    ggplot(data=subset(figure6_data,.Gamma==1&id==1),
           aes(x=.Beta,y=rho,color=factor(id),
               group=factor(id),
               ymin=cimin,ymax=cimax,fill="red",alpha=.3))+
      geom_point()+
      geom_line()+
      scale_x_continuous(trans='log10')+
      #scale_y_continuous(trans='log10')+
      geom_ribbon()+
      facet_wrap(~id)    
    
    
if(!is.null(pdfexportfolder)){    
    figbasename=paste0("figure6.pdf")
    output=file.path(pdfexportfolder,figbasename)
    SweaveLst::graph2pdffile("print(figureXX)",
                             widthe = 8.5,
                             height=2.88,
                             output = output)
    fs::file_show(file.path(getwd(),"latex/fig/figure2.pdf"))
    if(pngconvert){
      temppngfile<-tempfile(pattern = ".png")
      pdftools::pdf_convert(output,"png",dpi = 1200,filenames = temppngfile)
      file.copy(from=temppngfile,to=file.path(pdfexportfolder,"figure2.png"),overwrite = TRUE)
    }
  }
  figureXX
}

