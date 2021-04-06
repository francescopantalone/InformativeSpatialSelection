library(InformativeSpatialSelection)
set.seed(5)
Y<-Y3000
x1<-generatex_i(201)
x2=x1
X<-generatesquaregridU(x1=x1,x2=x1)
ell=sample(201^2,5)
Yc<-updateY(X,Y,ell=ell,y=c(2,1,0,-1,3))

plotriskmaps<-function(Yc,sel=1:nrow(Yc)){
  qq<-summary(apply(Yc,1,mean))
  Yc1<-1*Yc[sel,]>qq[3]
  Yc1sum<-apply(Yc1,2,sum)
  quantile(Yc1sum,probs=c(.25,.5,.75))
  Yc1<-Yc1[,order(Yc1sum)]
  XX<-cbind(X[sel,],data.frame(
    mean=apply(Yc1,1,mean),
    sd=apply(Yc1,1,sd),
    mean0.5=apply(Yc1[,1:60],1,mean),
    mean0.25=apply(Yc1[,1:750],1,mean),
    mean25.50=apply(Yc1[,751:1500],1,mean),
    mean50.75=apply(Yc1[,1501:2250],1,mean),
    mean75.100=apply(Yc1[,2251:3000],1,mean)),
    mean95.100=apply(Yc1[,2941:3000],1,mean),
    q1=Yc1[,1],
    q2=Yc1[,2],
    q750=Yc1[,750],
    q1500=Yc1[,1500],
    q1501=Yc1[,1501],
    q2250=Yc1[,2250],
    q2999=Yc1[,2999],
    q3000=Yc1[,3000])
  YY<-reshape2::melt(XX,id.vars=c("x1","x2"))
  
  ggplot2::ggplot(YY, ggplot2::aes(x1,x2)) + 
    ggplot2::geom_tile(ggplot2::aes(fill = value),size=2) + 
    ggplot2::xlab("")+
    ggplot2::ylab("")+
    ggplot2::scale_x_continuous(expand=c(0,0))+
    ggplot2::scale_y_continuous(expand=c(0,0))+
    ggplot2::theme(
      axis.text.x=ggplot2::element_blank(),
      axis.ticks.x=ggplot2::element_blank(),
      axis.text.y=ggplot2::element_blank(),
      axis.ticks.y=ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(), 
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank())+
    ggplot2::scale_fill_gradient(low = "white",     
                                 high = "black"#,guide=guide_colorbar(frame.colour=c("black"),frame.linewidth=2)
    )  +
    ggplot2::facet_wrap(~variable)+ 
    ggplot2::labs(fill = "signal")+ 
    ggplot2::theme(legend.position="bottom",
                   legend.justification="right",
                   legend.margin=ggplot2::margin(0,0,0,0),
                   legend.box.margin=ggplot2::margin(-20,#top
                                                     18,#right
                                                     10,#bottom
                                                     10))+#left
    geom_segment(data=data.frame(x   =c(0.5,0.5,0.5,  1),
                                 y   =c(0.5,0.5,  1,0.5),
                                 xend=c(  1,0.5,  1,  1),
                                 yend=c(0.5,  1,  1,  1)),
                                 aes(x=x,xend=xend,y=y,yend=yend),
                                 color="red",size=1)
}


plotriskmaps(Yc)
plotriskmaps(Yc,X$x1>.5&X$x2>.5)

qq[3]
ell
Yc[ell,1]