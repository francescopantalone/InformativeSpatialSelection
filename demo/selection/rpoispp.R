#This demo is to test the poisson process sampling spatstat::rpoispp

ulimit::memory_limit(12500)
#rpoispp
library("sp")
library("spatstat")
library("gstat")
library("plotly")
library("maps")
data(meuse,package="sp")
data(meuse.grid,package="sp")
meuse.grid2<-meuse.grid
sp::coordinates(meuse) = ~x+y
sp::coordinates(meuse.grid2) = ~x+y
lznr.vgm = variogram(log(zinc)~sqrt(dist), meuse)
lznr.fit = fit.variogram(lznr.vgm, model = vgm(1, "Exp", 300, 1))
lzn.vgm = variogram(log(zinc)~1, meuse)
lzn.fit = fit.variogram(lzn.vgm, model = vgm(1, "Sph", 900, 1))
lzn.kriged = krige(log(zinc)~1, meuse, meuse.grid, model = lzn.fit)
lzn.kriged = krige(log(zinc)~1, meuse, meuse.grid, model = lzn.fit)

lzn.fit = fit.variogram(lzn.vgm, model = vgm(1, "Gau", 900, 1))
lzn.kriged = krige(log(zinc)~1, meuse, meuse.grid2, model = lzn.fit)

plot_ly(cbind(meuse.grid[c("x","y")],z=lzn.kriged$var1.pred), x=~x,y=~y,z=~z,color=~z)

pred.fun<-function(x,y){
  A<-data.frame(x=x,y=y)
  coordinates(A)=~x+y
  krige(log(zinc)~1, meuse, A, model = lzn.fit)$var1.pred/600000}

     pp <- rpoispp(pred.fun, 
       lmax=max(lzn.kriged$var1.pred)/6000000, 
       win=as.owin(W=list(xrange=range(meuse$x),yrange=range(meuse$y))),nsim=1)

     
     
     plot(pp, cex = 0.5, pch = 19)
     rpoispp(pred.fun, 
             lmax=max(lzn.kriged$var1.pred), 
             win=owin(range(meuse$x),range(meuse$y)),nsim=100)
     
     rpoispp()



meuse_p <- ppp(meuse$x, meuse$y, window = owin(range(meuse$x), range(meuse$y)))


pp <- rpoispp(1, win=owin(c(0,10),c(0,10)))

pp <- rpoispp(1, win = owin(c(min(meuse$x), max(meuse$x)), c(min(meuse$y), max(meuse$y))), nsim = 100)


pp <- rpoispp(1, win = owin(range(meuse$x), range(meuse$y)), nsim = 100)


us <- map_data("usa")

# Simplified outline of the USA
us.simple <- function() {
  map <- data.frame(x = us$long, y = us$lat)
  map <- subset(map, !is.na(x))
  map <- map[!duplicated(paste(map$x, map$y)), ]
  points <- c(177, 563, 1296, 1683, 1857, 2420, 3649, 3887, 4034, 4241, 4326, 
              4422, 4624, 4861, 5010, 5071, 5331, 5559, 5744, 5833, 6097, 6315, 6608, 
              6748)
  map[points, ]
}
ow <- owin(range(us$long), range(us$lat), poly = us.simple())
plot(ow$xr, ow$yr, type = "n", xlab = "Long.", ylab = "Lat.")
polygon(ow$bdry[[1]], col = rgb(0.5, 0.5, 0.1, 0.4))