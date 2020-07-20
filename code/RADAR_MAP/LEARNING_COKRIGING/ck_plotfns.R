### R source code, file ck_plotfns.R
### author: DG Rossiter, d.g.rossiter@cornell.edu
###
### to accompany Technical Note "Co-kriging with the gstat package
###    of the R environment for statistical computing"
### plot Kriging results: predictions and variances
##     with postplot of sample/subsample points on predictions
##        colour: bpy, grey, red
##     and locations of sample/subsample points on prediction variance
##        colour: cm, green, red
##  arguments
##     kr.o     SpatialPointsDataFrame from kriging
##         coordinates must be named x, y
##     var1     prefix for field names *.pred, *.var, e.g. "var1"
##     samp.pts  SpatialPointsDataFrame sample points
##     subsamp.pts  SpatialPointsDataFrame subsample points
##     f       field name (quoted) or number for point data values
##     title
plot.kresults <- function(kr.o, var1, samp.pts, subsamp.pts,
                          f=1, title="") {
  to.eval <- paste("plot.kr <- levelplot(",
                   paste(var1,"pred",sep="."),
                   " ~ x+y, as.data.frame(kr.o),
            aspect='iso',
            col.regions=bpy.colors(64), cut=32,
            main=paste(title, 'Prediction'),
                panel = function(x, ...) {
                    panel.levelplot(x, ...);
                    panel.points(coordinates(samp.pts), col='grey',
                         pch=20,
                         # log scale, but still show minimum
                         cex=1.6 * (log10(samp.pts[[f]]) - 
                                0.9 * min(log10(samp.pts[[f]]))));
                    panel.points(coordinates(subsamp.pts), col='red',
                         pch=20,
                         cex=1.6 * (log10(subsamp.pts[[f]]) -
                                0.9 * min(log10(subsamp.pts[[f]]))));
                    panel.grid(h=-1, v=-1, col='')
             })"  );
  eval(parse(text=to.eval));
  to.eval <- paste("plot.kr.e <- levelplot(",
                   paste(var1,"var",sep="."),
                   "~ x+y, as.data.frame(kr.o),
      aspect='iso', 
        col.regions=cm.colors(64), cut=32,
        main=paste(title, 'Prediction variance'),
             panel = function(x, ...) {
                panel.levelplot(x, ...);
                panel.points(coordinates(samp.pts), col='green', 
                     pch=20, cex=.6);
                panel.points(coordinates(subsamp.pts), col='red',
                      pch=20, cex=.8); # subsample points larger
                panel.grid(h=-1, v=-1, col='darkgrey')
    })"  );
  eval(parse(text=to.eval));
  print(plot.kr, split = c(1,1,2,1), more=T);
  print(plot.kr.e, split = c(2,1,2,1), more=F)
}

### plot Kriging validation and cross-validation errors
##        colours: validation: bubble() default: palette()[2:3]
##        colours: x-valid: palette()[4:5]
##  arguments
##     kv.o     SpatialPointsDataFrame from kriging to validation points
##     var1     prefix for kriging field name *.pred, e.g. "var1"
##     valid.pts  SpatialPointsDataFrame with validation points
##     f        field name (quoted) or number for point data values
##     cv.o     SpatialPointsDataFrame from x-validation kriging
##     title
plot.valids <- function(kv.o, var1, valid.pts, f, cv.o, title="") {
  # validation errors
  to.eval <- paste("diff <- kv.o$", paste(var1,"pred",sep="."),
                   " - valid.pts[[f]]")
  eval(parse(text=to.eval))
  extreme <- max(abs(range(diff, as.data.frame(cv.o)$residual)))
  d <- SpatialPointsDataFrame(kv.o, data=as.data.frame(diff))
  b1 <- bubble(d,
               main=paste(title,"Validation errors"),
               maxsize = 3 * (max(abs(range(diff))))/extreme,
               panel = function(x, ...) {
                 panel.xyplot(x, ...);
                 panel.grid(h=-1, v=-1, col="darkgrey")}
  )
  b2 <- bubble(cv.o, z="residual",
               main=paste(title,"Cross-validation errors"), col=c(4,5),
               maxsize = 3 * (max(abs(range(cv.o$residual))))/extreme,
               panel = function(x, ...) {
                 panel.xyplot(x, ...);
                 panel.grid(h=-1, v=-1, col="darkgrey")}
  )
  print(b1, split=c(1, 1, 2, 1), more=T)
  print(b2, split=c(2, 1, 2, 1), more=F)
}

### compare Kriging results: predictions or prediction variance
##     with locations of sample/subsample points on predictions
##        colour: bpy, grey, red
##     and locations of sample/subsample points on prediction variance
##        colour: cm, green, red
##  arguments
##     k1.o, k2.o     SpatialPointsDataFrame'sfrom kriging
##         coordinates must be named x, y
##     var1.1, var1.2 prefix for field names *.pred, *.var, e.g. "var1"
##         in the two objects
##     samp.pts  SpatialPointsDataFrame sample points
##     subsamp.pts  SpatialPointsDataFrame subsample points
##     type   what to compare: "pred" or "var"
##     title.1, title.2   titles for the two kriging objects
plot.cf <- function(k1.o, k2.o, var1.1, var1.2, samp.pts, subsamp.pts,
                    type="pred", title.1="", title.2="") {
  # common scale
  to.eval <- paste("range <- range(k1.o$", paste(var1.1, type, sep="."),
                   ", k2.o$", paste(var1.2, type, sep="."), ")", sep="")
  eval(parse(text=to.eval));  # makes range
  breaks <- seq(range[1], range[2], length=32)
  to.eval <- paste("plot.k1 <- levelplot(",
                   paste(var1.1,type,sep="."),
                   " ~ x+y, as.data.frame(k1.o), aspect='iso', col.regions=",
                   ifelse(type=='pred', 'bpy.colors', 'cm.colors'),
                   "(64), at=breaks, main=title.1,
                   panel = function(x, ...) {
                          panel.levelplot(x, ...);
                          panel.grid(h=-1, v=-1, col='darkgrey')})"
  );
  eval(parse(text=to.eval));  # makes plot.k1
  to.eval <- paste("plot.k2 <- levelplot(",
                   paste(var1.2,type,sep="."),
                   " ~ x+y, as.data.frame(k2.o), aspect='iso', col.regions=",
                   ifelse(type=='pred', 'bpy.colors', 'cm.colors'),
                   "(64), at=breaks, main=title.2,
                   panel = function(x, ...) {
                          panel.levelplot(x, ...);
                          panel.grid(h=-1, v=-1, col='darkgrey')})"
  );
  eval(parse(text=to.eval));  # makes plot.k2
  to.eval <- paste("diff <- k2.o$", paste(var1.2, type, sep="."),
                   "- k1.o$", paste(var1.1, type, sep="."), sep="")
  print(to.eval)
  eval(parse(text=to.eval));  # makes diff
  tmp <- data.frame(x = coordinates(k1.o)[,"x"],
                    y = coordinates(k2.o)[,"y"], diff)
  plot.diff <-
    levelplot(diff ~ x+y, tmp, aspect="iso",
              col.regions=topo.colors(64), cut=32,
              main="Difference",
              panel = function(x, ...) {
                panel.levelplot(x, ...);
                panel.points(coordinates(samp.pts), col="white",
                             pch=20, cex=.6);
                panel.points(coordinates(subsamp.pts), col="red",
                             pch=20, cex=.8);
                panel.grid(h=-1, v=-1, col="darkgrey")
              })
  # display the plots in one figure
  print(plot.k1, split = c(1,1,2,2), more=T)
  print(plot.k2, split = c(2,1,2,2), more=T)
  print(plot.diff, split = c(1,2,2,2), more=F)
}
