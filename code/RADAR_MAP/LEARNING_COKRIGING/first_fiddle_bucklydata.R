# based on http://www.css.cornell.edu/faculty/dgr2/_static/files/R_PDF/CoKrigeR.pdf
#pp is like meuse
#gps is like meuse.pb (incomplete sample)

library(sp)
library(gstat)

#data(meuse)#
#str(meuse)

#load data
DEM = read.csv('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/buck_DEM.csv' ,header=T)
gps = read.csv('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/buck_allpoints.csv',header=T)
DEM_overpoints = read.csv('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/buck_DEM_overpoints.csv' ,header=T)

newDEM = DEM[,c("x","y")]

newDEM_op = DEM_overpoints[,c("x","y")]

colnames(DEM)[3] = "h"
#DEM$h <- NA
colnames(gps)[3] = "h"
#gps$h <- NA

pp <- rbind(DEM, gps)

require("lattice")

h1 = histogram(~ h, pp, xlab="gps points", col="thistle3", nint=12)
print(h1, split =c(1,1,2,1), more=T)
rm(h1)

#display stuff
h1 <-histogram(~ h, pp, xlab="OM", col="lightblue", nint=12)
h2 <-histogram(~ h , pp, xlab="Zn", col="red4", nint=12)
h3 <-histogram(~log10(h), pp, xlab="log10(OM)", col="lightblue", nint=12)
h4 <-histogram(~log10(h) , pp, xlab="log10(Zn)", col="red4", nint=12)
print(h1, split =c(1,1,2,2), more=T)
print(h3, split =c(2,1,2,2), more=T)
print(h2, split =c(1,2,2,2), more=T)
print(h4, split =c(2,2,2,2), more=F)
rm(h1, h2, h3, h4)

#Simulation of under-sampling
#meuse.pb <- meuse[seq(1,length(meuse$lead), by=3),c("x", "y", "lead", "om", "zinc")]str(meuse.pb)
#meuse.pb <-cbind(meuse.pb,ltpb =log10(meuse.pb$lead),ltom =log10(meuse.pb$om),ltzn =log10(meuse.pb$zinc))
#str(meuse.pb)

# set log values
pp <- cbind(pp,lth =log10(pp$h)) #lt stands for log10
pp <- cbind(pp,lth =log10(pp$h))

DEM <- cbind(DEM,lth =log10(DEM$h)) #lt stands for log10
gps <- cbind(gps,lth =log10(gps$h))

#Modelling the spatial structure of the target variable
class(pp)
coordinates(pp) <- ~ x + y
coordinates(gps) <- ~ x + y
coordinates(DEM) <- ~ x + y
coordinates(newDEM)  <- ~ x + y
coordinates(newDEM_op)  <- ~ x + y
# alternate command format: coordinates(meuse) <- c("x", "y")

xyplot(y ~ x,as.data.frame(pp), asp="iso",panel =function(x, ...) {
  panel.points(coordinates(gps),cex=1.8*(log10(gps$h) - 1.3),pch=1, col="red");
  panel.points(coordinates(DEM),cex=1.8*(log10(DEM$h) - 1.3),pch=1, col="blue");
  panel.grid(h=-1, v=-1, col="darkgrey")})

#variogram cloud
plot(v.lth.c <-variogram(lth  ~ 1, data=gps, cloud=T))
plot(v.lth.c <-variogram(lth  ~ 1, data=gps, cutoff=1800, cloud=T))

#empircal variogram showing number of point pairs in each bin
plot(v.lth <-variogram(lth  ~ 1, data=gps, cutoff=1800, width=5), pl=T)

# estimate variogram model form and params by eye
m.lth <-vgm(0.00004,"Pow",1,0)
plot(v.lth, pl=T, model=m.lth)

#fit model params by weighted leatsquares
(m.lth.f <-fit.variogram(v.lth, m.lth))
plot(v.lth, pl=T, model=m.lth.f)
rm(v.lth.c, v.lth, m.lth)

#Ordinary kridging
# interpolate
k.o <-krige(lth ~1, locations=gps, newdata=newDEM, model=m.lth.f)

summary(k.o)

source("/Users/home/whitefar/DATA/code/RADAR_MAP/LEARNING_COKRIGING/ck_plotfns.R")
plot.kresults(k.o, "var1", pp, gps, "h", "log10(h), OK")

# predict at the extra points
k <- krige(lth  ~ 1, gps, newDEM, m.lth.f)
# compute and summarize evaluation errors
summary(k)
diff <- k$var1.pred - DEM$lth
summary(diff)
sqrt(mean(sum(diff^2)))   # RMSE (precision)
mean(diff)           # mean error (bias)
median(DEM$lth)         # median error

#GOT TO HERE

#Modelling the covariable
xyplot(gps$lth ~ DEM$lth,   pch=20, cex=1.2,col="blue", ylab="log10(h_gps)", xlab="log10(h_DEM)")

with(DEM@data,cor(DEM$h, gps$h))

sum(is.na(meuse.pb$om))

with(meuse.pb@data,cor(ltom, ltpb, use = "complete"))

#model the omnidirectional spatal variogram
# all valid covariable observations, with coordinates
meuse.co <-subset(as.data.frame(meuse), !is.na(om),c(x, y, om))
# add log10-transformed variables for convenience
meuse.co <-cbind(meuse.co, ltom =log10(meuse.co$om))
str(meuse.co)

# convert to spatial object
coordinates(meuse.co) <- ~ x + y
# experimental variogram
v.ltom <-variogram(ltom ~ 1, meuse.co, cutoff=1800)
plot(v.ltom, pl=T)
# model by eye
m.ltom <-vgm(.035, "Sph", 800, .015)
# fit(m.ltom.f <-fit.variogram(v.ltom, m.ltom))

plot(v.ltom, pl=T, model=m.ltom.f)
# compare variogram structure to target variable
m.ltom.f$range[2]; m.ltpb.f$range[2].
round(m.ltom.f$psill[1]/sum(m.ltom.f$psill),2)
round(m.ltpb.f$psill[1]/sum(m.ltpb.f$psill),2)

#Building a data structure to model co-regionalisation

(g <-gstat(NULL, id = "ltpb", form = ltpb ~ 1, data=meuse.pb))

(g <-gstat(g, id = "ltom", form = ltom ~ 1, data=meuse.co))

v.cross <-variogram(g)
str(v.cross)

plot(v.cross, pl=T)

# Fitting a linear model of co-regionalisation
(g <-gstat(g, id = "ltpb", model = m.ltpb.f, fill.all=T))


(g <-fit.lmc(v.cross, g, fit.method=6, correct.diagonal=1.01))

plot(variogram(g), model=g$model)

# Comparing models of regionalisation and co-regionalisation
str(m.ltom.f)
str(g$data$ltpb, max.level = 1)
str(g$model$ltpb, max.level = 1)


g$model$ltom$psill - m.ltom.f$psill

sum(g$model$ltom$psill) -sum(m.ltom.f$psill)

sum(g$model$ltom$psill)

g$model$ltpb$psill - m.ltpb.f$psill

sum(g$model$ltpb$psill) -sum(m.ltpb.f$psill)

sum(g$model$ltpb$psill)

# Co-kriging with one co-variable

# interpolate
k.c <-predict(g, meuse.grid)

str(k.c)

# summarize predictions and their errors
summary(k.c$ltpb.pred);summary(k.c$ltpb.var)

# predict at the extra points
k <-predict(g, meuse.extra)

# compute and summarize prediction errors
diff <- k$ltpb.pred - meuse.extra$ltpb
summary(diff)

sqrt(sum(diff^2)/length(diff))

sum(diff)/length(diff)

diff <-as.data.frame(diff)
coordinates(diff) <-coordinates(meuse.extra)
bubble(diff, zcol="diff", pch=1,main="CK evaluation errors at undersampled points, log10(Pb)")
