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
pts = read.csv('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/buck_allpoints_DEMsampled.csv' ,header=T)
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



pts = read.csv('/Users/home/whitefar/DATA/FIELD_ANT_19/POST_FIELD/RADAR_MAP/PRACTICE_INTERP/buck_allpoints_DEMsampled.csv' ,header=T)
coordinates(pts)  <- ~ x + y
colnames(pts)[4] = "h_DEM"


#Modelling the covariable
xyplot(h_gps ~ h_DEM, data=pts,   pch=20, cex=1.2,col="blue", ylab="log10(h_gps)", xlab="log10(h_DEM)")

with(pts,cor(h_gps, h_DEM, use = "complete"))

#model the omnidirectional spatal variogram
# all valid covariable observations, with coordinates
pts.co <-subset(as.data.frame(pts), !is.na(h_DEM),c(x, y, h_DEM))
pts.com <-subset(as.data.frame(pts), !is.na(h_gps),c(x, y, h_gps))
# add log10-transformed variables for convenience
#meuse.co <-cbind(meuse.co, ltom =log10(meuse.co$om))
str(pts.co)

# convert to spatial object
coordinates(pts.co) <- ~ x + y
coordinates(pts.com) <- ~ x + y


#from above
#plot(v.lth <-variogram(lth  ~ 1, data=gps, cutoff=1800, width=5), pl=T)
# estimate variogram model form and params by eye
#m.lth <-vgm(0.00004,"Pow",1,0)
#plot(v.lth, pl=T, model=m.lth)

#fit model params by weighted leatsquares
#(m.lth.f <-fit.variogram(v.lth, m.lth))
#plot(v.lth, pl=T, model=m.lth.f)


# experimental variogram DEM
v.h_DEM <-variogram(h_DEM ~ 1, data = pts.co, cutoff=1800,width=5)
plot(v.h_DEM, pl=T)
# model by eye
m.h_DEM <-vgm(0.08, "Pow", 2, 0) #vgm(slope,model,exponent,x_inercept)
plot(v.h_DEM, pl=T, model=m.h_DEM)

(m.h_DEM.f <-fit.variogram(v.h_DEM, m.h_DEM))

plot(v.h_DEM, pl=T, model=m.h_DEM.f)
# experimental variogram DEM
v.h_gps <-variogram(h_gps ~ 1, pts.com, cutoff=1800,width=5)
plot(v.h_gps, pl=T)
# model by eye
m.h_gps <-vgm(0.08, "Pow", 2, 0) #vgm(slope,model,exponent,x_inercept)
plot(v.h_gps, pl=T, model=m.h_gps)

(m.h_gps.f <-fit.variogram(v.h_gps, m.h_gps))

plot(v.h_gps, pl=T, model=m.h_gps.f)
# compare variogram structure to target variable
m.h_DEM.f$range[2]; m.h_gps.f$range[2]
round(m.h_DEM.f$psill[1]/sum(m.h_DEM.f$psill),2)
round(m.h_gps.f$psill[1]/sum(m.h_gps.f$psill),2)

#Building a data structure to model co-regionalisation

(g <-gstat(NULL, id = "h_gps", form = h_gps ~ 1, data=pts))

(g <-gstat(g, id = "h_DEM", form = h_DEM ~ 1, data=pts))

v.cross <-variogram(g)
str(v.cross)

plot(v.cross, pl=T)

# Fitting a linear model of co-regionalisation
(g <-gstat(g, id = "h_gps", model = m.h_gps.f, fill.all=T))


(g <-fit.lmc(v.cross, g, fit.method=6, correct.diagonal=1.01))

plot(variogram(g), model=g$model)

# Comparing models of regionalisation and co-regionalisation
str(m.h_gps.f)
str(g$data$h_gps, max.level = 1)
str(g$model$h_gps, max.level = 1)


g$model$h_DEM$psill - m.h_gps.f$psill

sum(g$model$h_DEM$psill) -sum(m.h_DEM.f$psill)

sum(g$model$h_DEM$psill)

g$model$h_gps$psill - m.h_gps.f$psill

sum(g$model$h_gps$psill) -sum(m.h_gps.f$psill)

sum(g$model$h_gps$psill)

# Co-kriging with one co-variable

# interpolate
k.c <-predict(g, newDEM_op)


str(k.c)

# summarize predictions and their errors
summary(k.c$h_gps.pred);summary(k.c$h_gps.var)

# predict at the extra points
k <-predict(g, meuse.extra)

# compute and summarize prediction errors
diff <- k$h_gps.pred - meuse.extra$h_gps
summary(diff)

sqrt(sum(diff^2)/length(diff))

sum(diff)/length(diff)

diff <-as.data.frame(diff)
coordinates(diff) <-coordinates(meuse.extra)
bubble(diff, zcol="diff", pch=1,main="CK evaluation errors at undersampled points, log10(Pb)")
