#-------------------------------------------------------------------
# Optimization of Spatial Sampling with the R packageSamplingStrata"
# R script for case study 1: the "meuse" dataset
# Univariate case
#-------------------------------------------------------------------

library(lattice) # for meuse and meuse.grid
library(sp)
library(gstat)
library(automap)
library(magrittr)
library(ggplot2)
library(SamplingStrata)

# locations (155 observed points)
data("meuse")
summary(meuse)
# grid of points (3103)
data("meuse.grid")
meuse.grid$id <- c(1:nrow(meuse.grid))
cor(meuse[,c(3:6)])
par(mfrow=c(1,2))
png("./images/meuse.png")
par(mfrow=c(1,2))
plot(meuse.grid$x,meuse.grid$y,
     xlab="Longitude",
     ylab="Latitude")
title("meuse.grid dataset")
plot(meuse$x,meuse$y,
     xlab="Longitude",
     ylab="Latitude")
title("meuse dataset")
dev.off()

coordinates(meuse)<-c("x","y")
coordinates(meuse.grid)<-c("x","y")


png("./images/lead.png")
spplot(meuse,"lead",do.log=F,ylab="Lead localisation and concentration")
dev.off()
spplot(meuse.grid,"soil")
levels(meuse.grid$soil)
bubble(meuse,"lead",do.log=F,key.space="bottom")
meuse.grid <- as(meuse.grid, "SpatialPixelsDataFrame")
meuse.grid@bbox
ha <- (meuse.grid@bbox[1,2]-meuse.grid@bbox[1,1]) * (meuse.grid@bbox[2,2]-meuse.grid@bbox[2,1]) / 20000 
ha
v <- variogram(lead ~ 1, data=meuse, cutoff=2000, width=2000/20)
fit.vgm = autofitVariogram(lead ~ 1, meuse, model = "Exp")
plot(v, fit.vgm$var_model)
fit.vgm$var_model
# model     psill    range
# 1   Nug  2590.284   0.0000
# 2   Exp 14803.573 520.6309

# univariate case: lead as target variable
g <- gstat(NULL, "Pb", lead ~ 1, meuse)
# g <- gstat(g, "Cu", copper ~ 1, meuse)
# g <- gstat(g, "Pb", lead ~ 1, meuse)
# g <- gstat(g, "Zn", zinc ~ 1, meuse)
g
vm <- variogram(g)
vm.fit <- fit.lmc(vm, g, 
                  vgm(psill=fit.vgm$var_model$psill[2], 
                      model="Exp", 
                      range=fit.vgm$var_model$range[2], 
                      nugget=fit.vgm$var_model$psill[1]))


preds <- predict(vm.fit, meuse.grid)
names(preds)
summary(preds$Pb.pred)
preds$Pb.pred <- ifelse(preds$Pb.pred < 0,0,preds$Pb.pred)
summary(preds$Pb.pred)
df <- NULL
# df$Cd.pred <- preds@data$Cd.pred
# df$Cd.var <- preds@data$Cd.var
# df$Cu.pred <- preds@data$Cu.pred
# df$Cu.var <- preds@data$Cu.var
df$Pb.pred <- preds@data$Pb.pred
df$Pb.var <- preds@data$Pb.var
# df$Zn.pred <- preds@data$Zn.pred
# df$Zn.var <- preds@data$Zn.var
# df$dom <- meuse.grid@data$soil
df$dom1 <- 1
df <- as.data.frame(df)
df$id <- meuse.grid$id



#-----------------------------------------------------------
# Solution with spatial optimization on continuous variables
#-----------------------------------------------------------
frame <- buildFrameDF(df=df,
                      id="id",
                      X=c("Pb.pred"),
                      Y=c("Pb.pred"),
                      domainvalue = "dom1")
frame$var1 <- df$Pb.var
frame$lon <- meuse.grid$x
frame$lat <- meuse.grid$y

cv <- as.data.frame(list(DOM=rep("DOM1",1),
                         CV1=rep(0.05,1),
                         domainvalue=c(1:1) ))
cv
# DOM  CV1 domainvalue
# 1 DOM1 0.05           1

set.seed(1234)
solution <- optimizeStrataSpatial (
  errors=cv, 
  framesamp=frame,
  iter = 50,
  pops = 10,
  nStrata = 3,
  fitting = 1,
  range = fit.vgm$var_model$range[2],
  gamma=3,
  writeFiles = FALSE,
  showPlot = TRUE,
  parallel = FALSE
)

sum(round((solution$aggr_strata$SOLUZ)))

# expected_CV(solution$aggr_strata)
ss <- summaryStrata(x=solution$framenew,outstrata=solution$aggr_strata)
sink("./output/ss.txt")
ss
sink()
framenew <- solution$framenew
outstrata <- solution$aggr_strata
framenew$lon <- meuse.grid@coords[,1]
framenew$lat <- meuse.grid@coords[,2]

plotStrata2d(framenew,outstrata,domain=1,vars=c("X1","X1"))

s <- selectSample(framenew,outstrata)

expected_CV(outstrata)
val1 <- evalSolution(framenew,outstrata,nsampl=1000,progress=F)

val1$coeff_var 

frameres2 <- SpatialPointsDataFrame(data=framenew, coords=cbind(framenew$LON,framenew$LAT) )
frameres3 <- SpatialPixelsDataFrame(points=frameres2[c("lon","lat")], data=framenew)
frameres3$LABEL <- as.factor(frameres2$LABEL)
spplot(frameres3,c("LABEL"), col.regions=bpy.colors(5))

#-----------------------------------------------------------
# Solution with ospats
#-----------------------------------------------------------

#------------------ Ospats processing
setwd("./ospats")
# devtools::install_github("Non-Contradiction/JuliaCall")
library(JuliaCall)
julia <- julia_setup()
# julia <- julia_setup(JULIA_HOME = "C:\\Users\\Giulio\\AppData\\Local\\Julia-0.6.4\\bin")
#-------------------------------------
add_variance <- 0
file <- NULL
file$x <- frame$lon
file$y <- frame$lat
file$id <- c(1:nrow(frame))
file$z <- frame$Y1
file$var <- frame$var1+add_variance
# file$stratum <- frame$LABEL
file <- as.data.frame(file)
# file <- file[order(paste(file$x,file$y,sep="")),]
write.table(file,"example.txt",row.names=F,col.names=F,sep=",",quote=F)
#--------------------------------------------------------
julia_console()
include("main.jl")
exit
#--------------------------------------------------------
# fr <- read.table("stratification",sep="\t",header=FALSE)
fr <- read.table("stratification",sep=",",header=TRUE,dec=".",stringsAsFactors = FALSE)
colnames(fr) <- c("lon","lat","stratum")
table(fr$stratum)
frameone <- merge(frame,fr,by=c("lon","lat"))
strata <- aggrStrataSpatial(dataset=frameone,
                            fitting=1,
                            range = fit.vgm$var_model$range[2], 
                            vett = frameone$stratum,
                            gamma = 1,
                            dominio = 1)
strata$allocation <- bethel(strata,cv)
sink("./output/ospats_strata.txt")
strata
sink()
# stratum    N COST CENS DOM1        M1       S1 allocation
# 1       1 1065    1    0    1 143.00822 73.59579         43
# 2       2  677    1    0    1 254.62742 87.53205         33
# 3       3 1360    1    0    1  61.96497 71.99991         53
solution$aggr_strata
sum(bethel(strata,cv))
sum(bethel(solution$aggr_strata,cv))
table(framenew$LABEL)
table(fr$stratum)
spfr <- SpatialPointsDataFrame(data=fr, coords=cbind(fr$lon,fr$lat) )
spfr2 <- SpatialPixelsDataFrame(points=spfr[c("lon","lat")], data=fr)
spfr2$Ospats <- as.factor(fr$stratum)
spplot(spfr2,c("Ospats"), col.regions=bpy.colors(5))



save.image(file="lead_4strata.RData")

