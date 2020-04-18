#-------------------------------------------------------------------
# Optimization of Spatial Sampling with the R packageSamplingStrata"
# R script for case study 1: the "meuse" dataset
# Univariate case
#-------------------------------------------------------------------

library(sp)
library(gstat)
library(automap)
library(SamplingStrata)

# locations (155 observed points)
data("meuse")
summary(meuse)
# grid of points (3103)
data("meuse.grid")
meuse.grid$id <- c(1:nrow(meuse.grid))
cor(meuse[,c(3:6)])
# png("./images/meuse.png")
par(mfrow=c(1,2))
plot(meuse.grid$x,meuse.grid$y,
     xlab="Longitude",
     ylab="Latitude")
title("meuse.grid dataset")
plot(meuse$x,meuse$y,
     xlab="Longitude",
     ylab="Latitude")
title("meuse dataset")
par(mfrow=c(1,1))
# dev.off()

coordinates(meuse)<-c("x","y")
coordinates(meuse.grid)<-c("x","y")


# png("./images/lead.png")
spplot(meuse,"lead",do.log=F,ylab="Lead localisation and concentration")
# dev.off()
# spplot(meuse.grid,"soil")
# levels(meuse.grid$soil)
# bubble(meuse,"lead",do.log=F,key.space="bottom")
# meuse.grid <- as(meuse.grid, "SpatialPixelsDataFrame")
# meuse.grid@bbox
# ha <- (meuse.grid@bbox[1,2]-meuse.grid@bbox[1,1]) * (meuse.grid@bbox[2,2]-meuse.grid@bbox[2,1]) / 20000 
# ha
#################
# kriging
#################

v <- variogram(lead ~ dist + soil, data=meuse)
fit.vgm <- autofitVariogram(lead ~ elev + soil, meuse, 
                            model = c("Exp"))
plot(v, fit.vgm$var_model)
fit.vgm$var_model
# model    psill    range
# 1   Nug 1524.895   0.0000
# 2   Exp 8275.431 458.3303
g <- gstat(NULL,"Lead", lead ~ dist + soil, meuse)
g
vm <- variogram(g)
vm.fit <- fit.lmc(vm, g, vgm(psill=fit.vgm$var_model$psill[2], 
                             model="Exp", range=fit.vgm$var_model$range[2], 
                             nugget=fit.vgm$var_model$psill[1]))

# Evaluation of fitting not possible, then fitting = 1


# Prediction on the whole grid
preds <- predict(vm.fit, meuse.grid)
names(preds)
# [1] [1] "Lead.pred" "Lead.var"
summary(preds$Lead.pred)
preds$Lead.pred <- ifelse(preds$Lead.pred < 0,0,preds$Lead.pred)
summary(preds$Lead.pred)
df <- NULL
df$Lead.pred <- preds@data$Lead.pred
df$Lead.var <- preds@data$Lead.var
df$dom1 <- 1
df <- as.data.frame(df)
df$id <- meuse.grid$id
df$lon <- meuse.grid$x
df$lat <- meuse.grid$y
head(df)



#------------------------
# Solution SamplingStrata
#------------------------
frame <- buildFrameSpatial(df=df, 
                  id="id", 
                  X=c("Lead.pred"), 
                  Y=c("Lead.pred"), 
                  variance=c("Lead.var"), 
                  lon="lon", 
                  lat="lat", 
                  domainvalue="dom1")
# frame <- buildFrameDF(df=df,
#                       id="id",
#                       X=c("Lead.pred"),
#                       Y=c("Lead.pred"),
#                       domainvalue = "dom1")
# frame$var1 <- df$Lead.var
# frame$lon <- meuse.grid$x
# frame$lat <- meuse.grid$y

cv <- as.data.frame(list(DOM=rep("DOM1",1),
                         CV1=rep(0.05,1),
                         domainvalue=c(1:1) ))
cv
# DOM  CV1 domainvalue
# 1 DOM1 0.05           1

set.seed(1234)
solution <- optimStrata (
  method="spatial",
  errors=cv, 
  framesamp=frame,
  iter = 50,
  pops = 10,
  nStrata = 5,
  fitting = 1,
  range = fit.vgm$var_model$range[2],
  kappa = 1,
  writeFiles = FALSE,
  showPlot = TRUE,
  parallel = FALSE
)

sum(round((solution$aggr_strata$SOLUZ)))

# expected_CV(solution$aggr_strata)
ss <- summaryStrata(x=solution$framenew,outstrata=solution$aggr_strata)
# sink("./output/ss.txt")
ss
# sink()
framenew <- solution$framenew
outstrata <- solution$aggr_strata
framenew$lon <- meuse.grid@coords[,1]
framenew$lat <- meuse.grid@coords[,2]

# plotStrata2d(framenew,outstrata,domain=1,vars=c("X1","X1"))

s <- selectSample(framenew,outstrata)

expected_CV(outstrata)
val1 <- evalSolution(framenew,outstrata,nsampl=1000,progress=F)
val1$coeff_var 

frameres2 <- SpatialPointsDataFrame(data=framenew, coords=cbind(framenew$LON,framenew$LAT) )
frameres3 <- SpatialPixelsDataFrame(points=frameres2[c("lon","lat")], data=framenew)
frameres3$LABEL <- as.factor(frameres2$LABEL)
spplot(frameres3,c("LABEL"), col.regions=bpy.colors(5))

#---------------------
# Solution with Ospats
#---------------------

#------------------ Ospats processing -------------------------
# N.B. : (1) use only Julia 0.6.4 
# (2) Modify the complete path of ./ospats in main.jl at line 9
#--------------------------------------------------------------
setwd("./ospats")
#-------------------------------------
add_variance <- 0
file <- NULL
file$x <- frame$lon
file$y <- frame$lat
file$id <- c(1:nrow(frame))
file$z <- frame$Y1
file$var <- frame$var1+add_variance
file <- as.data.frame(file)
write.table(file,"meuse.txt",row.names=F,col.names=F,sep=",",quote=F)
#--------------------------------------------------------
# devtools::install_github("Non-Contradiction/JuliaCall")
# library(JuliaCall)
# julia <- julia_setup()
# julia_console()
# include("main.jl")
# exit
############################################
### If the above 5 lines do not work then 
### execute "main.jl" in Julia 0.6.4 console
############################################
#--------------------------------------------------------
# fr <- read.table("stratification",sep=",",header=TRUE,dec=".",stringsAsFactors = FALSE)
# colnames(fr) <- c("lon","lat","stratum")
fr <- read.table("stratification",header=FALSE)
colnames(fr) <- "stratum"
table(fr$stratum)
# frameone <- merge(frame,fr,by=c("lon","lat"))
frameone <- cbind(frame,fr)
strata <- aggrStrataSpatial(dataset=frameone,
                            fitting=1,
                            range = fit.vgm$var_model$range[2], 
                            vett = frameone$stratum,
                            kappa = 1,
                            dominio = 1)
strata$allocation <- bethel(strata,cv)
# sink("ospats_strata.txt")
strata
# stratum    N COST CENS DOM1        M1       S1 allocation
# 1       1  146    1    0    1 338.05849 77.96705          6
# 2       2  804    1    0    1 106.22683 68.26663         28
# 3       3  534    1    0    1 231.61554 74.00451         21
# 4       4  589    1    0    1 163.32563 71.54368         22
# 5       5 1030    1    0    1  52.98693 70.85171         38
# sink()
# Plot
spfr <- SpatialPointsDataFrame(data=frameone, coords=cbind(frameone$lon,frameone$lat) )
spfr2 <- SpatialPixelsDataFrame(points=spfr[c("lon","lat")], data=fr)
spfr2$Ospats <- as.factor(fr$stratum)
spplot(spfr2,c("Ospats"), col.regions=bpy.colors(5))
#---------------------------------
# Comparison of results
#---------------------------------
# Ospats
# sink("./output/ospats_strata.txt")
strata[order(strata$M1),]
# stratum    N COST CENS DOM1        M1       S1 allocation
# 5       5 1030    1    0    1  52.98693 70.85171         38
# 2       2  804    1    0    1 106.22683 68.26663         28
# 4       4  589    1    0    1 163.32563 71.54368         22
# 3       3  534    1    0    1 231.61554 74.00451         21
# 1       1  146    1    0    1 338.05849 77.96705          6
sum(bethel(strata,cv))
table(fr$stratum)
#---------------------------------
# SamplingStrata
# STRATO    N        M1       S1 COST CENS DOM1 X1     SOLUZ
# 1      1  279  27.55464 66.50251    1    0    1  1  9.441627
# 2      2 1134  72.36470 68.94567    1    0    1  2 39.785483
# 3      3  862 137.71894 71.28072    1    0    1  3 31.266836
# 4      4  659 220.02812 75.03269    1    0    1  4 25.161731
# 5      5  169 329.93712 79.53935    1    0    1  5  6.840271
sum(bethel(solution$aggr_strata,cv))
table(framenew$LABEL)

# save.image(file="lead_5strata.RData")

