#-------------------------------------------------------------------
# Optimization of Spatial Sampling with the R packageSamplingStrata"
# R script for case study 1: the "meuse" dataset
# Multivariate case 
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
coordinates(meuse)<-c("x","y")
coordinates(meuse.grid)<-c("x","y")
meuse.grid <- as(meuse.grid, "SpatialPixelsDataFrame")

#---------------------------------------------------
# Estimation of psill, nugget and range with automap
# Cadmium
v <- variogram(cadmium ~ dist + soil, data=meuse, cutoff=3000, width=3000/30)
plot(v)
fit.vgm = autofitVariogram(cadmium ~ dist + soil, meuse, 
                           model = c("Sph","Exp"))
plot(v, fit.vgm$var_model)
fit.vgm$var_model
# model    psill    range
# 1   Nug 2.821884   0.0000
# 2   Exp 5.121428 240.6374

# Copper
v <- variogram(copper ~ dist + soil, data=meuse, cutoff=3000, width=3000/30)
plot(v)
fit.vgm = autofitVariogram(copper ~ dist + soil, meuse,  
                           model = c("Sph","Exp"))
plot(v, fit.vgm$var_model)
fit.vgm$var_model
# model    psill   range
# 1   Nug   0.0000  0.0000
# 2   Exp 339.6246 98.4095

# Lead
v <- variogram(lead ~ dist + soil, data=meuse, cutoff=3000, width=3000/30)
plot(v)
fit.vgm = autofitVariogram(lead ~ dist + soil, meuse,   
                           model = c("Sph","Exp"))
plot(v, fit.vgm$var_model)
fit.vgm$var_model
# model    psill    range
# 1   Nug 3738.507    0.000
# 2   Sph 6243.891 1020.692

# Zinc
v <- variogram(zinc ~ dist + soil, data=meuse, cutoff=3000, width=3000/30)
plot(v)
fit.vgm = autofitVariogram(zinc ~ dist + soil, meuse,   
                           model = c("Sph","Exp"))
plot(v, fit.vgm$var_model)
fit.vgm$var_model
# model    psill    range
# 1   Nug 25328.92   0.0000
# 2   Exp 67450.62 400.5704

# multivariate case:cadmium, copper, lead, zinc as target variables
g <- NULL
g <- gstat(NULL, "Cd", cadmium ~ dist + soil, meuse)
g <- gstat(g, "Cu", copper ~ dist + soil, meuse)
g <- gstat(g, "Pb", lead ~ dist + soil, meuse)
g <- gstat(g, "Zn", zinc ~ dist + soil, meuse)
vm <- variogram(g)
# Using "Sph" because "Exp" does not work
vm.fit <- fit.lmc(vm, g, vgm(psill=1, model=c("Sph"), range=800, nugget=1))
vm.fit
plot(vm, vm.fit)
preds <- predict(vm.fit, meuse.grid)
names(preds)
# [1] "Cd.pred"   "Cd.var"    "Cu.pred"   "Cu.var"    "Pb.pred"  
# [6] "Pb.var"    "Zn.pred"   "Zn.var"    "cov.Cd.Cu" "cov.Cd.Pb"
# [11] "cov.Cu.Pb" "cov.Cd.Zn" "cov.Cu.Zn" "cov.Pb.Zn"

preds$Cd.pred <- ifelse(preds$Cd.pred < 0,0,preds$Cd.pred)
preds$Cu.pred <- ifelse(preds$Cu.pred < 0,0,preds$Cu.pred)
preds$Pb.pred <- ifelse(preds$Pb.pred < 0,0,preds$Pb.pred)
preds$Zn.pred <- ifelse(preds$Zn.pred < 0,0,preds$Zn.pred)
summary(preds)
# Optimization with SamplingStrata

df <- NULL
df$Cd.pred <- preds@data$Cd.pred
df$Cd.var <- preds@data$Cd.var
df$Cu.pred <- preds@data$Cu.pred
df$Cu.var <- preds@data$Cu.var
df$Pb.pred <- preds@data$Pb.pred
df$Pb.var <- preds@data$Pb.var
df$Zn.pred <- preds@data$Zn.pred
df$Zn.var <- preds@data$Zn.var
df$dom <- meuse.grid@data$soil
df$dom1 <- 1
df <- as.data.frame(df)
df$id <- meuse.grid$id

cv <- as.data.frame(list(DOM=rep("DOM1",1),
                         CV1=rep(0.05,1),
                         CV2=rep(0.05,1),
                         CV3=rep(0.05,1),
                         CV4=rep(0.05,1),
                         domainvalue=c(1:1) ))
cv
# DOM  CV1  CV2  CV3  CV4 domainvalue
# 1 DOM1 0.05 0.05 0.05 0.05           1
frame <- buildFrameDF(df=df,
                      id="id",
                      X=c(
                        "Cd.pred",
                        "Cu.pred",
                        "Pb.pred",
                        "Zn.pred"),
                      Y=c(
                        "Cd.pred",
                        "Cu.pred",
                        "Pb.pred",
                        "Zn.pred"),
                      domainvalue = "dom1")
frame$lon <- meuse.grid$x
frame$lat <- meuse.grid$y
frame$var1 <- df$Cd.var
frame$var2 <- df$Cu.var
frame$var3 <- df$Pb.var
frame$var4 <- df$Zn.var


kmeans <- KmeansSolutionSpatial(frame,
                                fitting=1,
                                range=800,
                                kappa=1,
                                gamma=0,
                                errors=cv,
                                maxclusters = 10)
# The difference in sample size after nStrata = 5 is negligible
# so we choose nStrata = 5

#---------------------------------------------------
# Solution with spatial optimization on continuous variables
#---------------------------------------------------

set.seed(4321)
solution <- optimizeStrataSpatial (
  errors=cv, 
  framesamp=frame,
  iter = 50,
  pops = 10,
  nStrata = 5,
  fitting = 1,
  range = 800,
  kappa = 1,
  gamma = 0,
  writeFiles = FALSE,
  showPlot = TRUE,
  parallel = FALSE
)
sum(round(solution$aggr_strata$SOLUZ))


ss <- summaryStrata(x=solution$framenew,outstrata=solution$aggr_strata)
# sink("./output/ssmulti.txt")
ss
# sink()
sum(ss$Allocation)
framenew <- solution$framenew
framenew$lon <- meuse.grid@coords[,1]
framenew$lat <- meuse.grid@coords[,2]
# framenew$strat <- paste(framenew$DOMAINVALUE,framenew$LABEL,sep="")
table(framenew$LABEL)

# Analysis of CVs
expected_CV(solution$aggr_strata)
# cv(Y1) cv(Y2) cv(Y3) cv(Y4)
# DOM1   0.05  0.028  0.033  0.033

# Simulation
val <- evalSolution(solution$framenew,solution$aggr_strata,nsampl=30)
val$coeff_var     
# CV1  CV2    CV3    CV4  dom
# 1 0.0173 0.01 0.0146 0.0134 DOM1

# selection of sample
sample <- selectSample(frame=framenew,outstrata=solution$aggr_strata)
head(sample)

# Plot
frameres2 <- SpatialPointsDataFrame(data=framenew, coords=cbind(framenew$LON,framenew$LAT) )
frameres3 <- SpatialPixelsDataFrame(points=frameres2[c("LON","LAT")], data=framenew)
frameres3$LABEL <- as.factor(frameres2$LABEL)
spplot(frameres3,c("LABEL"), col.regions=bpy.colors(5))

# save(solution,file="solution_meuse_multi_5strata.RData")
