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
fit.vgm.cadmium = autofitVariogram(cadmium ~ dist + soil, meuse, 
                           model = c("Sph","Exp"))
plot(v, fit.vgm.cadmium$var_model)
fit.vgm.cadmium$var_model
# model    psill    range
# 1   Nug 2.821884   0.0000
# 2   Exp 5.121428 240.6374

g <- NULL
g <- gstat(g,"cadmium", cadmium ~ dist + soil, meuse)
g
vm <- variogram(g)
vm.fit <- fit.lmc(vm, g, vgm(psill=fit.vgm.cadmium$var_model$psill[2], 
                             model=fit.vgm.cadmium$var_model$model[2], 
                             range=fit.vgm.cadmium$var_model$range[2], 
                             nugget=fit.vgm.cadmium$var_model$psill[1]))
# Prediction on the whole grid
preds.cadmium <- predict(vm.fit, meuse.grid)
names(preds.cadmium)

# Copper
v <- variogram(copper ~ dist + soil, data=meuse, cutoff=3000, width=3000/30)
plot(v)
fit.vgm.copper = autofitVariogram(copper ~ dist + soil, meuse,  
                           model = c("Sph","Exp"))
plot(v, fit.vgm.copper$var_model)
fit.vgm.copper$var_model
# model    psill   range
# 1   Nug   0.0000  0.0000
# 2   Exp 339.6246 98.4095

g <- NULL
g <- gstat(g,"copper", copper ~ dist + soil, meuse)
g
vm <- variogram(g)
vm.fit <- fit.lmc(vm, g, vgm(psill=fit.vgm.copper$var_model$psill[2], 
                             model=fit.vgm.copper$var_model$model[2], 
                             range=fit.vgm.copper$var_model$range[2], 
                             nugget=fit.vgm.copper$var_model$psill[1]))
# Prediction on the whole grid
preds.copper <- predict(vm.fit, meuse.grid)
names(preds.copper)

# Lead
v <- variogram(lead ~ dist + soil, data=meuse, cutoff=3000, width=3000/30)
plot(v)
fit.vgm.lead = autofitVariogram(lead ~ dist + soil, meuse,   
                           model = c("Sph","Exp"))
plot(v, fit.vgm.lead$var_model)
fit.vgm.lead$var_model
# model    psill    range
# 1   Nug 3738.507    0.000
# 2   Sph 6243.891 1020.692

g <- NULL
g <- gstat(g,"lead", lead ~ dist + soil, meuse)
g
vm <- variogram(g)
vm.fit <- fit.lmc(vm, g, vgm(psill=fit.vgm.lead$var_model$psill[2], 
                             model=fit.vgm.lead$var_model$model[2], 
                             range=fit.vgm.lead$var_model$range[2], 
                             nugget=fit.vgm.lead$var_model$psill[1]))
# Prediction on the whole grid
preds.lead <- predict(vm.fit, meuse.grid)
names(preds.lead)

# Zinc
v <- variogram(zinc ~ dist + soil, data=meuse, cutoff=3000, width=3000/30)
plot(v)
fit.vgm.zinc = autofitVariogram(zinc ~ dist + soil, meuse,   
                           model = c("Sph","Exp"))
plot(v, fit.vgm.zinc$var_model)
fit.vgm.zinc$var_model
# model    psill    range
# 1   Nug 25328.92   0.0000
# 2   Exp 67450.62 400.5704

g <- NULL
g <- gstat(g,"zinc", zinc ~ dist + soil, meuse)
g
vm <- variogram(g)
vm.fit <- fit.lmc(vm, g, vgm(psill=fit.vgm.zinc$var_model$psill[2], 
                             model=fit.vgm.zinc$var_model$model[2], 
                             range=fit.vgm.zinc$var_model$range[2], 
                             nugget=fit.vgm.zinc$var_model$psill[1]))
# Prediction on the whole grid
preds.zinc <- predict(vm.fit, meuse.grid)
names(preds.zinc)


# multivariate case:cadmium, copper, lead, zinc as target variables
# g <- NULL
# g <- gstat(NULL, "Cd", cadmium ~ dist + soil, meuse)
# g <- gstat(g, "Cu", copper ~ dist + soil, meuse)
# g <- gstat(g, "Pb", lead ~ dist + soil, meuse)
# g <- gstat(g, "Zn", zinc ~ dist + soil, meuse)
# vm <- variogram(g)
# # Using "Sph" because "Exp" does not work
# vm.fit <- fit.lmc(vm, g, vgm(psill=1, model=c("Sph"), range=800, nugget=1))
# vm.fit
# plot(vm, vm.fit)
# preds <- predict(vm.fit, meuse.grid)
# names(preds)
# # [1] "Cd.pred"   "Cd.var"    "Cu.pred"   "Cu.var"    "Pb.pred"  
# # [6] "Pb.var"    "Zn.pred"   "Zn.var"    "cov.Cd.Cu" "cov.Cd.Pb"
# # [11] "cov.Cu.Pb" "cov.Cd.Zn" "cov.Cu.Zn" "cov.Pb.Zn"

preds.cadmium$cadmium.pred <- ifelse(preds.cadmium$cadmium.pred < 0,0,preds.cadmium$cadmium.pred)
preds.copper$copper.pred <- ifelse(preds.copper$copper.pred < 0,0,preds.copper$copper.pred)
preds.lead$lead.pred <- ifelse(preds.lead$lead.pred < 0,0,preds.lead$lead.pred)
preds.zinc$zinc.pred <- ifelse(preds.zinc$zinc.pred < 0,0,preds.zinc$zinc.pred)
preds.cadmium$cadmium.var <- ifelse(preds.cadmium$cadmium.var < 0,0,preds.cadmium$cadmium.var)
preds.copper$copper.var <- ifelse(preds.copper$copper.var < 0,0,preds.copper$copper.var)
preds.lead$lead.var <- ifelse(preds.lead$lead.var < 0,0,preds.lead$lead.var)
preds.zinc$zinc.var <- ifelse(preds.zinc$zinc.var < 0,0,preds.zinc$zinc.var)

# Optimization with SamplingStrata

df <- NULL
df$cadmium.pred <- preds.cadmium$cadmium.pred
df$cadmium.var <- preds.cadmium$cadmium.var
df$copper.pred <- preds.copper$copper.pred
df$copper.var <- preds.copper$copper.var
df$lead.pred <- preds.lead$lead.pred
df$lead.var <- preds.lead$lead.var
df$zinc.pred <- preds.zinc$zinc.pred
df$zinc.var <- preds.zinc$zinc.var
<<<<<<< HEAD
df$lon <- meuse.grid$x
df$lat <- meuse.grid$y
=======
# df$dom <- meuse.grid@data$soil
>>>>>>> b72ad61ae26897bf65c1d7468f9d6846dcb31bed
df$dom <- 1
df <- as.data.frame(df)
df$id <- meuse.grid$id

head(df)
# cadmium.pred cadmium.var copper.pred copper.var lead.pred lead.var zinc.pred zinc.var dom id
# 7.917987    7.007012    65.64426   331.6897  265.9674 7043.944  898.3354 67350.56   1  1
# 8.308433    6.482721    69.01170   308.4855  267.2160 6403.047  918.6911 59349.98   1  2
# 7.898212    6.661317    66.93037   319.6262  260.5370 6574.794  877.8272 61860.84   1  3
# 7.281627    6.828875    63.61440   326.8032  248.9653 6770.737  818.0511 64442.07   1  4
# 8.824362    5.748581    74.34104   245.6440  269.2219 5744.181  949.7147 49761.03   1  5
# 8.333077    6.009188    71.27280   278.8897  261.5889 5909.198  899.6006 52859.14   1  6

cv <- as.data.frame(list(DOM=rep("DOM1",1),
                         CV1=rep(0.05,1),
                         CV2=rep(0.05,1),
                         CV3=rep(0.05,1),
                         CV4=rep(0.05,1),
                         domainvalue=c(1:1) ))
cv
# DOM  CV1  CV2  CV3  CV4 domainvalue
# 1 DOM1 0.05 0.05 0.05 0.05           1
<<<<<<< HEAD
frame <- buildFrameSpatial(df=df, 
                           id="id", 
                           X=c("cadmium.pred","copper.pred","lead.pred","zinc.pred"),
                           Y=c("cadmium.pred","copper.pred","lead.pred","zinc.pred"),
                           variance=c("cadmium.var","copper.var","lead.var","zinc.var"), 
                           lon="lon", 
                           lat="lat", 
                           domainvalue="dom")
=======
frame <- buildFrameDF(df=df,
                      id="id",
                      X=c(
                        "cadmium.pred",
                        "copper.pred",
                        "lead.pred",
                        "zinc.pred"),
                      Y=c(
                        "cadmium.pred",
                        "copper.pred",
                        "lead.pred",
                        "zinc.pred"),
                      domainvalue = "dom")
frame$lon <- meuse.grid$x
frame$lat <- meuse.grid$y
frame$var1 <- df$cadmium.var
frame$var2 <- df$copper.var
frame$var3 <- df$lead.var
frame$var4 <- df$zinc.var
>>>>>>> b72ad61ae26897bf65c1d7468f9d6846dcb31bed

range <- c(fit.vgm.cadmium$var_model$range[2],
           fit.vgm.copper$var_model$range[2],
           fit.vgm.lead$var_model$range[2],
           fit.vgm.zinc$var_model$range[2])

kmeans <- KmeansSolutionSpatial(frame,
                                fitting=c(1,1,1,1),
                                range=range,
                                kappa=1,
                                errors=cv,
                                maxclusters = 10)
# The difference in sample size after nStrata = 5 is negligible
# so we choose nStrata = 5

#---------------------------------------------------
# Solution with spatial optimization on continuous variables
#---------------------------------------------------

set.seed(4321)
solution <- optimStrata (
  method="spatial",
  errors=cv, 
  framesamp=frame,
  iter = 50,
  pops = 10,
  nStrata = 5,
  fitting = c(1,1,1,1),
  range = range,
  kappa = 1,
  writeFiles = FALSE,
  showPlot = TRUE,
  parallel = FALSE
)
sum(round(solution$aggr_strata$SOLUZ))


ss <- summaryStrata(x=solution$framenew,outstrata=solution$aggr_strata)
# sink("ssmulti.txt")
ss
# sink()
sum(ss$Allocation)
framenew <- solution$framenew
# framenew$lon <- meuse.grid@coords[,1]
# framenew$lat <- meuse.grid@coords[,2]
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
