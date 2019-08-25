#--------------------------------------------------------------------
# "Optimization of Spatial Sampling with the R packageSamplingStrata"
# Case study 2
# R script for simulated population
#--------------------------------------------------------------------
library(SamplingStrata)
library(rgdal)
library(spdep)
library(gstat)
library(automap)
library(rgeos)
library(geoR)
library(spatstat)
library(automap)
library(spatialreg)

#--------------------------------------------------------------
# PREPARATION OF DATA

# Read geographic data
dsn <- "./Bologna_EA_data/Bologna_EA_2011.shp"
Comune_BO_geo_orig<-readOGR(dsn=dsn,layer="Bologna_EA_2011")

# Subset of variables
Comune_BO_geo<-Comune_BO_geo_orig[,c("SEZ","P1","ST1")]
Comune_BO_geo@data$SEZ<-as.character(Comune_BO_geo@data$SEZ)
Comune_BO_geo@data$ST1<-ifelse(is.na(Comune_BO_geo@data$ST1),0,Comune_BO_geo@data$ST1)
Comune_BO_geo@data$P1<-ifelse(is.na(Comune_BO_geo@data$P1),0,Comune_BO_geo@data$P1)

# Centroids
grd<-gCentroid(Comune_BO_geo,byid=TRUE)

# Coordinates
cc<-as.data.frame(coordinates(grd))
coordinates(cc)<-cc[,c("x","y")]
# Spatial points dataframe
spoints<-SpatialPointsDataFrame(cc,Comune_BO_geo@data[,c("SEZ","P1","ST1")])
spoints@data$id <-c(1:nrow(spoints@data)) 

# Build contiguity matrix W and distance matrix d
contnb<-poly2nb(Comune_BO_geo, "queen"=TRUE)
# Two formats of contiguity matrix: 
#  - WM to transform variables 
#  - WL to interpolate
WM<-nb2mat(contnb)
WL<-nb2listw(contnb)
# Distance matrix
d<-spDists(x=Comune_BO_geo, y = Comune_BO_geo, longlat = FALSE, segments = FALSE, diagonal = FALSE)

# Build spatialized variable
P1W<-WM%*%Comune_BO_geo@data$P1
Comune_BO_geo@data$P1W<-P1W[,1]

# Errors generation
# (with range = 0 errors are not spatially correlated)
var_eps=2000
range_eps=0
beta1=1
beta2=1
set.seed(1234)
eps1 <- grf(nrow(spoints@data),  
            grid=spoints@coords, 
            cov.pars=c(var_eps, range_eps),           
            cov.model="exponential") 
head(eps1$data)
# if range = 0 use this:
Comune_BO_geo@data$eps1<-eps1$data[,1]
# if range > 0 use this:
# Comune_BO_geo@data$eps1<-eps1$data

# Create the variable of interest (target)
Comune_BO_geo@data$target <- beta1*Comune_BO_geo@data$P1+
                             beta2*Comune_BO_geo@data$P1W+
                             Comune_BO_geo@data$eps1
summary(Comune_BO_geo@data$target)

# Plot 
# spplot(Comune_BO_geo,"P1")
# spplot(Comune_BO_geo,"P1W")
# spplot(Comune_BO_geo,"Y1")
# spplot(Comune_BO_geo,"eps1")

# Moran test
lm_target <-lm(Comune_BO_geo@data$target~Comune_BO_geo@data$P1)
lm_target$coefficients
summary(lm_target)
moran.test(lm_target$residuals,listw = WL)

######################################################
# Selection of a sample of EAs
sample_rate <- 0.2
samplesize <- round(nrow(Comune_BO_geo)*sample_rate)
set.seed(1234)
camp <- sample(c(1:nrow(Comune_BO_geo)),samplesize)
# camp <- spoints
spoints_samp<-Comune_BO_geo[camp,]

# Frame 
df <- NULL
df$id <- Comune_BO_geo@data$SEZ
df$P1 <- Comune_BO_geo@data$P1
df <- as.data.frame(df)
df$dom <- 1

frame <- buildFrameDF(df=df,
                      id="id",
                      X="P1",
                      Y="P1",
                      domainvalue = "dom")
frame$lon <- cc@coords[,1]
frame$lat <- cc@coords[,2]
frame$target <- Comune_BO_geo$target
frame$P1 <- Comune_BO_geo$P1 
frame$W1 <- Comune_BO_geo$P1W

cv <- NULL
cv$DOM <- "DOM1"
cv$CV1 <- 0.03
cv$domainvalue <- 1
cv <- as.data.frame(cv)
cv

################################################################################################
# SOLUTION 1
# Linear model
lm_1 <- lm(target~P1,data=spoints_samp)
summary(lm_1)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 126.78536    6.65885   19.04   <2e-16 ***
#   P1            1.32323    0.02863   46.22   <2e-16 ***
#   ---
# Residual standard error: 99.83 on 465 degrees of freedom
# Multiple R-squared:  0.8212,	Adjusted R-squared:  0.8208 
# F-statistic:  2136 on 1 and 465 DF,  p-value: < 2.2e-16

plot(lm_1)

model <- NULL
model$type[1] <- "linear"
model$beta[1] <- summary(lm_1)$coefficients[2]
model$sig2[1] <- summary(lm_1)$sigma
model$gamma[1] <- 0
model <- as.data.frame(model)
model
# type     beta     sig2         gamma
# linear   1.323232 99.83204     0
set.seed(1234)
solution1 <- optimizeStrata2 (
  errors=cv, 
  framesamp=frame,
  model=model,
  iter = 50,
  pops = 10,
  nStrata = 5,
  writeFiles = FALSE,
  showPlot = TRUE,
  parallel = FALSE
)
sum(solution1$aggr_strata$SOLUZ)
framenew <- solution1$framenew
outstrata <- solution1$aggr_strata
outstrata
s1 <- summaryStrata(framenew,outstrata)
s1
framenew$Y2 <- framenew$TARGET

expected_CV(outstrata)
unlink("./simulation",recursive=TRUE)
val1 <- evalSolution(framenew,outstrata,nsampl=1000,progress=F)
val1$rel_bias
val1$coeff_var 

# Plot
colnames(framenew)[1] <- "SEZ"
bologna <- Comune_BO_geo
bologna@data <- merge(bologna@data,framenew[,c("SEZ","LABEL")])
bologna@data$LABEL <- as.factor(bologna@data$LABEL)
spplot(bologna,"LABEL")

################################################################################################
# SOLUTION 2 
# Universal kriging

# variogram
v <- variogram(target ~ P1, data=spoints_samp, cutoff=3000, width=3000/30)
plot(v)
# Estimation of psill, range and nugget with automap
fit.vgm = autofitVariogram(target ~ P1, spoints_samp, 
                           model = c("Exp", "Sph" ))
plot(v, fit.vgm$var_model)
fit.vgm$var_model
# model    psill    range
# 1   Nug 5753.814   0.0000
# 2   Exp 4041.294 958.0816

# prediction with gstat

g <- gstat(NULL, "v", target ~ P1, spoints_samp)
v <- variogram(g)
v.fit <- fit.lmc(v, g, 
                 vgm(psill=fit.vgm$var_model$psill[2], 
                     model="Exp", 
                     range=fit.vgm$var_model$range[2], 
                     nugget=fit.vgm$var_model$psill[1]))
preds <- predict(v.fit, Comune_BO_geo)

# Add estimated model variance to frame
frame$var1 <- preds$v.var

set.seed(1234)
solution2 <- optimizeStrataSpatial (
  errors=cv, 
  framesamp=frame,
  iter = 50,
  pops = 10,
  nStrata = 5,
  fitting = summary(lm_1)$r.squared, # 0.8212337 
  range = fit.vgm$var_model$range[2],
  gamma = 3,
  writeFiles = FALSE,
  showPlot = TRUE,
  parallel = FALSE
)
sum(round(solution2$aggr_strata$SOLUZ))
framenew <- solution2$framenew
outstrata <- solution2$aggr_strata
outstrata
s2 <- summaryStrata(framenew,outstrata)
s2
framenew$Y2 <- framenew$TARGET
expected_CV(outstrata)
unlink("./simulation",recursive=TRUE)
val2 <- evalSolution(framenew,outstrata,nsampl = 1000,progress = F)
val2$rel_bias
val2$coeff_var 

# Plot
colnames(framenew)[1] <- "SEZ"
bologna <- Comune_BO_geo
bologna@data <- merge(bologna@data,framenew[,c("SEZ","LABEL")])
bologna@data$LABEL <- as.factor(bologna@data$LABEL)
spplot(bologna,"LABEL")

################################################################################################
# SOLUTION 3 
# Spatial Linear Model

lm_2 <- lm(target ~ P1 + W1, data=frame[camp,])
summary(lm_2)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -5.24287    4.28806  -1.223    0.222    
# P1           0.99936    0.01490  67.067   <2e-16 ***
#   W1         1.03155    0.02402  42.940   <2e-16 ***
#   ---
# Residual standard error: 44.81 on 464 degrees of freedom
# Multiple R-squared:  0.9641,	Adjusted R-squared:  0.9639 



# plot(lm_2)

model <- NULL
model$type[1] <- "spatial"
model$beta[1] <- summary(lm_2)$coefficients[2]
model$beta2[1] <- summary(lm_2)$coefficients[3]
model$sig2[1] <- fit.vgm$var_model$psill[2]
model$range[1] <- fit.vgm$var_model$range[2]
model$gamma[1] <- NA
model <- as.data.frame(model)
model
# type      beta    beta2     sig2    range       gamma
# spatial 0.9993553 1.031549 4041.294 958.0816    NA
set.seed(1234)
solution3 <- optimizeStrata2 (
  errors=cv, 
  framesamp=frame,
  model=model,
  iter = 50,
  pops = 10,
  nStrata = 5,
  writeFiles = FALSE,
  showPlot = TRUE,
  parallel = FALSE
)
sum(solution3$aggr_strata$SOLUZ)
framenew <- solution3$framenew
outstrata <- solution3$aggr_strata
s3 <- summaryStrata(framenew,outstrata)
s3
framenew$Y2 <- framenew$TARGET

expected_CV(outstrata)
unlink("./simulation",recursive=TRUE)
val3 <- evalSolution(framenew,outstrata,nsampl=1000,progress=F)
val1$rel_bias
val3$coeff_var 

# Plot
colnames(framenew)[1] <- "SEZ"
bologna <- Comune_BO_geo
bologna@data <- merge(bologna@data,framenew[,c("SEZ","LABEL")])
bologna@data$LABEL <- as.factor(bologna@data$LABEL)
spplot(bologna,"LABEL")

####################################################################
sink("report_simulated_population.txt")
cat("\n ------------------------------\n")
cat("\n Report on simulated population\n")
cat("\n ------------------------------\n")
cat("\n *** Linear model ***")
cat("\nSample size",sum(solution1$aggr_strata$SOLUZ))
cat("\nStrata structure\n")
s1
cat("\n  CV(P1)  CV(target)\n")
val1$coeff_var
cat("\n")
cat("\n")
cat("\n *** Kriging ***")
cat("\nSample size",sum(solution2$aggr_strata$SOLUZ))
cat("\nStrata structure\n")
s2
cat("\n  CV(P1)  CV(target)\n")
val2$coeff_var
cat("\n")
cat("\n")
cat("\n *** Spatial model ***")
cat("\nSample size",sum(solution3$aggr_strata$SOLUZ))
cat("\nStrata structure\n")
s3
cat("\n  CV(P1)  CV(target)\n")
val3$coeff_var
cat("\n")
cat("\n")
sink()
####################################################################

save.image(file="Bologna simulated population.RData")

