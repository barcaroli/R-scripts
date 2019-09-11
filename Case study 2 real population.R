#--------------------------------------------------------------------
# "Optimization of Spatial Sampling with the R packageSamplingStrata"
# Case study 2
# R script for real population (foreign residents : variable ST1)
#--------------------------------------------------------------------
library(SamplingStrata)
library(rgdal) # readOGR
library(spdep) # poly2nb
library(gstat)
library(automap)
library(rgeos) # gCentroid

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




# Moran test 
lm_ST1 <-lm(Comune_BO_geo@data$ST1 ~ Comune_BO_geo@data$P1)
lm_ST1$coefficients
summary(lm_ST1)
moran.test(lm_ST1$residuals,listw = WL)

df <- NULL
df$id <- Comune_BO_geo@data$SEZ
df$P1 <- Comune_BO_geo@data$P1
df$ST1 <- Comune_BO_geo@data$ST1
df <- as.data.frame(df)
df$dom <- 1



######################################################
# Selection of a sample of EAs
sample_rate <- 0.2
samplesize <- round(nrow(Comune_BO_geo)*sample_rate)
set.seed(1234)
camp <- sample(c(1:nrow(Comune_BO_geo)),samplesize)
# camp <- spoints
spoints_samp<-Comune_BO_geo[camp,]


cv <- NULL
cv$DOM <- "DOM1"
cv$CV1 <- 0.03
cv$domainvalue <- 1
cv <- as.data.frame(cv)
cv

################################################################################################
# SOLUTION 0
# Linear model

frame <- buildFrameDF(df=df,
                      id="id",
                      X="ST1",
                      Y="ST1",
                      domainvalue = "dom")

set.seed(1234)
solution0 <- optimizeStrata2 (
  errors=cv, 
  framesamp=frame,
  model=NULL,
  iter = 50,
  pops = 10,
  nStrata = 5,
  writeFiles = FALSE,
  showPlot = TRUE,
  parallel = FALSE
)
sum(round(solution0$aggr_strata$SOLUZ))
framenew <- solution0$framenew
outstrata <- solution0$aggr_strata
outstrata
s0 <- summaryStrata(framenew,outstrata)
s0


expected_CV(outstrata)
unlink("./simulation",recursive=TRUE)
val0 <- evalSolution(framenew,outstrata,nsampl = 1000,progress=F)
val0$rel_bias
val0$coeff_var 

# Plot
colnames(framenew)[1] <- "SEZ"
bologna <- Comune_BO_geo
bologna@data <- merge(bologna@data,framenew[,c("SEZ","LABEL")])
bologna@data$LABEL <- as.factor(bologna@data$LABEL)
spplot(bologna,"LABEL")


################################################################################################
# SOLUTION 1
# Linear model

frame <- buildFrameDF(df=df,
                      id="id",
                      X="P1",
                      Y="P1",
                      domainvalue = "dom")
frame$lon <- cc@coords[,1]
frame$lat <- cc@coords[,2]
frame$ST1 <- Comune_BO_geo$ST1
frame$P1 <- Comune_BO_geo$P1 
frame$W1 <- Comune_BO_geo$P1W


lm_1 <- lm(ST1 ~ P1,data=spoints_samp)
summary(lm_1)

# Coefficients:
#   Estimate Std. Error t value            Pr(>|t|)    
# (Intercept)  0.40386    1.08621   0.372                0.71    
# P1           0.11371    0.00467  24.349 <0.0000000000000002 ***
#   ---
# Residual standard error: 16.28 on 465 degrees of freedom
# Multiple R-squared:  0.5604,	Adjusted R-squared:  0.5595 


# plot(lm_1)

# Heteroscedasticity index
df1 <- NULL
df1$x <- spoints_samp@data$P1
df1$p <- predict(lm_1,spoints_samp@data)
df1$e <- summary(lm_1)$residuals
df1 <- as.data.frame(df1)
source("computeGamma.R")
gamma_sigma_1 <- computeGamma(dataset=df1)
gamma_sigma_1

model <- NULL
model$type[1] <- "linear"
model$beta[1] <- summary(lm_1)$coefficients[2]
model$sig2[1] <- gamma_sigma_1[2]^2
model$gamma[1] <- gamma_sigma_1[1]*2
# model$gamma <- gamma_est
model <- as.data.frame(model)
model

set.seed(1234)
solution1 <- optimizeStrata2 (
  errors=cv, 
  framesamp=frame,
  model=model,
  iter = 75,
  pops = 10,
  nStrata = 5,
  writeFiles = FALSE,
  showPlot = TRUE,
  parallel = FALSE
)
sum(round(solution1$aggr_strata$SOLUZ))
framenew <- solution1$framenew
outstrata <- solution1$aggr_strata
outstrata
s1 <- summaryStrata(framenew,outstrata)
s1
framenew$Y2 <- framenew$ST1

expected_CV(outstrata)
unlink("./simulation",recursive=TRUE)
val1 <- evalSolution(framenew,outstrata,nsampl = 1000,progress=F)
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
v <- variogram(ST1 ~ P1, data=spoints_samp, cutoff=3000, width=3000/30)
plot(v)
# Estimation of psill, range and nugget with automap
fit.vgm = autofitVariogram(ST1 ~ P1, spoints_samp, 
                           model = c("Exp", "Sph" ))
plot(v, fit.vgm$var_model)
fit.vgm$var_model
# model     psill    range
# 1   Nug 219.16587   0.0000
# 2   Exp  52.69629 683.9035

# prediction with gstat

g <- gstat(NULL, "v", ST1 ~ P1, spoints_samp)
v <- variogram(g)
v.fit <- fit.lmc(v, g, 
                 vgm(psill=fit.vgm$var_model$psill[2], 
                     model="Exp", 
                     range=fit.vgm$var_model$range[2], 
                     nugget=fit.vgm$var_model$psill[1]))
preds <- predict(v.fit, Comune_BO_geo)
names(preds)

# Add estimated model variance to frame
frame$pred <- preds$v.pred
frame$var1 <- preds$v.var

# Compute fitting
plot(frame$ST1,frame$pred)
lm_pred <- lm(ST1 ~ pred,data=frame)
summary(lm_pred)
summary(lm_pred)$r.squared

# Compute heteroscedasticity index
df1 <- NULL
df1$x <- spoints_samp@data$P1
df1$e <- frame$ST1[camp] - frame$pred[camp]
df1 <- as.data.frame(df1)
source("computeGamma.R")
gamma_sigma_2 <- computeGamma(dataset=df1)
gamma_sigma_2


set.seed(1234)
solution2 <- optimizeStrataSpatial (
  errors=cv, 
  framesamp=frame,
  iter = 75,
  pops = 10,
  nStrata = 5,
  fitting = summary(lm_pred)$r.squared, # 0.6240108
  range = fit.vgm$var_model$range[2],
  kappa = 3,
  gamma = gamma_sigma_2[1]*2,
  writeFiles = FALSE,
  showPlot = TRUE,
  parallel = FALSE
)
sum(solution2$aggr_strata$SOLUZ)
framenew <- solution2$framenew
outstrata <- solution2$aggr_strata
outstrata
s2 <- summaryStrata(framenew,outstrata)
s2
framenew$Y2 <- framenew$ST1
expected_CV(outstrata)
unlink("./simulation",recursive=TRUE)
val2 <- evalSolution(framenew,outstrata,nsampl = 1000,progress = F)
val2$rel_bias
val2$coeff_var 

# Comparison with same sample size than Solution 1
size <- sum(solution1$aggr_strata$SOLUZ)
newstrata <- adjustSize(size,outstrata)
sum(newstrata$SOLUZ)
newstrata
unlink("./simulation",recursive=TRUE)
val2a <- evalSolution(framenew,newstrata,nsampl = 1000,progress = F)
val2a$rel_bias
val2a$coeff_var 

# Plot
colnames(framenew)[1] <- "SEZ"
bologna <- Comune_BO_geo
bologna@data <- merge(bologna@data,framenew[,c("SEZ","LABEL")])
bologna@data$LABEL <- as.factor(bologna@data$LABEL)
spplot(bologna,"LABEL")

################################################################################################
# SOLUTION 3 
# Spatial Linear Model

lm_2 <- lm(ST1 ~ P1 + P1W, data=spoints_samp)
summary(lm_2)
# plot(lm_2)

# Compute heteroscedasticity index
df1 <- NULL
df1$x <- spoints_samp@data$P1
df1$p <- predict(lm_2,data=frame[camp,])
df1$e <- frame$ST1[camp] - df1$p
df1 <- as.data.frame(df1)
source("computeGamma2_rev3.R")
gamma_sigma_3 <- computeGamma(dataset=df1)
gamma_sigma_3

# Estimate psill and range on residuals of lm_2
spoints_samp@data$fit_spatial <- predict(lm_2,spoints_samp@data)
spoints_samp@data$res_spatial <- summary(lm_2)$residuals

v2 <- variogram(res_spatial  ~ 1, data=spoints_samp, cutoff=3000, width=3000/30)
plot(v2)
fit.vgm2 = autofitVariogram(res_spatial  ~ 1, spoints_samp, model = c("Exp","Sph"))
plot(v2, fit.vgm2$var_model)
fit.vgm2$var_model
# fit.vgm$var_model

v3 <- variogram(fit_spatial  ~ res_spatial, data=spoints_samp, cutoff=3000, width=3000/30)
plot(v3)
fit.vgm3 = autofitVariogram(fit_spatial  ~ res_spatial, spoints_samp, model = c("Exp","Sph"))
plot(v3, fit.vgm3$var_model)
fit.vgm3$var_model



model <- NULL
model$type[1] <- "spatial"
model$beta[1] <- summary(lm_2)$coefficients[2]
model$beta2[1] <- summary(lm_2)$coefficients[3]
model$sig2[1] <- fit.vgm2$var_model$psill[2]
# model$sig2[1] <- gamma_sigma_3[2]^2
model$range[1] <- fit.vgm2$var_model$range[2]
# model$sig2_2[1] <- fit.vgm3$var_model$psill[2]
# model$range_2[1] <- fit.vgm3$var_model$range[2]
model$gamma[1] <- gamma_sigma_3[1]*2
# model$gamma[1] <- 0
model$fitting[1] <- summary(lm_2)$r.square
model <- as.data.frame(model)
model

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
framenew$Y2 <- framenew$ST1

expected_CV(outstrata)
unlink("./simulation",recursive=TRUE)
val3 <- evalSolution(framenew,outstrata,nsampl=1000,progress=F)
val3$rel_bias
val3$coeff_var 

# Comparison with same sample size of Solution 2
size <- sum(solution2$aggr_strata$SOLUZ)
newstrata <- adjustSize(size,outstrata)
sum(newstrata$SOLUZ)
newstrata
unlink("./simulation",recursive=TRUE)
val2a <- evalSolution(framenew,newstrata,nsampl = 1000,progress = F)
val2a$rel_bias
val2a$coeff_var 


# Plot
colnames(framenew)[1] <- "SEZ"
bologna <- Comune_BO_geo
bologna@data <- merge(bologna@data,framenew[,c("SEZ","LABEL")])
bologna@data$LABEL <- as.factor(bologna@data$LABEL)
spplot(bologna,"LABEL")

####################################################################
sink("report_real_population.txt")
cat("\n ---------------------------------------------\n")
cat("\n Report on real population (foreign residents)\n")
cat("\n ---------------------------------------------\n")
cat("\n")
cat("\n *** No model (Y1 = ST1) ***")
cat("\nSample size",sum(round(solution0$aggr_strata$SOLUZ)))
cat("\nStrata structure\n")
s0
cat("\n  CV(ST1)\n")
val0$coeff_var
cat("\n")
cat("\n")
cat("\n *** Linear model (gamma/sigma = ",gamma_sigma_1,") ***")
cat("\nSample size",sum(round(solution1$aggr_strata$SOLUZ)))
cat("\nStrata structure\n")
s1
cat("\n  CV(P1)  CV(ST1)\n")
val1$coeff_var
cat("\n")
cat("\n")
cat("\n *** Kriging (gamma/sigma = ",gamma_sigma_2,") ***")
cat("\nSample size",sum(round(solution2$aggr_strata$SOLUZ)))
cat("\nStrata structure\n")
s2
cat("\n  CV(P1)  CV(ST1)\n")
val2$coeff_var
cat("\n")
cat("\n")
cat("\n *** Spatial model (gamma/sigma = ",gamma_sigma_3,") ***")
cat("\nSample size",sum(round(solution3$aggr_strata$SOLUZ)))
cat("\nStrata structure\n")
s3
cat("\n  CV(P1)  CV(ST1)\n")
val3$coeff_var
cat("\n")
cat("\n")
sink()
####################################################################

save.image(file="Bologna real population.RData")

