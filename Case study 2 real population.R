#--------------------------------------------------------------------
# "Optimization of Spatial Sampling with the R packageSamplingStrata"
# Case study 2
# R script for real population (foreign residents : variable ST1)
#--------------------------------------------------------------------
library(SamplingStrata)
# library(MASS) # ginv
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
spoints <- comune_BO_geo
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
solution0 <- optimStrata (
  method="continuous",
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
# (Intercept) -0.053350   0.982107  -0.054               0.957    
# P1           0.118006   0.004238  27.844 <0.0000000000000002 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 15.41 on 465 degrees of freedom
# Multiple R-squared:  0.6251,	Adjusted R-squared:  0.6243 



# plot(lm_1)


# Breusch-Pagan test on heteroscedasticity
library(lmtest)
bptest(ST1 ~ P1,data=spoints_samp)
# Heteroscedasticity index
gamma_sigma_1 <- computeGamma(e=summary(lm_1)$residuals,
                              x=spoints_samp@data$P1,
                              nbins=6)
gamma_sigma_1
# gamma     sigma  r.square 
# 0.6146480 0.6210120 0.9662238 

model <- NULL
model$type[1] <- "linear"
model$beta[1] <- summary(lm_1)$coefficients[2]
model$gamma[1] <- gamma_sigma_1[1]
model$sig2[1] <- gamma_sigma_1[2]^2
model <- as.data.frame(model)
model
# type      beta    gamma      sig2
# 1 linear 0.1180064 0.614648 0.3856559

set.seed(1234)
solution1 <- optimStrata (
  method="continuous",
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
# CV1    CV2  dom
# 1 0.0094 0.029 DOM1

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
                           model = c("Exp", "Sph", "Mat" ))
                           # model = c("Exp", "Sph"))
                           
plot(v, fit.vgm$var_model)
fit.vgm$var_model
# model    psill    range kappa
# 1   Nug 114.3363    0.000   0.0
# 2   Mat 135.2903 1232.796   0.2

# prediction with gstat

g <- gstat(NULL, "v", ST1 ~ P1, spoints_samp)
v <- variogram(g)
v.fit <- fit.lmc(v, g, 
                 vgm(psill=fit.vgm$var_model$psill[2], 
                     model=fit.vgm$var_model$model[2], 
                     range=fit.vgm$var_model$range[2], 
                     nugget=fit.vgm$var_model$psill[1]))
preds <- predict(v.fit, Comune_BO_geo)
names(preds)
# plot(preds$v.pred,preds$v.var)

# Add predicted values and residuals variance to the frame
frame1 <- frame
frame1$Y1 <- preds$v.pred
frame1$Y1 <- ifelse(frame1$Y1 < 0, 0, frame1$Y1)
frame1$X1 <- frame1$Y1
frame1$var1 <- preds$v.var
head(frame1)

# preds = krige(ST1 ~ P1, spoints_samp, Comune_BO_geo, model = fit.vgm$var_model)
# frame1 <- frame
# frame1$Y1 <- preds$var1.pred
# frame1$Y1 <- ifelse(frame1$Y1 < 0, 0, frame1$Y1)
# frame1$X1 <- frame1$Y1
# frame1$var1 <- preds$var1.var
# head(frame1)

# Compute fitting
plot(frame1$Y1[camp],(frame$Y1[camp]-frame$ST1[camp]))
lm_pred <- lm(ST1 ~ Y1,data=frame1)
summary(lm_pred)
summary(lm_pred)$sigma^2
summary(lm_pred)$r.squared
# [1] 0.6351977

# Compute heteroscedasticity index
gamma_sigma_2 <- computeGamma(e=(frame1$Y1[camp]-frame1$ST1[camp]),
                              x=frame1$P1[camp],
                              nbins=6)
gamma_sigma_2
# gamma     sigma  r.square 
# 0.5006069 0.9888317 0.9372327 

# head(frame1)

# frame1$var2 <- NULL
frame1$var1 <- gamma_sigma_2[2]^2 * frame1$P1 ^ (gamma_sigma_2[1] * 2)
head(frame1)
set.seed(4321)
solution2 <- optimStrata (
  method="spatial",
  errors=cv, 
  framesamp=frame1,
  iter = 75,
  pops = 10,
  nStrata = 5,
  fitting = summary(lm_pred)$r.squared,
  # fitting = summary(lm_1)$r.squared, 
  range = fit.vgm$var_model$range[2],
  kappa = 1,
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
framenew$Y2 <- framenew$ST1
expected_CV(outstrata)
unlink("./simulation",recursive=TRUE)
val2 <- evalSolution(framenew,outstrata,nsampl = 1000,progress = F)
val2$coeff_var 
# CV1    CV2  dom
# 1 0.0105 0.0321 DOM1

# Comparison with same sample size of Solution 1
size <- sum(solution1$aggr_strata$SOLUZ)
newstrata <- adjustSize(size,outstrata)
sum(newstrata$SOLUZ)
newstrata
unlink("./simulation",recursive=TRUE)
val2b <- evalSolution(framenew,newstrata,nsampl = 1000,progress = F)
val2b$coeff_var 

# Plot
colnames(framenew)[1] <- "SEZ"
bologna <- Comune_BO_geo
bologna@data <- merge(bologna@data,framenew[,c("SEZ","LABEL")])
bologna@data$LABEL <- as.factor(bologna@data$LABEL)
spplot(bologna,"LABEL")

s <- selectSample(framenew,outstrata)
bologna <- merge(bologna,s[,c("SEZ","LABEL","WEIGHTS")],all.x=TRUE)
head(bologna@data)
bologna@data$sampled <- as.factor(ifelse(is.na(bologna@data$WEIGHTS),2,1))
table(bologna@data$sampled)
spplot(bologna,"sampled")

bologna@data[bologna@data$SEZ==2139,]

################################################################################################
# SOLUTION 3 
# Spatial Linear Model

lm_2 <- lm(ST1 ~ P1 + P1W, data=spoints_samp)
summary(lm_2)
# Coefficients:
#   Estimate Std. Error t value            Pr(>|t|)    
# (Intercept)  0.231107   1.457822   0.159               0.874    
# P1           0.118487   0.004616  25.667 <0.0000000000000002 ***
#   P1W         -0.002018   0.007638  -0.264               0.792    
# Residual standard error: 15.42 on 464 degrees of freedom
# Multiple R-squared:  0.6251,	Adjusted R-squared:  0.6235 

# plot(lm_2)

# Compute heteroscedasticity index
gamma_sigma_3 <- computeGamma(e=(frame$ST1[camp] - predict(lm_2,data=frame[camp,])),
                              x=spoints_samp@data$P1,
                              nbins=6)
gamma_sigma_3
# gamma     sigma  r.square 
# 0.6132473 0.6257764 0.9654370 
# Estimate psill and range on residuals of lm_2
spoints_samp@data$fit_spatial <- predict(lm_2,spoints_samp@data)
spoints_samp@data$res_spatial <- summary(lm_2)$residuals

v2 <- variogram(res_spatial  ~ 1, data=spoints_samp, cutoff=3000, width=3000/30)
plot(v2)
fit.vgm2 = autofitVariogram(res_spatial  ~ 1, spoints_samp, model = c("Exp","Sph","Mat"))
plot(v2, fit.vgm2$var_model)
fit.vgm2$var_model
# model    psill    range kappa
# 1   Nug 114.3825    0.000   0.0
# 2   Mat 135.2430 1232.444   0.2
# fit.vgm$var_model

# v3 <- variogram(fit_spatial  ~ res_spatial, data=spoints_samp, cutoff=3000, width=3000/30)
# plot(v3)
# fit.vgm3 = autofitVariogram(fit_spatial  ~ res_spatial, spoints_samp, model = c("Exp","Sph"))
# plot(v3, fit.vgm3$var_model)
# fit.vgm3$var_model



model <- NULL
model$type[1] <- "spatial"
model$beta[1] <- summary(lm_2)$coefficients[2]
model$beta2[1] <- summary(lm_2)$coefficients[3]
# model$sig2[1] <- fit.vgm2$var_model$psill[2]
model$sig2[1] <- gamma_sigma_3[2]^2
model$range[1] <- fit.vgm2$var_model$range[2]
# model$sig2_2[1] <- fit.vgm3$var_model$psill[2]
# model$range_2[1] <- fit.vgm3$var_model$range[2]
model$gamma[1] <- gamma_sigma_3[1]
# model$gamma[1] <- 0
model$fitting[1] <- summary(lm_2)$r.square
model <- as.data.frame(model)
model
# type      beta        beta2      sig2    range     gamma   fitting
# 1 spatial 0.1184873 -0.002018293 0.3915961 1232.444 0.6132473 0.6251384

set.seed(1234)
solution3 <- optimizeStrata2 (
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
# CV1    CV2  dom
# 1 0.0094 0.0302 DOM1

# Comparison with same sample size of Solution 1
size <- sum(solution1$aggr_strata$SOLUZ)
newstrata <- adjustSize(size,outstrata)
sum(newstrata$SOLUZ)
newstrata
unlink("./simulation",recursive=TRUE)
val3a <- evalSolution(framenew,newstrata,nsampl = 1000,progress = F)
val3a$rel_bias
val3a$coeff_var 


# Plot
colnames(framenew)[1] <- "SEZ"
bologna <- Comune_BO_geo
bologna@data <- merge(bologna@data,framenew[,c("SEZ","LABEL")])
bologna@data$LABEL <- as.factor(bologna@data$LABEL)
spplot(bologna,"LABEL")

####################################################################
sink("reportRealPop.txt")
# cat("\n ---------------------------------------------\n")
# cat("\n Report on real population (foreign residents)\n")
# cat("\n ---------------------------------------------\n")
# cat("\n")
# cat("\n *** No model (Y1 = ST1) ***")
# cat("\nSample size",sum(round(solution0$aggr_strata$SOLUZ)))
# cat("\nStrata structure\n")
# s0
# cat("\n  CV(ST1)\n")
# val0$coeff_var
# cat("\n")
# cat("\n")
cat("\n *** Linear model (gamma/sigma = ",gamma_sigma_1,") ***")
cat("\nR squared: ",summary(lm_1)$r.squared)
cat("\nSample size",sum(round(solution1$aggr_strata$SOLUZ)))
cat("\nStrata structure\n")
s1
cat("\n  CV(P1)  CV(ST1)\n")
val1$coeff_var
# cat("\n")
# cat("\n")
# cat("\n *** Kriging (gamma/sigma = ",gamma_sigma_2,") ***")
# cat("\nR squared: ",summary(lm_pred)$r.squared)
# cat("\nSample size",sum(round(solution2$aggr_strata$SOLUZ)))
# cat("\nStrata structure\n")
# s2
# cat("\n  CV(P1)  CV(ST1)\n")
# val2$coeff_var
# cat("... with the same sample size than Linear Model")
# cat("\n  CV(P1)  CV(ST1)\n")
# val2a$coeff_var
# cat("\n")
# cat("\n")
cat("\n *** Kriging (gamma/sigma = ",gamma_sigma_2,") ***")
cat("\nR squared: ",summary(lm_pred)$r.squared)
cat("\nSample size",sum(round(solution2$aggr_strata$SOLUZ)))
cat("\nStrata structure\n")
s2
cat("\n  CV(P1)  CV(ST1)\n")
val2$coeff_var
cat("... with the same sample size than Linear Model")
cat("\n  CV(P1)  CV(ST1)\n")
val2b$coeff_var
# cat("\n")
# cat("\n")
cat("\n *** Spatial model (gamma/sigma = ",gamma_sigma_3,") ***")
cat("\nR squared: ",summary(lm_2)$r.squared)
cat("\nSample size",sum(round(solution3$aggr_strata$SOLUZ)))
cat("\nStrata structure\n")
s3
cat("\n  CV(P1)  CV(ST1)\n")
val3$coeff_var
cat("... with the same sample size than Linear Model")
cat("\n  CV(P1)  CV(ST1)\n")
val3a$coeff_var
# cat("\n")
# cat("\n")
sink()
####################################################################

save.image(file="Bologna real population.RData")

