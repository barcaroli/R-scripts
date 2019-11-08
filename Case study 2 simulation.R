#--------------------------------------------------------------------
# "Optimization of Spatial Sampling with the R packageSamplingStrata"
# Case study 2
# R script for evaluating performance of the three models
# varying conditions of generation of a simulated population
# (values of parameters in file "simulation.csv")
# Input: simulation.csv (values of parameters for each iteration)
# Output: simul_results.csv
#--------------------------------------------------------------------

filename <- "simulation_gamma.csv"

library(SamplingStrata)
library(rgdal) # readOGR
library(spdep) # poly2nb
library(gstat)
library(automap)
library(rgeos) # gCentroid
library(geoR) # grf



simula <- function(itera,
                   Comune_BO_geo,
                   spoints,
                   cc,
                   range_eps,
                   var_eps,
                   gamma,
                   beta1,
                   beta2) {
  cat("\nIteration",itera,"\n")
  set.seed(4321)
  eps1 <- grf(nrow(spoints@data),  
              grid=spoints@coords, 
              cov.pars=c(var_eps, range_eps),           
              cov.model="exponential") 
  head(eps1$data)

  Comune_BO_geo@data$eps1<-eps1$data
  Comune_BO_geo@data$target <- beta1*Comune_BO_geo@data$P1+
    beta2*Comune_BO_geo@data$P1W+
    Comune_BO_geo@data$eps1 * Comune_BO_geo@data$P1^gamma
  Comune_BO_geo@data$target <- ifelse(Comune_BO_geo@data$target < 0, 0, Comune_BO_geo@data$target)
  summary(Comune_BO_geo@data$target)


  ######################################################
  # Selection of a sample of EAs
  sample_rate <- 0.2
  samplesize <- round(nrow(Comune_BO_geo)*sample_rate)
  set.seed(4321)
  camp <- sample(c(1:nrow(Comune_BO_geo)),samplesize)
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
  
  # Heteroscedasticity index
  gamma_sigma_1 <- computeGamma(e=summary(lm_1)$residuals,
                                x=spoints_samp@data$P1,
                                nbins=6)
  gamma_sigma_1
  gamma_sigma_1[1] <- ifelse(gamma_sigma_1[1]<0,0,gamma_sigma_1[1])
  
  # plot(lm_1)
  
  model_linear <- NULL
  model_linear$type[1] <- "linear"
  model_linear$beta[1] <- summary(lm_1)$coefficients[2]
  model_linear$sig2[1] <- gamma_sigma_1[2]^2
  model_linear$gamma[1] <- gamma_sigma_1[1]
  model_linear <- as.data.frame(model_linear)
  model_linear

  set.seed(4321)
  solution1 <- optimizeStrata2 (
    errors=cv, 
    framesamp=frame,
    model=model_linear,
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

  ################################################################################################
  # SOLUTION 2 
  # Universal kriging
  
  # variogram
  v <- variogram(target ~ P1, data=spoints_samp)
  plot(v)
  # Estimation of psill, range and nugget with automap
  fit.vgm = autofitVariogram(target ~ P1, spoints_samp, 
                             # model = c("Exp","Sph","Gau","Mat","Log","Exc","Ste", "Cir", "Lin", "Bes", "Pen", "Per", "Wav", "Hol", "Log", "Pow", "Spl", "Leg", "Err", "Int"))
                             model = c("Exp","Sph","Mat"))
  plot(v, fit.vgm$var_model)
  fit.vgm$var_model

  
  # prediction with gstat
  
  g <- gstat(NULL, "v", target ~ P1, spoints_samp)
  v <- variogram(g)
  v.fit <- fit.lmc(v, g, 
                   vgm(psill=fit.vgm$var_model$psill[2], 
                       model=fit.vgm$var_model$model[2], 
                       range=fit.vgm$var_model$range[2], 
                       nugget=fit.vgm$var_model$psill[1]))
  preds <- predict(v.fit, Comune_BO_geo)
  # Add predicted values and residuals variance to the frame
  frame1 <- frame
  frame1$Y1 <- preds$v.pred
  frame1$Y1 <- ifelse(frame1$Y1 < 0, 0, frame1$Y1)
  frame1$X1 <- frame1$Y1
  frame1$var1 <- preds$v.var
  
  
  gamma_sigma_2 <- computeGamma(e=(frame1$Y1[camp]-frame1$target[camp]),
                                x=frame1$P1[camp],
                                nbins=6)
  gamma_sigma_2
  
  # Compute fitting
  # plot(frame1$Y1[camp],(frame1$Y1[camp]-frame1$target[camp]))
  lm_pred <- lm(target ~ Y1,data=frame1)
  summary(lm_pred)
  summary(lm_pred)$sigma^2
  summary(lm_pred)$r.squared

  frame1$var1 <- gamma_sigma_2[2]^2 * frame1$P1 ^ (gamma_sigma_2[1] * 2)
  
  
  set.seed(4321)
  solution2 <- optimizeStrataSpatial (
    errors=cv, 
    framesamp=frame1,
    iter = 50,
    pops = 10,
    nStrata = 5,
    # fitting = summary(lm_pred)$r.squared,
    fitting = summary(lm_1)$r.squared, 
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
  framenew$Y2 <- framenew$TARGET
  expected_CV(outstrata)
  unlink("./simulation",recursive=TRUE)
  val2 <- evalSolution(framenew,outstrata,nsampl = 1000,progress = F)
  val2$rel_bias
  val2$coeff_var 
  
  # Comparison with same sample size of Solution 1
  size <- sum(solution1$aggr_strata$SOLUZ)
  newstrata <- adjustSize(size,outstrata)
  sum(newstrata$SOLUZ)
  newstrata
  unlink("./simulation",recursive=TRUE)
  val2a <- evalSolution(framenew,newstrata,nsampl = 1000,progress = F)
  val2a$coeff_var 
  
  ################################################################################################
  # SOLUTION 3 
  # Spatial Linear Model
  
  lm_2 <- lm(target ~ P1 + P1W, data=spoints_samp)
  summary(lm_2)
  
  # plot(lm_2)
  
  gamma_sigma_3 <- computeGamma(e=(frame$target[camp] - predict(lm_2,data=frame[camp,])),
                                x=spoints_samp@data$P1,
                                nbins=6)
  gamma_sigma_3
  # Estimate psill and range on residuals of lm_2
  spoints_samp@data$fit_spatial <- predict(lm_2,spoints_samp@data)
  spoints_samp@data$res_spatial <- summary(lm_2)$residuals
  
  v2 <- variogram(res_spatial  ~ 1, data=spoints_samp)
  plot(v2)
  fit.vgm2 = autofitVariogram(res_spatial  ~ 1, spoints_samp, model = c("Exp","Sph","Mat"))
  plot(v2, fit.vgm2$var_model)
  fit.vgm2$var_model
  # fit.vgm$var_model
  
  # v3 <- variogram(fit_spatial  ~ res_spatial, data=spoints_samp, cutoff=3000, width=3000/30)
  # plot(v3)
  # fit.vgm3 = autofitVariogram(fit_spatial  ~ res_spatial, spoints_samp, model = c("Exp","Sph"))
  # plot(v3, fit.vgm3$var_model)
  # fit.vgm3$var_model
  
  
  
  model_spatial <- NULL
  model_spatial$type[1] <- "spatial"
  model_spatial$beta[1] <- summary(lm_2)$coefficients[2]
  model_spatial$beta2[1] <- summary(lm_2)$coefficients[3]
  # model_spatial$sig2[1] <- summary(lm_2)$sigma^2
  # model_spatial$sig2[1] <- fit.vgm2$var_model_spatial$psill[2]
  model_spatial$sig2[1] <- gamma_sigma_3[2]^2
  model_spatial$range[1] <- fit.vgm2$var_model$range[2]
  # model_spatial$sig2_2[1] <- fit.vgm3$var_model_spatial$psill[2]
  # model_spatial$range_2[1] <- fit.vgm3$var_model_spatial$range[2]
  model_spatial$gamma[1] <- gamma_sigma_3[1]
  model_spatial$fitting[1] <- summary(lm_2)$r.square
  model_spatial <- as.data.frame(model_spatial)
  model_spatial
  
  set.seed(1234)
  solution3 <- optimizeStrata2 (
    errors=cv, 
    framesamp=frame,
    model=model_spatial,
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
  
  # Comparison with same sample size of Solution 1
  size <- sum(solution1$aggr_strata$SOLUZ)
  newstrata <- adjustSize(size,outstrata)
  sum(newstrata$SOLUZ)
  newstrata
  unlink("./simulation",recursive=TRUE)
  val3a <- evalSolution(framenew,newstrata,nsampl = 1000,progress = F)
  val3a$rel_bias
  val3a$coeff_var 
  
  result <- list(n1 = sum(round(solution1$aggr_strata$SOLUZ)),
                 cv1 = as.numeric(val1$coeff_var[2]),
                 gamma1 = as.numeric(gamma_sigma_1[1]*2),
                 var1 = as.numeric(gamma_sigma_1[1]^2),
                 n2 = sum(round(solution2$aggr_strata$SOLUZ)),
                 cv2 = as.numeric(val2$coeff_var[2]),
                 gamma2 <- as.numeric(gamma_sigma_2[1]*2),
                 var2 <- as.numeric(gamma_sigma_2[1]^2),
                 fitting2 = summary(lm_pred)$r.squared, 
                 range2 = fit.vgm$var_model$range[2],
                 cv2a = as.numeric(val2a$coeff_var[2]),
                 # n2b = sum(round(solution2b$aggr_strata$SOLUZ)),
                 # cv2b = as.numeric(val2b$coeff_var[2]),
                 # cv2c = as.numeric(val2c$coeff_var[2]),
                 n3 = sum(round(solution3$aggr_strata$SOLUZ)),
                 cv3 = as.numeric(val3$coeff_var[2]),
                 gamma3 <- as.numeric(gamma_sigma_3[1]*2),
                 var3 <- as.numeric(gamma_sigma_3[1]^2),
                 fitting3 = summary(lm_2)$r.squared, 
                 range3 = fit.vgm2$var_model$range[2],
                 cv3a = as.numeric(val3a$coeff_var[2])
                 )
  return(result)
}

setwd("C:/Users/Giulio/Google Drive/Paper Spatial Sampling/script per paper/Case study 2/simulations")
load("Bologna.RData")

iters <- read.csv2(filename,dec=".")
res <- iters
res$n1 <- NA
res$cv1 <- NA
res$gamma1 <- NA
res$var1 <- NA
res$n2 <- NA
res$cv2 <- NA
res$gamma2 <- NA
res$var2 <- NA
res$fitting2 <- NA
res$range2 <- NA
res$cv2a <- NA
# res$n2b <- NA
# res$cv2b <- NA
# res$cv2c <- NA
res$n3 <- NA
res$cv3 <- NA
res$gamma3 <- NA
res$var3 <- NA
res$fitting3 <- NA
res$range3 <- NA
res$cv3a <- NA
res <- as.data.frame(res)
for (i in (1:nrow(iters))) {
# for (i in (1:1)) {
  results <- simula(itera = i,
                     Comune_BO_geo,
                     spoints,
                     cc,
                     range_eps = iters$range[i],
                     var_eps = iters$var_eps[i],
                     gamma = iters$gamma[i],
                     beta1 = 1,
                     beta2 = 1)
  res$n1[i] <- results$n1
  res$cv1[i] <- as.numeric(results$cv1)
  res$gamma1[i] <- results$gamma1
  res$var1[i] <- results$var1
  res$n2[i] <- results$n2
  res$cv2[i] <- as.numeric(results$cv2)
  res$gamma2[i] <- results[[7]]
  res$var2[i] <- results[[8]]
  res$fitting2[i] <- results$fitting2
  res$range2[i] <- results$range2
  res$cv2a[i] <- as.numeric(results$cv2a)
  # res$n2b[i] = results$n2b
  # res$cv2b[i] = results$cv2b
  # res$cv2c[i] = results$cv2c
  res$n3[i] <- results$n3
  res$cv3[i] <- as.numeric(results$cv3)
  res$gamma3[i] <- results[[14]]
  res$var3[i] <- results[[15]]
  res$fitting3[i] <- results$fitting3
  res$range3[i] <- results$range3 
  res$cv3a[i] <- as.numeric(results$cv3a)
}
write.table(res,"simul_results_gamma_4_var_eps_2000_4.csv",row.names=F,col.names=T,dec=".",quote=F,sep=";")
