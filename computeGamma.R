# Estimating the coefficient of heteroscedasticity 

computeGamma <- function (y, x, dataset) {
  # fitting of the model 
  st <- paste("df <- ",dataset,sep="")
  eval(parse(text=st))
  ris_gamma<-NULL
  stmt <- paste("mod1 <- lm(",y," ~ ",x,", data=",dataset,")",sep="")
  eval(parse(text=stmt))
  summary(mod1)
  df$e <- mod1$residuals
  df$p <- mod1$fitted.values
  # clustering of predicted in 20 bins
  # df$p_bins <- var.bin(df$p,20)
  df$p_bins <- var.bin(df$P1,20)
  std_eps <- sqrt(tapply(df$e,df$p_bins,var))
  x_eps<-tapply(df$P1,df$p_bins,mean)
  plot(x_eps,std_eps)
  plot(log(x_eps),log(std_eps))
  lm_gamma<-lm(log(std_eps)~log(x_eps))
  lm_gamma$coefficients
  gamma_est<-lm_gamma$coefficients[2]
  si<-exp(lm_gamma$coefficients[1])
  si
  gamma_est
  ris_gamma[1]<-gamma_est
  ris_gamma[2]<-si
  
  return(ris_gamma)
}



