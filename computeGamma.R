# PDQ method of estimating the coefficient of heteroscedasticity 
# (from J.R.Knaub "Estimating the Coefficient of Heteroscedasticity" June 2019) 


computeGamma <- function (y, x, dataset) {
  # fitting of the model 
  stmt <- paste("mod1 <- lm(",y," ~ ",x,", data=",dataset,")",sep="")
  eval(parse(text=stmt))
  # residuals
  e <- mod1$residuals
  # fitted values
  ystar <- mod1$fitted.values
  # y values
  stmt <- paste("y <- ", dataset,"$",y,sep="")
  eval(parse(text=stmt))
  # tranformation of fitted values (e0)
  e0 <- abs(e)/sqrt(abs(ystar))
  # tranformation of fitted values (e0)
  ynew <- log(abs(y-ystar))
  # OLS on transformed values
  mod2 <- lm(ynew[ynew > 0] ~ log(ystar[ynew > 0]))
  gamma <- as.numeric(mod2$coefficients[2])
  return(gamma)
}



