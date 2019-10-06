est <- read.csv2("simul_results_gamma_4_var_eps_200.csv",dec=".")
est <- est[1:15,]

maxi = (max(c(est$n1),(est$n2b),(est$n3)))
plot(est$n1, col = "red", lty = 4, type = "p", pch = "&", cex = 1,cex.lab = 1.3,
     xlab = "Heteroscedasticity", 
     ylab = "Sample size",
     ylim = c(0,maxi), axes = F)
points(est$n2b, col = "black", lty = 4, type = "p", pch = "$", cex = 1)
points(est$n3, col = "blue", lty = 4, type = "p", pch = "£", cex = 1)

abline(0,0)
axis(1, at=1:nrow(est), labels = est$gamma,las=2)
axis(2, at=seq(0,maxi+10,50), labels=seq(0,maxi+10,50), las=1)
legend("topleft",
       legend = c("Linear model","Kriging","Spatial Linear Model"), 
       pch = c("&","$","£"),
       col = c("red","black","blue"),
       ncol = 1, cex = 1.2, text.font = 1)

# title("Performance in terms of sampling size",cex.main=1.5)

maxi = (max(c(est$cv1),(est$cv2b),(est$cv3)))
plot(est$cv1, col = "red", lty = 4, type = "p", pch = "&", cex = 1,cex.lab = 1.3,
     xlab = "Heteroscedasticity", 
     ylab = "CV",
     ylim = c(0,maxi+0.01), axes = F)
points(est$cv2b, col = "black", lty = 4, type = "p", pch = "$", cex = 1)
points(est$cv3, col = "blue", lty = 4, type = "p", pch = "£", cex = 1)

abline(0,0)
abline(a=0.03, b=0.0,col="blue")
axis(1, at=1:nrow(est), labels = est$gamma,las=2)
axis(2, at=seq(0,maxi+10,0.01), labels=seq(0,maxi+10,0.01), las=1)
legend("topleft",
       legend = c("Linear model","Kriging","Spatial Linear Model"), 
       pch = c("&","$","£"),
       col = c("red","black","blue"),
       ncol = 1, cex = 1.2, text.font = 1)

# title("Performance in terms of CV on target variable",cex.main=1.5)

maxi = (max(c(est$cv1),(est$cv2c),(est$cv3a)))
plot(est$cv1, col = "red", lty = 4, type = "p", pch = "&", cex = 1,cex.lab = 1.3,
     xlab = "Heteroscedasticity", 
     ylab = "Equalized CV",
     ylim = c(0,maxi+0.01), axes = F)
points(est$cv2c, col = "black", lty = 4, type = "p", pch = "$", cex = 1)
points(est$cv3a, col = "blue", lty = 4, type = "p", pch = "£", cex = 1)

abline(0,0)
abline(a=0.03, b=0.0,col="blue")
axis(1, at=1:nrow(est), labels = est$gamma,las=2)
axis(2, at=seq(0,maxi+10,0.01), labels=seq(0,maxi+10,0.01), las=1)
legend("topleft",
       legend = c("Linear model","Kriging","Spatial Linear Model"), 
       pch = c("&","$","£"),
       col = c("red","black","blue"),
       ncol = 1, cex = 1.2, text.font = 1)

# title("Performance in terms of equalized CV on target variable",cex.main=1.2)
