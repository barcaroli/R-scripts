library(lattice)
library(sp)
library(gstat)
library(dplyr) # for "glimpse"
library(ggplot2)
library(scales) # for "comma"
library(magrittr)

library(SamplingStrata)

# locations (155 observed points)
data("meuse")
summary(meuse)
# grid of points (3103)
data("meuse.grid")
meuse.grid$id <- c(1:nrow(meuse.grid))
cor(meuse[,c(3:6)])
plot(meuse$x,meuse$y)
plot(meuse.grid$x,meuse.grid$y)
coordinates(meuse)<-c("x","y")
coordinates(meuse.grid)<-c("x","y")
spplot(meuse,"zinc",do.=T)
spplot(meuse.grid,"soil")
bubble(meuse,"zinc",do.=T,key.space="bottom")
meuse.grid <- as(meuse.grid, "SpatialPixelsDataFrame")
# multivariate case:cadmium, copper, lead, zinc as target variables
g <- gstat(NULL, "Cd", (cadmium) ~ 1, meuse)
g <- gstat(g, "Cu", (copper) ~ 1, meuse)
g <- gstat(g, "Pb", (lead) ~ 1, meuse)
g <- gstat(g, "Zn", (zinc) ~ 1, meuse)
g
vm <- variogram(g)
vm.fit <- fit.lmc(vm, g, vgm(psill=1, model="Sph", range=800, nugget=1))
vm.fit
plot(vm, vm.fit)

cok.maps <- predict(vm.fit, meuse.grid)
names(cok.maps)


coeffvar <- function(x) {
  cv <- sd(x) / mean(x)
  return(cv)
}
tapply(frame$Y1,frame$domainvalue,coeffvar)
tapply(frame$Y2,frame$domainvalue,coeffvar)
tapply(frame$Y3,frame$domainvalue,coeffvar)
tapply(frame$Y4,frame$domainvalue,coeffvar)

df <- NULL
df$Cd.pred <- cok.maps@data$Cd.pred
df$Cd.var <- cok.maps@data$Cd.var
df$Cu.pred <- cok.maps@data$Cu.pred
df$Cu.var <- cok.maps@data$Cu.var
df$Pb.pred <- cok.maps@data$Pb.pred
df$Pb.var <- cok.maps@data$Pb.var
df$Zn.pred <- cok.maps@data$Zn.pred
df$Zn.var <- cok.maps@data$Zn.var
df$dom <- meuse.grid@data$soil
df$dom1 <- 1
df <- as.data.frame(df)
df$id <- meuse.grid$id
df$soil <- meuse.grid$soil

library(SamplingStrata)

cv <- NULL
cv$DOM <- rep("DOM1",3)
cv$CV1 <- rep(0.10,3)
cv$CV2 <- rep(0.05,3)
cv$CV3 <- rep(0.05,3)
cv$CV4 <- rep(0.05,3)
cv$domainvalue <- 1:3
cv <- as.data.frame(cv)
cv


#---------------------------------------------------
# kmean spatial solution 
#---------------------------------------------------
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
                      domainvalue = "soil")
frame$lon <- meuse.grid$x
frame$lat <- meuse.grid$y
frame$var1 <- df$Cd.var
frame$var2 <- df$Cu.var
frame$var3 <- df$Pb.var
frame$var4 <- df$Zn.var


kmeans <- KmeansSolutionSpatial(frame,
                                fitting=1,
                                range=800,
                                gamma=1,
                                errors=cv,
                                maxclusters = 5)
table(kmeans$domainvalue,kmeans$suggestions)

strataKm <- NULL
for (i in unique(frame$domainvalue)) {
  strata <- aggrStrataSpatial(dataset=frame,
                              fitting=1,
                              range=800,
                              gamma=1,
                              vett=kmeans$suggestions[frame$domainvalue==i],
                              dominio=i)
  strata$SOLUZ <- bethel(strata,cv[i,])
  strataKm <- rbind(strataKm,strata)
}
strataKm
sum(strataKm$SOLUZ)
strataKm$STRATO <- strataKm$stratum
framenew <- frame
framenew$LABEL <- kmeans$suggestions
ssKm <- summaryStrata(framenew,strataKm)
ssKm


sugg <- prepareSuggestionSpatial(ssKm)
sugg

frameres1 <- SpatialPointsDataFrame(data=framenew, coords=cbind(framenew$lon,framenew$lat) )
frameres1$LABEL <- as.factor(frameres1$LABEL)
spplot(frameres1,"LABEL")

#---------------------------------------------------
# Solution with spatial optimization on continuous variables
#---------------------------------------------------

set.seed(4321)
solution <- optimizeStrataSpatial2 (
  errors=cv, 
  framesamp=frame,
  iter = 50,
  pops = 10,
  nStrata = 5,
  mut_chance = NA,
  # suggestions = sugg,
  fitting = 1,
  range = 800,
  gamma = 1,
  writeFiles = FALSE,
  showPlot = TRUE,
  parallel = TRUE
)
save(solution,file="solution_meuse_dom.RData")
sum(solution$aggr_strata$SOLUZ)

outstrata <- solution$aggr_strata
cv1 <- cv[1,]
cv1$CV1 <- 010
sol <- bethel(outstrata,cv1,printa=TRUE)
sum(sol)
sol

expected_CV(solution$aggr_strata)
ss <- summaryStrata(x=solution$framenew,outstrata=solution$aggr_strata)
ss
framenew <- solution$framenew
framenew$lon <- meuse.grid@coords[,1]
framenew$lat <- meuse.grid@coords[,2]
# framenew$strat <- paste(framenew$DOMAINVALUE,framenew$LABEL,sep="")
table(framenew$LABEL)
# selection of sample
sample <- selectSample(frame=framenew,outstrata=solution$aggr_strata)
head(sample)


frameres2 <- SpatialPointsDataFrame(data=framenew, coords=cbind(framenew$LON,framenew$LAT) )
frameres2$STR <- paste(frameres2$DOMAINVALUE,frameres2$LABEL,sep="")
frameres2$STR <- as.factor(frameres2$STR)
spplot(frameres2,"STR")

samp <- SpatialPointsDataFrame(data=sample, coords=cbind(sample$LON,sample$LAT) )
samp$STR <- paste(samp$DOMAINVALUE,samp$LABEL,sep="")
samp$STR <- as.factor(samp$STR)
samp$STR <- as.factor(samp$STR)
spplot(samp,"STR")

val <- evalSolution(solution$framenew,solution$aggr_strata,nsampl=30)
val$coeff_var                  
