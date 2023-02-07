## ----setup, include=FALSE--------------------------------------------------------
require(knitr)
knitr::opts_chunk$set(echo = TRUE)
r <- getOption("repos")
r["CRAN"] <- "https://ftp.osuosl.org/pub/cran/"
options(repos = r)


## ----load packages, include=FALSE------------------------------------------------

#function to install and load required packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#load or install these packages:
packages <- c("ROCR", "tidyverse", "ks", "mapview", "colorRamps","rgeos", "VGAM", "AICcmodavg", "MuMIn", "corrgram", "GGally","caret", "DescTools", "car", "sf", "terra", "stars", "tmap")

#run function to install packages
ipak(packages)


## --------------------------------------------------------------------------------
wolfkde <- read.csv(here::here("Data","wolfkde5.csv"), header=TRUE, sep = ",", na.strings="NA", dec=".")
wolfkde3 <-na.omit(wolfkde)
wolfkde3$usedFactor <-as.factor(wolfkde3$usedFactor) ## make sure usedFactor is a factor

kernelHR <- st_read(here::here("Data","homerangeALL.shp"))

tmap_mode("plot")
tm_shape(kernelHR) + tm_polygons()

plot(kernelHR)

ext(kernelHR)
kernels <- rast()
ext(kernels) <- c(xmin=546836, xmax=612093, ymin=5662036, ymax=5748911) 


## --------------------------------------------------------------------------------
deer_w<-rast("Data/deer_w2.tiff")
moose_w<-rast("Data/moose_w2.tif")
elk_w<-rast("Data/elk_w2.tif") # already brought in above
sheep_w<-rast("Data/sheep_w2.tif")
goat_w<-rast("Data/goat_w2.tif")
wolf_w<-rast("Data/wolf_w2.tif")
elevation2<-rast("Data/Elevation2.tif") #resampled
disthumanaccess2<-rast("Data/DistFromHumanAccess2.tif") #resampled in lab 4
disthhu2<-rast("Data/DistFromHighHumanAccess2.tif") #resampled in lab 4
landcover2 <- rast("Data/landcover16.tif") ## resampled to same extent as lab 4
ext(landcover2)
plot(landcover2)
plot(kernelHR, add=TRUE, col = NA) #note - col = NA maps just the polygon border


## --------------------------------------------------------------------------------
cats(landcover2)


## --------------------------------------------------------------------------------

elevation2 <- resample(elevation2, elk_w)
disthhu2 <- resample(disthhu2, elk_w)
disthumanaccess2 <- resample(disthumanaccess2, elk_w)
landcover2 <- resample(landcover2, elk_w)



## --------------------------------------------------------------------------------


alpine <- ifel(landcover2 == 15 | landcover2 == 16, 1,ifel(is.na(landcover2),NA,0))
burn <- ifel(landcover2 == 12 | landcover2 == 13 |landcover2 == 14, 1, ifel(is.na(landcover2),NA,0))
closedConif <- ifel(landcover2 == 3, 1, ifel(is.na(landcover2),NA,0))
herb <- ifel(landcover2 == 7, 1, ifel(is.na(landcover2),NA,0))
mixed <- ifel(landcover2 == 5, 1, ifel(is.na(landcover2),NA,0))
rockIce <- ifel(landcover2 == 10, 1, ifel(is.na(landcover2),NA,0))
water <- ifel(landcover2 == 9, 1, ifel(is.na(landcover2),NA,0))
modConif <- ifel(landcover2 == 2, 1, ifel(is.na(landcover2),NA,0))
decid <- ifel(landcover2 == 10, 1, ifel(is.na(landcover2),NA,0))
plot(rockIce)
plot(closedConif)

# note that open conifer is implicitly being defined as the intercept


## --------------------------------------------------------------------------------

all_rasters<-c(deer_w, moose_w, elk_w, sheep_w, goat_w, wolf_w,elevation2, disthumanaccess2, disthhu2, landcover2, alpine, burn, closedConif, modConif, herb, mixed, rockIce, water, decid)

plot(all_rasters) ## note limit of plotting 9 layers


## --------------------------------------------------------------------------------
names = c("deer_w", "moose_w", "elk_w", "sheep_w", "goat_w", "wolf_w","elevation2", "disthumanaccess2", "disthhu2", "landcover2", "alpine", "burn", "closedConif", "modConif", "herb", "mixed", "rockIce", "water", "decid")

writeRaster(all_rasters,"Output/lab6Stack.tif",names = 'names', overwrite = TRUE)


## --------------------------------------------------------------------------------

## top Biotic model was model 41
top.biotic <- glm(used ~ DistFromHumanAccess2+deer_w2 + goat_w2, family=binomial(logit), data=wolfkde3)
summary(top.biotic)
#double check the VIF for each final model, just to be sure
vif(top.biotic)


## --------------------------------------------------------------------------------
top.env <- glm(used ~ Elevation2 + DistFromHighHumanAccess2 + openConif+modConif+closedConif+mixed+herb+shrub+water+burn, family=binomial(logit), data=wolfkde3)
summary(top.env)
vif(top.env)


## --------------------------------------------------------------------------------
models = list(top.biotic, top.env)
modnames = c("top biotic", "top env")
aictab(cand.set = models, modnames = modnames)
# aictab(models, modnames) ## short form. 


## --------------------------------------------------------------------------------
par(mfrow = c(2,2))
plot(top.env)


## --------------------------------------------------------------------------------
##### Saving predictions manually, an example with the environment model
wolfkde3$fitted.top.env <- fitted(top.env)
#### this is the predicted probability from the model

wolfkde3$residuals.top.env <- residuals(top.env)
## these are the deviations from the predictions for each row (data point)

wolfkde3$rstudent.top.env <- rstudent(top.env)
## This is a standardized residual - the studentized residual

wolfkde3$hatvalues.top.env <- hatvalues(top.env)
#### this is the first of the leverage statistics, the larger hat value is, the bigger the influence on the fitted value

wolfkde3$cooks.distance.top.env <- cooks.distance(top.env)
#### this is the Cooks leverage statistic, the larger hat value is, the bigger the influence on the fitted value

wolfkde3$obsNumber <- 1:nrow(wolfkde3) ## just added a row number for plotting


## --------------------------------------------------------------------------------
ggplot(wolfkde3, aes(fitted.top.env, residuals.top.env)) + geom_point() + geom_text(aes(label = obsNumber, colour = used))
## Manual residual plot

ggplot(wolfkde3, aes(wolfkde3$residuals.top.env, wolfkde3$cooks.distance.top.env)) + geom_point() + geom_text(aes(label = obsNumber, colour = used))
## shows us some points at high cooks values that might be having a big influence

ggplot(wolfkde3, aes(wolfkde3$cooks.distance.top.env, wolfkde3$hatvalues.top.env)) + geom_point() + geom_text(aes(label = obsNumber, colour = used))


## --------------------------------------------------------------------------------
## e.g., datapoint 16
wolfkde3[16,]


## ---- warning = FALSE------------------------------------------------------------
scatterplot(fitted.top.env~Elevation2, reg.line=lm, smooth=TRUE, spread=TRUE, boxplots='xy', span=0.5, xlab="elevation", ylab="residual", cex=1.5, cex.axis=1.4, cex.lab=1.4, data=wolfkde3)

hist(wolfkde3$fitted.top.env, scale="frequency", breaks="Sturges", col="darkgray")


## --------------------------------------------------------------------------------
ggplot(wolfkde3, aes(x=wolfkde3$fitted.top.env, fill=usedFactor)) + geom_histogram(binwidth=0.05, position="identity", alpha=0.7) + xlab("Predicted Probability of Wolf Use") + theme(axis.title.x=element_text(size=16)) #+ facet_grid(pack ~ ., scales="free")


## --------------------------------------------------------------------------------
ggplot(wolfkde3, aes(x=fitted.top.env, y=..density.., fill=usedFactor)) + geom_histogram(binwidth=0.05, position="identity", alpha=0.7) + xlab("Predicted Probability of Wolf Use") + theme(axis.title.x=element_text(size=16)) + facet_grid(pack ~ ., scales="free")



## --------------------------------------------------------------------------------
par(mfrow = c(2,2))
plot(top.biotic)

##### Saving predictions manually, an example with the environment model
wolfkde3$fitted.top.biotic <- fitted(top.biotic)
#### this is the predicted probability from the model

wolfkde3$residuals.top.biotic <- residuals(top.biotic)
## these are the deviations from the predictions for each row (data point)

wolfkde3$rstudent.top.biotic <- rstudent(top.biotic)
## This is a standardized residual - the studentized residual

wolfkde3$hatvalues.top.biotic <- hatvalues(top.biotic)
#### this is the first of the leverage statistics, the larger hat value is, the bigger the influence on the fitted value

wolfkde3$cooks.distance.top.biotic <- cooks.distance(top.biotic)
#### This isthe Cooks leverage statistic


## Making manual residual verus predicted plots
ggplot(wolfkde3, aes(wolfkde3$cooks.distance.top.biotic, wolfkde3$hatvalues.top.biotic)) + geom_point() + geom_text(aes(label = obsNumber, colour = used))


## --------------------------------------------------------------------------------
wolfkde3[30,]


## ---- warning = FALSE------------------------------------------------------------
scatterplot(fitted.top.biotic~Elevation2, reg.line=lm, smooth=TRUE, spread=TRUE, boxplots='xy', span=0.5, xlab="elevation", ylab="residual", cex=1.5, cex.axis=1.4, cex.lab=1.4, data=wolfkde3)

## next, the histogram of predicted probaiblities
hist(wolfkde3$fitted.top.biotic, scale="frequency", breaks="Sturges", col="darkgray")


## --------------------------------------------------------------------------------
ggplot(wolfkde3, aes(x=wolfkde3$fitted.top.biotic, fill=usedFactor)) + geom_histogram(binwidth=0.05, position="identity", alpha=0.7) + xlab("Predicted Probability of Wolf Use") + theme(axis.title.x=element_text(size=16)) #+ facet_grid(pack ~ ., scales="free")

#### This plot shows that somewhere around 0.25 - 0.40 eyeballing it it looks like we could 'cut' used and available points? 

ggplot(wolfkde3, aes(x=fitted.top.biotic, y=..density.., fill=usedFactor)) + geom_histogram(binwidth=0.05, position="identity", alpha=0.7) + xlab("Predicted Probability of Wolf Use") + theme(axis.title.x=element_text(size=16)) + facet_grid(pack ~ ., scales="free")
#### But note this 'cut' point looks different for both wolf packs?



## --------------------------------------------------------------------------------
ggplot(wolfkde3, aes(x=fitted.top.biotic, y=fitted.top.env)) + geom_point() + stat_smooth(method="lm")


## --------------------------------------------------------------------------------
ggplot(wolfkde3, aes(x=fitted.top.biotic, y=fitted.top.env, fill = pack)) + geom_point() + stat_smooth(method="lm") 


## --------------------------------------------------------------------------------
ggplot(wolfkde3, aes(x=fitted.top.biotic, y=fitted.top.env, fill = pack)) + geom_point() + stat_smooth()


## --------------------------------------------------------------------------------
#require(DescTools)
PseudoR2(top.biotic, c("McFadden", "CoxSnell", "Nagel"))


## --------------------------------------------------------------------------------
# First we will arbitrarily define the cutpoint between 1 and 0's using p = 0.5
ppused = wolfkde3$fitted.top.biotic>0.5
table(ppused,wolfkde3$used)


## --------------------------------------------------------------------------------
167/(167+229)


## --------------------------------------------------------------------------------
1655/(1655+67)


## --------------------------------------------------------------------------------
ppused = wolfkde3$fitted.top.biotic>0.25
table(ppused,wolfkde3$used)


## --------------------------------------------------------------------------------
304/(304+92)


## --------------------------------------------------------------------------------
1376 / (1376+346)


## --------------------------------------------------------------------------------
ppused = wolfkde3$fitted.top.biotic>0.10
table(ppused,wolfkde3$used)


## --------------------------------------------------------------------------------
357/(357+39)


## --------------------------------------------------------------------------------
1001 / (1001+721)


## --------------------------------------------------------------------------------
ppused = wolfkde3$fitted.top.biotic>0.70
table(ppused,wolfkde3$used)


## --------------------------------------------------------------------------------
66/(330+66)


## --------------------------------------------------------------------------------
17 / (1705+17)


## --------------------------------------------------------------------------------
#require(caret)
wolfkde3$pr.top.biotic.used <- ifelse(wolfkde3$fitted.top.biotic>0.5, 1, 0)
xtab1<-table(wolfkde3$pr.top.biotic.used, wolfkde3$used)
xtab1

#?confusionMatrix
confusionMatrix(xtab1) ## Note this is incorrectly specifying what the 1's are. 

confusionMatrix(xtab1, positive = "1") ## Thanks to spring 2019 student Forest HAyes for finding this mistake :)


## --------------------------------------------------------------------------------
#require(ROCR)
pp = predict(top.biotic,type="response")
pred = prediction(pp, wolfkde3$used)

perf3 <- performance(pred, "sens", x.measure = "cutoff")
plot(perf3)


## --------------------------------------------------------------------------------
perf4 <- performance(pred, "spec", x.measure = "cutoff")
plot(perf4)


## --------------------------------------------------------------------------------
perfClass <- performance(pred, "tpr","fpr") # change 2nd and/or 3rd arguments for other metrics
fpr <- perfClass@x.values[[1]]
tpr <- perfClass@y.values[[1]]
sum <- tpr + (1-fpr)
index <- which.max(sum)
cutoff <- perfClass@alpha.values[[1]][[index]]
cutoff


## --------------------------------------------------------------------------------
table(wolfkde3$used)
396/(1722+396)


## --------------------------------------------------------------------------------
plot(perf3, col="blue") # Sensitivity
plot(perf4, add = TRUE) # Specificity
abline(v=cutoff, col="red") ## optimal cutpoint


## --------------------------------------------------------------------------------
plot(perfClass)
abline(a=0, b= 1)


## --------------------------------------------------------------------------------
BMauc <- performance(pred, measure="auc") 
str(BMauc)
auc <- as.numeric(BMauc@y.values)
auc


## --------------------------------------------------------------------------------
plot(perfClass, colorize = T, lwd = 5, print.cutoffs.at=seq(0,1,by=0.1),
     text.adj=c(1.2,1.2),
     main = "ROC Curve")
text(0.5, 0.5, "AUC = 0.867")
abline(v=cutoff, col = "red", lwd = 3)


## --------------------------------------------------------------------------------
acc.perf = performance(pred, measure = "acc")
plot(acc.perf)


## --------------------------------------------------------------------------------
### now lets try a p = of our cutoff
ppused = wolfkde3$fitted.top.biotic>cutoff
table(ppused,wolfkde3$used)
#### Now, what is specificity? (i.e., the probability of classifying the 1's correctly?)
320/(320+76)
#### about 80% - Great! But - what happened to our sensitivity (i.e., the probability of classifying the 0's correctly?)
1344 / (1344+378)
#### so our probability of classifying 0's correctly is about 78%


## --------------------------------------------------------------------------------
wolfkde3$pr.top.biotic.used2 <- ifelse(wolfkde3$fitted.top.biotic>cutoff, 1, 0)
xtab2<-table(wolfkde3$pr.top.biotic.used2, wolfkde3$used)
xtab2

#?confusionMatrix
confusionMatrix(xtab2)


## --------------------------------------------------------------------------------
## this is our best model classifying used and avail locations into 1 and 0's. 
ggplot(wolfkde3, aes(x=wolfkde3$fitted.top.biotic, fill=usedFactor)) + geom_histogram(binwidth=0.05, position="identity", alpha=0.7) + xlab("Predicted Probability of Wolf Use") + theme(axis.title.x=element_text(size=16)) + geom_vline(xintercept = cutoff, col="red")


## --------------------------------------------------------------------------------
source("Rcode/kxv.R", verbose = FALSE)


## ---- warning = FALSE------------------------------------------------------------
# Kfolds with a 'fuzz' factor
kxvPrintFlag=FALSE
kxvPlotFlag=TRUE
kxvFuzzFactor = 0.01
kfolds = kxvglm(top.env$formula, data=wolfkde3, k=5, nbin=10) ## note we get a lot of ties here and some error messages, this is because of all the categories. Read more about this in the clumpy data.pdf in the Vignette folder.  
kfolds


## --------------------------------------------------------------------------------
# Kfolds by each pack with a 'fuzz' factor
kxvPrintFlag=FALSE
kxvPlotFlag=TRUE
kxvFuzzFactor = 0.01
kfolds2 = kxvglm(top.env$formula, data=wolfkde3, k=5, nbin=10, partition="pack")
kfolds2


## --------------------------------------------------------------------------------
# 1. Create a vector of random "folds" in this case 5, 1:5
wolfkde3$rand.vec = sample(1:5,nrow(wolfkde3),replace=TRUE)

#2. Run the model for 1 single random subset of the data == 1
top.env.1= glm(used ~ Elevation2 + DistFromHighHumanAccess2 + openConif+modConif+closedConif+mixed+herb+shrub+water+burn, family=binomial(logit), data=wolfkde3, subset=rand.vec==1) ## note this last subset = rand.vec==1 is just 1 of the 5 subsets. 

# 3. Make predictions for points not used in this random subset (2:5) to fit the model.
pred.prob = predict(top.env.1,newdata=wolfkde3[wolfkde3$rand.vec!=1,],type="response") ## note that != means everything but subset 1. 

# 4. Make quantiles for the predictions - this calculates the 'bin's of the categories of habitat availability
q.pp = quantile(pred.prob,probs=seq(0,1,.1))

# 5. Then for each of 10 bins, put each row of data into a bin
bin = rep(NA,length(pred.prob))
for (i in 1:10){
	bin[pred.prob>=q.pp[i]&pred.prob<q.pp[i+1]] = i
}

## 5. This then makes a count of just the used locations for all other K folds 2:5 
used1 = wolfkde3$used[wolfkde3$rand.vec!=1]

## We then make a table of them
rand.vec.1.table <- table(used1,bin)
rand.vec.1.table


## --------------------------------------------------------------------------------
cor.test(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), c(0,2,0,8,6,15,24,50,99,110), method="spearman") 


## --------------------------------------------------------------------------------
rand.vec.1.table <- as.data.frame(rand.vec.1.table)
ggplot(rand.vec.1.table, aes(as.numeric(bin), Freq, col = used1)) + geom_point(size=5) + geom_line()


## --------------------------------------------------------------------------------
par(mfrow = c(1,1)) # reset graphical parameters

ggplot(wolfkde3, aes(EASTING, NORTHING, col = fitted.top.biotic)) + geom_point(size=5) + coord_equal() +  scale_colour_gradient(low = 'yellow', high = 'red')
ggplot(wolfkde3, aes(EASTING, NORTHING, col = fitted.top.env)) + geom_point(size=5) + coord_equal() +  scale_colour_gradient(low = 'yellow', high = 'red')


## --------------------------------------------------------------------------------
par(mfrow = c(1,1))
top.biotic$coefficients
biotic.coefs <- top.biotic$coefficients[c(1:4)]


## --------------------------------------------------------------------------------
rast.top.biotic <- exp(biotic.coefs[1] + biotic.coefs[2]*disthumanaccess2 + biotic.coefs[3]*deer_w + biotic.coefs[4]*goat_w) / (1 +exp(biotic.coefs[1] + biotic.coefs[2]*disthumanaccess2 + biotic.coefs[3]*deer_w + biotic.coefs[4]*goat_w ))


## --------------------------------------------------------------------------------
wolfyht<-st_read("Data/wolfyht.shp")
# plot predicted raster within extent of kernels raster
plot(rast.top.biotic, col=colorRampPalette(c("yellow", "orange", "red"))(255), ext=kernels)
plot(kernelHR, add=TRUE, col = NA)
plot(wolfyht, col='blue', pch = 16, add=TRUE)


## --------------------------------------------------------------------------------
#hist(rast.top.biotic@data@values) ## weird, this worked in R but not in R markdown. Try it on your own. 


## --------------------------------------------------------------------------------
##
bv.raster<-rast()
ext(bv.raster) <- c(xmin=570000, xmax=600000, ymin=5665000, ymax=5685000) 
plot(rast.top.biotic, col=colorRampPalette(c("yellow", "orange", "red"))(255), ext=bv.raster)
plot(kernelHR, add=TRUE, col = NA)
plot(wolfyht, col='blue', pch = 16, add=TRUE)


## --------------------------------------------------------------------------------
##
rd.raster<-rast()
ext(rd.raster) <- c(xmin=540000, xmax=600000, ymin=5700000, ymax=5730000) 
plot(rast.top.biotic, col=colorRampPalette(c("yellow", "orange", "red"))(255), ext=rd.raster)
plot(kernelHR, add=TRUE, col = NA)
plot(wolfyht, col='blue', pch = 16, add=TRUE)

