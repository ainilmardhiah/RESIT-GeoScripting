library(raster)
library(sf)
library(rgeos)
## Download, unzip and load the data
tempdir = "temp"
if (!dir.exists(tempdir)) dir.create(tempdir)
download.file(url = 'https://github.com/GeoScripting-WUR/AdvancedRasterAnalysis/archive/gh-pages.zip', destfile = file.path(tempdir, 'data.zip'), method = 'auto')
unzip(file.path(tempdir, 'data.zip'), exdir=tempdir)
# Move from the temp dir to your working directory and clean up
file.rename(file.path(tempdir, 'AdvancedRasterAnalysis-gh-pages', 'data'), 'data')
unlink(file.path(tempdir, 'AdvancedRasterAnalysis-gh-pages'), TRUE)
file.remove(file.path(tempdir, 'data.zip'))
file.remove(tempdir)

# load the data
load("data/GewataB2.rda")
load("data/GewataB3.rda")
load("data/GewataB4.rda")

# execute and check out the attributes
GewataB2

# some basic statistic of the data
cellStats(GewataB2, stat = mean)
cellStats(GewataB2, stat = max)
cellStats(GewataB2, stat = min)

# or if you want to execute statistik all of band between 3 different raster
max(c(maxValue(GewataB2), maxValue(GewataB3), maxValue(GewataB4)))
summary(GewataB2)

# if you have three different band and you want to put into want raster

gewata <- brick(GewataB2, GewataB3, GewataB4)

# execute the histogram in window
hist(gewata, maxpixels = 100)

# if you want to show all of in one
par(mfrow = c(1,1))
hist(gewata, xlim = c(0,5000), ylim = c(0,750000), breaks = seq(0,5000, by = 100))
pairs(gewata)


# calculate the NDVI forgewata
NDVI <- overlay(GewataB4, GewataB3, fun = function(x,y){(x-y)/(x+y)})
plot(NDVI)

# load dthe data and check it out
load("./data/vcfGewata.rda")

plot(vcfGewata)

# execute the vcfgewata
summary(vcfGewata)
hist(vcfGewata)

# because the maximum number of tree is 100, so all outside and greater than 100 must be a flags for water or noise, so all pixel value larger than 1000 is NA
vcfGewata[vcfGewata > 100] <- NA
plot(vcfGewata)
summary(vcfGewata)

# now we need to rescale the raster vcf, it is not significantly required for RF, but only will help for the interpreattion

gewata <- calc(gewata, fun = function(x)x/10000)

# now create the coovariate include the raster gewata, NDVI and vcf
covs <- stack(gewata, NDVI, vcfGewata)
covs

# now change the name to make youeasy to identify
names(covs) <- c("band2","band3", "band4", "NDVI", "VCF")
plot(covs)

# now create training data prepartion
trainingPoly <- readRDS("data/trainingPoly.rds")

plot(NDVI)
plot(trainingPoly, add = TRUE)

str(trainingPoly$Class)

# now create new field

trainingPoly$Code <- as.numeric(trainingPoly$Class)


# asssign code value to raster cell
classes <- rasterize(trainingPoly, NDVI, field = "Code")

# now plot into the nice visualization
cols <- c("orange", "dark green", "light blue")
plot(classes, col = cols, legend = FALSE, main = " Landuse Calssification")
legend("bottomleft", legend = c("cropland", "forest", " wetland"), fill = cols, bg = "white")

## next
covmasked <- mask(covs, classes)
names(classes) <- "class"
trainingstack <- addLayer(covmasked, classes)
plot(trainingstack)


# after that we need to convert it into dataframe data.frame to representing all training data
#before it, we extract the value from training stack()

valuetable <- getValues(trainingstack)
# erase the NA value with na.omit

valuetable <- na.omit(valuetable)

# now we create a matrix as dataframe
valuetable <- as.data.frame(valuetable)
head(valuetable)


# now we are already have the data frame of training data, now we need to convert as factor for the class
valuetable$class <- factor(valuetable$class, levels = c(1:3))


# now we create 3 different dataframe
val_crop <- subset(valuetable, class ==1)
val_forest <- subset(valuetable, class ==2)
val_wetland <- subset(valuetable, class ==3)

# now we create the histogram for each class
hist(val_crop$NDVI, main = "cropland", xlab = "NDVI", xlim = c(0,1), ylim = c(0, 4000), col="orange")
hist(val_crop$NDVI, main = "cropland", xlab = "NDVI", xlim = c(0, 1), ylim = c(0, 4000), col = "orange")
hist(val_forest$NDVI, main = "forest", xlab = "NDVI", xlim = c(0, 1), ylim = c(0, 4000), col = "dark green")
hist(val_wetland$NDVI, main = "wetland", xlab = "NDVI", xlim = c(0, 1), ylim = c(0, 4000), col = "light blue")

# now we want to create scatterplots
plot(band4 ~ band3, data = val_crop, pch = ".", col = "orange", xlim = c(0, 0.2), ylim = c(0, 0.5))
points(band4 ~ band3, data = val_forest, pch = ".", col = "dark green")
points(band4 ~ band3, data = val_wetland, pch = ".", col = "light blue")
legend("topright", legend=c("cropland", "forest", "wetland"), fill=c("orange", "dark green", "light blue"), bg="white")

# now we try random forest 
if(!"randomForest" %in% rownames(installed.packages())){install.packages("randomForest")}

library(randomForest)

modelRF <- randomForest(x=valuetable[ ,c(1:5)], y=valuetable$class,
                        importance = TRUE)

class(modelRF)
str(modelRF)
names(modelRF)


# inspect the confusion matrix
modelRF$confusion

# now to make it mor redeable
colnames(modelRF$confusion) <- c("cropland", "forest", "wetland", "class.error")
rownames(modelRF$confusion) <- c("cropland", "forest", "wetland")
modelRF$confusion


# now if you want to see it in the plot
varImpPlot(modelRF)



names(covs)

# predict landcover using model rf
predLC <- predict(covs, model=modelRF, na.rm=TRUE)
cols <- c("orange", "dark green", "light blue")
plot(predLC, col=cols, legend=FALSE)
legend("bottomright",
       legend=c("cropland", "forest", "wetland"),
       fill=cols, bg="white")




formask <- setValues(raster(predLC), NA)
## Assign 1 to formask to all cells corresponding to the forest class
formask[predLC==2] <- 1
plot(formask, col="dark green", legend = FALSE)


if(!file.exists(fn <- "data/clumformask.grd")) {
  forestclumps <- clump(formask, directions=8, filename=fn)
} else {
  forestclumps <- raster(fn)
}


plot(forestclumps, col=topo.colors(nrow(forestclumps)))


clumpFreq <- freq(forestclumps)
head(clumpFreq)
tail(clumpFreq)







# working with thematic raster, we dont use random forest
load("data/lulcGewata.rda")
## Check out the distribution of the values
freq(lulcGewata)

load("data/LUTGewata.rda")
LUTGewata


lulc <- as.factor(lulcGewata)
# assign a raster attribute table (RAT)
levels(lulc) <- LUTGewata
lulc


classes <- layerize(lulc)
# Layer names follow the order of classes in the LUT
names(classes) <- LUTGewata$Class
plot(classes, legend=FALSE)



forest <- raster(classes, 5)
# is equivalent to:
forest <- classes[[5]]
# or (since the layers are named):
forest <- classes$forest
## Replace 0's (non-forest) with NA's
forest[forest==0] <- NA
plot(forest, col="dark green", legend=FALSE)





# I am trying to crop a raster to given vector boundaries, so that anything not covered by a polygon becomes a missing value.
# How can I do that? A minimal working example is:
library(raster)
library(sf)

# Download a vector file
FileName = "Wageningen-boundary.csv"
download.file("https://raw.githubusercontent.com/GeoScripting-WUR/IntroToVector/gh-pages/data/Wageningen-boundary.csv", FileName)
WagBoundary = st_read(FileName)
st_crs(WagBoundary) = 4326

# Download a raster file
ExampleRaster = getData('worldclim', var='tmin', res=0.5, lon=st_bbox(WagBoundary)$xmin, lat=st_bbox(WagBoundary)$ymin)[[1]]

# Crop
Result = crop(ExampleRaster, WagBoundary)
plot(Result)








# alhamdulilah, tinggal bahas soal ujian, practice exam, sama assignmentttttt, omg banyak yaa
#lah kagak guna sign makanya nggk update di github guee