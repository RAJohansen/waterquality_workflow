#Title: waterquality_workflow: A case study and workflow for detecting and quantifying cyanobacterial harmful algal blooms (CHABs) from Sentinel-2 Imagery
# Author: Richard A. Johansen
# Date: December 6th 2018
# Source: https://github.com/RAJohansen/waterquality_workflow
# Citation: Richard A. Johansen. (2018). RAJohansen/waterquality_workflow: waterquality_workflow: A case study and workflow for detecting and quantifying cyanobacterial harmful algal blooms (CHABs) from Sentinel-2 Imagery (Version v.0.2). Zenodo. http://doi.org/10.5281/zenodo.2003619

### Initial Requirements and R Packages ----------------------------------------
#Packages must be installed using install.packages("Package Name") or 
# devtools::install_github("Package Name"), if this is the first time you are using these packages.
library(tidyverse)
library(raster)
library(waterquality)
library(sen2r)
library(sf)
library(gdalUtils)
library(magrittr)
library(rgdal)
library(caret)

#### Sen2Cor Atmospheric Correction---------------------------------------------
#Make Sure Dependencies are installed
check_sen2r_deps()
#Once all are installed close dependencies GUI

#Run sen2cor
#default output will be in a new folder at the same level as the input level 1C
#Calculate time elapsed using system.time Process takes ~30 minutes
sen2cor(".../S2A_MSIL1C_20160808T162342_N0204_R040_T16SGJ_20160808T162611.SAFE")

### Preprocessing imagery ------------------------------------------------------
#Recommended*** Option 1: Preprocess  Clipped & Masked Reflectance Imagery

#Image Directory Folder
Image_Directory = "C:/temp/S2A_MSIL2A_20160808T162342_N0204_R040_T16SGJ_20160808T162611.SAFE/GRANULE/L2A_T16SGJ_A005900_20160808T162611/IMG_DATA/R20m"
# Extracts all raster files from Image Directory with extention .jp2
Rasters = dir(Image_Directory, pattern = "*.jp2$", full.names = TRUE)

#Import Shapefile of Area of Interest
AOI = st_read("C:/temp/waterquality Vignette/Harsha_Lake.gpkg")

#Reproject AOI if needed (Example - UTM Zone 16)
Projection = "+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
AOI = st_transform(AOI, Projection)
AOI = as(AOI, "Spatial")

# Create Template for Final Raster Stack
#For Sentinel-2 we want our final output to be in 20m so use Band 5
raster_B5 = raster(Rasters[[5]])
raster_template = raster(raster_B5)

# Resample, Crop, & Stack All Images
raster_stack = Rasters %>% 
  lapply(raster) %>% 
  lapply(resample, raster_template) %>%
  lapply(crop, AOI) %>%
  stack()

# Mask cropped image for further reduction of the stacked image
S2_Ref_Image <- mask(raster_stack,AOI)

#Band Subset
S2_Harsha <- S2_Ref_Image[[1:9]] #Sentinel-2 Algorithms Only Use Bands 1-8A

#Save Final Stacked Image as Tiff
writeRaster(x = S2_Harsha,
            filename= "C:/temp/S2_Harsha_08082016.tif", # save as a tif
            datatype ="FLT4S", # save as a float 4 significant digits
            overwrite = FALSE) #Overwrites same named file


#Option 2: Preprocess Full Image Reflectance Imagery
#Image Directory Folder
Image_Directory = "C:/temp/S2A_MSIL2A_20180918T154911_N0206_R054_T18STD_20180918T205506.SAFE/GRANULE/L2A_T18STD_A016925_20180918T160118/IMG_DATA/R20m"
# Extacts all raster files from Image Directory with extention .jp2
Rasters = dir(Image_Directory, pattern = "*.jp2$", full.names = TRUE)
S2_Ref_Image = stack(Rasters) 

#Band Subset
S2_Harsha <- S2_Ref_Image[[1:9]] #Sentinel-2 Algorithms Only Use Bands 1-8A

#Save Final Stacked Image as Tiff
writeRaster(x = S2_Harsha,
            filename= "C:/temp/S2_Harsha.tif", # save as a tif
            datatype ="FLT4S", # save as a float 4 significant digits
            overwrite = FALSE) #Overwrites same named file


###Calculate Water Quality indices using waterquality---------------------------
S2_raster <- stack("C:/R_Packages/USACE_WQ/Data/S2_Harsha_08082016.tif") 
S2_wq <- wq_calc(raster_stack = S2_raster, alg = "all", sat = "sentinel2")
writeRaster(x = S2_wq,
            filename= "C:/temp/S2_Harsha_WQ_08082016.tif", # save as a tif
            datatype ="FLT4S", # save as a float 4 significant digits
            overwrite = FALSE) #Overwrites same named file
### Extract Values from Raster imagery from Shapefile---------------------------
#Option 1: Import Shapefile
wq_points <- shapefile('C:/temp/samples.shp')

#Option 2: Create spatial file from spreadsheet
#Read water quality data from csv
library(sp)
wq_df <- read.csv("C:/R_Packages/USACE_WQ/Data/Harsha_WQ_Measurements_08082016.csv")

# Get long and lat from your data.frame. Make sure that the order is in long/lat.
xy <- wq_df[,c(4,3)]
# Create spatial objects (points) using xy and rest of water quality data
#define projection! Harsha Lake is in utm zone 16
wq_points <- SpatialPointsDataFrame(coords = xy, data = wq_df, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
##View points on map
#require(mapview)
#mapview(wq_points)
writeOGR(obj=wq_points, dsn="C:/R_Packages/USACE_WQ/Data/Harsha_wq_points.gpkg", layer="wq_points", driver="GPKG") # this is in geographical projection

#Extract values from tiff using points
#Input raster image
raster <- stack("C:/R_Packages/USACE_WQ/Data/S2_Harsha_WQ_08082016.tif")
##View raster on map
#require(mapview)
#mapview(raster[[1]])
waterquality_data <- data.frame(wq_points, extract(S2_wq, wq_points))
write.csv(waterquality_data, file = "C:/R_Packages/USACE_WQ/Data/Harsha_WQ_Algs_08082016.csv")

### Inspect and Remove mix pixels-----------------------------------------------
#Manually inspected image for mixed pixels and clouds.
#One sample site (H03) removed because pixel appears mixed land covers due to beach
#Read Data frame 
df <- read.csv("C:/R_Packages/USACE_WQ/Data/Harsha_WQ_Algs_08082016.csv")

#Removed H03
df <- df[c(1:2,4:42),]

### Calculate Descriptive Stats and Create Scatterplot matrix-------------------
#Calculate Summary statistics
df_stats <- pastecs::stat.desc(df[,c(5:7,9,15:18)])

#Scatter Plot of highest performing algorithms 
#Optional: Add Text to graph
# Run fit model
fit <- lm(Chl_ugL~BGA_PC_RFU, data = df)
rmse <- round(sqrt(mean(resid(fit)^2)), 2)
coefs <- coef(fit)
b0 <- round(coefs[1], 2)
b1 <- round(coefs[2],2)
r2 <- round(summary(fit)$r.squared, 2)
eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                r^2 == .(r2) * "," ~~ RMSE == .(rmse))

#Chl & PC & Turbidity Scatterplot Matrix
panel.lm <- function (x, y,  pch = par("pch"), col.lm = "red",  ...) {   
  ymin <- min(y)
  ymax <- max(y)
  xmin <- min(x)
  xmax <- max(x)
  ylim <- c(min(ymin,xmin),max(ymax,xmax))
  xlim <- ylim
  points(x, y, pch = pch,ylim = ylim, xlim= xlim,...)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    abline(lm(y[ok]~ x[ok]), 
           col = col.lm, ...)
}

pairs(df[c(15,16,18)], panel=panel.lm)


#Chl & Al10SABI  
ggplot(df, aes(Chl_ugL, Al10SABI)) +
  geom_point() +
  geom_smooth(method = "lm",se = FALSE) +
  theme_classic() + 
  labs(x = "Chlorophyll (µg/L)", y = "Al10SABI", subtitle = eqn) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

#BGA_PC_RFU & Go04MCI  
ggplot(df, aes(BGA_PC_RFU, Go04MCI)) +
  geom_point() +
  geom_smooth(method = "lm",se = FALSE) +
  theme_classic() +
  labs(x = "BGA/PC (RFU)", y = "Go04MCI", subtitle = eqn) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

#Turbidity & Go04MCI  
ggplot(df, aes(Turbid_NTU, Go04MCI)) +
  geom_point() +
  geom_smooth(method = "lm",se = FALSE) +
  theme_classic() +
  labs(x = "Turbidity (NTU)", y = "Go04MCI", subtitle = eqn) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

###Conduct Cross-Validated linear regression Analysis--------------------------
# Set seed using date created for reproducibility
set.seed(2018-11-09)

# Run LM and Cross-Validation Function
# DONT NOT ALTER
extract_lm_caret = function(y, x, df){
  my_formula = as.formula(paste(y, "~", x))
  caret_model = train(form = my_formula,
                      data = df,
                      method = "lm",
                      na.action = na.exclude,
                      #repeated k-fold validation
                      trControl = trainControl(method = "repeatedcv",
                                               number = 3, repeats = 5)) 
  my_lm = caret_model$finalModel
  CV_R_Squared = getTrainPerf(caret_model)[, "TrainRsquared"]
  RMSE = getTrainPerf(caret_model)[, "TrainRMSE"]
  MAE = getTrainPerf(caret_model)[, "TrainMAE"]
  R_Squared = summary(my_lm)$r.squared
  P_Value = summary(my_lm)$coefficients[8]
  Slope = summary(my_lm)$coefficients[2]
  Intercept = summary(my_lm)$coefficients[1]
  data_frame(R_Squared = R_Squared, Slope = Slope, Intercept = Intercept, P_Value = P_Value,
             CV_R_Squared = CV_R_Squared, RMSE = RMSE, MAE = MAE)
}

#Define Parameters
#Column numbers of water quality indices from df
indices <- 24:49
#Water quality parameter to evaluate against indices
WQ_parameter <- "Turbid_NTU" #In situ water quality parameter (Must match name from column in df)

#Run LM and Cross-Validation
Algorithms = names(df)[indices]
names(Algorithms) = Algorithms
Harsha_Turbidity_08082016 = Algorithms %>% 
  map_dfr(~extract_lm_caret(y = WQ_parameter, x = ., df = df), .id="Algorithms")

### Export Results -------------------------------------------------------------
write_csv(Harsha_Turbidity_08082016, path=".../Harsha_08082016_Turbidity_Results.csv")
