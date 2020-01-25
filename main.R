# Metadata

# Import packages
library(raster)

# Load functions
source("R/raster_thresholding.R")

# Create data and output folders and download data from URL
data_folder   <- "./data"
output_folder <- "./output"

if (!dir.exists(data_folder)){
  dir.create(data_folder)}

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Load MNDWI data
#mndwi_URL = "https://drive.google.com/open?id=1H2LJKpcRhppegSBjTBoG8E8O677zgu4T"
#download.file(url=mndwi_URL, destfile='./data/mndwiSeries.zip', method='auto', overwrite = TRUE)
unzip('./data/mndwiSeries.zip', exdir = data_folder, overwrite = TRUE)
tiffile = list.files(data_folder, pattern = glob2rx("*.tif"), full.names = TRUE)
MNDWI = stack(tiffile)



