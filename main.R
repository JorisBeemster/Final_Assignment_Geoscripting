# Metadata

# Import packages
library(raster)
library(sp)

# Create data and output folders and download data from URL
data_folder   <- "./data"
output_folder <- "./output"

if (!dir.exists(data_folder)){
  dir.create(data_folder)}

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Load MNDWI data
mndwi_URL = "https://drive.google.com/open?id=1vXGFJuFcCepCvKT50IEDq-9PcvBwQFHo"
download.file(url=mndwi_URL, destfile='data/mndwiSeries.tif', method='auto', overwrite = TRUE)