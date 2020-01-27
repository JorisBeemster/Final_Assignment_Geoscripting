# Metadata

# Import packages
library(raster)
library(rtiff)
library(sp)
library(rgdal)
library(gdalUtils)

# Load functions
source("R/raster_thresholding.R")
source("R/boundary_detection.R")
source("R/heightFromCL.R")

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
tiffile = list.files(data_folder, pattern = glob2rx("mndwi*.tif"), full.names = TRUE)
MNDWI = stack(tiffile)
MNDWI_info = gdalinfo(tiffile)

# Get dates of bands
MNDWI_dates = MNDWI_info[grepl(glob2rx('*LC08*'), MNDWI_info)]
MNDWI_dates = substring(MNDWI_dates, 45)
  
# classify as water or land
water_classification = raster_thresholding(MNDWI)

# boundary detection
coastlines = boundary_detection(water_classification)

# download AHN data
AHN_urls = c('https://geodata.nationaalgeoregister.nl/ahn3/extract/ahn3_5m_dtm/M5_49CZ2.ZIP',
            'https://geodata.nationaalgeoregister.nl/ahn3/extract/ahn3_5m_dtm/M5_49DZ1.ZIP',
            'https://geodata.nationaalgeoregister.nl/ahn3/extract/ahn3_5m_dtm/M5_55AN2.ZIP',
            'https://geodata.nationaalgeoregister.nl/ahn3/extract/ahn3_5m_dtm/M5_55BN1.ZIP',
            'https://geodata.nationaalgeoregister.nl/ahn3/extract/ahn3_5m_dtm/M5_49CZ1.ZIP',
            'https://geodata.nationaalgeoregister.nl/ahn3/extract/ahn3_5m_dtm/M5_55AN1.ZIP')
AHN_mapsheets =  c('49CZ2', '49DZ1', '55AN2', '55BN1','49CZ1','55AN1')
for (i in 1:length(AHN_urls)){
  filename = paste(data_folder, '/', AHN_mapsheets[i], sep = '')
  download.file(url = AHN_urls[i], destfile = filename, method = 'auto', overwrite = TRUE)
  unzip(filename, exdir = data_folder, overwrite = TRUE)
  file.remove(filename)
}
AHN_files = list.files(data_folder, pattern = glob2rx('m5*'), full.names = TRUE)
AHN_5m = merge(raster(AHN_files[1]), 
               raster(AHN_files[2]), 
               raster(AHN_files[3]), 
               raster(AHN_files[4]),
               raster(AHN_files[5]),
               raster(AHN_files[6]))

# Reproject computed coastlines
coastlines = spTransform(coastlines, CRS(proj4string(AHN_5m)))

# Crop AHN to coastlines
AHN_5m = crop(AHN_5m, coastlines)

# Load RWS water table heights
unzip('./data/waterstanden.zip', exdir = data_folder, overwrite = TRUE)
csvfiles = list.files(data_folder, pattern = glob2rx("*.csv"), full.names = TRUE)
waterstanden = read.csv(csvfiles, sep=";")

# RWS measurement stations
RWS_station_coords = matrix(c(576756.115680625, 5691113.19915419, 581080.037560033, 5692490.64235438), nrow = 2, byrow = TRUE)
RWS_station_names = c("baalhoek", "schaar_vd_noord")
RWS_stations = SpatialPointsDataFrame(RWS_station_coords, data.frame(RWS_station_names), proj4string = CRS("+init=epsg:32631"))
RWS_stations = spTransform(RWS_stations, CRS(proj4string(AHN_5m)))
baalhoek_buffer = buffer(RWS_stations[1,], 1500)
schaar_vd_noord_buffer = buffer(RWS_stations[2,], 1500)

# Get coastlines in proximity of measurement stations
cl_baalhoek = crop(coastlines, baalhoek_buffer)
cl_schaar_vd_noord = crop(coastlines, schaar_vd_noord_buffer)

# Get water levels and standard deviations
H_baalhoek =        heightFromCL(cl_baalhoek       , AHN_5m)
H_schaar_vd_noord = heightFromCL(cl_schaar_vd_noord, AHN_5m)

