#Metadata

# Import packages
library(raster)
library(rtiff)
library(sp)
library(rgdal)
library(gdalUtils)
library(lubridate)
library(ggplot2)
library(hydroGOF)

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
data_URL = "https://drive.google.com/uc?export=download&id=1hktrL56zI6murwzCqg_m1TXxI-dr9-kB"
download.file(url=data_URL, destfile='./data/projectdata.zip', method = 'auto', overwrite = TRUE)
unzip('./data/projectdata.zip', exdir = data_folder, overwrite = TRUE)
tiffile = list.files('./data/geoscripting_data/', pattern = glob2rx("mndwi*.tif"), full.names = TRUE)
MNDWI = stack(tiffile)
MNDWI_info = gdalinfo(tiffile)

# Get dates of bands
MNDWI_dates = MNDWI_info[grepl(glob2rx('*LC08*'), MNDWI_info)]
MNDWI_dates = substring(MNDWI_dates, 45)
MNDWI_dates = as.POSIXct(MNDWI_dates, format="%Y-%m-%dT%H-%M-%S")
MNDWI_dates = round_date(MNDWI_dates, "10 mins")

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
AHN= merge(raster(AHN_files[1]), 
           raster(AHN_files[2]), 
           raster(AHN_files[3]), 
           raster(AHN_files[4]),
           raster(AHN_files[5]),
           raster(AHN_files[6]))

# Reproject according to AHN
coastlines = spTransform(coastlines, CRS(proj4string(AHN)))
MNDWI = projectRaster(MNDWI, crs = CRS(proj4string(AHN)))
water_classification = projectRaster(water_classification, crs = CRS(proj4string(AHN)))

# Crop AHN to coastlines
AHN = crop(AHN, coastlines)

# Load RWS water table heights
csvfiles = list.files('./data/projectdata.zip', pattern = glob2rx("*.csv"), full.names = TRUE)
waterlevels = read.csv(csvfiles, sep=";", dec = ",")
waterlevels = waterlevels[,c(2,20,21,23,36,37,38)]

# Merge date and time columns
waterlevels$DATUMENTIJD = as.POSIXct(paste(waterlevels$WAARNEMINGDATUM, waterlevels$WAARNEMINGTIJD),
                                      format="%d-%m-%Y %H:%M:%S")
waterlevels = waterlevels[,c(1,4:8)]

# To GMT
waterlevels$DATUMENTIJD = waterlevels$DATUMENTIJD - hours(1)

# Convert from cm to m20200125_011
waterlevels$NUMERIEKEWAARDE[waterlevels$NUMERIEKEWAARDE == 999999999] = NA
waterlevels$NUMERIEKEWAARDE = waterlevels$NUMERIEKEWAARDE/100

# Keep relevant dates
relevant_dates = waterlevels$DATUMENTIJD %in% MNDWI_dates
waterlevels = waterlevels[relevant_dates,]

# Split water levels
H_baalhoek_rws = waterlevels[waterlevels$MEETPUNT_IDENTIFICATIE == "Baalhoek",]
H_schaar_vd_noord_rws = waterlevels[waterlevels$MEETPUNT_IDENTIFICATIE == "Schaar van de Noord",]

# RWS measurement stations
RWS_station_coords = matrix(c(H_baalhoek_rws$X[1]       , H_baalhoek_rws$Y[1],
                              H_schaar_vd_noord_rws$X[1], H_schaar_vd_noord_rws$Y[1]),
                            nrow = 2, byrow = TRUE)
RWS_station_names = c("baalhoek", "schaar_vd_noord")
RWS_stations = SpatialPointsDataFrame(RWS_station_coords,
                                      data.frame(RWS_station_names), 
                                      proj4string = CRS(paste("+init=epsg:", H_baalhoek_rws$EPSG[1], sep = "")))
RWS_stations = spTransform(RWS_stations, CRS(proj4string(AHN)))
baalhoek_buffer = buffer(RWS_stations[1,], 1000)
schaar_vd_noord_buffer = buffer(RWS_stations[2,], 1000)

# Get coastlines in proximity of measurement stations
cl_baalhoek = crop(coastlines, baalhoek_buffer)
cl_schaar_vd_noord = crop(coastlines, schaar_vd_noord_buffer)

# Get water levels and standard deviations
H_baalhoek =        heightFromCL(cl_baalhoek       , AHN)
H_schaar_vd_noord = heightFromCL(cl_schaar_vd_noord, AHN)

# Create dataframe with results 
stations        = c(as.character(H_baalhoek_rws$MEETPUNT_IDENTIFICATIE), as.character(H_schaar_vd_noord_rws$MEETPUNT_IDENTIFICATIE))
H_RWS           = c(H_baalhoek_rws$NUMERIEKEWAARDE, H_schaar_vd_noord_rws$NUMERIEKEWAARDE)
H_Landsat       = c(H_baalhoek$Mean, H_schaar_vd_noord$Mean[MNDWI_dates %in% H_schaar_vd_noord_rws$DATUMENTIJD])
H_Landsat_STDEV = c(H_baalhoek$STDEV, H_schaar_vd_noord$STDEV[MNDWI_dates %in% H_schaar_vd_noord_rws$DATUMENTIJD])
H = data.frame(stations, H_RWS, H_Landsat, H_Landsat_STDEV)

limits = c(min(c(H$H_RWS, H$H_Landsat), na.rm = TRUE), 
           max(c(H$H_RWS, H$H_Landsat), na.rm = TRUE))
NSE_overall         = NSE(H$H_Landsat, H$H_RWS)
NSE_baalhoek        = round(NSE(H$H_Landsat[H$stations == "Baalhoek"], H$H_RWS[H$stations == "Baalhoek"]), digits = 2)
NSE_schaar_vd_noord = round(NSE(H$H_Landsat[H$stations == "Schaar van de Noord"], H$H_RWS[H$stations == "Schaar van de Noord"]), digits = 2)

dat_text <- data.frame(
  label = c(paste("NSE =", NSE_baalhoek), paste("NSE =", NSE_schaar_vd_noord)),
  stations   = c("Baalhoek", "Schaar van de Noord")
)

# Plotting results and intermediate steps
ggplot(H, aes(x = H_RWS, y = H_Landsat)) + 
  geom_point() +
  geom_errorbar(aes(ymin = H_Landsat - H_Landsat_STDEV, ymax = H_Landsat + H_Landsat_STDEV)) +
  geom_abline(intercept = 0) +
  coord_equal(xlim = limits, ylim = limits) +
  facet_grid(. ~ stations) + 
  labs(x = "Water Level Rijkswaterstaat (m)", y = "Water Level Landsat (m)") +
  ggtitle("Model performance") +
  geom_text(dat_text, mapping = aes(x = -2.5, y = 2.25, label = label), size = 3.5) +
  theme_linedraw() 
ggsave("output/model_performance.png", plot = last_plot())

orangetoblue = colorRampPalette(c("Orange", "White", "Blue"))
png('output/MNDWI_Saeftinghe.png')
plot(MNDWI[[1]], xlab = "x (m)", ylab = "y (m)", main = "MNDWI of Saeftinghe", 
     col = orangetoblue(100), zlim = c(-1,1))
text(66500, 377900, MNDWI_dates[1], col = "white")
dev.off()

png('output/Water_Saeftinghe.png')
plot(water_classification[[1]], xlab = "x (m)", ylab = "y (m)", main = "Water classification Saeftinghe" , 
   col = orangetoblue(100), legend = FALSE)
text(66500, 377900, MNDWI_dates[1], col = "white")
legend("bottomright", legend=c("land", "water"), fill=orangetoblue(2), bg="white")
dev.off()


png('output/coastline.png')
plot(coastlines[1,], xlab = "x (m)", ylab = "y (m)", main = "Coastline Saeftinghe (2013-09-30)", axes = TRUE)
dev.off()

png('output/coastline_with_buffer.png')
plot(coastlines[1,], xlab = "x (m)", ylab = "y (m)", main = "Coastline Saeftinghe (2013-09-30)", axes = TRUE)
plot(RWS_stations, add = TRUE)
plot(cl_baalhoek[1,], col = 'red', add = TRUE)
plot(baalhoek_buffer, add= TRUE)
plot(cl_schaar_vd_noord[1,], col = 'red', add = TRUE)
plot(schaar_vd_noord_buffer, add= TRUE)
dev.off()

png('output/coastline_with_AHN.png')
plot(AHN, xlab = "x (m)", ylab = "y (m)", main = "AHN with coastline (2013-09-30)")
plot(RWS_stations, add = TRUE)
plot(coastlines[1,], add = TRUE)
plot(cl_baalhoek[1,], col = 'red', add = TRUE)
plot(baalhoek_buffer, add= TRUE)
plot(cl_schaar_vd_noord[1,], col = 'red', add = TRUE)
plot(schaar_vd_noord_buffer, add= TRUE)
dev.off()
