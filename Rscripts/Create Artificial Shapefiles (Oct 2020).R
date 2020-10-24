# Create a dummy shapefile
#
# Author: Earl Duncan
# Created: 19/05/2017
# Updated: 13/03/2019 Updated points2polygons for row.names
#          28/02/2020 Updated for 20x20 and 50x50
#==========================================================================

# Set working directory and filepaths
fp.wd <- "C:/Users/Duncanew/Creating New Shapefiles"
setwd(fp.wd)
fp.coords <- function(Shapefile) paste0(fp.wd, "/Coordinates - ", Shapefile, ".csv")
fp.shp <- function(Shapefile) paste0(fp.wd, "/Shapefile - ", Shapefile, ".shp")

# Load packages
library(maptools)	# For coordinates()
library(rgdal)		# For writeOGR()

# UDF
points2polygons <- function(df, data){
	get.grpPoly <- function(group, ID, df){
		Polygon(coordinates(df[df$id == ID & df$group == group,]))
	}
	get.spPoly <- function(ID, df){
		Polygons(lapply(unique(df[df$id == ID,]$group), get.grpPoly, ID, df), ID)
	}
	spPolygons <- SpatialPolygons(lapply(unique(df$id),get.spPoly, df))
	row.names(data) <- data$id
	SpatialPolygonsDataFrame(spPolygons, match.ID = TRUE, data = data)
}

#--------------------------------------------------------------------------

# Set shapefile
Shapefile <- "15 x 15"
# Shapefile <- "20 x 20"
# Shapefile <- "50 x 50"

# Read in coordinate data
coords <- read.csv(fp.coords(Shapefile), header = TRUE)

# Convert coordinates from data.frame to SpatialPointsDataFrame
coordinates(coords) <- ~x + y

# Convert SpatialPointsDataFrame to a SpatialPolygonDataFrame
ID <- data.frame(id = unique(coords$id))
map <- points2polygons(coords, ID)

# Save the new shapefile
writeOGR(map, dsn = fp.wd, layer = paste0("Shapefile - ", Shapefile), 
	driver = "ESRI Shapefile", overwrite_layer = TRUE)

# EOF