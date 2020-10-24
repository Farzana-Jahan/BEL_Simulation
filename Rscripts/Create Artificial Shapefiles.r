# Create a regular grid shapefile
#
# Author: Earl Duncan
# Created: 24/10/2020
# Updated: 
#==========================================================================

# Set working directory and filepaths
fp.wd <- getwd()


# Load packages
library(maptools)	# For coordinates()
library(rgdal)		# For writeOGR()

# UDF for converting points to polygons
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

# UDF for creating coordinates for a region of size N by N
get.coords <- function(N){
    id <- rep(1:(N^2), each = 5)
    x <- c(0, 1, 1, 0, 0) + (id-1) %% N
    y <- c(0, 0, 1, 1, 0) + (id-1) %% N
    order <- 1:5
    d <- data.frame(id, x, y, order, hole = FALSE, piece = id, group = id)
    return(d)
}


#--------------------------------------------------------------------------

# Choose grid size (square root of number of areas)
N <- 50 #(for 100 areas)

# Create coordinate data
coords <- get.coords(N)

# Convert coordinates from data.frame to SpatialPointsDataFrame
coordinates(coords) <- ~x + y

# Convert SpatialPointsDataFrame to a SpatialPolygonDataFrame
ID <- data.frame(id = unique(coords$id))
map <- points2polygons(coords, ID)

# Save the new shapefile
writeOGR(map, dsn = fp.wd, layer = "My New Shapefile_2500areas", 
	driver = "ESRI Shapefile", overwrite_layer = TRUE)

# EOF