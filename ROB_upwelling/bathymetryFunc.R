###############################################################################
## DESCRIPTION: Given coordinates, this functions finds bathy data for further use
## USAGE: The data are then used to discern the width and steepness of the continental shelf
## ARGUMENTS:
## DETAILS:
## VALUE:
## AUTHORS(S): Robert Schlegel
## REFERENCE(S):
## EXAMPLE(S):
## DEPENDS: library(RNetCDF)
## USED BY:
##############################################################################
require(FNN)
require(fossil)
require(rgeos)
require(maptools)
require(geosphere)
require(marmap)
source("func/earthdist.R")

##############################################################################
# Load the bathy data for SA
#load("data/bathy/bathy.RData")
# Load the coastline for Africa
#load("graph/africa_coast.RData")
# Load the coastline for southern Africa
#load("graph/southern_africa_coast.RData")
# Load the coastline for southern Africa
#load("graph/south_africa_coast.RData")
#siteList <- read.csv("setupParams/site_list_v3.4.csv")
#site <- siteList[1, ]
#site <- southern_africa_coast2[1,]

##############################################################################
## This function finds shore normal transects for a given area 
  ## Then finds the bathy along those lines so that it may find the shelf width and slope
shelf.area.width.slope <- function(xlim, ylim, bathymetry){
  # Expand range to avoid edge constraints from coastal shape file edges
  #range_lon <- xlim[2]-xlim[1]; range_lat <- ylim[2]-ylim[1]
  xlim2 <- c(xlim[1]-2, xlim[2]+2); ylim2 <- c(ylim[1]-2, ylim[2]+2)
  if(xlim[1] < -20) {
    xlim[1] <- xlim[1]+360
  } else {
    xlim[1] <- xlim[1]
  }
  if(xlim[2] < -20) {
    xlim[2] <- xlim[2]+360
  } else {
    xlim[2] <- xlim[2]
  }
  # Possible coastline groups
  coast_group <- c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1) # Each continent has a value
  # Extract coastal shapefile section for use guiding transects
  coast_i <- fortify(getRgshhsMap("coords_to_extract/gshhs_i.b", xlim = xlim, ylim = ylim))
  coast_i <- droplevels(subset(coast_i, group %in% coast_group)) # Remove islands
  coast_i <- coast_i[1:(length(coast_i$group)-3), ] # Remove last three rows which are created to make a complete path for plotting
  coast_c <- fortify(getRgshhsMap("coords_to_extract/gshhs_c.b", xlim = xlim, ylim = ylim2))
  coast_c <- droplevels(subset(coast_c, group %in% coast_group)) 
  coast_c <- coast_c[1:(length(coast_c$group)-3), ]
  # Find the site on the coastline and it's nearest neighbour points
  shelf_stats <- data.frame()
  for(i in 1:length(coast_i$lat)){
    coords <- coast_i[i, 1:2]
    coords2 <- knnx.index(coast_c[,1:2], as.matrix(coords), k = 1)
    coords3 <- data.frame(coast_c[c(coords2-1, coords2+1),]) # This line controls how wide the spread of points used to determine the shore-normal transect is. A wider spread gives a more even, less sensitive transect
    coords3 <- coords3[2:1,] # Flip so that shore is pointing the correct direction... not certain this will work the same all over the world. May need to come back to this.
    # Define the shore normal transect bearing
    heading <- earth.bear(coords3[1,1], coords3[1,2], coords3[2,1], coords3[2,2]) + 90
    if(heading >= 360){
      heading <- heading-360
    } else {
      heading <- heading
    }
    # Find width using get.transect function... not currently working
    #distance <- as.data.frame(destPoint(p = coords, b = heading, d = 400000)) # The coordinates for a point 400km from the shore along the correct shore normal transect
    #transect <- get.transect(as.bathy(bathy), 18, -30, 15, -30, distance = TRUE)
    #transect <- get.transect(bathy, distance[1], distance[2], coords[1], coords[2], distance = TRUE)
    # Find depth along bearing for 400km using fastest nearest neighbour search
    distances <- seq(from = 0, to = 400000, by = 1000)
    distances2 <- as.data.frame(destPoint(p = coords, b = heading, d = distances))
    sitesIdx <- knnx.index(bathymetry[,1:2], as.matrix(distances2), k = 1)
    bathy2 <- bathymetry[sitesIdx,]
    bathy2 <- bathy2[complete.cases(bathy2[,3]),]
    if(length(bathy2$depth) < 2){
      info <- data.frame(coords, heading, width_200 = NA, depth_200 = NA, angle_200 = NA,
                         width_2000 = NA, depth_2000 = NA, angle_2000 = NA,
                         width_4000 = NA, depth_4000 = NA, angle_4000 = NA)
    } else {
      bathy3 <- data.frame(sapply(bathy2[1:2], deg2rad), depth = bathy2$depth)
      # Find the point at which the depth first hits/ exceeds 200m
      shelf_end <- bathy3[bathy3$depth <= -200,][1,]
      # Calculate how far it takes to reach the 200m isobath
      shelf_width <- gcd.hf(bathy3$lon[1], bathy3$lat[1], shelf_end$lon, shelf_end$lat)
      # From this calculate the shelf slope angle
      shelf_depth <- abs(shelf_end$depth)
      shelf_angle <- tan(shelf_width/shelf_depth)
      # Find the point at which the depth first hits/ exceeds 2000m
      rise_2000_end <- bathy3[bathy3$depth <= -2000,][1,]
      # Calculate how far it takes to reach the 2000m isobath from the 200m isobath
      rise_2000_width <- gcd.hf(shelf_end$lon, shelf_end$lat, rise_2000_end$lon, rise_2000_end$lat)
      # From this calculate the shelf slope angle
      rise_2000_depth <- abs(rise_2000_end$depth - shelf_end$depth)
      rise_2000_angle <- tan(rise_2000_width/rise_2000_depth)
      # Find the point at which the depth first hits/ exceeds 4000m
      rise_4000_end <- bathy3[bathy3$depth <= -4000,][1,]
      # Calculate how far it takes to reach the 4000m isobath from the 200m isobath
      rise_4000_width <- gcd.hf(shelf_end$lon, shelf_end$lat, rise_4000_end$lon, rise_4000_end$lat)
      # From this calculate the shelf slope angle
      rise_4000_depth <- abs(rise_4000_end$depth - shelf_end$depth)
      rise_4000_angle <- tan(rise_4000_width/rise_4000_depth)
      # Put it all together
      info <- data.frame(coords, heading, width_200 = shelf_width, depth_200 = shelf_depth, angle_200 = shelf_angle,
                         width_2000 = rise_2000_width, depth_2000 = rise_2000_depth, angle_2000 = rise_2000_angle,
                         width_4000 = rise_4000_width, depth_4000 = rise_4000_depth, angle_4000 = rise_4000_angle)
    }
    shelf_stats <- rbind(shelf_stats, info)
  }
  return(shelf_stats)
}

##############################################################################
## This function finds the shore normal transect, then finds the bathy along that line
  ## So that it may find the shelf width and slope
    ## NB: "site" must have columns: site, lon and lat in no particular order
shelf.width.slope <- function(site){
  # Find the site on the coastline and it's nearest neighbour points
  coords <- data.frame(lon = site$lon, lat = site$lat)
  coords2 <- knnx.index(africa_coast[,1:2], as.matrix(coords), k = 1)
  coords3 <- data.frame(site = site$site, africa_coast[c(coords2-1, coords2+1),]) # This line controls how wide the spread of points used to determine the shore-normal transect is. A wider spread gives a more even, less sensitive transect
  coords3 <- coords3[2:1,1:3]
  # Define the shore normal transect bearing
  heading <- earth.bear(coords3[1,2], coords3[1,3], coords3[2,2], coords3[2,3]) + 90
  if(heading >= 360){
    heading <- heading-360
  } else {
    heading <- heading
  }
  # Find depth along bearing for 400km
  distances <- seq(from = 0, to = 400000, by = 1000)
  distances2 <- as.data.frame(destPoint(p = coords, b = heading, d = distances))
  sitesIdx <- knnx.index(bathy[,1:2], as.matrix(distances2), k = 1)
  bathy2 <- bathy[sitesIdx,]
  bathy2 <- bathy2[complete.cases(bathy2$depth),]
  if(length(bathy2$depth) < 2){
    info <- data.frame(site = site$site, coords, width_200 = NA, depth_200 = NA, angle_200 = NA,
                       width_2000 = NA, depth_2000 = NA, angle_2000 = NA,
                       width_4000 = NA, depth_4000 = NA, angle_4000 = NA)
  } else {
    bathy3 <- data.frame(sapply(bathy2[1:2], deg2rad), depth = bathy2$depth)
    # Find the point at which the depth first hits/ exceeds 200m
    shelf_end <- bathy3[bathy3$depth <= -200,][1,]
    # Calculate how far it takes to reach the 200m isobath
    shelf_width <- gcd.hf(bathy3$lon[1], bathy3$lat[1], shelf_end$lon, shelf_end$lat)
    # From this calculate the shelf slope angle
    shelf_depth <- abs(shelf_end$depth)
    shelf_angle <- tan(shelf_width/shelf_depth)
    # Find the point at which the depth first hits/ exceeds 2000m
    rise_2000_end <- bathy3[bathy3$depth <= -2000,][1,]
    # Calculate how far it takes to reach the 2000m isobath from the 200m isobath
    rise_2000_width <- gcd.hf(shelf_end$lon, shelf_end$lat, rise_2000_end$lon, rise_2000_end$lat)
    # From this calculate the shelf slope angle
    rise_2000_depth <- abs(rise_2000_end$depth - shelf_end$depth)
    rise_2000_angle <- tan(rise_2000_width/rise_2000_depth)
    # Find the point at which the depth first hits/ exceeds 4000m
    rise_4000_end <- bathy3[bathy3$depth <= -4000,][1,]
    # Calculate how far it takes to reach the 4000m isobath from the 200m isobath
    rise_4000_width <- gcd.hf(shelf_end$lon, shelf_end$lat, rise_4000_end$lon, rise_4000_end$lat)
    # From this calculate the shelf slope angle
    rise_4000_depth <- abs(rise_4000_end$depth - shelf_end$depth)
    rise_4000_angle <- tan(rise_4000_width/rise_4000_depth)
    # Put it all together
    info <- data.frame(site = site$site, coords, width_200 = shelf_width, depth_200 = shelf_depth, angle_200 = shelf_angle,
                       width_2000 = rise_2000_width, depth_2000 = rise_2000_depth, angle_2000 = rise_2000_angle,
                       width_4000 = rise_4000_width, depth_4000 = rise_4000_depth, angle_4000 = rise_4000_angle)
    # Retun
  }
  return(info)
}

##############################################################################
## This function finds the shore normal transect only
  ## NB: "siteList" mus t have columns: site, lon and lat in no particular order
shore.normal.transect <- function(site){
  # Find the site on the coastline and it's nearest neighbour points
  coords <- data.frame(lon = site$lon, lat = site$lat)
  coords2 <- knnx.index(africa_coast[,1:2], as.matrix(coords), k = 1)
  coords3 <- data.frame(site = site$site, africa_coast[c(coords2-1, coords2+1),]) # This line controls how wide the spread of points used to determine the shore-normal transect is. A wider spread gives a more even, less sensitive transect
  coords3 <- coords3[2:1,1:3]
  # Define the shore normal transect bearing
  heading <- earth.bear(coords3[1,2], coords3[1,3], coords3[2,2], coords3[2,3]) + 90
  if(heading >= 360){
    heading <- heading-360
  } else {
    heading <- heading
  }
  distances <- seq(from = 0, to = 300000, by = 1000)
  distances2 <- as.data.frame(destPoint(p = coords, b = heading, d = distances))
  # Put it all together
  info <- data.frame(site = site$site, lon = distances2$lon, lat = distances2$lat)
  return(info)
}

##############################################################################
### Testing

## Test loops
#transects_site <- data.frame()
#for(i in 1:length(siteList$site)){
#  site <- siteList[i,]
#  transect <- shore.normal.transect(site)
#  transects_site  <- rbind(transects_site, transect)
#}

#africa_coast2 <- africa_coast
#africa_coast2 <- subset(africa_coast2, group == "1.1") # Keep only coastline
#africa_coast2$site <- seq(1, length(africa_coast2$long))
#africa_coast2$lon <- africa_coast2$long
#africa_coast2 <- africa_coast2[3:length(africa_coast2$long)-1, ] # Manually cut off the last bits
#transects_coast <- data.frame()
#for(i in 1:length(africa_coast2$site)){
#  site <- africa_coast2[i,]
#  transect <- shore.normal.transect(site)
#  transects_coast  <- rbind(transects_coast, transect)
#}

## Test graph...
#ggplot(data = coastline, aes(x = long, y = lat)) +
#  stat_contour(data = bathy, aes(x = lon, y = lat, z = depth, alpha = ..level..), 
#               col = "black", size = 0.2, binwidth = 200) +
#  geom_path(aes(x = long, y = lat, group = group), 
#            fill = "grey80", colour = "black", size = 0.2, show_guide = FALSE) +
#  geom_path(data = sa_provinces_new, aes(x = long, y = lat, group = group), size = 0.2) +
#  #geom_point(data = transects_site, aes(x = lon, y = lat), size = 0.2) +
#  geom_point(data = transects_coast, aes(x = lon, y = lat), size = 0.2) +
#  #geom_point(data = shelf_stats_site, aes(x = lon, y = lat)) +
#  guides(font = "bold") + coord_equal() +
#  scale_alpha_continuous(breaks = c(-200, -500, -1000, -1500, -2000, -3000, -4000),
#                         guide_legend(title = "Depth (m)")) +
#  coord_map(xlim = c(12, 36), ylim = c(-38, -26.8), projection = "mercator")
