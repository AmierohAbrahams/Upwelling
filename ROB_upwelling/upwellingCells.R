#############################################################################
## This script does:
# 1. Load in data;
  # https://cran.r-project.org/web/packages/marmap/vignettes/marmap-DataAnalysis.pdf
  # Then save bathymetry data frames for future use;
# 2. Calculate shelf stats;
  # Then save stats for future use;
# 3. Prep for plotting;
# 4. Plot the bathymetry map and ribbon;
# 5. Saves as...
#############################################################################

#############################################################################
## DEPENDS ON:
require(ggplot2); require(grid); require(plyr); require(FNN); require(marmap); require(PBSmapping)
#require(maps); require(mapdata)
source("func/bathymetryFunc.R"); source("func/statsFunc.R")
#"setupParams/site_list_v3.4.csv"
#"graph/southern_africa_coast.RData"
#############################################################################

#############################################################################
## USED BY:
# Nothing
#############################################################################

#############################################################################
## CREATES:
#"graph/upwelling/Benguela.pdf"
#"graph/upwelling/Humboldt.pdf"
#"graph/upwelling/Canary.pdf"
#"graph/upwelling/California.pdf"
#############################################################################

##########################################################################
## Load data

# Benguela
b_lat <- c(-32.5, -22.5); b_lon <- c(5, 20)
#b_bathy <- as.xyz(getNOAA.bathy(lon1 = b_lon[1], lon2 = b_lon[2], lat1 = b_lat[1], lat2 = b_lat[2], resolution = 4))
#colnames(b_bathy) <- c("lon", "lat", "depth")
#save(b_bathy, file = "data/bathy/b_bathy.RData")
load("data/bathy/b_bathy.RData")
#plot(as.bathy(b_bathy), image = T)

# Humboldt
h_lat <- c(-17.5, -7.5); h_lon <- c(-85, -70)
#h_bathy <- as.xyz(getNOAA.bathy(lon1 = h_lon[1], lon2 = h_lon[2], lat1 = h_lat[1], lat2 = h_lat[2], resolution =  4))
#colnames(h_bathy) <- c("lon", "lat", "depth")
#save(h_bathy, file = "data/bathy/h_bathy.RData")
load("data/bathy/h_bathy.RData")
#plot(as.bathy(h_bathy), image = T)

# Canary
cn_lat <- c(25, 35); cn_lon <- c(-20, -5)
#cn_bathy <- as.xyz(getNOAA.bathy(lon1 = cn_lon[1], lon2 = cn_lon[2], lat1 = cn_lat[1], lat2 = cn_lat[2], resolution =  4))
#colnames(cn_bathy) <- c("lon", "lat", "depth")
#save(cn_bathy, file = "data/bathy/cn_bathy.RData")
load("data/bathy/cn_bathy.RData")
#plot(as.bathy(cn_bathy), image = T)

# California
cl_lat <- c(35, 45); cl_lon <- c(-135, -120)
#cl_bathy <- as.xyz(getNOAA.bathy(lon1 = cl_lon[1], lon2 = cl_lon[2], lat1 = cl_lat[1], lat2 = cl_lat[2], resolution =  4))
#colnames(cl_bathy) <- c("lon", "lat", "depth")
#save(cl_bathy, file = "data/bathy/cl_bathy.RData")
load("data/bathy/cl_bathy.RData")
#plot(as.bathy(cl_bathy), image = T)


bbox <- data.frame(BC = c(-35, -25, 15, 20), # Benguela Current
                   CC = c(25, 35, -20, -5), # Canary Current
                   CalC = c(30, 50, -115, -130), # California Current
                   HC = c(-17.5, -7.5, -85, -70), # Humboldt Current
                   row.names = c("latmin", "latmax", "lonmin", "lonmax"))




##########################################################################
## Prep for bathy plot

# Extract the coastline and border shape files
BC_coastline <- fortify(getRgshhsMap("/home/amieroh/Documents/Masters_2019/Upwelling/gshhs_i.b", xlim = b_lon, ylim = b_lat))
BC_borders <- fortify(importGSHHS("/home/amieroh/Documents/Masters_2019/Upwelling/wdb_borders_i.b", xlim = b_lon, ylim = b_lat)) # Doesn't plot well...
# BC_borders <- data.frame(lon = map("worldHires", xlim = b_lon, ylim = b_lat, boundary = FALSE, plot = FALSE)$x,
#                   lat = map("worldHires", xlim = b_lon, ylim = b_lat, boundary = FALSE, plot = FALSE)$y) # Also not...

HC_coastline <- fortify(getRgshhsMap("/home/amieroh/Documents/Masters_2019/Upwelling/gshhs_i.b", xlim = c(h_lon[1]+360, h_lon[2]+360), ylim = h_lat))
HC_borders <- fortify(importGSHHS("/home/amieroh/Documents/Masters_2019/Upwelling/wdb_borders_i.b", xlim = h_lon, ylim = h_lat))

CC_coastline <- fortify(getRgshhsMap("/home/amieroh/Documents/Masters_2019/Upwelling/gshhs_i.b", xlim = cn_lon, ylim = cn_lat))
#cn_borders <- fortify(importGSHHS("/home/amieroh/Documents/Masters_2019/Upwelling/gshhs_i.b", xlim = cn_lon, ylim = cn_lat))

CalC_coastline <- fortify(getRgshhsMap("/home/amieroh/Documents/Masters_2019/Upwelling/gshhs_i.b", xlim = c(cl_lon[1]+360, cl_lon[2]+360), ylim = cl_lat))
#cl_borders <- fortify(importGSHHS("/home/amieroh/Documents/Masters_2019/Upwelling/gshhs_i.b", xlim = cl_lon, ylim = cl_lat))

# Functions to bin depths for better plotting
round_200 <- function(x){
  if(x > -200){
    x <- -200
  } else{
    x <- x
  }
}

cleanBathy <- function(x){
  x2 <- x[complete.cases(x$depth),]
  x2$depth <- sapply(x2$depth, round_200)
  x2$depth <- round_any(x2$depth, 200)
  return(x2)
}

b_bathy2 <- cleanBathy(b_bathy)
h_bathy2 <- cleanBathy(h_bathy)
cn_bathy2 <- cleanBathy(cn_bathy)
cl_bathy2 <- cleanBathy(cl_bathy)

##########################################################################
## Bathy plot function
bathyPlot <- function(bathymetry, coastline, xlim, ylim){
  ggplot(data = coastline, aes(x = long, y = lat)) + theme_bw() +
    geom_raster(data = bathymetry, aes(x = lon, y = lat, fill = depth)) +
    stat_contour(data = bathymetry, aes(x = lon, y = lat, z = depth, alpha = ..level..), 
                 colour = "black", size = 0.2, binwidth = 200, na.rm = TRUE, show_guide = FALSE) +
    geom_polygon(aes(group = group), fill = "grey80", colour = "black", size = 0.2, show_guide = FALSE) +
    #geom_path(data = borders, aes(x = lon, y = lat)) +
    guides(font = "bold", colour = "White") + coord_equal() +
    scale_alpha_continuous(breaks = c(-200, -1000, -2000, -3000, -4000, -5000),
                           guide_legend(title = "Depth (m)")) +
    scale_fill_gradient(low = "steelblue4", high = "steelblue1", na.value = "steelblue4", 
                        breaks = c(-1000, -2000, -3000, -4000, -5000),
                        guide_legend(title = "Depth (m)")) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    coord_equal() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          axis.title = element_blank(),
          #axis.text = element_blank(),
          #axis.ticks = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.title = element_blank(),
          legend.text = element_text(size = 7, colour = "White"),
          legend.title = element_text(size = 7, colour = "White"),
          legend.key = element_rect(colour = NA, size = 0.2),
          legend.key.height = unit(0.4, "cm"),
          legend.background = element_blank(),
          legend.justification = c(1,0), 
          legend.position = c(0.17, 0.07))
}

##########################################################################
# Currently unused

## The scale bar and north arrow
# Create parameters for a scale bar:
scales <-  c(0, 100000, 200000) # Sets initial distances for scale in metres
scaleBar <- data.frame("DISTANCE" = scales/1000, # Changes the distance to look like km
                       "b.lon" = 29.5, # Set beginning longitude point for scale bar
                       "b.lat" = -34.8, # Set beginning latitude point for scale bar
                       destPoint(p = c(29.5, -34.8), b = 90, d = scales)) # Set start point, direction and length of scale bar
scaleLength <- scaleBar[3,] # The ending point

# Create parameters for a North arrow:
nArrow <- data.frame("b.lon" = scaleBar[3,4],
                     "b.lat" = -34.4,
                     destPoint(p = c(scaleBar[3,4], -34.4), b = 0, d = 50000))

### Insert city locations and labels
geom_point(data = sub_dist, 
           aes(lon, lat), shape = 16, size = 2.0, colour = "White", alpha = 0.5) + # Insert major cities
  geom_text(data = westSites, aes(lon, lat, label = site),
            hjust = 1.2, vjust = 0.3, size = 2.2, colour = "White") + # West coast cities
  geom_text(data = southSites, aes(lon, lat, label = site),
            hjust = 0.3, vjust = 2, size = 2.2, colour = "White") + # South coast cities
  geom_text(data = eastSites, aes(lon, lat, label = site),
            hjust = -0.2, vjust = 0.4, size = 2.2, colour = "White") + # East coast cities
  ### Scale bar and N arrow
  # Insert scale bar bottom:
  geom_segment(data = scaleLength, 
               aes(x = b.lon, y = b.lat, xend = lon, yend = b.lat), colour = "White", size = 0.3) +
  # Insert scale bar tips:
  geom_segment(data = scaleBar, 
               aes(x = lon, y = b.lat, xend = lon, yend = b.lat + 0.1), colour = "White", size = 0.3) +
  # Label distances:
  annotate("text", label = c("0", "100", "200"), 
           x = scaleBar[,4], y = scaleBar[,3] + 0.2, colour = "White", size = 2) +
  annotate("text", label = "km", x = scaleBar[3,4] + .2, y = scaleBar[1,3], colour = "White", size = 2) + 
  # Insert N arrow:
  geom_segment(data = nArrow, aes(x = b.lon, y = b.lat, xend = lon, yend = lat), 
               arrow = arrow(length = unit(0.2, "cm")), colour = "White") +
  annotate("text", label = "N", x = nArrow$b.lon + .15, y = nArrow$b.lat + .1, size = 2, colour = "White") + # Label N arrow
  ### End insert scale bar and arrow.

map("world", xlim = b_lon, ylim = b_lat, interior = FALSE)
map("worldHires", xlim = b_lon, ylim = b_lat, interior = FALSE)
map2 <- data.frame(lon = map("world", xlim = b_lon, ylim = b_lat, interior = FALSE, plot = FALSE)$x, 
                   lat = map("world", xlim = b_lon, ylim = b_lat, interior = FALSE, plot = FALSE)$y)
coastOrder <- order(coastline_inter$lat, coastline_inter$lon)
coastline_inter <- coastline_inter[coastOrder,]

ggplot(data = coastline_inter, aes(x = lon, y = lat)) + theme_bw() + geom_path()
h_coastline2 <- droplevels(subset(h_coastline, group == 3.1))
ggplot(data = h_coastline2, aes(x = long, y = lat)) + theme_bw() + geom_path()

ne10 <- read.oce("coords_to_extract/ne_10m/ne_10m_coastline.shp")

map3 <- data.frame(lon = ne10@data$longitude, lat = ne10@data$longitude)

##########################################################################
## Create a bathy plot for each cell

# Benguela Plot
b_plot <- bathyPlot(b_bathy2, b_coastline, b_lon, b_lat)
b_plot

# Humboldt Plot
h_plot <- bathyPlot(h_bathy2, h_coastline, h_lon, h_lat)
h_plot

# Canary Plot
cn_plot <- bathyPlot(cn_bathy2, cn_coastline, cn_lon, cn_lat)
#cn_plot + xlim(cn_lon)+ylim(cn_lat) # Testing
cn_plot

# California Plot
cl_plot <- bathyPlot(cl_bathy2, cl_coastline, cl_lon, cl_lat)
cl_plot

##########################################################################
### Prep for bathy ribbon

# Calculate shelf stats
source("func/bathymetryFunc.R")
b_stats <- cbind(site = rep("Benguela"), shelf.area.width.slope(b_lon, b_lat, b_bathy))
h_stats <- cbind(site = rep("Humboldt"), shelf.area.width.slope(h_lon, h_lat, h_bathy))
cn_stats <- cbind(site = rep("Canary"), shelf.area.width.slope(cn_lon, cn_lat, cn_bathy))
cl_stats <- cbind(site = rep("California"), shelf.area.width.slope(cl_lon, cl_lat, cl_bathy))
cl_stats <- cl_stats[c(1:75, 164:291), ] # Remove spurious bay data


# Function that calculates distance from each coastal point
coastalDist <- function(x){
  x2 <- data.frame(x$site, x$long, x$lat)
  x2[2:3] <- sapply(x2[2:3], deg2rad)
  x2$dist <- PairsDists(x2)
  a <- rownames(x2)
  b <- round(x2$dist, 1)
  c <- b[1:length(b)-1,]
  d <- c(0, c)
  e <- cumsum(d)
  x2$dist <- d # Corrected distances
  x2$cum_dist <- e # Cummulative distances
  x <- cbind(x, x2[,4:5])
}

# Calculate distances between points along coast
b_stats2 <- coastalDist(b_stats)
h_stats2 <- coastalDist(h_stats)
cn_stats2 <- coastalDist(cn_stats)
cl_stats2 <- coastalDist(cl_stats)

##########################################################################
## Bathy ribbon function
bathyRibbon <- function(x_stats){
  ggplot() + theme_bw() +
    geom_path(data = x_stats, aes(x = width_200, y = cum_dist, colour = angle_2000), 
              lineend = "round", linejoin = "round") +
    #geom_point(aes(colour = shelf_stats_coast2$angle)) + # Makes plot line show temp as colours
    geom_vline(data = x_stats, linetype = "dashed", size = 0.5, 
               aes(xintercept = mean(width_200), colour = mean(angle_2000))) +
    scale_colour_gradient2(low = "blue", mid = "yellow", high = "red", na.value = "grey",
                           midpoint = 0.07,
                           guide_legend(title = "Rise Angle"), 
                           breaks = c(0.01, 0.05, 0.09, 0.13),
                           limits = c(0.0, 0.14)) +
    #guides(colour = "colorbar") +
    #guides(colour = FALSE) +
    xlab(expression(paste("Shelf Width (km)"))) + # y axis label
    scale_x_reverse(breaks = seq(0,200,50),
                    limits = c(250, 0)) +
    #scale_y_continuous(breaks = x_stats$cum_dist) +#, labels = sub_dist$site, expand = c(0, 0)) +
    theme(panel.border = element_blank(), # Creates border
          plot.background = element_blank(), # Makes background transparent
          #panel.grid.major.x = element_line(colour = "black", linetype = "dotted", size = 0.2), # x axis grid lines
          #panel.grid.major.y = element_line(colour = NA, linetype = "dotted",size = 0.2), # y axis grid lines
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), # Removes minor grid lines
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          axis.text.y = element_blank(), # Change y axis text size
          axis.title.y = element_blank(), # Removes y axis label
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = 8),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.key = element_rect(colour = NA, size = 0.2),
          legend.key.height = unit(0.4, "cm"),
          legend.background = element_blank(),
          legend.justification = c(1,0),
          legend.position = c(0.17, 0.07))
}

##########################################################################
## Create a bathy ribbon for each cell

# Benguela Ribbon
b_ribbon <- bathyRibbon(b_stats2)
b_ribbon

# Humboldt Ribbon
h_ribbon <- bathyRibbon(h_stats2)
h_ribbon

# Canary Ribbon
cn_ribbon <- bathyRibbon(cn_stats2)
cn_ribbon

# California Ribbon
cl_ribbon <- bathyRibbon(cl_stats2)
cl_ribbon

##########################################################################
## Stitch them all together

# Benguela Figure
pdf("graph/upwelling/Benguela.pdf", width = 12, height = 6, pointsize = 10) # Set PDF dimensions
vp1 <- viewport(x = 0.16, y = 0.0, w = 1.0, h = 1.0, just = c("left", "bottom")) # Africa
vp2 <- viewport(x = 0.01, y = 0.0, w = 0.3, h = 1.0, just = c("left", "bottom"))  # Ribbon
print(b_plot, vp = vp1)
print(b_ribbon, vp = vp2)
dev.off()

# Humboldt figure
pdf("graph/upwelling/Humboldt.pdf", width = 12, height = 6, pointsize = 10) # Set PDF dimensions
vp1 <- viewport(x = 0.16, y = 0.0, w = 1.0, h = 1.0, just = c("left", "bottom")) # Africa
vp2 <- viewport(x = 0.01, y = 0.0, w = 0.3, h = 1.0, just = c("left", "bottom"))  # Ribbon
print(h_plot, vp = vp1)
print(h_ribbon, vp = vp2)
dev.off()

# Canary figure
pdf("graph/upwelling/Canary.pdf", width = 12, height = 6, pointsize = 10) # Set PDF dimensions
vp1 <- viewport(x = 0.16, y = 0.0, w = 1.0, h = 1.0, just = c("left", "bottom")) # Africa
vp2 <- viewport(x = 0.01, y = 0.0, w = 0.3, h = 1.0, just = c("left", "bottom"))  # Ribbon
print(cn_plot, vp = vp1)
print(cn_ribbon, vp = vp2)
dev.off()

# California figure
pdf("graph/upwelling/California.pdf", width = 12, height = 6, pointsize = 10) # Set PDF dimensions
vp1 <- viewport(x = 0.16, y = 0.0, w = 1.0, h = 1.0, just = c("left", "bottom")) # Africa
vp2 <- viewport(x = 0.01, y = 0.0, w = 0.3, h = 1.0, just = c("left", "bottom"))  # Ribbon
print(cl_plot, vp = vp1)
print(cl_ribbon, vp = vp2)
dev.off()

