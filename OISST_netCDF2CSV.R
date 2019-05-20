# netCDF2csv
# OISST 

# Loading libraries
library(ncdf4) # library for processing netCDFs
library(ncdf4)
library(data.table)
library(tidyverse)
library(reshape2)
library(plyr)
library(lubridate)
library(stringr)
library(doMC); doMC::registerDoMC(cores = 4)

# bbox were determined using the parameters in the following paper
# Reduced Nearshore Warming Associated With Eastern Boundary Upwelling Systems 
# Rui Seabra 1 , Rubén Varela 2 , António M. Santos 1,3 , Moncho Gómez-Gesteira 2 ,
# Claudia Meneghesso 1,3 , David S. Wethey 4 and Fernando P. Lima 1 *

# Under Pressure: Climate Change, Upwelling, and Eastern Boundary  Upwelling Ecosystems
# Marisol García-Reyes 1 *, William J. Sydeman 1 , David S. Schoeman 2 ,
# Ryan R. Rykaczewski 3 , Bryan A. Black 4 , Albertus J. Smit 5 and Steven J. Bograd 6

bbox <- data.frame(BC =  c(-35, -15, 15, 20), # Benguela Current
                   CC = c(18, 45, -23, -10), # Canary Current
                   CalC = c(25, 48, -130, -110), # California Current
                   HC = c(-40, -10 , -80, -72), # Humboldt Current
                   row.names = c("latmin", "latmax", "lonmin", "lonmax"))

# Setting path
nc.dir <- "/home/amieroh/Documents/spatial/data/netCDF"
csv.dir <- "/home/amieroh/Documents/spatial/data/csv"

read_nc <- function(ncFile, region = region, csvDir = csvDir) {
  coords <- bbox[, region]
  nc <- nc_open(ncFile)
  pathLen <- nchar(nc.dir) + 1 
  fNameStem <-
    substr(ncFile, pathLen + 1, pathLen + 13)
  fDate <- substr(ncFile, pathLen + 15, pathLen + 22)
  LatIdx <- which(nc$dim$lat$vals > coords[1] & nc$dim$lat$vals < coords[2])
  LonIdx <- which(nc$dim$lon$vals > coords[3] & nc$dim$lon$vals < coords[4])
  sst <- ncvar_get(nc,
                   varid = "sst",
                   start = c(LonIdx[1], LatIdx[1], 1, 1),
                   count = c(length(LonIdx), length(LatIdx), 1, 1)) %>%
    round(4)
  dimnames(sst) <- list(lon = nc$dim$lon$vals[LonIdx],
                        lat = nc$dim$lat$vals[LatIdx])
  nc_close(nc)
  sst <-
    as.data.table(melt(sst, value.name = "temp"), row.names = NULL) %>%
    mutate(t = ymd(fDate)) %>%
    na.omit()
  fwrite(sst,
         file = paste(csvDir, "/", region, "-", fNameStem, ".", strtDate, "-", endDate, ".csv", sep = ""),
         append = TRUE, col.names = FALSE)
  rm(sst)
}

# the list of files
ncList <- list.files(path = nc.dir, pattern = "*.nc", full.names = TRUE, include.dirs = TRUE)
strtDate <- str_sub(ncList[1], start = 15, end = 22)
endDate <- str_sub(ncList[length(ncList)], start = 15, end = 22)

# apply the function
system.time(llply(ncList, read_nc, region = "BC", csvDir = csv.dir, .parallel = TRUE))
system.time(llply(ncList, read_nc, region = "CC", csvDir = csv.dir, .parallel = TRUE))
system.time(llply(ncList, read_nc, region = "CalC", csvDir = csv.dir, .parallel = TRUE))
system.time(llply(ncList, read_nc, region = "HC", csvDir = csv.dir, .parallel = TRUE))
