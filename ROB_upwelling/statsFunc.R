# Function that calculates the Root Mean Squared Error for multisite data
      ## NB: expects these column labels: site, src, temp, date (in no particular order)
multisiteRMSE <- function(data){
  data <- data[complete.cases(data),]
  RMSE <- data.frame()
  for (i in 1:length(levels(data$site))){
    data2 <- droplevels(subset(data, site == levels(data$site)[i]))
    w <- as.character(levels(as.factor(data2$src))) # Create index of data source names
    m <- combn(length(levels(as.factor(w))), 2)# Create matrix by which to compare all sources
    for(j in 1:ncol(m)){
      x <- subset(data2, src == w[m[1,j]]) # Keep only one source of data
      y <- subset(data2, src == w[m[2,j]])
      x <- subset(x, (date %in% y$date)) # Match overlapping times
      y <- subset(y, (date %in% x$date))
      if (length(x$date) == 0) {
        error <- 0
      } else {
        error <- x$temp-y$temp # Create error value, same as residuals if a linear model were run
        error <- round(rmse(error), 4) # Calculate RMSE
        if (is.na(error) == TRUE){
          error <- 0
        }# Correct NaN values
      }
      overlap <- length(x$date) # Show overlap in days as an integer
      src <- paste(w[m[1,j]], w[m[2,j]], sep = "/") # Show sources
      site <- levels(data$site)[i] # Show site name
      data3 <- data.frame(site, src, error, overlap) # Bind it together
      RMSE <- rbind(RMSE, data3) # Add it to the data.frame
    }
  }
  names(RMSE)[3] <- "RMSE"
  return(RMSE)
}

# Function that returns Root Mean Squared Error
rmse <- function(error){
  sqrt(mean(error^2, na.rm = TRUE))
}

# Function that returns Mean Absolute Error
mae <- function(error){
  mean(abs(error), na.rm = TRUE)
}

# Function that calculates significance of a linear model by site and source
      ## NB: expects these column labels: site, src, temp, date (in no particular order)
significance <- function(data){
  stats <- data.frame()
  for(i in 1:length(levels(data$site))){
    data2 <- subset(data, site == levels(data$site)[i])
    for(j in 1:length(levels(as.factor(data2$src)))){
      data3 <- subset(data2, src == levels(as.factor(data2$src))[j])
      x <- lm(data3$temp ~ seq(1:length(data3$date)))
      y <- round(coef(summary(x))[2,4],3)
      w <- round(summary(x)$adj.r.squared,2)
      v <- coef(summary(x))[2,1] # round(coef(summary(x))[2,1] *365*10, 4) # not used; apply conversion to per decade afterwards for daily, monthly and annual data, as the multiplication will differ...
      z <- data.frame(as.character(data3$site[1]),data3$src[1],y,w,v)
      stats <- rbind(stats, z)
    }
  }
  colnames(stats) <- c("site", "src", "p", "R2", "DT")
  return(stats)
}

# Function that calculates t.test or Wilcoxon depending on homogeneity of variances, also returns R2 and temperature difference for the time compared
      ## NB: expects these column labels: site, src, temp, date (in no particular order)
variance <- function(data){
  data <- data[complete.cases(data),]
  stats <- data.frame()
  for (i in 1:length(levels(data$site))){
    data2 <- droplevels(subset(data, site == levels(data$site)[i]))
    w <- as.character(levels(as.factor(data2$src))) # Create index of data source names
    m <- combn(length(levels(as.factor(w))), 2)# Create matrix by which to compare all sources
    for(j in 1:ncol(m)){
      x <- subset(data2, src == w[m[1,j]]) # Keep only one source of data
      y <- subset(data2, src == w[m[2,j]])
      x <- subset(x, (date %in% y$date)) # Match overlapping times
      y <- subset(y, (date %in% x$date))
      if (length(x$date) == 0) {
        data3 <- data.frame(as.character(data2$site[1]), paste(w[m[1,j]], w[m[2,j]], sep = "/"), 0, 0, 0, "none", 0) 
        colnames(data3) <- c("site", "src", "t", "p", "R2", "test", "temp diff")
        stats <- rbind(stats, data3)
      } else {
        tempx <- (ddply(x, .(site, src), summarize, temp = mean(temp, na.rm = TRUE)))
        tempy <- (ddply(y, .(site, src), summarize, temp = mean(temp, na.rm = TRUE)))
        diff <- tempx$temp - tempy$temp # Calculate temperature difference
        var <- var.test(x$temp, y$temp) # Calculate homogeneity of variances
        if (var$p.value > 0.05) {
          ttest <- t.test(x$temp, y$temp, paired = F, var.equal = T)
        } else if (var$p.value <= 0.05) {
          ttest <- t.test(x$temp, y$temp, paired = F, var.equal = F)
        }
        lmodel <- lm(x$temp ~ y$temp)
        R2 <- round(summary(lmodel)$adj.r.squared,2)
        data3 <- data.frame(as.character(x$site[1]), paste(w[m[1,j]], w[m[2,j]], sep = "/"), 
                            round(ttest$statistic, 3), round(ttest$p.value, 3), R2, "t-test", round(diff, 2))
        colnames(data3) <- c("site", "src", "t", "p", "R2", "test", "temp diff")
        stats <- rbind(stats, data3)
        #} else if (var$p.value <= 0.05) { # Wilcox tests assume the data points being compared are randomly sampled. Meaning this wouldn't be an appropriate test
        #  wilcox <- wilcox.test(x$temp, y$temp)
        #  lmodel <- lm(x$temp ~ y$temp)
        #  R2 <- round(summary(lmodel)$adj.r.squared,2)
        #  data3 <- data.frame(as.character(x$site[1]), paste(w[m[1,j]], w[m[2,j]], sep = "/"), 
        #                      round(wilcox$statistic, 3), round(wilcox$p.value, 3), R2, "Wilcoxon", round(diff, 2))
        #  colnames(data3) <- c("site", "src", "t", "p", "R2", "test", "temp diff")
        #  stats <- rbind(stats, data3)
      }
    }
  }
  row.names(stats) <- NULL
  return(stats)
}

# Function that calculates all of the standard deviations and mean temperatures per site in a long list of data
      ## NB: expects these column labels: site, src, temp (in no particular order)
standev <- function(data){
  data <- data[complete.cases(data),]
  stats <- data.frame()
  for(i in 1:length(levels(data$site))){
    x <- droplevels(subset(data, site == levels(data$site)[i]))
    for(j in 1:length(levels(as.factor((x$src))))){
      y <- droplevels(subset(x, src == levels(as.factor(x$src))[j]))
      #y1 <- y[complete.cases(y),]
      #if(length(y1$site) == 0){
      #  stats <- data.frame(site = y$site[1], src = y$src[1], sd = NA, temp = NA)  
      #} else{} # An idea to allow this function to process time series where entire months are missing
      stdev <- sd(y$temp)
      temp <- mean(y$temp)
      site <- y$site[1]
      src <- y$src[1]
      data2 <- data.frame(site,src,stdev, temp)
      stats <- rbind(stats, data2)
    }
  }
  colnames(stats) <- c("site", "src", "sd", "temp")
  return(stats)
}

# Function that creates random data sets
randomDataSets <- function(years){
  load("prep/insituDailyMultisite_v3.3.RData") # Used for calculating the range of SD's
  range <- standev(insituDailyMultisite_v3.3) # Greyed out during testing because it is slow
  three_sd <- c(min(range$sd), mean(range$sd), max(range$sd))
  six_change <- (seq(0.15, 0.90, length = 6)/3650) # Rather than using changes observed from data, use these
  sd_labels <- paste("SD", as.character(c(round(min(range$sd),2), round(mean(range$sd),2), round(max(range$sd),2))), sep = " = ")
  change_labels <- paste("DT", as.character(seq(0.15, 0.90, length = 6)), sep = " = ")
  resolution <- c("0.01", "0.05", "0.1", "0.5", "1.0")
  dataset <- data.frame()
  ## Populate the data set
  ds1 <- data.frame()
  for(i in 1:3){
    #set.seed(204) # Use if you want consistent "random" results
    data <- rnorm(365*years, 19, three_sd[i]) # The three different standard deviations
    label <- rep(sd_labels[i], 365*years)
    index <-  seq(1:(365*years))
    data1 <- data.frame(index, data, label)
    ds1 <- rbind(ds1, data1)
  }
  ds2 <- data.frame()
  for(i in 1:3){
    data <- droplevels(subset(ds1, label == sd_labels[i]))
    for(j in 1:6){
      data2 <- data
      data2[,2] <- data2[,2]+(data2[,1]*six_change[j])
      data2[,3] <- paste(as.character(data2[,3]), change_labels[j], sep = ", ")
      ds2 <- rbind(ds2, data2)
    } 
  }
  ds3 <- data.frame()
  for(i in 1:18){
    data <- droplevels(subset(ds2, label == levels(as.factor(ds2$label))[i]))
    data1 <- data.frame(data[1], round(data[2], 2), data[3], resolution = resolution[1])
    data2 <- data.frame(data[1], round((data[2]*2),1)/2, data[3], resolution = resolution[2])
    data3 <- data.frame(data[1], round(data[2], 1), data[3], resolution = resolution[3])
    data4 <- data.frame(data[1], round((data[2]*2),0)/2, data[3], resolution = resolution[4])
    data5 <- data.frame(data[1], round(data[2], 0), data[3], resolution = resolution[5])
    ds3 <- rbind(data1, data2, data3, data4, data5)
    dataset <- rbind(dataset, ds3)
  }
  colnames(dataset) <- c("index", "temp", "stats", "resolution")
  dataset$stats <- as.factor(dataset$stats)
  dataset$resolution <- as.character(dataset$resolution)
  row.names(dataset) <- NULL
  return(dataset)
}

# Function that calculates quantile regressions by site and source for daily data more than 10 years long
      ## NB: You may change which quantiles you want, but you have to enter five quantiles, as below
qreg <- function(data, quantiles = c(0.95,0.75,0.5,0.25,0.05)){
  library(quantreg)
  index <- data.frame()
  data <- data[complete.cases(data),]
  for(i in 1:length(levels(data$site))){
    data2 <- subset(data, site == levels(data$site)[i])
    for(j in 1:length(levels(as.factor(data2$src)))){
      data3 <- subset(data2, src == levels(as.factor(data2$src))[j])
      x <- seq(1:length(data3$temp))
      y <- data3$temp
      z <- rq(y~x, tau = quantiles)
      #summary(z)
      coef <- data.frame(round(z$coefficients,2))
      #resid <- round(z$residuals,2)
      #rho  <- z$rho
      #fit <- z$fitted.values
      data4 <- data.frame(data3$site[1], data3$src[1], coef[1,])
      index <- rbind(index, data4)
    }
  }
  row.names(index) <- NULL
  colnames(index) <- c("site", "src", quantiles[1], quantiles[2], quantiles[3], quantiles[4], quantiles[5])
  return(index)
}

# Function that calculates the mean DT, SD, temp, and quantiles for the three coastal regions
coastStats <- function(data){
  sites <- read.csv("setupParams/site_list_v3.3.csv")
  #wc <- droplevels(subset(data, site %in% levels(droplevels(subset(sites, coast == "wc"))$site)))
  wc <- droplevels(subset(sites, coast == "wc"))
  sc <- droplevels(subset(sites, coast == "sc"))
  ec <- droplevels(subset(sites, coast == "ec"))
  delta <- significance(data)
  sd <- standev(data)
  qr <- qreg(data, c(0.95,0.75,0.50,0.25,0.05))
  stats <- data.frame(delta[,c(1,2,5)], sd[,3:4], qr[,3:7])
  #stats$coast <- as.factor(c(rep("wc", length(wc[,1])),
  #                           rep("sc", length(sc[,1])),
  #                           rep("ec", length(ec[,1])))) # Doesn't work due to multiple sources
  means <- data.frame(rbind(wc = sapply(droplevels(subset(stats, site %in% levels(wc$site)))[,3:10], mean),
                            sc = sapply(droplevels(subset(stats, site %in% levels(sc$site)))[,3:10], mean),
                            ec = sapply(droplevels(subset(stats, site %in% levels(ec$site)))[,3:10], mean)))
  means <- data.frame(coast = row.names(means), means)
  row.names(means) <- NULL
  return(means)
}

# Function that compares a time series to the Nino 3.4 Index
nino <- function(data, lag = 4, maMonths = 5){
  index <- data.frame()
  data <- data[complete.cases(data),]
  ninoIndex <- read.csv("setupParams/Nino_3.4_Index.csv")
  ninoIndex$date <- parse_date_time(ninoIndex$date, "y")
  ninoIndex2 <- melt(ninoIndex, id.vars = "date", variable.name = "month", value.name = "temp")
  ninoIndex2 <- ninoIndex2[order(ninoIndex2$date),]
  ninoIndex2$index <- paste(ninoIndex2$date, ninoIndex2$month, sep = "/")
  for(i in 1:length(levels(data$site))){
    data2 <- subset(data, site == levels(data$site)[i])
    for(j in 1:length(levels(as.factor(data2$src)))){
      data3 <- subset(data2, src == levels(as.factor(data2$src))[j])
      data3$month <- month(data3$date, label = TRUE, abbr = TRUE)
      monthly <- data3
      monthly$date <- floor_date(monthly$date, "month")
      monthly <- ddply(monthly, .(site, src, date, month), summarize, temp = mean(temp, na.rm = TRUE))
      clim <- data3
      clim <- ddply(clim, .(site, src, month), summarize, temp = mean(temp, na.rm = TRUE))
      anomalies <- data.frame()
      if(length(monthly$temp) <= 12){
        data4 <- data.frame(monthly[1, 1:2], R2 = NA, length = length(monthly$temp))
        row.names(data4) <- NULL
        index <- rbind(index, data4)
      } else {
        for(k in 1:length(levels(monthly$month))){
          monthly2 <- droplevels(subset(monthly, month == levels(monthly$month)[k]))
          clim2 <- subset(clim, month == levels(clim$month)[k])
          monthly2$temp <- monthly2$temp-clim2$temp
          anomalies <- rbind(anomalies, monthly2)
        }
        anomalies <- anomalies[order(anomalies$date),]
        anomalies$temp <- ma(anomalies$temp, maMonths) # x month running mean ##NB: 5 month running mean cuts off first and last two months of time series
        anomalies <- anomalies[complete.cases(anomalies$temp),]
        anomalies$index <- paste(floor_date(anomalies$date, "year"), anomalies$month, sep = "/")
        ninoIndex3 <- subset(ninoIndex2, index %in% anomalies$index)
        anomalies <- subset(anomalies, index %in% ninoIndex3$index)
        if(length(anomalies$temp >= 12)){ # Requires 12 months of data to calculate R2
          anomalies <- anomalies[(lag+1):length(anomalies$index), ] # Lag time series by chosen amount (default 4 months)
          ninoIndex3 <- ninoIndex3[1:length(anomalies$index), ]
          lmodel <- lm(anomalies$temp ~ ninoIndex3$temp)
          R2 <- round(summary(lmodel)$adj.r.squared,2)
          data4 <- data.frame(anomalies[1,1:2], R2, length = length(anomalies$index))
          index <- rbind(index, data4)
        } else {
          data4 <- data.frame(data3[1, c(1,4)], R2 = NA, length = length(anomalies$index))
          row.names(data4) <- NULL
          index <- rbind(index, data4)
        }
      }
    }
  }
  row.names(index) <- NULL
  return(index)
}

# A quick little function that streamlines the testing of other functions, ss stands for "Select Site"
ss <- function(data, siteName){
  data2 <- droplevels(subset(data, site == siteName))
  return(data2)
}

# Function that calculates moving averages (running means)
ma <- function(x, n = 5){
  filter(x, rep(1/n, n), sides = 2) # This allows for leading data to also have an effect on the mean
}

# Function for fitting a sine curve to monthly means (one mean for each of the 12 months). It simply expects a vector of twelve values, and for now its only purpose is to calculate the "lag" (see below). Currently it is used in "proc/ordinations.R"
seasFit <- function(x) {
  A <- (max(x) - min(x)) / 2 # amplitude
  t <- seq(1:12)
  F <- 2 # this is the annual mean
  lag <- 1 # a factor that indicates how far from the start of the year the warmest seawater temp is found
  m <- 0 # zero slope in linear model
  c <- 16 # the intercept (with zero slope this is the annual mean)
  x_fit1 <- nls(x ~ A * (sin((t + lag) * (pi / 6))) + F, start = list(F = F, lag = lag), data = list(x))
  #x_fit2 <- nls(x ~ (m * t) + c, start = list(c = c)) # a linear function --- not used
  return(coef(x_fit1)[[2]])
}

# Function that calculates the min max and range of a series of a time series
      ## NB: expects these column labels: site, src, temp (in no particular order)
mmr <- function(data){
  index <- data.frame()
  data <- data[complete.cases(data),]
  for(i in 1:length(levels(data$site))){
    data2 <- subset(data, site == levels(data$site)[i])
    for(j in 1:length(levels(as.factor(data2$src)))){
      data3 <- subset(data2, src == levels(as.factor(data2$src))[j])
      x <- min(data3$temp)
      y <- max(data3$temp)
      z <- y-x
      data4 <- data.frame(site = data3$site[1], src = data3$src[1], min = x, max = y, range = z)
      index <- rbind(index, data4)
    }
  }
  row.names(index) <- NULL
  return(index)
}

# Function that extracts Feb and Aug temps from a clim time series, also calculates the overall mean
  # These are then used for plotting
      ## NB: expects these column labels: site, src, lon, lat, temp, date (in no particular order)
        ## Also looks for a distance column, if one does not exist it creates one with the value "0km" in every row
climPlotData <- function(data){
  if("distance" %in% colnames(data) == F){
    data <- data[,c("site", "src", "lon", "lat", "temp", "date")]
    mean_temp <- ddply(data, .(site, src, lon, lat), summarize, temp = mean(temp, na.rm = TRUE))
    mean_temp$date <- rep("mean")
    feb_temp <- droplevels(subset(data, date == "Feb"))
    aug_temp <- droplevels(subset(data, date == "Aug"))
    temps <- rbind(feb_temp, aug_temp, mean_temp)
    temps$distance <- rep("0km")
    temps <- temps[,c(1:4,7,5:6)]
    return(temps)
  } else if ("distance" %in% colnames(data) == T){
    data <- data[,c("site", "src", "lon", "lat", "distance", "temp", "date")]
    mean_temp <- ddply(data, .(site, src, distance, lon, lat), summarize, temp = mean(temp, na.rm = TRUE))
    mean_temp$date <- rep("mean")
    feb_temp <- droplevels(subset(data, date == "Feb"))
    aug_temp <- droplevels(subset(data, date == "Aug"))
    temps <- rbind(feb_temp, aug_temp, mean_temp)
    return(temps)
  }
}