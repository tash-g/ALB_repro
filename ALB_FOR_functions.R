
# PROCESS_INTERACTION_ESTIMATES: CALCULATE ESTIMATES FOR EACH INTE --------

process_interaction_estimates <- function(var_name, fixef_data, posterior_samples.df, colony_interaction) {
  # Calculate estimates and confidence intervals
  est <- fixef_data$Estimate[which(rownames(fixef_data) == gsub("b_", "", var_name))] +
    fixef_data$Estimate[which(rownames(fixef_data) == gsub("b_", "", colony_interaction))]
  LCI <- fixef_data$`l-95% CI`[which(rownames(fixef_data) == gsub("b_", "", var_name))] +
    fixef_data$`l-95% CI`[which(rownames(fixef_data) == gsub("b_", "", colony_interaction))]
  UCI <- fixef_data$`u-95% CI`[which(rownames(fixef_data) == gsub("b_", "", var_name))] +
    fixef_data$`u-95% CI`[which(rownames(fixef_data) == gsub("b_", "", colony_interaction))]
  
  # Compute probabilities
  samples <- posterior_samples.df[[var_name]] + posterior_samples.df[[colony_interaction]]
  above_zero <- round(above_zero.func(samples), digits = 2)
  below_zero <- round(below_zero.func(samples), digits = 2)
  
  # Create a row to append
  result <- cbind(matrix(NA, 1, 12), est, LCI, UCI, below_zero, above_zero)
  rownames(result) <- paste("col", var_name, sep = ".")
  
  return(result)
}


# ABOVE AND BELOW ZERO: CALCULATE PROPORTION OF SAMPLES -------------------

above_zero.func <- function (variable) { mean(variable > 0) }
below_zero.func <- function (variable) { mean(variable < 0) }

# SUMMARYSE ---------------------------------------------------------------

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# ANGLE_DIFF: CALCULATE ANGULAR DIFFERENCES -------------------------------
angle_diff <- function(theta1, theta2){
  theta <- abs(theta1 - theta2) %% 360 
  return(ifelse(theta > 180, 360 - theta, theta))
}


# CALCULATE BEARING -------------------------------------------------------

calculate_bearing <- function(lon, lat) {
  require(geosphere)
  if (length(lon) > 1) {
    return(c(bearing(cbind(lon[-length(lon)], lat[-length(lat)]), cbind(lon[-1], lat[-1])), NA))
  } else {
    return(NA)
  }
}

# CALCULATE WIND DIRECTION -----------------------------------------------------
calc_windDir <-function(u,v){
  (270-atan2(v,u)*180/pi)%%360 
}

# LOAD RDATA TO SPECIFIC NAME -----------------------------------------------
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# CLEAN ZEROES FROM FILENAMES ---------------------------------------------

clean_filenames <- function(filenames) {
  filenames <- gsub("(_0+)(\\.csv)$", "\\2", filenames)
  filenames <- gsub(" +0+(\\.csv)$", "\\1", filenames)
  return(filenames)
}

# FIND_MEDIAN_DATE -------------------------------------------------------------

## Find median of dates

find_median_date <- function (dates) {
  
  dates <- as.POSIXlt(dates, format = "%Y-%m-%d")
  dates.julian <- median(dates$yday)
  median.julian <- as.Date(dates.julian, origin=as.Date("1960-01-01"))
  median.julian <- format(as.Date(median.julian), "%m-%d")
  
  return(median.julian)
  
}


# FILE_CLEANER.KRG -------------------------------------------------------------

## Used to match filenames to metadata

file_cleaner.krg <- function(myfile) {
  
  myfile <- gsub("_?0*\\.csv$", "", myfile)
  myfile <- gsub(" ", "", myfile)
  
  return(myfile)
  
}

# SUMMARY SE --------------------------------------------------------------

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# FIND BEHAVIOUR FROM IMMERSION ------------------------------------------------

find_behaviour <- function(immersionValues) {
  
  behaviour <- rep("NA", length = length(immersionValues))
  
  next_immersion <- c(immersionValues[2:length(immersionValues)], NA)
  prev_immersion <- c(NA, immersionValues[1:length(immersionValues) - 1])
  
  behaviour[immersionValues <= max(immersionValues*0.05)] <- "flight/col"
  behaviour[immersionValues >= max(immersionValues*0.95) &
              next_immersion == max(immersionValues) |
              prev_immersion == max(immersionValues)] <- "rest"
  behaviour[behaviour == "NA"] <- "foraging"
  
  return(behaviour)
}

# FIND TRIP DURATION FROM PROCESSED IMMERSION DATA ------------------------

# For troubleshooting
# mydf <- mydf
# rolling.immersion <- mydf$rolling.percentage
# datetime <- mydf$datetime
# immersion <- mydf$immersion
# thresh <- thresh.immersion


tripFinder <- function(datetime, immersion, rolling.immersion, thresh) {
  
  stts <- c()
  ends <- c()
  
  i <- 1
  r <- 1
  
  while (r < length(immersion)) {
    
    if (rolling.immersion[r] > thresh) {  # Point at which cross threshold
      
      pointR <- r   
      
      while (immersion[r] > 0 & r > 1) { # Move backwards to point where raw immersion = 0
        
        r <- r - 1 }  
      
      # Bug in R means POSIXct at midnight excludes time component
      posix.startDate <- format(as.POSIXct(datetime[r], tz="UTC"), format="%Y-%m-%d %H:%M:%S")
      
      stts[i] <- as.character(posix.startDate) #  Label start of trip
      
      r <- pointR + 1     # Advance r by 1 
      
      
      while (rolling.immersion[r] > thresh & r < length(immersion)-3) {   # While immersion > threshold, move forwards to find end of trip
        
        r <- r + 1 }    
      
      r <- r + 1
      
      while (sum(immersion[r:(r+3)]) > 1 & r < length(immersion)-3) {  # Move forward until raw immersion = 0 for 30 mins
        
        r <- r + 1 }
      
      r <- r + 1  # Add 10 min commuting time penalty for return 
      
      # Bug in R means POSIXct at midnight excludes time component
      posix.endDate <- format(as.POSIXct(datetime[r], tz="UTC"), format="%Y-%m-%d %H:%M:%S")
      
      ends[i] <- as.character(posix.endDate) # Set end time 
      
      r <- r + 1
      
      i <- i + 1 
      
    } else { r <- r + 1} } 
  
  trips <- data.frame(starttime = stts, endtime = ends)
  
  return(trips)
  
}


# GCD.HF - Calculate great circle distances using haversine formula ------------

gcd.hf <- function(long1, lat1, long2, lat2) { 
  R <- 6371 # Earth mean radius [km]
  deg2rad <- function(deg) return(deg*pi/180)
  long1 <- deg2rad(long1)
  long2 <- deg2rad(long2)
  lat1 <- deg2rad(lat1)
  lat2 <- deg2rad(lat2)
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  
  coo2 <- function(x){2 * asin(min(1,sqrt(x)))}
  c <- lapply(a, coo2)
  c <- do.call("c", c)
  d = R * c
  return(d) # Distance in km
}

# INTERPOLATE GPS - interpolate GPS to given interval --------------------------

interpolate_gps <- function(lat, lon, datetime, int) {
  
  require(geosphere); require(pracma)
  
  datetime <- datetime + 1
  df <- data.frame(cbind(as.character(datetime), lon, lat))
  colnames(df)[1] <- "datetime"
  df$datetime <- as.POSIXct(df$datetime, format = "%Y-%m-%d %H:%M:%S")
  df[,2:3] <- lapply(df[,2:3], as.numeric)
  
  ## Resample GPS to exact fixes
  df$x1 <-  c(abs(df$lon[1:(nrow(df)-1)]-df$lon[2:(nrow(df))]),1)
  df$x2 <-  c(abs(df$lat[1:(nrow(df)-1)]-df$lat[2:(nrow(df))]),1)
  df$x3 <- c(as.numeric(df$datetime)[1:(length(as.numeric(df$datetime))-1)]- as.numeric(df$datetime)[2:(length(as.numeric(df$datetime)))],-1)
  df <- df[which(df$x1 != 0 & df$x2 != 0 & df$x3 < 0),]
  
  df$distance_next_metres <- c(0,distHaversine(data.frame(df$lon[1:(nrow(df)-1)],
                                                          df$lat[1:(nrow(df)-1)]),
                                               data.frame(df$lon[2:(nrow(df))],
                                                          df$lat[2:(nrow(df))]))) 
  
  df$cumdist <- cumsum(df$distance_next_metres)
  df <- df[order(df$datetime),]
  
  nrt <- nrow(df)
  idt <- seq(from=df$datetime[1], to=df$datetime[nrt], by=int)
  ix = pchip(as.numeric(df$datetime), df$lon, as.numeric(idt))
  iy = pchip(as.numeric(df$datetime), df$lat, as.numeric(idt))
  icumdist = pchip(as.numeric(df$datetime), df$cumdist, as.numeric(idt))
  
  df <- data.frame(datetime = idt,
                   Longitude = ix,
                   Latitude = iy)
  
  df$datetime <- df$datetime-1
  
  return(df)
  
}


## GET_TIME_RANGE : Get time range for timestamping ----------------------------
get_time_range <- function(df) {
  if (df$datetime[1] < df$datetime[nrow(df)]) {
    sttime <- df$datetime[1]
    endtime <- df$datetime[nrow(df)] + df$duration[nrow(df)]
  } else {
    sttime <- df$datetime[1] - df$duration[1]
    endtime <- df$datetime[nrow(df)]
  }
  return(data.frame(sttime = sttime, endtime = endtime))
}



# ## CONVERT HIGHRES - Converts high res GLS to 10-minute bins  ----------------
convert_highres <- function(df){
  
  df$datetime <- time_sorter(df$datetime)
  
  expansion <- rep(df$state, df$duration)
  
  # Create 10 minute chunks and sum wet events within them 
  chunks <- split(expansion, ceiling(seq_along(expansion)/600))
  wet_freq <- plyr::ldply(chunks, function(c) sum(c=="wet"))
  
  time_range <- get_time_range(df)
  
  ## Make the dataframe 
  outlength <- ceiling(sum(df$duration)/600)
  df2 <- data.frame(datetime = seq.POSIXt(time_range$sttime, by = "600 secs", length.out = outlength),
                                 immersion = wet_freq$V1/3)
  
  
  return(df2)
  
}



# ## TASH'S HACKY DATE TIME FORMATTER -------------------------------------

# testdates <- c("23/05/2024 15:56",
#                "2022-04-11 15:11:11",
#                "12/06/06 16:12",
#                "2021-03-17 17:03",
#                "04-Jan-16 16:17",
#                "19-Aug-23 17:14:24",
#                "09-02-2012 18:45:12")
# 
# time_sorter(testdates)

time_sorter <- function(x){
  
  date1 <- sapply(strsplit(x, " "), "[[", 1)
  time1 <- sapply(strsplit(x, " "), "[[", 2)
  seconds <- sapply(strsplit(time1, ":"), function(x) length(x))
  
  date_format <- character(length(date1))
  
  for (i in seq_along(date1)) {
    
    # Find where four digit year is (if present)
    split_date <- unlist(strsplit(date1[i], "[-/]"))
    yr_pos <- grep("^\\d{4}$", split_date)
    yr_pos <- ifelse(length(yr_pos) == 0, 0, yr_pos)
    
    if (grepl("/", date1[i]) & nchar(date1[i]) == 8) {
      date_format[i] <- as.character(as.Date(date1[i], tryFormats = c("%Y-%m-%d", "%d/%m/%y", "%d-%b-%y", "%d-%b-%Y")))
    } else if (yr_pos == 3) {
      date_format[i] <- as.character(as.Date(date1[i], tryFormats = c("%d-%m-%Y", "%d/%m/%Y", "%d-%b-%y", "%d-%b-%Y", "%d/%m/%Y")))
    } else  {
      date_format[i] <- as.character(as.Date(date1[i], tryFormats = c("%Y-%m-%d", "%d-%m-%Y", "%d/%m/%Y", "%d-%b-%y", "%d-%b-%Y", "%d/%m/%Y")))
    }
  }
  
  time_format <- ifelse(seconds == 2, paste0(time1, ":00"), time1)
  
  datetime <- paste0(date_format, " ", time_format)
  datetime <- as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  
  return(datetime)
}



# 
# time_sorter <- function(x){
#   date1 <- sapply(strsplit(x, " "), "[[", 1)[1]
#   time1 <- sapply(strsplit(x, " "), "[[", 2)[1]
#   seconds <- length(strsplit(time1, ":")[[1]])
#   
#   if(grepl("[A-Z]", date1) & seconds == 2) {
#     return(as.POSIXct(x, format = "%d-%b-%y %H:%M", tz = "UTC")) 
#     
#   } else {
#     
#     if(grepl("[A-Z]", date1) & seconds == 3) {
#       return(as.POSIXct(x, format = "%d-%b-%y %H:%M:%S", tz = "UTC")) 
#       
#     } else {
#       
#       if(grepl("-", date1) == TRUE) { 
#         return(as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) 
#         
#       } else {
#         
#         if(nchar(date1) == 10 & seconds == 3 & nchar(sapply(strsplit(date1, "/"), "[[", 1)) == 4){
#           return(as.POSIXct(x, format = "%Y/%m/%d %H:%M:%S", tz = "UTC")) 
#         } else {
#           
#           if(nchar(date1) == 10 & seconds == 3 & nchar(sapply(strsplit(date1, "/"), "[[", 1)) == 2){
#             return(as.POSIXct(x, format = "%d/%m/%Y %H:%M:%S", tz = "UTC"))
#           } else {
#             
#             if(nchar(date1) == 8 & seconds == 3){
#               return(as.POSIXct(x, format = "%d/%m/%y %H:%M:%S", tz = "UTC"))
#             } else {
#               
#               if(nchar(date1) == 10 & seconds == 2){
#                 return(as.POSIXct(x, format = "%d/%m/%Y %H:%M", tz = "UTC"))
#               } else {
#                 
#                 if(nchar(date1) == 8 & seconds == 2){
#                   return(as.POSIXct(x, format = "%d/%m/%y %H:%M", tz = "UTC"))
#                 }
#               }
#             }
#           }
#         }
#         
#       }
#     }
#   }
# }
