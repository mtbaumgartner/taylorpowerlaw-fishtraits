#' Script to organize fish data from RivFishTIME database
#' It takes raw files from database and exports a list with selected time series
#' based on specified criteria
#' Author: Matheus T. Baumgartner

rm(list = ls())

# Load original raw data
TimeSeries <- read.csv("1873_2_RivFishTIME_TimeseriesTable.csv", header = T)
SurveyTable <- read.csv("1873_2_RivFishTIME_SurveyTable.csv", header = T)

# Select only USA sites
TimeSeries <- subset(TimeSeries, Country == "USA")
plot(Latitude ~ Longitude, TimeSeries)

# Retrieve the 1) Number of species, 2) Number of sampling years, 3) Abundance unit
for(i in 1:nrow(TimeSeries)){
  tab <- subset(SurveyTable, TimeSeriesID == TimeSeries$TimeSeriesID[i])
  
  TimeSeries$NumSpecies[i] <- length(unique(tab$Species))
  TimeSeries$TSlength[i] <- length(unique(tab$Year))
  TimeSeries$UnitAbundance[i] <- unique(tab$UnitAbundance)
  cat(i, "\r")
}

# Select only Count data and electrofishing sampling
TimeSeriesSelected <- subset(TimeSeries, UnitAbundance == "Count")

TimeSeriesSelected <- subset(TimeSeriesSelected, Protocol %in% c("Electrofishing",
                                                                 "Electrofishing_backpack",
                                                                 "Electrofishing_boat"))
plot(Latitude ~ Longitude, TimeSeriesSelected)

# Export as .csv
write.csv(TimeSeriesSelected, file = "TimeSeriesSelected.csv", row.names = F)
