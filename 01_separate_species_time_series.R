#' Script to organize selected fish species data into a single time series
#' for each species to facilitate selection and organization
#' Author: Matheus T. Baumgartner

rm(list = ls())

# Load table with selected time series and original raw data
TimeSeriesSelected <- read.csv("TimeSeriesSelected.csv", header = T, stringsAsFactors = T)
SurveyTable <- SurveyTable <- read.csv("1873_2_RivFishTIME_SurveyTable.csv", header = T)

library(reshape2)

### Separate time series for each species as site (columns) and year (rows)

# Subset selected time series
SurveyTable <- subset(SurveyTable, TimeSeriesID %in% TimeSeriesSelected$TimeSeriesID)

# List of Species Names
SpecNames <- sort(unique(SurveyTable$Species))

TimeSeriesSpecies <- list()
for(i in 1:length(SpecNames)){
  
  df <- subset(SurveyTable, Species == SpecNames[i])
  
  tab <- dcast(df, Year ~ TimeSeriesID, value.var = "Abundance", fun.aggregate = sum)
  rownames(tab) <- tab[,1]
  tab <- data.frame(tab[,-1])
  
  TimeSeriesSpecies[[i]] <- tab
  cat("Time series organized for species", i, "\r")
}
names(TimeSeriesSpecies) <- SpecNames


## Filter data according to criteria
# Remove species with less than 5 subpopulations
TimeSeriesSpecies <- TimeSeriesSpecies[sapply(TimeSeriesSpecies, ncol) >= 5]

# Remove species with less than 15 years of data
TimeSeriesSpecies <- TimeSeriesSpecies[sapply(TimeSeriesSpecies, nrow) >= 15]

# Export object with times series
save(TimeSeriesSpecies, file = "TimeSeriesSpecies.RData")
