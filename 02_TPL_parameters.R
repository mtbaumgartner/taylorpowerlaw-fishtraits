#' Script to estimate Taylor's power law parameters (alpha and beta_mu)
#' using the number of species and evenness as covariates, with spatial structure
#' Author: Matheus T. Baumgartner

rm(list = ls())

library(reshape2)
library(vegan)
library(nlme)
library(MuMIn)

# Load time series of species
load("TimeSeriesSpecies.RData")

# Load original raw data
TimeSeries <- read.csv("1873_2_RivFishTIME_TimeseriesTable.csv", header = T)

# Extract point coordinates (to use here and later in QGIS)
TimeSeriesID <- lapply(TimeSeriesSpecies, colnames)
TimeSeriesID <- unique(unlist(TimeSeriesID))

PointCoord <- data.frame(SiteID = TimeSeriesID)

for(i in 1:nrow(PointCoord)){
  
  st <- subset(TimeSeries, TimeSeriesID == PointCoord$SiteID[i])
  
  PointCoord$Lat[i] <- st$Latitude
  PointCoord$Long[i] <- st$Longitude
  
  cat(i, "\r")
}

# Export point coordinates
write.csv(PointCoord, file = "PointCoordinates.csv", row.names = F)

# --- Organize a single main dataset
# Load survey data
SurveyTable <- read.csv("1873_2_RivFishTIME_SurveyTable.csv", header = T)

# Subset selected time series
SurveyTable <- subset(SurveyTable, TimeSeriesID %in% PointCoord$SiteID)

# Organize a single dataset
TS <- dcast(SurveyTable, TimeSeriesID ~ Species, value.var = "Abundance", fun.aggregate = sum)
rownames(TS) <- TS[,1]
TS <- TS[,-1]
TS <- TS[PointCoord$SiteID,] # Reorder based on Point Coordinates

# --- Calculate species richness (S) and Pielou's evenness (J)
sites <- PointCoord
sites$S <- specnumber(TS)
sites$J <- diversity(TS, index = "shannon") / log(sites$S)
sites$J[is.na(sites$J)] <- 0

# --- Calculate TPL parameters

TPL_outputs <- data.frame(Species = names(TimeSeriesSpecies))

for(i in 1:nrow(TPL_outputs)){
  
  dat <- TimeSeriesSpecies[[i]]
  
  TPL_outputs[i,2] <- ncol(dat) # number of sites
  TPL_outputs[i,3] <- nrow(dat) # number of years
  
  df <- data.frame(mu = log(apply(dat, 2, mean)),
                   sigma2 = log(apply(dat, 2, var)),
                   sites[sites$SiteID %in% colnames(dat),-1])
  
  # Add some very small noise to coordinates to avoid zero distances
  df$Lat <- df$Lat + runif(nrow(df), -1e5, 1e5)
  df$Long <- df$Long + runif(nrow(df), -1e5, 1e5)
  
  # OLS
  fit_lm <- lm(sigma2 ~ mu + S + J, data = df)
  
  TPL_outputs[i,4:7] <- fit_lm$coefficients # lm coefficients
  TPL_outputs[i,8] <- summary(fit_lm)$r.squared # lm R2
  
  # GLS
  fit_gls <- tryCatch({
    gls(sigma2 ~ mu + S + J, data = df,
        correlation = corRatio(form = ~ Lat + Long))
  }, error = function(x) NULL)
  
  if(is.null(fit_gls)){
    next
  }else{
    TPL_outputs[i,9:12] <- fit_gls$coefficients # gls coefficients
    
    mod_aic <- MuMIn::AICc(fit_lm, fit_gls)
    TPL_outputs[i,13] <- mod_aic[1,2] - mod_aic[2,2] # AIC (lm - gls)
  }
  
  cat("Models fitted for species", i, "\r")
}
names(TPL_outputs) <- c("Species", "Nsites", "Nyears",
                        paste0("OLS_", names(coef(fit_lm))), "OLS_R2",
                        paste0("GLS_", names(coef(fit_gls))),
                        "deltaAIC")
names(TPL_outputs)[c(4,9)] <- c("OLS_a", "GLS_a")

# Export results as .csv
write.csv(TPL_outputs, file = "TPL_outputs.csv", row.names = F)
