#' Script to obtain the environmental covariates and build PCA to use in
#' regressions of TPL parameters
#' Author: Matheus T. Baumgartner

############################################################################
#### Run the code below after retrieving all environmental data in QGIS ####

rm(list = ls())

library(readxl)
library(corrplot)
library(vegan)
library(ggplot2)
library(factoextra)

# Load environmental data
env <- data.frame(read_xls("tpl_env_data.xls", sheet = "tpl_env_data"))
rownames(env) <- env[,1] # pass row names
env <- env[,-1]
env$J <- as.numeric(env$J) # declare as numeric and replace NAs with 0
env$J[is.na(env$J)] <- 0

# Visualize histograms for each environmental variable
for(i in 1:ncol(env)){ hist(env[,i], main = names(env)[i]) }

# Log-transform some variables
# ELEV, snw_pc_c11, slp_dg_uav, sgr_dk_rav, inu_pc_cmx, DIST_UP_KM, UPLAND_SKM, dis_m3_pyr, HFP2009
env[,c(5,13,17,19,20,21,22,23,25)] <- log10(env[,c(5,13,17,19,20,21,22,23,25)] + 1)

# Remove unused variables
names(env)
env_sel <- env[,-c(3,4,5,6,7,8,10,11,18,21)]
names(env_sel)

# Rename variables
names(env_sel) <- c("lat", "long", "rain", "raincv", "snw", "tmp",
                    "aet", "cmi", "slp", "sgr", "inu", "dup",
                    "dis", "sord", "hfp")

# Plot to show collinearity
corrplot(cor(env_sel), order = "hclust", addrect = 3)

# PCA for multicollinearity
pca_env <- prcomp(env_sel[,-c(1,2)], scale = T) # remove lat and long
summary(pca_env)
biplot(pca_env)
pca_env$rotation[,1:3] # loadings

eigenvals(pca_env)
bstick(pca_env) # Selected: 2 axes

# Build data frame
covars <- data.frame(scale(env[,1:4]), scale(pca_env$x[,1:2]))
names(covars)[5:6] <- c("PC1_env","PC2_env")
cor(covars)

# Export as .csv
covars <- data.frame(SiteID = rownames(covars),
                     covars)

write.csv(covars, file = "covariates.csv", row.names = F)

# Graph
tiff("Fig_PCA_env.tiff",
     width = 16, height = 14, units = "cm",
     res = 300, compression = "lzw")
fviz_pca_biplot(pca_env,
                geom = c("point")) +
  labs(title = "PCA - Environmental") +
  theme_minimal()
dev.off()
