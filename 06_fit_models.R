#' Script to fit models and produce figures for the article
#' Author: Matheus T. Baumgartner

rm(list = ls())

library(ape)
library(nlme)
library(phytools)
library(phylolm)
library(rr2)
library(GGally)
library(ggplot2)
library(tidyr)
library(readxl)

# Load time series data
load("TimeSeriesSpecies.RData")

# Load TPL parameters
TPL <- read.csv("TPL_outputs.csv", sep = ",")

# Load PCA for environmental covariates
env <- read.csv("covariates.csv", sep = ",")

# Load PCA for traits
traits <- read.csv("portfolio.csv", sep = ",")

# Filter time series, TPL parameters and environmental data to match trait data
TimeSeriesSpecies <- TimeSeriesSpecies[names(TimeSeriesSpecies) %in% traits$Species]
TPL <- subset(TPL, Species %in% traits$Species)
env <- env[env$SiteID %in% unique(unlist(lapply(TimeSeriesSpecies, colnames))),]

# Get environmental centroids for each species
df <- data.frame(traits, PC1_env = NA, PC2_env = NA)

for(i in 1:nrow(df)){
  dat <- subset(env, SiteID %in% colnames(TimeSeriesSpecies[[df$Species[i]]]))
  df[i,22:23] <- colMeans(dat[,6:7])
}

# Load phylogeny
phylo <- read.nexus("phylogeny.nexus")

# Calculate phylogenetic distances
phylo_dist <- as.dist(cophenetic.phylo(phylo))
plot(hclust(phylo_dist))

# Match names in df with those in phylogenetic tree
df$Species[df$Species == 'Ericymba buccata'] <- 'Notropis buccatus'
df$Species <- gsub('[ ]', '_', df$Species)

# Check phylogenetic signal in alhpa and beta using Pagel's lambda
alphas <- df$GLS_a
names(alphas) <- df$Species
phylosig(phylo, alphas, method = "lambda", test = T, nsim = 999) # no signal

betas <- df$GLS_mu
names(betas) <- df$Species
phylosig(phylo, betas, method = "lambda", test = T, nsim = 999) # signal


# ---- Fit linear models and PGLS ----
# Check distribution of TPL parameters
par(mfrow = c(1, 2))
hist(df$GLS_a)
hist(df$GLS_mu)
layout(1)

# Recalculate phylogenetic distances
phylo_dist <- as.dist(cophenetic.phylo(phylo))


## Fit GLS to select the best phylogenetic correlation structure based on AIC
gls_b <- gls(GLS_mu ~ (PC1_trophic + PC2_trophic +
                         PC1_life + PC2_life +
                         PC1_habitat + PC2_habitat +
                         PC1_env + PC2_env + Nsites + Nyears) ^ 2,
             data = df,
             correlation = corPagel(1, phylo, form = ~Species), method = "REML")
             #correlation = corBrownian(1, phylo, form = ~Species), method = "REML")
             #correlation = corBlomberg(1, phylo, form = ~Species, fixed = T), method = "REML")
             #correlation = corMartins(1, phylo, form = ~Species), method = "REML")
             #correlation = corGrafen(1, phylo, form = ~Species, fixed = T), method = "REML")
plot(gls_b)
AIC(gls_b) # used to compare correlation structures
summary(gls_b)
intervals(gls_b)

qqnorm(resid(gls_b))
qqline(resid(gls_b))

# ---- Backward model selection using stepwise with phylolm package ----
# Remove outliers from each model while checking for normality of residuals and R2

# alpha
rownames(df) <- df$Species
bestmod_a <- phylostep(GLS_a ~ (PC1_trophic + PC2_trophic +
                                  PC1_life + PC2_life +
                                  PC1_habitat + PC2_habitat +
                                  PC1_env + PC2_env + Nsites + Nyears) ^ 2,
                       #data = df, phy = phylo,
                       data = subset(df, !Species %in% c("Fundulus_catenatus",
                                                         "Campostoma_pauciradii",
                                                         "Hypentelium_etowanum",
                                                         "Luxilus_chrysocephalus")),
                       phy = drop.tip(phylo, c("Fundulus_catenatus",
                                               "Campostoma_pauciradii",
                                               "Hypentelium_etowanum",
                                               "Luxilus_chrysocephalus")),
                       model = "lambda", direction = "backward")
summary(bestmod_a)
plot(bestmod_a)
abline(0, 1)

qqnorm(resid(bestmod_a))
qqline(resid(bestmod_a))

hist(resid(bestmod_a))
shapiro.test(resid(bestmod_a))
ks.test(resid(bestmod_a), "pnorm")


# beta
bestmod_b <- phylostep(GLS_mu ~ (PC1_trophic + PC2_trophic +
                                   PC1_life + PC2_life +
                                   PC1_habitat + PC2_habitat +
                                   PC1_env + PC2_env + Nsites + Nyears) ^ 2,
                       #data = df, phy = phylo,
                       data = subset(df, !Species %in% c("Lepisosteus_oculatus",
                                                         "Lepisosteus_platostomus")),
                       phy = drop.tip(phylo, c("Lepisosteus_oculatus",
                                               "Lepisosteus_platostomus")),
                       model = "lambda", direction = "backward")
summary(bestmod_b)
plot(bestmod_b)
abline(0, 1)

qqnorm(resid(bestmod_b))
qqline(resid(bestmod_b))

hist(resid(bestmod_b))
shapiro.test(resid(bestmod_b))
ks.test(resid(bestmod_b), "pnorm")


# ---- Fit final models using selected model structure (PGLS) ----

# Remove species that were outliers in either alpha or beta

# alpha
fit_a <- gls(update(formula(gls_b), formula(bestmod_a)),
             data = subset(df, !Species %in% c("Fundulus_catenatus",
                                               "Campostoma_pauciradii",
                                               "Hypentelium_etowanum",
                                               "Lepisosteus_oculatus",
                                               "Lepisosteus_platostomus",
                                               "Luxilus_chrysocephalus")),
             correlation = corPagel(1, drop.tip(phylo, c("Fundulus_catenatus",
                                                         "Campostoma_pauciradii",
                                                         "Hypentelium_etowanum",
                                                         "Lepisosteus_oculatus",
                                                         "Lepisosteus_platostomus",
                                                         "Luxilus_chrysocephalus")),
                                    form = ~Species))
summary(fit_a)
intervals(fit_a)
plot(fit_a)
rcompanion::nagelkerke(fit_a) # pseudo R-squared
car::Anova(fit_a, type = "III")

qqnorm(resid(fit_a))
qqline(resid(fit_a))

hist(resid(fit_a))
shapiro.test(resid(fit_a))


# beta
fit_b <- gls(update(formula(gls_b), formula(bestmod_b)),
             data = subset(df, !Species %in% c("Fundulus_catenatus",
                                               "Campostoma_pauciradii",
                                               "Hypentelium_etowanum",
                                               "Lepisosteus_oculatus",
                                               "Lepisosteus_platostomus",
                                               "Luxilus_chrysocephalus")),
             correlation = corPagel(1, drop.tip(phylo, c("Fundulus_catenatus",
                                                         "Campostoma_pauciradii",
                                                         "Hypentelium_etowanum",
                                                         "Lepisosteus_oculatus",
                                                         "Lepisosteus_platostomus",
                                                         "Luxilus_chrysocephalus")),
                                    form = ~Species))
summary(fit_b)
intervals(fit_b)
plot(fit_b)
rcompanion::nagelkerke(fit_b) # pseudo R-squared
car::Anova(fit_b, type = "III")

qqnorm(resid(fit_b))
qqline(resid(fit_b))

hist(resid(fit_b))
shapiro.test(resid(fit_b))


# Calculate partial r-squared values based on Ives (2018) - Sys. Biol.
R2_lik(fit_a);R2_pred(fit_a)
R2_lik(fit_b);R2_pred(fit_b)


# Remove species that were outliers in either alpha or beta
df <- subset(df, !Species %in% c("Fundulus_catenatus",
                                 "Campostoma_pauciradii",
                                 "Hypentelium_etowanum",
                                 "Lepisosteus_oculatus",
                                 "Lepisosteus_platostomus",
                                 "Luxilus_chrysocephalus"))


# Histograms of GLS coefficients for S and J
hist(df$GLS_S, xlab = expression(beta[S]), main = "", freq = F, breaks = 10)
hist(df$GLS_J, xlab = expression(beta[J]), main = "", freq = F, breaks = 10)
plot(GLS_S ~ GLS_J, df,
     xlab = expression(beta[J]),
     ylab = expression(beta[S]))

# Test if gls was better than lm for estimating TPL parameters based on AIC
df$deltaAIC[is.infinite(df$deltaAIC)] <- 0
hist(df$deltaAIC)
mean(df$deltaAIC)
t.test(df$deltaAIC, mu = 0)
wilcox.test(df$deltaAIC, mu = 0, alternative = "less")


# Test for differences in the values of beta between OLS and GLS
plot(OLS_mu ~ GLS_mu, df)
abline(a = 0, b = 1)
t.test(df$OLS_mu,df$GLS_mu, paired = T)


# Plot alpha vs beta
plot(GLS_mu ~ GLS_a, df,
     xlab = expression(alpha),
     ylab = expression(beta[mu]))
cor.test(df$GLS_mu, df$GLS_a)

# Save
save.image("fit_models.RData")


# Figures
# For SM
ggpairs(df[,c(9,10,2,3)])


# For MS: boxplot of TPL parameters (Figure 3)
df_TPL <- df[,(9:12)]
df_TPL <- data.frame(pivot_longer(df_TPL, cols = 1:4, names_to = "pred", values_to = "values"))
df_TPL$pred <- factor(df_TPL$pred, levels = rev(c("GLS_a", "GLS_mu", "GLS_S", "GLS_J")),
                      ordered = T)

tiff("Fig_boxplot.tiff",
     width = 12, height = 12, units = "cm",
     res = 300, compression = "lzw")
ggplot(df_TPL, aes(x = pred, y = values)) +
  #ylim(-5.5, 5.5) +
  xlab(expression("Predictor of"~ln(sigma^2))) +
  ylab("Parameter value") +
  scale_x_discrete(labels = c(expression(beta[J]),
                              expression(beta[S]),
                              expression(beta[mu]),
                              expression(alpha))) +
  geom_abline(slope = 0, intercept = c(0, 1, 2),
              linetype = "dashed", color = c("black", "red", "red")) +
  geom_boxplot(fill = c("yellow", "yellow", "red", "red")) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(color = "black", size = 14),
        axis.text.x = element_text(color = "black", size = 12),
        axis.title = element_text(size = 14))
dev.off()



# Retrieve statistics for manuscript
MS <- df
MS$Species <- gsub('[ ]', '_', MS$Species)

mean(MS$GLS_a);sd(MS$GLS_a)
mean(MS$GLS_mu);sd(MS$GLS_mu)

plot(MS$GLS_a, MS$GLS_mu)


## Load fish traits and environmental variables to correlate with TPL parameters
fishtraits <- data.frame(read_xls("fishtraits.xls", sheet = "traits_binary2"))
fishtraits$Names[fishtraits$Names == 'Ericymba buccata'] <- 'Notropis buccatus'
fishtraits$Names <- gsub('[ ]', '_', fishtraits$Names)
fishtraits <- subset(fishtraits, Names %in% MS$Species)

# alpha
tiff("Fig_scatterplots.tiff",
     width = 23, height = 28, units = "cm",
     res = 300, compression = "lzw")
par(mfrow = c(4,3))
plot(MS$GLS_a ~ fishtraits$maxtl,
     xlab = "maxtl", ylab = expression(alpha))

plot(MS$GLS_a ~ fishtraits$matuage,
     xlab = "matuage", ylab = expression(alpha))

plot(MS$GLS_a ~ fishtraits$longevity,
     xlab = "longevity", ylab = expression(alpha))

plot(MS$GLS_mu ~ fishtraits$maxtl,
     xlab = "maxtl", ylab = expression(beta[mu]))

plot(MS$GLS_mu ~ fishtraits$matuage,
     xlab = "matuage", ylab = expression(beta[mu]))

plot(MS$GLS_mu ~ fishtraits$longevity,
     xlab = "longevity", ylab = expression(beta[mu]))

plot(MS$GLS_a ~ df$PC1_env,
     xlab = expression(PC1[env]), ylab = expression(alpha))

plot(MS$GLS_mu ~ df$PC2_env,
     xlab = expression(PC2[env]), ylab = expression(beta[mu]))

plot(MS$GLS_a ~ log(MS$Nsites),
     xlab = expression(log(N[sites])), ylab = expression(alpha))

plot(MS$GLS_a ~ log(MS$Nyears),
     xlab = expression(log(N[years])), ylab = expression(alpha))

plot(MS$GLS_mu ~ log(MS$Nsites),
     xlab = expression(log(N[sites])), ylab = expression(beta[mu]))

plot(MS$GLS_mu ~ log(MS$Nyears),
     xlab = expression(log(N[years])), ylab = expression(beta[mu]))
dev.off()


# Check phylogenetic signal in traits using Pagel's lambda
phylo$tip.label <- gsub('[ ]', '_', phylo$tip.label)

for(i in 6:ncol(fishtraits)){
  tr <- fishtraits[,i]
  names(tr) <- fishtraits$Names
  cat(names(fishtraits)[i], "\n")
  print(phylosig(drop.tip(phylo, phylo$tip.label[!phylo$tip.label %in% fishtraits$Names]),
           tr, method = "lambda", test = T, nsim = 999))
}

