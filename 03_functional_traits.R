#' Script to build niche dimension from functional traits
#' Author: Matheus T. Baumgartner

rm(list = ls())

library(readxl)
library(vegan)

# Load TPL outputs
TPL_outputs <- read.csv("TPL_outputs.csv", sep = ",")

# Remove species whose GLS did not converge
TPL_outputs <- na.omit(TPL_outputs)

# Load trait data
fishtraits <- data.frame(read_xls("fishtraits.xls", sheet = "traits_binary2"))

# Assign species as row names
rownames(fishtraits) <- fishtraits[,1]
fishtraits <- fishtraits[,-1]

# Remove species without full trait data
fishtraits <- na.omit(fishtraits)

# Remove species whose GLS did not converge
fishtraits <- subset(fishtraits, Species %in% TPL_outputs$Species)

# Match data between TPL and traits
TPL_outputs <- subset(TPL_outputs, Species %in% fishtraits$Species)

# Split data according to niche dimensions
trophic <- fishtraits[,5:13]
life    <- fishtraits[,14:20]
habitat <- fishtraits[,21:35]


# --- PCoA with Gower dissimilarities for each niche dimension ----

## Trophic
pc_trophic <- cmdscale(vegdist(trophic, method = "gower"), eig = T)
plot(pc_trophic$points)
pc_trophic$eig[1:2]/sum(pc_trophic$eig)*100
cor(trophic, pc_trophic$points[,1:2]) # Pearson is equivalent to point-bisserial coefficient

## Life History
pc_life <- cmdscale(vegdist(life[,-7], method = "gower"), eig = T)
plot(pc_life$points)
pc_life$eig[1:2]/sum(pc_life$eig)*100
cor(life, pc_life$points[,1:2])

## Habitat Use
pc_habitat <- cmdscale(vegdist(habitat, method = "gower"), eig = T)
plot(pc_habitat$points)
pc_habitat$eig[1:2]/sum(pc_habitat$eig)*100
cor(habitat, pc_habitat$points[,1:2])

## Substrate
pc_substrate <- cmdscale(vegdist(habitat[,1:11], method = "gower"), eig = T)
plot(pc_substrate$points)
pc_substrate$eig[1:2]/sum(pc_substrate$eig)*100
cor(habitat[,1:11], pc_substrate$points[,1:2])

# Use PCs from substrate in PCoA for habitat
pc_habsub <- cmdscale(vegdist(data.frame(pc_substrate$points[,1:2], habitat[,12:15]), method = "can"), eig = T)
plot(pc_habsub$points)
pc_habsub$eig[1:2]/sum(pc_habsub$eig)*100
cor(data.frame(pc_substrate$points[,1:2], habitat[,12:15]), pc_habsub$points[,1:2])


# Combine results
df_portfolio <- data.frame(TPL_outputs,
                           pc_trophic$points[,1:2],
                           pc_life$points[,1:2],
                           pc_habitat$points[,1:2],
                           pc_habsub$points[,1:2])
names(df_portfolio)[14:21] <- c("PC1_trophic", "PC2_trophic",
                                "PC1_life", "PC2_life",
                                "PC1_habitat", "PC2_habitat",
                                "PC1_habsub", "PC2_habsub")

# Export as .csv
write.csv(df_portfolio, file = "portfolio.csv", row.names = F)

# ---- Figures ----
library(ggplot2)

# Remove outliers (identified later)
df_portfolio <- subset(df_portfolio, !Species %in% c("Fundulus_catenatus",
                                                     "Campostoma_pauciradii",
                                                     "Hypentelium_etowanum",
                                                     "Lepisosteus_oculatus",
                                                     "Lepisosteus_platostomus",
                                                     "Luxilus_chrysocephalus"))

## PCoA Trophic
arrows_trophic <- data.frame(cor(trophic, pc_trophic$points[,1:2], method = "pearson"))
names(arrows_trophic) <- c("PCoA1", "PCoA2")

tiff("Fig_PCoA_trophic.tiff",
     width = 16, height = 14, units = "cm",
     res = 300, compression = "lzw")
ggplot(df_portfolio,
       aes(x = PC1_trophic, y = PC2_trophic)) +
  labs(title = "(a) Trophic") +
  xlab("PCoA1 (55.95%)") +
  ylab("PCoA2 (33.81%)") +
  xlim(-0.5, 0.5) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(size = 3, aes(colour = GLS_mu)) +
  scale_colour_gradient2(low = 'white', high = 'darkgreen') +
  #scale_colour_gradient2(low = 'blue', mid = 'white', high = 'darkgreen') +
  geom_segment(data = arrows_trophic,
               aes(x = 0, y = 0, xend = PCoA1/2, yend = PCoA2/2),
               arrow = arrow(length = unit(4, "mm")),
               color = "red") +
  geom_text(data = arrows_trophic,
            aes(x = PCoA1/1.8, y = PCoA2/1.8, label = rownames(arrows_trophic)),
            size = 5, color = "red") +
  theme_minimal() +
  theme(legend.position = "none",
        title = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))
dev.off()


## PCoA Life History
arrows_life <- data.frame(cor(life, pc_life$points[,1:2], method = "pearson"))
names(arrows_life) <- c("PCoA1", "PCoA2")

tiff("Fig_PCoA_life.tiff",
     width = 16, height = 14, units = "cm",
     res = 300, compression = "lzw")
ggplot(df_portfolio,
       aes(x = PC1_life, y = PC2_life)) +
  labs(title = "(b) Life-history") +
  xlab("PCoA1 (63.03%)") +
  ylab("PCoA2 (17.40%)") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(size = 3, aes(colour = GLS_mu)) +
  scale_colour_gradient2(low = 'white', high = 'darkgreen') +
  guides(colour = guide_colourbar(title = expression(beta[mu]))) +
  geom_segment(data = arrows_life,
               aes(x = 0, y = 0, xend = PCoA1/3, yend = PCoA2/3),
               arrow = arrow(length = unit(4, "mm")),
               color = "red") +
  geom_text(data = arrows_life,
            aes(x = PCoA1/2.8, y = PCoA2/2.8, label = rownames(arrows_life)),
            size = 5, color = "red") +
  theme_minimal() +
  theme(legend.position = c(0.10, 0.20),
        legend.background = element_rect(),
        legend.box.background = element_blank(),
        title = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))
dev.off()


## PCoA Habitat Use
arrows_habitat <- data.frame(cor(habitat, pc_habitat$points[,1:2], method = "pearson"))
names(arrows_habitat) <- c("PCoA1", "PCoA2")

tiff("Fig_PCoA_habitat.tiff",
     width = 16, height = 14, units = "cm",
     res = 300, compression = "lzw")
ggplot(df_portfolio,
       aes(x = PC1_habitat, y = PC2_habitat)) +
  labs(title = "(c) Habitat use") +
  xlab("PCoA1 (42.44%)") +
  ylab("PCoA2 (31.18%)") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(size = 3, aes(colour = GLS_mu)) +
  scale_colour_gradient2(low = 'white', high = 'darkgreen') +
  guides(colour = guide_colourbar(title = expression(beta))) +
  geom_segment(data = arrows_habitat,
               aes(x = 0, y = 0, xend = PCoA1/3, yend = PCoA2/3),
               arrow = arrow(length = unit(4, "mm")),
               color = "red") +
  geom_text(data = arrows_habitat,
            aes(x = PCoA1/2.8, y = PCoA2/2.8, label = rownames(arrows_habitat)),
            size = 5, color = "red") +
  theme_minimal() +
  theme(legend.position = "none",
        title = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))
dev.off()


# Joint figure
library(ggpubr)
tiff("Fig_PCoAs.tiff",
     width = 16, height = 40, units = "cm",
     res = 300, compression = "lzw")
ggarrange(
  ggplot(df_portfolio,
         aes(x = PC1_trophic, y = PC2_trophic)) +
    labs(title = "(a) Trophic") +
    xlab("PCoA1 (55.95%)") +
    ylab("PCoA2 (33.81%)") +
    xlim(-0.5, 0.5) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point(size = 3, aes(colour = GLS_mu)) +
    scale_colour_gradient2(low = 'white', high = 'darkgreen') +
    #scale_colour_gradient2(low = 'blue', mid = 'white', high = 'darkgreen') +
    geom_segment(data = arrows_trophic,
                 aes(x = 0, y = 0, xend = PCoA1/2, yend = PCoA2/2),
                 arrow = arrow(length = unit(4, "mm")),
                 color = "red") +
    geom_text(data = arrows_trophic,
              aes(x = PCoA1/1.8, y = PCoA2/1.8, label = rownames(arrows_trophic)),
              size = 5, color = "red") +
    theme_minimal() +
    theme(legend.position = "none",
          title = element_text(size = 15, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12, color = "black")),
  
  ggplot(df_portfolio,
         aes(x = PC1_life, y = PC2_life)) +
    labs(title = "(b) Life-history") +
    xlab("PCoA1 (63.03%)") +
    ylab("PCoA2 (17.40%)") +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point(size = 3, aes(colour = GLS_mu)) +
    scale_colour_gradient2(low = 'white', high = 'darkgreen') +
    guides(colour = guide_colourbar(title = expression(beta[mu]))) +
    geom_segment(data = arrows_life,
                 aes(x = 0, y = 0, xend = PCoA1/3, yend = PCoA2/3),
                 arrow = arrow(length = unit(4, "mm")),
                 color = "red") +
    geom_text(data = arrows_life,
              aes(x = PCoA1/2.8, y = PCoA2/2.8, label = rownames(arrows_life)),
              size = 5, color = "red") +
    theme_minimal() +
    theme(legend.position = c(0.10, 0.20),
          legend.background = element_rect(),
          legend.box.background = element_blank(),
          title = element_text(size = 15, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12, color = "black")),
  
  ggplot(df_portfolio,
         aes(x = PC1_habitat, y = PC2_habitat)) +
    labs(title = "(c) Habitat use") +
    xlab("PCoA1 (42.44%)") +
    ylab("PCoA2 (31.18%)") +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point(size = 3, aes(colour = GLS_mu)) +
    scale_colour_gradient2(low = 'white', high = 'darkgreen') +
    guides(colour = guide_colourbar(title = expression(beta))) +
    geom_segment(data = arrows_habitat,
                 aes(x = 0, y = 0, xend = PCoA1/3, yend = PCoA2/3),
                 arrow = arrow(length = unit(4, "mm")),
                 color = "red") +
    geom_text(data = arrows_habitat,
              aes(x = PCoA1/2.8, y = PCoA2/2.8, label = rownames(arrows_habitat)),
              size = 5, color = "red") +
    theme_minimal() +
    theme(legend.position = "none",
          title = element_text(size = 15, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12, color = "black")),
  ncol = 1
)
dev.off()
