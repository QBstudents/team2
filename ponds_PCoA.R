###### Inputting Data
# relative abundance matrix
Ponds.rel <- as.matrix(read.csv("/Users/ericanadolski/GitHub/TeamProject/Ponds.rel.csv", row.names=1))
#sites 54-58 are missing some environmental data

# subsetting total (DNA) and active (cDNA) libraries for each site

keep_DNA <- c("BC001-DNA","BC002-DNA","BC003-DNA","BC004-DNA","BC005-DNA","BC010-DNA","BC015-DNA","BC016-DNA","BC018-DNA","BC020-DNA","BC048-DNA","BC049-DNA","BC051-DNA","BC105-DNA","BC108-DNA","BC262-DNA","BCL01-DNA","BCL03-DNA","HNF132-DNA","HNF133-DNA","HNF134-DNA","HNF144-DNA","HNF168-DNA","HNF185-DNA","HNF187-DNA","HNF189-DNA","HNF190-DNA","HNF191-DNA","HNF216-DNA","HNF217-DNA","HNF221-DNA","HNF224-DNA","HNF225-DNA","HNF229-DNA","HNF236-DNA","HNF237-DNA","HNF242-DNA","HNF250-DNA","HNF267-DNA","HNF269-DNA","HNF279-DNA","YSF004-DNA","YSF117-DNA","YSF295-DNA","YSF296-DNA","YSF298-DNA","YSF300-DNA","YSF44-DNA","YSF45-DNA","YSF46-DNA","YSF47-DNA","YSF65-DNA","YSF66-DNA","YSF67-DNA","YSF69-DNA","YSF70-DNA","YSF71-DNA","YSF74-DNA")

keep_cDNA <- c("BC001-cDNA","BC002-cDNA","BC003-cDNA","BC004-cDNA","BC005-cDNA","BC010-cDNA","BC015-cDNA","BC016-cDNA","BC018-cDNA","BC020-cDNA","BC048-cDNA","BC049-cDNA","BC051-cDNA","BC105-cDNA","BC108-cDNA","BC262-cDNA","BCL01-cDNA","BCL03-cDNA","HNF132-cDNA","HNF133-cDNA","HNF134-cDNA","HNF144-cDNA","HNF168-cDNA","HNF185-cDNA","HNF187-cDNA","HNF189-cDNA","HNF190-cDNA","HNF191-cDNA","HNF216-cDNA","HNF217-cDNA","HNF221-cDNA","HNF224-cDNA","HNF225-cDNA","HNF229-cDNA","HNF236-cDNA","HNF237-cDNA","HNF242-cDNA","HNF250-cDNA","HNF267-cDNA","HNF269-cDNA","HNF279-cDNA","YSF004-cDNA","YSF117-cDNA","YSF295-cDNA","YSF296-cDNA","YSF298-cDNA","YSF300-cDNA","YSF44-cDNA","YSF45-cDNA","YSF46-cDNA","YSF47-cDNA","YSF65-cDNA","YSF66-cDNA","YSF67-cDNA","YSF69-cDNA","YSF70-cDNA","YSF71-cDNA","YSF74-cDNA")

# Active libraries - only 55
active <- Ponds.rel[rownames(Ponds.rel) %in% keep_cDNA, ]      # Extract rows from matrix

# total libraries - 58
total <- Ponds.rel[rownames(Ponds.rel) %in% keep_DNA, ]      # Extract rows from matrix

# modifying row names
row.names(total) <- c("B001","B002","B003","B004","B005","B010","B015","B016","B018","B020","B048","B049","B051","B105","B108","B262","BL01","BL03","H132","H133","H134","H144","H168","H185","H187","H189","H190","H191","H216","H217","H221","H224","H225","H229","H236","H237","H242","H250","H267","H269","H279","Y004","Y117","Y295","Y296","Y298","Y300","Y44","Y45","Y46","Y47","Y65","Y66","Y67","Y69","Y70","Y71","Y74")

row.names(active) <- c("B001","B003","B004","B005","B015","B016","B018","B020","B048","B049","B051","B105","B108","B262","BL01","BL03","H132","H133","H134","H144","H168","H185","H187","H189","H190","H191","H216","H217","H224","H225","H229","H236","H237","H242","H250","H267","H269","H279","Y004","Y117","Y295","Y296","Y298","Y300","Y44","Y45","Y46","Y47","Y65","Y66","Y67","Y69","Y70","Y71","Y74")


############ Visualization Analyses
# Bray Curtis resemblance matrix
total.db <- vegdist(total, method="bray")

####### Principal Component Analysis - TOTAL
total.pcoa <- cmdscale(total.db, eig=TRUE, k=3)

exvar1 <- round(total.pcoa$eig[1] / sum(total.pcoa$eig), 3) * 100
exvar2 <- round(total.pcoa$eig[2] / sum(total.pcoa$eig), 3) * 100
exvar3 <- round(total.pcoa$eig[3] / sum(total.pcoa$eig), 3) * 100
total.sum.eig <- sum(exvar1, exvar2, exvar3)

# PCoA Plot PC1 x PC2
plot(total.pcoa$points[ ,1], total.pcoa$points[ ,2], #ylim = c(-0.2, 0.7),
     xlab= paste("PCoA 1 (", exvar1, "%)", sep = ""),
     ylab= paste("PCoA 2 (", exvar2, "%)", sep = ""),
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5,
     cex.axis=1.2, axes=FALSE);
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1);
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1);
abline(h = 0, v = 0, lty = 3);
box(lwd = 2);                                        
points(total.pcoa$points[ ,1], total.pcoa$points[ ,2],
       pch = 1, cex = 2, bg = "red", col = "red");
text(total.pcoa$points[ ,1], total.pcoa$points[ ,2],
     labels = row.names(total.pcoa$points), adj=1)
## would be great to get the points colored by state park

######## Principal Component Analysis - ACTIVE
# bray curtis
active.db <- vegdist(active, method="bray")
active.pcoa <- cmdscale(active.db, eig=TRUE, k=3)

a.exvar1 <- round(active.pcoa$eig[1] / sum(active.pcoa$eig), 3) * 100
a.exvar2 <- round(active.pcoa$eig[2] / sum(active.pcoa$eig), 3) * 100
a.exvar3 <- round(active.pcoa$eig[3] / sum(active.pcoa$eig), 3) * 100
active.sum.eig <- sum(a.exvar1, a.exvar2, a.exvar3)

# PCoA Plot PC1 x PC2
plot(active.pcoa$points[ ,1], active.pcoa$points[ ,2], #ylim = c(-0.2, 0.7),
     xlab= paste("PCoA 1 (", a.exvar1, "%)", sep = ""),
     ylab= paste("PCoA 2 (", a.exvar2, "%)", sep = ""),
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5,
     cex.axis=1.2, axes=FALSE);
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1);
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1);
abline(h = 0, v = 0, lty = 3);
box(lwd = 2);                                        
points(active.pcoa$points[ ,1], active.pcoa$points[ ,2],
       pch = 1, cex = 2, bg = "red", col = "red");
text(active.pcoa$points[ ,1], active.pcoa$points[ ,2],
     labels = row.names(active.pcoa$points), adj=1)

# Chunck of code from synthesis question (Biodiversity worksheet 2)

library(vegan)
library(dplyr)
library(tidyverse)
library(ggplot2)

data <- load("/Users/joyobrien/GitHub/team2/INPond_Initial.RData")

DormDecay_env <- readRDS("/Users/joyobrien/GitHub/team2/DormDecay_env.rds")

# DNA
total_matrix <- Pond97[grep('-DNA', rownames(Pond97)),]

# Permanova 
adonis <- adonis2(total_matrix ~ DormDecay_env$Location, method = "bray", permutations = 999)
print(adonis)

library(indicspecies)
indval_pond <- multipatt(total_matrix, cluster = DormDecay_env$Location, func = "IndVal.g",
                         control = how(nperm = 999))
summary(indval_pond)
