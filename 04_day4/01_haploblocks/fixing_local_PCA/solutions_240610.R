setwd("/Users/maffe217/Documents/2024/Courses/adaptation_genomics/resources")
setwd("/Users/maffe217/Documents/2024/Courses/adaptation_genomics/04_day4/01_haploblocks/") 

# open library
library(lostruct)
options(datatable.fread.input.cmd.message = FALSE)  # disable a useless message
snps <- vcf_windower("00_localPCA/capelin_NWA_sorted.bcf", size = 5, type = "snp", 
                     sites = vcf_positions("00_localPCA/capelin_NWA_sorted.bcf"))

snps(5)
region(snps) (5)

pcs <- eigen_windows(snps, k = 2)
dim(pcs) # check dimension
head (pcs[, 1:10]) # look at the first 10 columns
pcs_noNA <- pcs[-which(is.na(pcs[, 1])), ] # because of NA, some windows were not computed by pca. we will remove them

# retrieve positions
window_pos <- region(snps)()
head(window_pos)

# keep windows without NA
window_pos_noNA <- window_pos[-which(is.na(pcs[, 1])), ]
# merge
pca_matrix_noNA <- cbind(window_pos_noNA, pcs_noNA)
head(pca_matrix_noNA[, 1:10])

# save the file
write.table(pca_matrix_noNA, "00_localPCA/pca_matrix_5SNPs.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# load libraries
library(data.table)
library(lostruct)

# load matrix
pca_matrix <- read.table("00_localPCA/pca_matrix_5SNPs.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
head(pca_matrix)
# split columns with positions information and PC
window_pos <- pca_matrix[, 1:3]
pcs <- as.matrix(pca_matrix[, 4:dim(pca_matrix)[2]])

pcdist <- pc_dist(pcs,npc = 2)
mds_axe <- cmdscale(pcdist, k = 10)
head(mds_axe)

mds_matrix <- cbind(window_pos, mds_axe)
write.table(mds_matrix, "00_localPCA/mds_matrix_5SNPs.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Locally:
pca_matrix <- read.table("00_localPCA/pca_matrix_5SNPs.txt", header = TRUE)
pca_matrix[1:10, 1:10]
n_windows <- dim(pca_matrix)[1] # the number of windows we have

Nind <- 240
i = 15 # for the 15th window

pc1_i <- t(pca_matrix[i, 7:(Nind + 6)]) # scores along PC1
pc2_i <- t(pca_matrix[i, (Nind + 7):(2 * Nind + 6)]) # scores along PC2
var1 <- round(pca_matrix[i, 5] / pca_matrix[i, 4], 2) * 100 # % of variance explained by PC1
var2 <- round(pca_matrix[i, 6] / pca_matrix[i, 4], 2) * 100 # % of variance explained by PC2
midpos_i <- (pca_matrix[i, 2] + pca_matrix[i, 3]) / 2 # average position of the window
window_i <- paste(pca_matrix[i, 1], midpos_i , sep = "_") # paste the name of CHR and the midposition

plot(pc1_i, pc2_i, pch = 20, xlab = paste("PC1", var1, "%"), ylab = paste("PC2", var2, "%"), main = window_i)

geno <- read.table("00_localPCA/canada.012")[, -1] # load geno
geno[1:6, 1:6] # check the geno matrix
global.pca <- prcomp(geno) # run the pca
plot(global.pca$x[, 1], global.pca$x[, 2]) # plot the pca
PC_of_interest <- global.pca$x[, 1] # if you want to look at correlation with PC1

# initialise the vector
corr_vector <- vector(length = n_windows)
# loop over windows to store correlation factor
for (i in 1 : n_windows) {
  pc1_i <- t(pca_matrix[i, 7:(Nind + 6)]) # scores along PC1
  corr_vector[i] <- abs(cor(PC_of_interest, pc1_i)[1, 1])
}

#
pca_correlation <- cbind(pca_matrix[, 1:3], corr_vector)
pca_correlation$midpos <- (pca_correlation$start + pca_correlation$end) / 2
head(pca_correlation)

plot<-ggplot(pca_correlation, aes(x = midpos, y = corr_vector, colour = chrom)) +
  geom_point() +
  theme_classic() +
  facet_grid(cols = vars(chrom), scales = "free_x", space = "free_x")

ggsave(plot, height=4, width=5,
       file="/Users/maffe217/Library/CloudStorage/Dropbox/2024_physalia_adaptation_course/physalia_adaptation_course-2024/04_day4/01_haploblocks/00_localPCA/images/local_pca_corr_PC1_240610.png")

# Linkage disequilibrium ####
library(ggplot2)
# load data
chr4.ld <- read.table("03_ld/maf0.05_chr4.ld", header = TRUE)
head(chr4.ld)

# plot (very simple solution with ggplot. maybe you can find something nicer..?
ggplot(chr4.ld, aes(x = BP_A, y = BP_B, col = R2)) + theme_classic() + geom_point(size = 1, shape = 15) + 
  scale_colour_gradientn(colours = c("lightgrey", "deepskyblue3", "blue", "blue3", "navyblue", "black"), limits = c(0,1), name = "R2")

# load data homozygotes
AA_chr4.ld <- read.table("03_ld/AA_maf0.05_chr4.ld", header = TRUE)
head(AA_chr4.ld)

ggplot(AA_chr4.ld, aes(x = BP_A, y = BP_B, col = R2)) + theme_classic() + geom_point(size = 1, shape = 15) + 
  scale_colour_gradientn(colours = c("lightgrey", "deepskyblue3", "blue", "blue3", "navyblue", "black"), limits = c(0,1),  name = "R2")

# plotting both heatmap on the same graph...
ggplot(chr4.ld, aes(x = BP_A, y = BP_B, col = R2)) + theme_classic() + 
  geom_point(shape = 15) + 
  geom_point(data = AA_chr4.ld, aes(x = BP_B, y = BP_A, col = R2), size = 1, shape = 15) +
  scale_colour_gradientn(colours = c("lightgrey", "deepskyblue3", "blue", "blue3", "navyblue", "black"), limits = c(0, 1),  name = "R2")

# load data homozygotes
BB_chr4.ld <- read.table("03_ld/BB_maf0.05_chr4.ld", header = TRUE)
head(BB_chr4.ld)

ggplot(BB_chr4.ld, aes(x = BP_A, y = BP_B, col = R2)) + theme_classic() + geom_point(size = 1, shape = 15) + 
  scale_colour_gradientn(colours = c("lightgrey", "deepskyblue3", "blue", "blue3", "navyblue", "black"), limits = c(0,1),  name = "R2")


# plotting both heatmap on the same graph...
ggplot(chr4.ld, aes(x = BP_A, y = BP_B, col = R2)) + theme_classic() + 
  geom_point(shape = 15) + 
  geom_point(data = BB_chr4.ld, aes(x = BP_B, y = BP_A, col = R2), size = 1, shape = 15) +
  scale_colour_gradientn(colours = c("lightgrey", "deepskyblue3", "blue", "blue3", "navyblue", "black"), limits = c(0, 1),  name = "R2")

# Fst between groups:
AA_BB <- read.table("04_divergence/AA_BB.weir.fst", header = TRUE)
head(AA_BB)

ggplot(AA_BB, aes(x = POS / 1000000, y = WEIR_AND_COCKERHAM_FST, col = CHROM)) +
  geom_point() + geom_smooth() +
  theme_classic() +
  facet_grid(cols = vars(CHROM), scales = "free_x", space = "free_x") +
  labs(x = "position (in MB)")

AA_BB.win <- read.table("04_divergence/AA_BB.windowed.weir.fst", header = TRUE)
head(AA_BB.win)

ggplot(AA_BB.win, aes(x = (BIN_START+BIN_END) / 2000000, y = WEIGHTED_FST, col = CHROM)) +
  geom_point() + geom_smooth() +
  theme_classic() +
  facet_grid(cols = vars(CHROM), scales = "free_x", space = "free_x") +
  labs(x = "position (in MB)")

AA_AB <- read.table("04_divergence/AA_AB.weir.fst", header = TRUE)
head(AA_AB)
outlier_Chr4 <- AA_AB[AA_AB$WEIR_AND_COCKERHAM_FST >= 0.1, ]
write.table(outlier_Chr4, "04_divergence/outlier_Chr4.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# HW equilibrium
AB.hwe <- read.table("05_heterozygosity/AB_formatted.hwe", header = TRUE)

# rename columns
colnames(AB.hwe) <- c("CHR", "POS", "Homo1_obs", "Het_Obs", "Homo2_Obs", "Homo1_Exp", "Het_Exp","Homo2_Exp", "Chisq_HWE", "P_HWE", "P_HET_DEFICIT", "P_HET_EXCESS")
head(AB.hwe)
# calculate the fraction of observed heterozygotes
AB.hwe$het_fraction <- AB.hwe$Het_Obs / (AB.hwe$Homo1_obs + AB.hwe$Het_Obs + AB.hwe$Homo2_Obs)

# plot
ggplot(AB.hwe, aes(x = POS / 1000000, y = het_fraction, col = CHR)) +
  geom_point(alpha = 0.5) + geom_smooth() +
  theme_classic() +
  facet_grid(cols = vars(CHR), scales = "free_x", space = "free_x") +
  labs(x = "position (in MB)")

BB.hwe <- read.table("05_heterozygosity/BB_formatted.hwe", header = TRUE)
AA.hwe <- read.table("05_heterozygosity/AA_formatted.hwe", header = TRUE)
# rename columns
colnames(AA.hwe) <- c("CHR", "POS", "Homo1_obs", "Het_Obs", "Homo2_Obs", "Homo1_Exp", "Het_Exp","Homo2_Exp", "Chisq_HWE", "P_HWE", "P_HET_DEFICIT", "P_HET_EXCESS")
colnames(BB.hwe) <- c("CHR", "POS", "Homo1_obs", "Het_Obs", "Homo2_Obs", "Homo1_Exp", "Het_Exp","Homo2_Exp", "Chisq_HWE", "P_HWE", "P_HET_DEFICIT", "P_HET_EXCESS")
# calculate the fraction of observed heterozygotes
AA.hwe$het_fraction <- AA.hwe$Het_Obs / (AA.hwe$Homo1_obs + AA.hwe$Het_Obs + AA.hwe$Homo2_Obs)
BB.hwe$het_fraction <- BB.hwe$Het_Obs / (BB.hwe$Homo1_obs + BB.hwe$Het_Obs + BB.hwe$Homo2_Obs)

BB.hwe$geno <- "BB"
AB.hwe$geno <- "AB"
AA.hwe$geno <- "AA"
all.hwe <- rbind(AA.hwe, AB.hwe, BB.hwe)
head(all.hwe)

# plot, colouring by genotype
ggplot(all.hwe, aes(x = POS / 1000000, y = het_fraction, group = geno, col = geno)) +
  geom_point(alpha = 0.5) + geom_smooth() +
  theme_classic() +
  facet_grid(cols = vars(CHR), scales = "free_x", space = "free_x") +
  labs(x = "position (in MB)")

ggplot(all.hwe, aes(x = geno, y = het_fraction, col = geno)) +
  geom_violin() +
  stat_summary(fun = median, geom = "point", size = 2) +
  facet_grid(cols = vars(CHR))+
  theme_classic()

# plots inside and outside the inversion
library(tidyverse)

left_breakpoint<-4.8*1e6
right_breakpoint<-16.6*1e6
chromosome_inversion<-"Chr4"

# create a data.frame with the SNPs inside the inversion
all.hwe.chr4.inversion<-all.hwe %>% filter(CHR==chromosome_inversion & POS >= left_breakpoint & POS <=right_breakpoint)

# find chromosome 4 positions that are not inside the inversion
# create a filtering vector with the positions inside the inversion:
positions_in_inversion<-all.hwe.chr4.inversion$POS
# filter data.frame:
all.hwe.chr4.outside<-all.hwe %>% 
  filter(CHR==chromosome_inversion) %>%
  filter(!POS %in% positions_in_inversion)

# add a classifier to column location:
all.hwe.chr4.inversion$location<-"inversion"
all.hwe.chr4.outside$location<-"outside"

#merge data.frames:
all.hwe.chr4<-rbind(all.hwe.chr4.inversion, all.hwe.chr4.outside)

# Plot. You can see that inside the inversion, heterozigosity is higher
ggplot(all.hwe.chr4, aes(x = geno, y = het_fraction, col = geno)) +
  geom_violin() +
  stat_summary(fun = median, geom = "point", size = 2) +
  facet_grid(cols = vars(location))+
  theme_classic()


ggsave("06_images/Hobs_violin_inside_vs_outside.png")

# Genotyping ####
geno <- read.table("02_data/canada.chr4inv.012")[, -1] # load geno
indv <- read.table("02_data/canada.chr4inv.012.indv") # load individuals info
rownames(geno) <- indv[, 1]
geno[1:6, 1:6] # check the geno matrix
geno.pca <- prcomp(geno) # run the pca
plot(geno.pca$x[, 1], geno.pca$x[, 2]) # plot the pca

geno_kmean <- kmeans(geno.pca$x[, 1], c(min(geno.pca$x[, 1]), (min(geno.pca$x[, 1]) + max(geno.pca$x[, 1])) / 2, max(geno.pca$x[, 1])))
geno_kmean
plot(geno.pca$x[, 1],geno.pca$x[, 2], col = geno_kmean$cluster, pch = 20)

info_cluster_inv <- cbind(indv, geno_kmean$cluster)
colnames(info_cluster_inv) <- c("id_inv", "cluster_inv")
head(info_cluster_inv)

# prepare the list
AA_ind <- info_cluster_inv[info_cluster_inv$cluster_inv == "1", 1]
AB_ind <- info_cluster_inv[info_cluster_inv$cluster_inv == "2", 1]
BB_ind <- info_cluster_inv[info_cluster_inv$cluster_inv == "3", 1]

