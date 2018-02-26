setwd("~/github/2018_Ecological_Genomics/data/")
library("DESeq2")

library("ggplot2")

countsTable <- read.delim('allcountsdataRN_noIT.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_noIT.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)


######################################

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ devstage + sex + population)
# In this typical model, the "sex effect" represents the overall effect controlling for differences due to population and devstage. page 26-27 manual DESeq2.pdf
# The last term in the model is what is tested.  In this case sex.
#  This is not the same as an interaction.


dim(dds)
# [1] 17483    48

dds <- dds[ rowSums(counts(dds)) > 1, ]
dim(dds)
# [1] 16851    48  # Nice, just lost about 600 genes with not more than one read across the 48 samples

dds <- DESeq(dds, modelMatrixType="standard")

resultsNames(dds)
# [1] "Intercept"           "devstage_L3L_vs_AD4" "devstage_PD1_vs_AD4" "devstage_PP1_vs_AD4" "sex_M_vs_F"          "population_WA_vs_NC"

# the following line extracts results for a main effects (here population)
res_pop <- results(dds, name="population_WA_vs_NC", alpha=0.05)
res_pop <- res_pop[order(res_pop$padj),]
summary(res_pop)
# out of 16851 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 159, 0.94% 
# LFC < 0 (down)   : 147, 0.87% 
# outliers [1]     : 133, 0.79% 
# low counts [2]   : 1939, 12% 
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
head(res_pop)
# log2 fold change (MAP): population WA vs NC 
# Wald test p-value: population WA vs NC 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat       pvalue         padj
# <numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
#   OTAU008667-RA 231.87111     -0.7803753 0.1132661 -6.889755 5.588875e-12 4.129899e-08
# OTAU012562-RA 251.77436     -0.7828048 0.1135355 -6.894800 5.394059e-12 4.129899e-08
# OTAU011160-RA  10.24152     -1.0863185 0.1836927 -5.913781 3.343429e-09 1.438959e-05
# OTAU012716-RA 188.87768      1.0137238 0.1721500  5.888609 3.894604e-09 1.438959e-05
# OTAU002976-RA 998.42315     -0.6159655 0.1075972 -5.724735 1.035950e-08 3.062062e-05
# OTAU014686-RA 603.66091     -0.8168345 0.1468336 -5.562995 2.651835e-08 6.531911e-05

length(which(res_pop$padj > 0.05))

############## output to save and for enrichment analyses.

write.csv(res_pop, file = "DGE_NCvsWA_pop.csv", row.names = T, quote = F)


################################################################
# Data quality assessment by sample clustering and visualization 

plotMA(res_pop, main="DESeq2", ylim=c(-2,2))
abline(h=c(-1,1), col="dodgerblue", lwd=2)

# How would this look if we pulled out diffs between sexes? or devstages?

res_sex <- results(dds, name="sex_M_vs_F", alpha=0.05)
plotMA(res_sex, main="DESeq2", ylim=c(-2,2))
abline(h=c(-1,1), col="dodgerblue", lwd=2)

res_sex <- results(dds, name="sex_M_vs_F", alpha=0.05)
plotMA(res_sex, main="DESeq2", ylim=c(-3,3))
abline(h=c(-1,1), col="dodgerblue", lwd=2)


################################  let's look at the top individual genes

d <- plotCounts(dds, gene="OTAU002976-RA", intgroup=(c("population","sex","devstage")), returnData=TRUE)
d
p <- ggplot(d, aes(x= devstage, y=count, shape = sex, colour = population)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) + scale_x_discrete(limits=c("L3L","PP1","PD1","AD4"))
p

ggsave("point_plot_DGE_top_pop_OTAU012562.png", p, width=8, height=4, dpi=300)



# 2.2 Data quality assessment by sample clustering and visualization 
￼
pdf(file="PCA_1v2_stage.sex.pdf", height=5.5, width=5.5)
data <- plotPCA(rld, intgroup=c("population", "devstage", "sex"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

data$devstage <- factor(data$devstage, levels=c("L3L","PP1","PD1","AD4"), labels=c("L3L","PP1","PD1","AD4"))
# puts our levels in order

# plots the data generated from the plotPCA function using the ggplot package
ggplot(data, aes(PC1, PC2, color=sex, shape=devstage)) +
  geom_point(size=4, alpha = 0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal() + theme(text = element_text(size=15), panel.grid.major = element_line(colour = "grey"))
dev.off()

# how does the plot look showing the population factor as color
ggplot(data, aes(PC1, PC2, color=population, shape=devstage)) +
  geom_point(size=4, alpha = 0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal() + theme(text = element_text(size=15), panel.grid.major = element_line(colour = "grey"))

# the next several lines make a cluster heatmap of all the samples
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$population, rld$devstage, rld$sex, sep="-")
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
library("pheatmap")

pdf(file="clustering_DistMatrix.pdf", height=5, width=5.5)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
dev.off()


#########







############################################################

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ population * devstage * sex)
# Interaction

dim(dds)
# [1] 17483    48

dds <- dds[ rowSums(counts(dds)) > 1, ]
dim(dds)
# [1] 16851    48  # Nice, just lost about 600 genes with not more than one read across the 48 samples

dds <- DESeq(dds, modelMatrixType="standard")

resultsNames(dds)
# [1] "Intercept"                     "population_WA_vs_NC"           "devstage_L3L_vs_AD4"           "devstage_PD1_vs_AD4"           "devstage_PP1_vs_AD4"           "sex_M_vs_F"                   
# [7] "populationWA.devstageL3L"      "populationWA.devstagePD1"      "populationWA.devstagePP1"      "populationWA.sexM"             "devstageL3L.sexM"              "devstagePD1.sexM"             
# [13] "devstagePP1.sexM"              "populationWA.devstageL3L.sexM" "populationWA.devstagePD1.sexM" "populationWA.devstagePP1.sexM"

# the following line extracts results for a main effects (here population)
res_pop <- results(dds, name="population_WA_vs_NC", alpha=0.01)
res_pop <- res_pop[order(res_pop$padj),]
head(res_pop)

# the following line extracts results for one of the 3-way interactions
res_3int1 <- results(dds, name="populationWA.devstageL3L.sexM", alpha=0.01)
res_3int1 <- res_3int1[order(res_3int1$padj),]
head(res_3int1,13)

sig_3int1 <- res_3int1[which(res_3int1$padj < 0.05), ]
dim(sig_3int1) # 12 significant here with this 3-way interaction

# the following line extracts results for the second of the 3-way interactions
res_3int2 <- results(dds, name="populationWA.devstagePD1.sexM", alpha=0.01)
res_3int2 <- res_3int2[order(res_3int2$padj),]
head(res_3int2)

sig_3int2 <- res_3int2[which(res_3int2$padj < 0.05), ]
dim(sig_3int2) # 9 significant here with this 3-way interaction

# the following line extracts results for the third of the 3-way interactions
res_3int3 <- results(dds, name="populationWA.devstagePP1.sexM", alpha=0.01)
res_3int3 <- res_3int3[order(res_3int3$padj),]
head(res_3int3)

sig_3int3 <- res_3int3[which(res_3int3$padj < 0.05), ]
dim(sig_3int3) # 7 significant here with this 3-way interaction

## pull out only results with padj <0.05 for all 3 3-way inte
sig_3int1_df<-as.data.frame(sig_3int1)
sig_3int1_df$Row.names<-rownames(sig_3int1_df)
dim(sig_3int1_df)

sig_3int2_df<-as.data.frame(sig_3int2)
sig_3int2_df$Row.names<-rownames(sig_3int2_df)
dim(sig_3int2_df)

sig_3int3_df<-as.data.frame(sig_3int3)
sig_3int3_df$Row.names<-rownames(sig_3int3_df)
dim(sig_3int3_df)

genesOfInterest<-c(sig_3int1_df$Row.names, sig_3int2_df$Row.names, sig_3int3_df$Row.names)
length(genesOfInterest) #check; yes, 28 total

m <- assay(vsd)[genesOfInterest, ]  # pulls out a matrix m of the normalized counts for each of the sig genes for all samples; can specifiy which columns/samples if wanting to plot just some
head(m)

head(rowMeans(m))
m2 <- m - rowMeans(m)  # subtracts the row mean from the value in each cell in a row; making a new matrix m2
head(m2)
sd <- apply(m,1,sd) # calculates the standard deviation for each row; making a list of sds for each gene
head(sd)
z <- m2/sd  # divides the value in each cell for each row of the m2 matrix by the sd of that row; makes a new matrix z to be used in heatmaps
dim(z)
head(z)

sers <- apply(z,1,sd) # double check that the standard deviation of each row in the new z matrix should be 1
head(sers) # check


#make heatmap
pheatmap(z, labels_col=colData(vsd)$pn, cluster_cols=F, cluster_rows=F)
pheatmap(z, labels_col=colData(vsd)$pn, cluster_cols=F, cluster_rows=T)
pheatmap(z, labels_col=colData(vsd)$pn, cluster_cols=T, cluster_rows=T)

# try heatmaps without normaliztion
pheatmap(m, labels_col=colData(vsd)$pn, cluster_cols=F, cluster_rows=F)
pheatmap(m, labels_col=colData(vsd)$pn, cluster_cols=F, cluster_rows=T)
pheatmap(m, labels_col=colData(vsd)$pn, cluster_cols=T, cluster_rows=T)


# Day 7 Transcriptomics (Feb 21, 2018)

colData$group <- factor(paste0(colData$population, "-", colData$devstage, "-", colData$sex))
head(colData)

dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~group)
dds <- dds[rowSums(counts(dds))>1,]
dim(dds)

dds <- DESeq(dds, parallel = T)

resultsNames(dds)

# to comment out all highlighted tx, use command-shift-C
# :)

# [1] "Intercept"     "groupNC.AD4.F" "groupNC.AD4.M" "groupNC.L3L.F"
# [5] "groupNC.L3L.M" "groupNC.PD1.F" "groupNC.PD1.M" "groupNC.PP1.F"
# [9] "groupNC.PP1.M" "groupWA.AD4.F" "groupWA.AD4.M" "groupWA.L3L.F"
# [13] "groupWA.L3L.M" "groupWA.PD1.F" "groupWA.PD1.M" "groupWA.PP1.F"
# [17] "groupWA.PP1.M"

# Want to generate contrasts between particular levels of the group variable to test custom hypotheses:

res_pop_PP1_F <- results(dds, contrast=list(c("groupNC.PP1.F"),c("groupWA.PP1.F")), listValues=c(1/2,-1/2), alpha=0.05)

res_pop_PP1_F <- res_pop_PP1_F[order(res_pop_PP1_F$padj),]

head(res_pop_PP1_F)

######## Pull out significant DEG's to make a heatmap
sig_pop_PP1_F <- res_pop_PP1_F[which(res_pop_PP1_F$padj <0.05),]
dim(sig_pop_PP1_F)

sig_pop_PP1_F_df <- as.data.frame(sig_pop_PP1_F)
sig_pop_PP1_F_df$Row.names <- rownames(sig_pop_PP1_F_df)

genesOfInterest_pop_PP1_F <- c(sig_pop_PP1_F_df$Row.names)

vsd <- vst(dds, blind=F)

dds$combined = factor(paste0(dds$population, "-", dds$devstage, "-", dds$sex))
dds$combined <- factor(dds$combined, levels=c("WA-L3L-F","WA-L3L-M","WA-PP1-F","WA-PP1-M","WA-PD1-F","WA-PD1-M","WA-AD4-F","WA-AD4-M","NC-L3L-F","NC-L3L-M","NC-PP1-F","NC-PP1-M","NC-PD1-F","NC-PD1-M","NC-AD4-F","NC-AD4-M"), labels=c("WA-L3L-F","WA-L3L-M","WA-PP1-F","WA-PP1-M","WA-PD1-F","WA-PD1-M","WA-AD4-F","WA-AD4-M","NC-L3L-F","NC-L3L-M","NC-PP1-F","NC-PP1-M","NC-PD1-F","NC-PD1-M","NC-AD4-F","NC-AD4-M"))

baseMeanPerGrp <- sapply( levels(dds$combined), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$combined == lvl] ) )

head(baseMeanPerGrp)
dim(baseMeanPerGrp)

# Pulls out normalized counts (average per 3 reps) for all of our significant genes
m <- baseMeanPerGrp[genesOfInterest_pop_PP1_F, c("WA-PP1-F","WA-PP1-M","NC-PP1-F","NC-PP1-M")]

mat_scaled <- t(apply(m, 1, scale)) # sets a mean of zero and a SD=1 for each gene

library(pheatmap)

pheatmap(mat_scaled, 
         labels_col=c("WA-PP1-F","WA-PP1-M","NC-PP1-F","NC-PP1-M"), 
         cluster_cols=F,
         cluster_rows=T)

# Now, let's export our counts to use as input to WCGNA

norm.counts <- counts(dds, normalized=T)
dim(norm.counts)

write.csv(norm.counts, file="beetle_norm_counts.csv", row.names=T, quote=F)

############### Side discussion:
# Apparentaly, there are 2 different version of DESeq2 that lead to different results for significant DEG's
# What's the degree of overlap?

ethan <- read.csv("~/Desktop/Gene_list", header=T)
overlap <- length(which(sig_pop_PP1_F_df$Row.names %in% ethan$X)==T)

# 58 in common; 19 unique to v18; 6 unique to v14
###############

