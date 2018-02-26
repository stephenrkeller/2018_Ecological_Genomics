# WGCNA analysis of beetle RNASeq data from normalized counts of NC and WA populations

source("http://bioconductor.org/biocLite.R") 
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
install.packages("WGCNA")

library("WGCNA")

beetData <- read.csv("beetle_norm_counts.csv")
dim(beetData)

# Our initial dataset; no filtering
bdatExpr0 <- as.data.frame(t(beetData[,-1]))
dim(bdatExpr0)
names(bdatExpr0) = beetData$X

rownames(bdatExpr0) = names(beetData)[-c(1:1)]

# Use their function to pull out the good genes
# Beetle good sample genes, 'bgsg'

bgsg = goodSamplesGenes(bdatExpr0, verbose=3)

# Make a cluster tree to look for outlier samples

sampleTree <- hclust(dist(bdatExpr0), method="average")
# plot
sizeGrWindow(12,9)
par(cex=0.6)
par(mar=c(0,4,2,0))

plot(sampleTree, main="Sample clustering to detect outliers", sub="", xlab="",cex.lab=1.5, cex.axis=1.5, cex.main=2)

btraitData <- read.table("cols_data_noIT_num.txt", header=T)
head(btraitData)

# Network construction and module detection
powers = c(c(1:10,seq(from=12,to=20, by=2)))

# Call the network topology function

sft = pickSoftThreshold(bdatExpr0, powerVector = powers, verbose=5)

sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold(power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     type="n",main=paste("Scale Independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=1, col="red")
abline(h=0.9, col="red")

# Mean connectivity as a function of the soft-thresholding power

plot(sft$fitIndices[,1],sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean connectivity",
     type="n",main="Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5],labels=powers, cex=cex1, col="red")

# Now onto bulding the dendrogram; TOM = "topology modules"

net = blockwiseModules(bdatExpr0, power=6, TOMType="unsigned",
                       minModuleSize=30,reassignThreshold=0, mergeCutHeight=0.25,
                       numericLabels=T, pamRespectsDendro=F, saveTOMs=T, saveTOMFileBase="beetleTOM",
                       verbose=3)

sizeGrWindow(12,9)
mergedColors = labels2colors(net$colors)

# Plot the dendrogram with the module colors underneath

plotDendroAndColors(net$dendrograms[[1]],mergedColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels=F, hang=0.03, addGuide=T, guideHang=0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
geneTree = net$dendrograms[[1]]

# Define numbers of genes and samples
nGenes = ncol(bdatExpr0)
nSamples = nrow(bdatExpr0)

MEs0 = moduleEigengenes(bdatExpr0, moduleColors)$eignengenes
MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, btraitData, use="p")
moduleTraitsPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix) = dim(modeuleTraitCor)
par(mar=c(6,8.5,3,3))

labeledHeatmap(Matrix=moduleTraitCor, xLabels=names(btraitData),yLabels=names(MEs), 
               ySymbols=names(MEs),colorLabels=F, colors=greenWhiteRed(50),
               textMatrix=textMatrix,setStdMargins=F, cex.text=0.5,zlim=c(-1,1),
               main="Module-trait relationships")