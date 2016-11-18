########################################
#
# Code to reproduce the specificity score analysis in  
# Javierre, Burren, Wilder, Kreuzhuber, Hill et al. Cell 167, 1369-1384 (2016), DOI: 10.1016/j.cell.2016.09.037.
#
########################################

#####
# load required libraries
####

library(gplots)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)



##########
# specificity score function
##########

specificityScore <- function(x,distMatrix=NULL){
  # function to calculate cell type specificity scores for a vector of values x
  # Inputs:
  # x is of length # cell types -i.e. one value for each cell type
  # distMatrix is a square, symmetric matrix of length # cell types giving the distance between cell types
  # Output:
  # A vector of length # cell types, giving the specificity score for each cell type
  
  nCellType <- length(x)
  specScore <- vector(length=nCellType)
  
  # if no distMatrix, defaults to equal distances
  if(is.null(distMatrix)) {
    distMatrix <- matrix(1,nrow=nCellType,ncol=nCellType)
    diag(distMatrix) <- 0
  }
  
  for (i in 1:nCellType) {
    # calculate specificity score for each cell type
    if(is.na(x[i])) {
      specScore[i] <- NA
    } else {
      idx <- !is.na(x)
      specScore[i] <- sum((x[i]-x[idx])*(distMatrix[i,idx])) / sum(distMatrix[i,])
    }
  }
  return(specScore)
}



##########
# processing of PCHi-C interaction data
##########

# Start with PCHi-C peak matrix containing CHiCAGO scores for all interactions that pass a cutoff of >=5 in at least one cell type.

# peak matrix file is available at https://osf.io/u8tzp/.
pks <- read.table("PCHiC_peak_matrix_cutoff5.txt",header=TRUE)
pksOrig <- pks

cellTypes <- colnames(pks)[12:28]
cellTypesOrig <- cellTypes

# remove cell types from peak matrix for which there are no annotations or expression data
pks <- pks[,!(colnames(pks) %in% c("EP","FoeT","tCD4","aCD4","naCD4","nCD8","tCD8","nB","tB"))]
cellTypes <- colnames(pks)[12:19]

# remove interactions with no interaction score >5 in any of the remaining cell types
pks <- pks[apply(pks[,cellTypes]>=5,1,any),] # reduces 728,838 interactions down to 527,726

# load in PIR active status annotations
activeIndicators = as.data.frame(readRDS("PIRactivity.Rds"))

# replace interaction scores with no active PIR by NA
pksNA <- pks[,cellTypes]
pksNA[!activeIndicators[,cellTypes]] <- NA
pksNA <- cbind(pks[,!(colnames(pks) %in% cellTypes)],pksNA)

# replace interaction scores with no active regulatory region at other end by zero
pksZero <- pks[,cellTypes]
pksZero[!activeIndicators[,cellTypes]] <- 0
pksZero <- cbind(pks[,!(colnames(pks) %in% cellTypes)],pksZero)

# remove interactions that do not have at least one cell type with an active PIR and with score >5
removedRowsNoneActiveAbove5 <- !apply(pksNA[,cellTypes]>=5,1,any,na.rm=TRUE)
pksZero <- pksZero[!removedRowsNoneActiveAbove5,] # leaves 139,835 interactions

# asinh transform
pksZero[,cellTypes] <- asinh(pksZero[,cellTypes])

# upper bound, any values >4.3 are set to 4.3
tmp <- pksZero[,cellTypes]
tmp[tmp>4.3] <- 4.3
pksZero[,cellTypes] <- tmp


# For each bait (gene), form a matrix consisting of the interaction scores for each PIR for that bait
# Rows are cell types, columns are PIRs.
# Matrices stored in a list
geneIDs <- unique(pksZero$baitID)        # 15,130 genes (baits)
gene <- by(pksZero,pksZero$baitID,`[`,,c("oeID",cellTypes))
gene <- lapply(gene,function(x) {
  row.names(x) <- x$oeID
  x[,1] <- NULL
  return(t(x))
})

# assign bait IDs and gene names to list
tmp <- as.vector(pksZero$baitName) 
geneNames <- sapply(1:length(gene),function(i){
  tmp[which(pksZero$baitID==geneIDs[i])[1]]
})
names(gene) <- paste(geneIDs,geneNames,sep=": ")

# remove baits that have 0 or more than 1 protein coding interactions or have a missing annotation
baitAnnotations = readRDS("baitAnnotations.Rds")
proteinCoding <- grepl("protein_coding",baitAnnotations$annotation)
proteinCoding[baitAnnotations$annotation==""] <- NA
useGene <- vector(length=length(geneIDs))
for(i in 1:length(geneIDs)) {
  nProtCod <- sum(proteinCoding[baitAnnotations$baitID==geneIDs[i]])
  if(!is.na(nProtCod) && nProtCod==1) {
    useGene[i] <- TRUE
  }
}
gene <- gene[useGene] # results in 9,805 genes (baits)
geneNames <- geneNames[useGene]
geneIDs <- geneIDs[useGene]

# remove baits that have more than one annotation
useGene <- !grepl(";",geneNames)
gene <- gene[useGene] # results in 7,004 genes (baits)
geneNames <- geneNames[useGene]
geneIDs <- geneIDs[useGene]



#####
# calculate gene specificity scores (based on PCHi-C data)
####

# calculate distances between cell types based on original full peak matrix
pksOrig[,cellTypesOrig] <- asinh(pksOrig[,cellTypesOrig]) # asinh transform
tmp <- pksOrig[,cellTypesOrig]
tmp[tmp>4.3] <- 4.3 # upper bound, any values >4.3 are set to 4.3
pksOrig[,cellTypesOrig] <- tmp
distMatrix <- dist(t(pksOrig[,cellTypesOrig]))
distMatrix <- as.matrix(distMatrix)[cellTypes,cellTypes]

# calculate specificity scores for each PIR, for each gene
nGenes <- length(gene)
geneSpecificityMatrices <- lapply(gene,function(g) {
  out <- apply(g,2,specificityScore,distMatrix)
  rownames(out) <- rownames(g)
  return(out)
})

# For each gene, average specificity scores across PIRs to obtain a single specificity score for each cell type
geneSpecificityScoresIntn <- sapply(geneSpecificityMatrices,rowMeans,na.rm=TRUE)
rownames(geneSpecificityScoresIntn) <- rownames(gene[[1]])
colnames(geneSpecificityScoresIntn) <- names(gene)



##########
# kmeans clustering of gene specificity scores (based on PCHi-C data)
#########

set.seed(79634308) # set seed to reproduce results in paper
km <- kmeans(t(geneSpecificityScoresIntn),centers=12,iter.max=10000,nstart=10000)



##########
# Figure 4B - heatmap of gene specificity scores (based on PCHi-C data) with genes grouped into kmeans clusters
##########

clusterCols <- c("#000000","#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6","#63FFAC")

# set colours for heatmap
my_palette1 <- colorRampPalette(c("blue","white","red"))(n = 256)
my_palette1 <- c("#CCCCCC",my_palette1)
m <- max(abs(geneSpecificityScoresIntn),na.rm=TRUE)+1e-16
col.breaks <- seq(from=-m-(2*m/256),to=m,length.out=258)

# sort genes into clusters and reorder data accordingly
sortedCluster <- sort(km$cluster,index.return=TRUE)
geneSpecificityScoresIntnSorted <- geneSpecificityScoresIntn[,sortedCluster$ix]
geneIDsSorted <- geneIDs[sortedCluster$ix]
geneNamesSorted <- geneNames[sortedCluster$ix]
geneSorted <- gene[sortedCluster$ix]
geneSpecificityMatricesSorted <- geneSpecificityMatrices[sortedCluster$ix]
  
# reorder clusters so similar clusters are together
tmp <- aggregate(t(geneSpecificityScoresIntnSorted),list(sortedCluster$x),mean,na.rm=TRUE)[,-1]
tmp <- hclust(dist(tmp))$order
tmp <- sapply(sortedCluster$x,function(x) which(tmp==x))
sortedCluster2 <- sort(tmp,index.return=TRUE)
geneSpecificityScoresIntnSorted <- geneSpecificityScoresIntnSorted[,sortedCluster2$ix]
geneIDsSorted <- geneIDsSorted[sortedCluster2$ix]
geneNamesSorted <- geneNamesSorted[sortedCluster2$ix]
geneSorted <- geneSorted[sortedCluster2$ix]
geneSpecificityMatricesSorted <- geneSpecificityMatricesSorted[sortedCluster2$ix]
clusterIDsorted <- sortedCluster2$x
mappingIdx <- sortedCluster$ix[sortedCluster2$ix]
  
clusterCols2 <- clusterCols[clusterIDsorted]
  
tmp <- geneSpecificityScoresIntnSorted
tmp[is.na(tmp)] <- -m-(2*m/300)
pdf("Figure4B.pdf",width=20, height=6)
hmap <- heatmap.2(tmp,scale="none",trace="none",dendrogram='row',Rowv=TRUE,Colv=FALSE,col=my_palette1,margins=c(5,10),
                  labCol=NA,breaks=col.breaks,ColSideColors=clusterCols2,tracecol=NA,key.title=NA,cexRow=2,
                    key.xlab="Gene specificity score (PCHi-C data)",
                    key.ytickfun=function(){return(list(labels=FALSE, tick=FALSE))},key.ylab=NA)
dev.off()
rowOrder <- rev(hmap$rowInd) # get row (cell type) order from heatmap
  
# summary table containing gene(bait) metadata, cluster IDs and specificity scores
idx <- sapply(geneIDsSorted,function(x) which(pksZero$baitID==x)[1])
geneData <- data.frame(pksZero[idx,c(4:5,1:3)],clusterID=as.factor(clusterIDsorted),t(geneSpecificityScoresIntnSorted[rowOrder,]))
geneData$baitName <- as.character(geneData$baitName)
rownames(geneData) <- NULL

geneSorted <- lapply(geneSorted,function(x) x[rowOrder,])
geneSpecificityMatricesSorted <- lapply(geneSpecificityMatricesSorted,function(x) x[rowOrder,])



##########
# load and process expression data
#########

# expression data file is available at https://osf.io/u8tzp/.
expnData <- read.table("GeneExpressionMatrix.txt",header=TRUE)

# currently log-transformed; change into asinh-transform
expnData[,2:39] <- asinh(exp(expnData[,2:39]))

# average values for each cell type
cn <- colnames(expnData)
expnData <- data.frame(ENSEMBL_GENEID=expnData$ENSEMBL_GENEID,
                       Mon=apply(expnData[,grep("Mon",cn)],1,mean),
                       Mac0=apply(expnData[,grep("Mac0",cn)],1,mean),
                       Mac1=apply(expnData[,grep("Mac1",cn)],1,mean),
                       Mac2=apply(expnData[,grep("Mac2",cn)],1,mean),
                       Neu=apply(expnData[,grep("Neu",cn)],1,mean),
                       MK=apply(expnData[,grep("MK",cn)],1,mean),
                       nCD4=apply(expnData[,grep("nCD4",cn)],1,mean),
                       Ery=apply(expnData[,grep("Ery",cn)],1,mean),
                       baitID=expnData$BaitID)

# remove baits that map to >1 ENSEMBL_GENEID
expnData <- expnData[!(duplicated(expnData$baitID)|duplicated(expnData$baitID,fromLast=TRUE)),]

# keep only 7,004 baits used in PCHi-C clustering
# 3 baits are missing and take NA values
expnData <- expnData[match(geneData$baitID,expnData$baitID),]



##########
# calculate gene specificity scores (based on expression data)
#########

# calculate distances between cell types based on expression data
distMatrixExpn <- as.matrix(dist(t(expnData[,cellTypes])))

# calculate specificity scores
geneSpecificityScoresExpn <- apply(expnData[,cellTypes],1,specificityScore,distMatrixExpn)
rownames(geneSpecificityScoresExpn) <- cellTypes
colnames(geneSpecificityScoresExpn) <- expnData$baitID
geneSpecificityScoresExpn <- t(geneSpecificityScoresExpn)



##########
# Figure 4C/S4B - comparison of expression-based gene specificity scores and PCHi-C-based gene specificity scores
##########

pchicMean <- aggregate(as.matrix(geneData[,cellTypes]),list(geneData$clusterID),mean)[,cellTypes]
expnMean <- aggregate(as.matrix(geneSpecificityScoresExpn),list(geneData$clusterID),mean,na.rm=TRUE)[,cellTypes]
pchicSD <- aggregate(as.matrix(geneData[,cellTypes]),list(geneData$clusterID),sd)[,cellTypes]
expnSD <- aggregate(as.matrix(geneSpecificityScoresExpn),list(geneData$clusterID),sd,na.rm=TRUE)[,cellTypes]

pdf("Figure4C_S4B.pdf",width=24,height=12)
pt <- list()
for(i in 1:length(cellTypes)) {
  df <- data.frame(x=pchicMean[,i],y=expnMean[,i],xmin=pchicMean[,i]-pchicSD[,i],xmax=pchicMean[,i]+pchicSD[,i],ymin=expnMean[,i]-expnSD[,i],ymax=expnMean[,i]+expnSD[,i],clusterID=as.factor(1:12))
  pt[[i]] <- ggplot(data = df,aes(x = x,y = y)) +  
    geom_errorbar(aes(ymin = ymin,ymax = ymax,colour=clusterID)) + 
    geom_errorbarh(aes(xmin = xmin,xmax = xmax,colour=clusterID)) + 
    geom_point(aes(colour=clusterID),size=3) +
    scale_colour_manual(values=clusterCols[1:12]) + ggtitle(cellTypes[i]) + 
    labs(x="Mean gene specificity score (PCHi-C data)",y="Mean gene specificity score (expression data)") + theme_bw()
}
grid.arrange(pt[[1]],pt[[2]],pt[[3]],pt[[4]],pt[[5]],pt[[6]],pt[[7]],pt[[8]],nrow=2)
dev.off()



##########
# Figure 4D/S4C - heatmap of PCHi-C-based gene specificity scores for 100 genes most specific to a given cell type
##########

# pages 1 and 2 of generated pdf are Figures 4D (nCD4) and S4C (Mon) respectively. Remaining 6 pages are for the other cell types.
nBait <- 100
my_palette1 <- colorRampPalette(c("blue","white","red"))(n = 256)
pdf(file="Figure4D_S4C.pdf",height=5,width=15)
for(i in c(8,1,2,3,4,5,6,7)) {
  tmp <- geneSpecificityScoresExpn
  tmp <- tmp[order(tmp[,cellTypes[i]],decreasing=TRUE),]
  baitIDs <- unique(as.numeric(rownames(tmp)[1:nBait]))
  geneDataSubset <- geneData[geneData$baitID %in% baitIDs,]
  
  m <- max(abs(geneData[,cellTypes]),na.rm=TRUE)
  col.breaks <- seq(from=-m,to=m,length.out=257)
  clusterCols2 <- clusterCols[geneDataSubset$clusterID]
  heatmap.2(t(geneDataSubset[,cellTypes[c(8,3,4,2,1,6,7,5)]]),trace="none",scale="none",Rowv=NULL,Colv=NULL,labCol=NA,
            col=my_palette1,margins=c(5,10),ColSideColors=clusterCols2,breaks=col.breaks,tracecol=NA,key.title=NA,cexRow=2,
            key.xlab="Gene specificity score (PCHi-C data)",
            key.ytickfun=function(){return(list(labels=FALSE, tick=FALSE))},key.ylab=NA,
            xlab=paste0("Top ",nBait," ",cellTypes[i],"-specific genes (based on expression)"))
}
dev.off()



##########
# Figure 4E - heatmap of cluster enrichment scores for the top 100 cell type-specifically expressed genes
##########

nClust <- table(geneData$clusterID) # number of genes in each cluster
nClustSubset <- matrix(nrow=length(cellTypes),ncol=12)
enrichmentScores <- matrix(nrow=length(cellTypes),ncol=12)
nBait <- 100
for(i in 1:length(cellTypes)) {
  tmp <- geneSpecificityScoresExpn
  tmp <- tmp[order(tmp[,cellTypes[i]],decreasing=TRUE),]
  baitIDs <- unique(as.numeric(rownames(tmp)[1:nBait]))
  
  # number of genes in each cluster out of the 100 genes that are most specific, based on expression, to cell type cellTypes[i]
  nClustSubset[i,] <- table(geneData[geneData$baitID %in% baitIDs,]$clusterID)
  enrichmentScores[i,] <- (nClustSubset[i,]/nBait)-(nClust/sum(nClust))
}

pdf(file="Figure4E.pdf",height=5,width=8)
m <- max(abs(enrichmentScores),na.rm=TRUE)
col.breaks <- seq(from=-m,to=m,length.out=256)
my_palette1 <- colorRampPalette(c("#00008B","white","#FFB90F"))(n = 256)
heatmap.2(enrichmentScores[c(8,3,4,2,1,6,7,5),],trace="none",scale="none",main=NA,Rowv=NULL,Colv=NULL,labRow=cellTypes[c(8,3,4,2,1,6,7,5)],labCol=NA,col=my_palette1,
          margins=c(5,10),ColSideColors=clusterCols,density.info="none",key.title=NA,
          key.ytickfun=function(){return(list(labels=FALSE, tick=FALSE))},key.ylab=NA,
          key.xlab="Cluster enrichment for the top 100\n cell type-specifically expressed genes")
dev.off()