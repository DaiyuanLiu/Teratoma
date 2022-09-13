#--------------1.bin--------------
projHeme2 <- addIterativeLSI(
  ArchRProj = projHeme2,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 1, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = 1, 
    sampleCells = 20000, 
    n.start = 10), 
  varFeatures = 200000, 
  dimsToUse = 1:50,force = TRUE)


projHeme2 <- addClusters(
  input = projHeme2,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.6,force = TRUE)

# UMAP
projHeme2 <- addUMAP(
  ArchRProj = projHeme2, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.3, 
  metric = "cosine",
  seed =18,force = TRUE)

plotEmbedding(ArchRProj = projHeme2, 
              colorBy = "cellColData", labelSize = 5,size = 0.5,
              name = "Clusters",baseSize = 15,legendSize = 8,
              embedding = "UMAP")

#---------------peak-------------
projHeme2 <- addGroupCoverages(ArchRProj = projHeme2, groupBy = "Clusters",threads = 1)
pathToMacs2 <- '/home/ggj/anaconda3/bin/macs2'

projHeme2 <- addReproduciblePeakSet(
  ArchRProj = projHeme2, 
  groupBy = "Clusters", 
  pathToMacs2 = pathToMacs2)

projHeme2 <- addPeakMatrix(projHeme2)
saveArchRProject(ArchRProj = projHeme2, outputDirectory = "Save-ProjHeme2", load = F)

getAvailableMatrices(projHeme2)

projHeme3 <- addIterativeLSI(
  ArchRProj = projHeme2,
  useMatrix = "PeakMatrix", 
  name = "IterativeLSI", 
  iterations = 1, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = 1, 
    sampleCells = 20000, 
    n.start = 10), 
  varFeatures = 150000, 
  dimsToUse = 1:30,force = TRUE) 

projHeme3 <- addClusters( 
  input = projHeme3,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,force = TRUE)

# UMAP
projHeme3 <- addUMAP(
  ArchRProj = projHeme3, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 10, 
  minDist = 0.3, 
  metric = "cosine",
  seed =111,force = TRUE)

plotEmbedding(ArchRProj = projHeme3, 
              colorBy = "cellColData", labelSize = 5,size = 0.8,
              name = "Clusters",baseSize = 15,legendSize = 8,
              embedding = "UMAP")

saveArchRProject(ArchRProj = projHeme3, outputDirectory = "Save-ProjHeme3", load = FALSE)

# marker gene
markersGS2 <- getMarkerFeatures(
  ArchRProj = projHeme3, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")

# markerList
markerList2 <- getMarkers(markersGS2, cutOff = "FDR <= 0.5 & Log2FC >= 0.25") %>% as.data.frame()

#--------transfer----------
library(Seurat)
load("dge_integrate.rda")
meta_integrate = read.table('meta_integrate.txt',sep = '\t')
seob = CreateSeuratObject(counts = dge_integrate,meta.data = meta_integrate)
seRNA =  as.SingleCellExperiment(seob) 

table(seob$cell_type)

projHeme3 <- addGeneIntegrationMatrix(
  ArchRProj = projHeme3, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "cell_type",
  nameCell = "predictedCell_Un", 
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un") 

cM <- as.matrix(confusionMatrix(projHeme3$Clusters, projHeme3$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM))
preClust

unique(unique(projHeme3$predictedGroup_Un))

cgut <- paste0(c('mid/hind gut','H9_ESC'),collapse = "|")
cgut
cNongut <- paste0(c('foregut',"neural prog","choroid plexus",'mesen','fibroblast_DLK1'),collapse = "|")
cNongut

clustgut <- rownames(cM)[grep(cgut, preClust)]
clustgut
clustNongut <- rownames(cM)[grep(cNongut, preClust)]
clustNongut

rnagut <- colnames(seRNA)[grep(cgut, colData(seRNA)$cell_type)]
head(rnagut)
rnaNongut <- colnames(seRNA)[grep(cNongut, colData(seRNA)$cell_type)]
head(rnaNongut)

groupList <- SimpleList(
  gut = SimpleList(
    ATAC = projHeme3$cellNames[projHeme3$Clusters %in% clustgut],
    RNA = rnagut),
  Nongut = SimpleList(
    ATAC = projHeme3$cellNames[projHeme3$Clusters %in% clustNongut],
    RNA = rnaNongut))

projHeme3 <- addGeneIntegrationMatrix(
  ArchRProj = projHeme3, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE, #不加到Arrow文件
  groupList = groupList,
  groupRNA = "cell_type",
  nameCell = "predictedCell_Co",
  nameGroup = "predictedGroup_Co",
  nameScore = "predictedScore_Co")

pal <- paletteDiscrete(values = colData(seRNA)$cell_type)
pal

p1 <- plotEmbedding(
  projHeme3, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un", 
  pal = pal,size =0.01,alpha=0.5)

p2 <- plotEmbedding(
  projHeme3, 
  colorBy = "cellColData", 
  name = "predictedGroup_Co", 
  pal = pal,size =0.01,alpha=0.5)
p1+p2
table(projHeme2$predictedGroup_Co)

cM <- confusionMatrix(projHeme3$Clusters, projHeme3$predictedGroup_Co)
labelOld <- rownames(cM)
labelOld

labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew

remapClust <- c(
  "neural prog" = "C6",
  "mesen" = "C7",
  "foregut" = "C2",
  "mid/hind gut" = "C3",
  "choroid plexus" = "C4",
  "H9_ESC" = "C1",
  "fibroblast_DLK1" = "C5")
remapClust <- remapClust[names(remapClust) %in% labelNew]

labelNew2 <- mapLabels(labelNew, oldLabels = remapClust, newLabels = names(remapClust))
labelNew2

projHeme3$Clusters2 <- mapLabels(projHeme3$Clusters, newLabels = labelNew2, oldLabels = labelOld)

pdf(file = 'umap_transfer.pdf',width = 5,height = 5)
plotEmbedding(projHeme3, colorBy = "cellColData", name = "Clusters2",rastr=F)
dev.off()

projHeme3 <- addGeneIntegrationMatrix(
  ArchRProj = projHeme3, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = TRUE,
  force= TRUE,
  groupList = groupList,
  groupRNA = "cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore")

saveArchRProject(ArchRProj = projHeme3, outputDirectory = "Save-ProjHeme3", load = FALSE)

# marker gene
getAvailableMatrices(projHeme2)

projHeme3 <- addImputeWeights(projHeme3)

markersGS2 <- getMarkerFeatures(
  ArchRProj = projHeme3, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")

# markerList
markerList2 <- getMarkers(markersGS2, cutOff = "FDR <= 0.5 & Log2FC >= 0.25") %>% as.data.frame()
save(markersGS2,markerList2,file = 'marker.rda')
write.table(markerList2,file='markers.txt',sep='\t',quote=F)

# marker plot
library(ArchR)
markerGenes = c('POU5F1','CDX2')
p <- plotEmbedding(
  ArchRProj = projHeme3, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
 imputeWeights = getImputeWeights(projHeme3))

library(RColorBrewer)
p2 <- lapply(p, function(x){
  x + guides(color = NULL, fill = FALSE,raster = TRUE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank())})
do.call(cowplot::plot_grid, c(list(ncol = 7),p2))

#--------------peak calling based on Clusters2  ----------------
projHeme4 <- addGroupCoverages(ArchRProj = projHeme3, force=T,
                               groupBy = "Clusters2",threads = 1)
pathToMacs2 = findMacs2()
pathToMacs2 <- '/home/ggj/anaconda3/bin/macs2'


projHeme4 <- addReproduciblePeakSet(
  ArchRProj = projHeme4, 
  groupBy = "Clusters2", 
  pathToMacs2 = pathToMacs2)

projHeme4 <- addPeakMatrix(projHeme4)
saveArchRProject(ArchRProj = projHeme4, outputDirectory = "Save-ProjHeme4", load = F)


# peak_ratio 
library(ggsci)
projHeme4 = `Save-ArchR-Project`
peak = projHeme4@peakSet@elementMetadata@listData[["peakType"]]
table(peak)
peak_ratio = table(peak) %>% as.data.frame()
peak_ratio$proportion = (peak_ratio$Freq/sum(peak_ratio$Freq)) * 100
peak_ratio$peak = factor(peak_ratio$peak,levels = c('Intronic','Distal',
                                                    'Promoter','Exonic'))
colnames(peak_ratio)
ggplot(peak_ratio,aes(peak,proportion,fill = peak))+
  geom_bar(width = 0.3,stat = 'identity')+
  coord_polar(theta = 'y')+
  geom_text(aes(label = paste0(round(proportion,2),'% (',Freq,")")),
            position = position_stack(vjust = .5))+
  scale_fill_brewer(palette = 'Set3')+
  theme_void()+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

