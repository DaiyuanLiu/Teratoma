# motif enrichment
projHeme4 <- addMotifAnnotations(ArchRProj = projHeme4, motifSet = "cisbp", name = "Motif")
motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projHeme4,
  peakAnnotation = "Motif",
  cutOff = "Pval <= 0.05 & Log2FC >= 0.5")
motifsUp

# ChromVAR deviation
if("Motif" %ni% names(projHeme4@peakAnnotation)){
  projHeme4 <- addMotifAnnotations(ArchRProj = projHeme4, motifSet = "cisbp", name = "Motif")}

projHeme4 <- addBgdPeaks(projHeme4)

projHeme4 <- addDeviationsMatrix(
  ArchRProj = projHeme4, 
  peakAnnotation = "Motif",
  force = TRUE)

plotVarDev <- getVarDeviations(projHeme4, name = "MotifMatrix", plot = TRUE)
plotVarDev
save(plotVarDev,file = 'plotVarDev.rda')

motifs <- c("CDX2",'CDX1')
markerMotifs <- getFeatures(projHeme4, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs

projHeme4 <- addImputeWeights(projHeme4)
p <- plotGroups(ArchRProj = projHeme4, 
                groupBy = "Clusters2", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,alpha=0.75,
                imputeWeights = getImputeWeights(projHeme4))

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 10) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 10) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))


# footprint
motifPositions <- getPositions(projHeme4)
motifPositions

motifs <- c("HNF1A", "HNF1B", "HNF4A", "HNF4G")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs

seFoot <- getFootprints(
  ArchRProj = projHeme4, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters2")

# Subtract Tn5 Bias
plotFootprints(
  seFoot = seFoot,
  ArchRProj = projHeme4, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias_HNF",
  addDOC = FALSE,
  smoothWindow = 5,
  baseSize = 10)

# TSS footprint
seTSS <- getFootprints(
  ArchRProj = projHeme4, 
  positions = GRangesList(TSS = getTSS(projHeme4)), 
  groupBy = "Clusters2",
  flank = 2000)

plotFootprints(
  seFoot = seTSS,
  ArchRProj = projHeme4, 
  normMethod = "None",
  plotName = "TSS-No-Normalization",
  addDOC = FALSE,
  flank = 2000,
  flankNorm = 100,
  baseSize = 10)

# positive TF
seGroupMotif <- getGroupSE(ArchRProj = projHeme4, useMatrix = "MotifMatrix", groupBy = "Clusters2")
seGroupMotif

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

# GeneIntegrationMatrix
corGIM_MM <- correlateMatrices(
  ArchRProj = projHeme4,
  useMatrix1 = "GeneIntegrationMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "IterativeLSI")

corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])

ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05))
a = as.data.frame(corGIM_MM)

saveArchRProject(ArchRProj = projHeme4, outputDirectory = "Save-ProjHeme4", load = F)
