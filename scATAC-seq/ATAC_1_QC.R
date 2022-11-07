library(ArchR)
addArchRThreads(threads = 30)
addArchRGenome("hg19")

inputFiles <- c("col13_human_fragment.bed.gz", "col14_human_fragment.bed.gz",
                "col25_human_fragment.bed.gz", "col26_human_fragment.bed.gz",
                "col28_human_fragment.bed.gz", "col29_human_fragment.bed.gz",
                "col33_human_fragment.bed.gz", "col34_human_fragment.bed.gz",
                "H9_ESC_fragment.bed.gz")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = c('col13','col14','col25','col26','col28','col29','col33','col34','H9_ESC'),
  minTSS = 4, 
  minFrags = 1000, 
  threads =1, 
  TileMatParams = list(tileSize = 2000),
  addTileMat = TRUE,
  addGeneScoreMat = TRUE)

# Doublet
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", 
  LSIMethod = 1)

# ArchRProject
projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE)
projHeme1

meta_atac = read.table('./scATAC_metadata.txt.gz',sep = '\t',comment.char = '')
ProjHeme1@cellColData@listData = meta_atac

# QC plot
summary(projHeme2$nFrags)
summary(projHeme2$TSSEnrichment)
df <- getCellColData(projHeme2, select = c("log10(nFrags)", "TSSEnrichment"))
pdf(file = 'QC.pdf',width = 6,height = 6)
ggPoint(
  x = df[,1],
  y = df[,2],rastr=T,
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(2.8,4.2),
  ylim = c(0, 38))+
  geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
dev.off()

# filter doublets
projHeme2 <- filterDoublets(projHeme1)

# save
saveArchRProject(ArchRProj = projHeme2, outputDirectory = "Save-ProjHeme2", load = F)
