table(projHeme4$Clusters2)
trajectory <- c("H9_ES", "neural prog", "choroid plexus")

projHeme4 <- addTrajectory(
  ArchRProj = projHeme4, 
  name = "neural", 
  groupBy = "Clusters2",
  trajectory = trajectory, 
  embedding = "UMAP", 
  force = TRUE)

head(projHeme4$neural[!is.na(projHeme4$neural)])

plotTrajectory(projHeme4, 
               trajectory = "neural", 
               colorBy = "cellColData", size = 0.4,
               name = "neural",plotAs='points',continuousSet = 'whiteBlue')


p1 <- plotTrajectory(projHeme4, 
                     trajectory = "neural", 
                     colorBy = "GeneScoreMatrix", baseSize =10,plotAs='points',
                     name = "TTR", 
                     continuousSet = "horizonExtra")
p1

p2 <- plotTrajectory(projHeme4, 
                     trajectory = "neural", threads=1,
                     colorBy = "GeneIntegrationMatrix", 
                     name = "TTR", plotAs='points',baseSize =10,
                     continuousSet = "blueYellow")
p2


# heatmap
trajMM  <- getTrajectory(ArchRProj = projHeme4, 
                         name = "neural", 
                         useMatrix = "MotifMatrix", 
                         log2Norm = FALSE)

plotTrajectoryHeatmap(trajMM, labelTop =50,limits =c(-0.5,2),
                      pal = paletteContinuous(set = "whiteBlue")) #labelTop = 50,varCutOff = 0.9
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))

# GeneScoreMatrix
trajGSM <- getTrajectory(ArchRProj = projHeme4, name = "neural", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- trajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))

# GeneIntegrationMatrix
trajGIM <- getTrajectory(ArchRProj = projHeme4, name = "neural", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))

p1+p2+p3

# PeakMatrix
trajPM  <- getTrajectory(ArchRProj = projHeme4, name = "neural", useMatrix = "PeakMatrix", log2Norm = TRUE)
plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))

corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]

trajCombined <- trajGIM2
assay(trajCombined,withDimnames=FALSE) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))

ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)

ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
