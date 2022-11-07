library(Seurat)
count = read.table('*.txt.gz',sep = '\t')
meta = read.table('./metadata_human_total.txt',sep = '\t')
seob = CreateSeuratObject(counts = count,meta.data = meta)

library(stringr)
str_subset(rownames(seob),'MT-')
seob[["percent.mt"]] <- PercentageFeatureSet(seob, pattern = "^MT-")

seob <- NormalizeData(seob, normalization.method = "LogNormalize")
seob <- CellCycleScoring(seob, s.features = cc.genes$s.genes, 
                         g2m.features = cc.genes$g2m.genes)

seob_list <- SplitObject(seob, split.by = "sample")

for(i in 1:length(seob_list)){
  seob <- seob_list[[i]]
  seob <- NormalizeData(seob, normalization.method = "LogNormalize")
  seob_list[[i]] <- seob 
}

seob <- merge(x = seob_list[[1]],
              y = seob_list[-1], 
              add.cell.ids = names(seob_list))

seob <- FindVariableFeatures(seob, selection.method = "vst",nfeatures = 3000)
seob <- ScaleData(seob, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))

library(Seurat)
library(SeuratDisk)
total <-CreateSeuratObject(counts = seob@assays[["RNA"]]@scale.data,meta.data = seob@meta.data) 
SaveH5Seurat(total,filename = 'total.h5Seurat')
Convert('total.h5Seurat',dest = 'h5ad')

# Considering of the memory and time, scanpy was used to run PCA, Cluster and UMAP.

leiden <- read.csv('./leiden_35_1.6.csv',sep=',',row.names = 1)
seob$leiden <- leiden
Idents(seob) <- 'leiden'
seob$leiden <- factor(seob$leiden,levels = c(as.character(0:28)))

markers <- FindAllMarkers(seob, only.pos = TRUE,logfc.threshold = 0.25)
markers <- markers[order(markers$cluster, -markers$avg_log2FC, markers$p_val),]
write.table(markers,file='markers.txt',sep='\t',quote=F)

