library(Seurat)
seuratObj<- subset(seuratObj, cells=colnames(seuratObj)[(seuratObj$tumorType %in% "WT")])

seuratObj$CD83 <-rep( "not-expressed", length(seuratObj$seurat_clusters))
seuratObj$CD83 [as.numeric(seuratObj@assays$RNA["Cd83",])>0] <- "expressed"

seuratObj$CD45 <-rep( "not-expressed", length(seuratObj$seurat_clusters))
seuratObj$CD45 [as.numeric(seuratObj@assays$RNA["Ptprc",])>0] <- "expressed"

seuratObj$markers <- rep( "", length(seuratObj$seurat_clusters))
seuratObj$markers[as.numeric(seuratObj@assays$RNA["Cd83",])>0] <- "CD83"
seuratObj$markers[as.numeric(seuratObj@assays$RNA["Ptprc",])>0] <- paste0(seuratObj$markers[as.numeric(seuratObj@assays$RNA["Ptprc",])>0], "_CD45")

seuratObj$markers[seuratObj$markers==""] <- "1"
seuratObj$markers[seuratObj$markers=="_CD45"]<- "2_CD45"
seuratObj$markers[seuratObj$markers=="CD83"]<- "4_CD83"
seuratObj$markers[seuratObj$markers=="CD83_CD45"]<- "3_CD83_CD45"
seuratObj$isTumorSimple <- "non_tumor"
seuratObj$isTumorSimple[as.numeric(seuratObj@assays$RNA["EGFP",])>0] <- "tumor"
seuratObj$isTumorSimple[as.numeric(seuratObj@assays$RNA["Pdgfra",])>thresholds[1]] <- "tumor"
seuratObj$isTumorSimple[as.numeric(seuratObj@assays$RNA["Sox2",])>thresholds[3]] <- "tumor"
seuratObj$isTumorSimple[as.numeric(seuratObj@assays$RNA["Egfr",])>thresholds[2]] <- "tumor"


seuratObj<- subset(seuratObj, cells=colnames(seuratObj)[which(seuratObj$isTumorSimple=="tumor")])

library(Seurat)
p1 <- DimPlot( seuratObj, group.by="markers")
df <- unique(data.frame(col=ggplot_build(p1)$data[1][[1]][,1], id= seuratObj$markers))
rownames(df) <- NULL
df$col[df$id=="2_CD45"] <- "green"
df$col[df$id=="4_CD83"] <- "red"
df$col[df$id=="3_CD83_CD45"] <- "yellow"
df$col[df$id=="1"] <- "gray"
df <- df[order(df$id),]


p<- FeaturePlot(seuratObj, split.by="tumorType",features="Cd83", cols=c("lightgrey", "red"), order = F, min.cutoff="q9")
