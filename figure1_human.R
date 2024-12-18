library(Seurat)
seuratObj$CD83 <-rep( "not-expressed", length(seuratObj$seurat_clusters))
seuratObj$CD83 [as.numeric(seuratObj@assays$RNA["CD83",])>0] <- "expressed"

seuratObj$CD45 <-rep( "not-expressed", length(seuratObj$seurat_clusters))
seuratObj$CD45 [as.numeric(seuratObj@assays$RNA["PTPRC",])>0] <- "expressed"

seuratObj$markers <- rep( "", length(seuratObj$seurat_clusters))
seuratObj$markers[as.numeric(seuratObj@assays$RNA["CD83",])>0] <- "CD83"
seuratObj$markers[as.numeric(seuratObj@assays$RNA["PTPRC",])>0] <- paste0(seuratObj$markers[as.numeric(seuratObj@assays$RNA["PTPRC",])>0], "_CD45")

seuratObj$markers[seuratObj$markers==""] <- "1"
seuratObj$markers[seuratObj$markers=="_CD45"]<- "2_CD45"
seuratObj$markers[seuratObj$markers=="CD83"]<- "4_CD83"
seuratObj$markers[seuratObj$markers=="CD83_CD45"]<- "3_CD83_CD45"

seuratObj<- subset(seuratObj, cells=colnames(seuratObj)[which(seuratObj$isTumorSimple=="tumor")])

seuratObj<- subset(seuratObj, cells=colnames(seuratObj)[!(seuratObj$tumorType %in% "Normal")])

p1 <- DimPlot( seuratObj, group.by="markers")
df <- unique(data.frame(col=ggplot_build(p1)$data[1][[1]][,1], id= seuratObj$markers))
rownames(df) <- NULL
df$col[df$id=="2_CD45"] <- "green"
df$col[df$id=="4_CD83"] <- "red"
df$col[df$id=="3_CD83_CD45"] <- "yellow"
df$col[df$id=="1"] <- "gray"
df <- df[order(df$id),]

p1 <- DimPlot( seuratObj, group.by="markers",   order=F, raster=F, cols=df$col, split.by = "tumorType")+ggtitle('IDH Mutant+WT')+ theme(legend.position="none") 
ggsave(paste0("DimPlot_CD45_CD83_splitted_bytumortype_TUMORCELLSONLY_orderF.pdf"), p1 , width = 10, height = 5)
