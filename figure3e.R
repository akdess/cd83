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
seuratObj$isTumorSimple <- "non_tumor"
seuratObj$isTumorSimple[which(tumors>=1)] <- "tumor"

seuratObj<- subset(seuratObj, cells=colnames(seuratObj)[which(seuratObj$isTumorSimple=="tumor")])
seuratObj<- subset(seuratObj, cells=colnames(seuratObj)[!(seuratObj$tumorType %in% "Normal")])


pan_genes <- read.csv("PanGenes.csv")

gene_sets_by_category <- split(pan_genes, as.character(pan_genes[, 1]))
filtered_gene_sets <- lapply(gene_sets_by_category, function(genes) intersect(genes$Gene, rownames(seuratObj)))
filtered_gene_sets <- filtered_gene_sets[which(sapply(filtered_gene_sets, length) > 0)]

cibersort_data <- openxlsx::read.xlsx("cibersort.xlsx", sheet = 1)

cibersort_data <- cibersort_data[-1, ]  # Remove the first row
colnames(cibersort_data) <- cibersort_data[1, ]  # Set column names from the first row

binary_gene_sets <- apply(cibersort_data, 2, function(column) rownames(cibersort_data)[which(column == 1)])

combined_gene_sets <- c(filtered_gene_sets, binary_gene_sets)
final_filtered_sets <- lapply(combined_gene_sets, function(genes) intersect(genes, rownames(seuratObj)))
final_filtered_sets <- final_filtered_sets[which(sapply(final_filtered_sets, length) > 0)]

# Load CD83 signature and append it to the final sets
load("cd83_signature.rda")
final_filtered_sets <- c(final_filtered_sets, cd83_signature = list(cd83_signature))

seuratObj<-Seurat::AddModuleScore(seuratObj, features = final_filtered_sets, name = "GeneSet")
metadata <- seuratObj@meta.data

ModuleScoreIndex<-which(colnames(metadata)=="GeneSet1"):dim(metadata)[2]
#Matrix of data
metadf<-metadata[,ModuleScoreIndex]
colnames(metadf)<-names(final_filtered_sets)
colnames(seuratObj@meta.data)[ModuleScoreIndex] <- names(final_filtered_sets)
p<- cor(metadf)

cor.pval <- apply(metadf, 1, function(x) cor.test(metadf[, 130], x)$p.val)

p<- pheatmap(cor(metadf), breaks=breaks_heatmap, col=rev(brewer.pal(length(breaks_heatmap)-1, "RdBu")))
