library(Seurat)
seuratObj<- subset(seuratObj, cells=colnames(seuratObj)[which(seuratObj$seurat_clusters %in% c(2,7, 11, 24, 26))])
seuratObj<- subset(seuratObj, cells=colnames(seuratObj)[!(seuratObj$tumorType %in% "Normal")])

load("cd83_features.rda")

DefaultAssay(seuratObj) <- "RNA"
data <- (seuratObj@assays$RNA@counts)
samples  <- names(which(table(seuratObj$orig.ident)>10))
data2 <- lapply(1:length(samples), function(x) apply(data[, which(seuratObj$orig.ident==samples[x]), drop = FALSE], 1, sum))
raw.data <- do.call(cbind, data2)
colnames(raw.data) <- samples
raw.data<- t(apply(raw.data, 1, function(x) as.numeric(x)))
colnames(raw.data) <- samples


data2 <- CreateSeuratObject(counts = raw.data)
data2 <- NormalizeData(data2, verbose = FALSE)
itgam <- data2@assays$RNA$data[ "ITGAM", rownames(dat)]

to.plot <- data.frame(ITGAM=itgam,cd83_score=dat$set, location=dat$location,tumorType=dat$tumorType, name=dat$name )
p<- ggplot(to.plot, aes(x=cd83_score, y=ITGAM, label=name))+ geom_point() +  geom_text_repel() + theme_cowplot()+stat_cor(method="pearson") +   geom_smooth(method="lm",se=FALSE)  +ylab("ITGAM expression in macrophage cells") + xlab("Cd83 gene set scoring in tumor cells")
ggsave(paste0("ITGAM_cd83_score.pdf"),p, width=10, height=5)
