
library(cowplot)
library(CellChat)

load("human_HIGHvsLOW.cellchat.rda")
object.list <- list(LOW = cellchat_control,HIGH = cellchat_cKO)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
library(cowplot)
gg1 <- plot_grid(
compareInteractions(cellchat, show.legend = T, group =  rep(factor(c("LOW","HIGH"), levels = c("LOW","HIGH")), 1)),
compareInteractions(cellchat, show.legend = T, group = rep(factor(c("LOW","HIGH"), levels = c("LOW","HIGH")), 1), measure = "weight"))


project <- "human_HIGHvsLOW"

pdf(paste0(project, "_comparisonNumInteractions.pdf"),height = 5, width = 10)
gg1
dev.off()


gg <- list()
for (i in 1) {
  gg[[i]] <- rankNet(cellchat, mode = "comparison", stacked = T, comparison = c(1,2), do.stat = T) + theme(legend.position="top")
}
g <- patchwork::wrap_plots(plots = gg)
g
cowplot::save_plot(filename=paste0(project, "_comparison_contributions_signalingPathways.pdf"), plot=g, base_width = 3, base_height = 3)


load("mouse_gofvsWT.cellchat.rda")
project <- "mouse_gofvsWT"

object.list <- list(WT = cellchat_control,GOF = cellchat_cKO)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

library(cowplot)
gg1 <- plot_grid(
compareInteractions(cellchat, show.legend = T, group =  rep(factor(c("WT","GOF"), levels = c("WT","GOF")), 1)),
compareInteractions(cellchat, show.legend = T, group = rep(factor(c("WT","GOF"), levels = c("WT","GOF")), 1), measure = "weight"))

pdf(paste0(project, "_comparisonNumInteractions.pdf"),height = 5, width = 10)

gg1

dev.off()


gg <- list()
for (i in 1) {
  gg[[i]] <- rankNet(cellchat, mode = "comparison", stacked = T, comparison = c(1,2), do.stat = T) + theme(legend.position="top")
}
g <- patchwork::wrap_plots(plots = gg)
g
cowplot::save_plot(filename=paste0(project, "_comparison_contributions_signalingPathways.pdf"), plot=g, base_width = 3, base_height = 3)

