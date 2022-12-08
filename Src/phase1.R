# AML MicroArray Analysis

setwd('~/Desktop/HW/introduction-to-bioinformatics/project/')

library(GEOquery)
library(limma)
library(umap)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)

gset <- getGEO("GSE48558", GSEMatrix =TRUE, getGPL=T, destdir='Data/')
gset <- gset[[1]]

# gsms <- paste0("00000000000001111221211211221111122121203120312010",
#                "01001001101101101141501131366415015501550415041507",
#                "22072207220211111111121111111222222223333333005000",
#                "44444445666635555555")
gsms <- paste0("0000000000000XXXXXXXXXXXXXXXXXXXXXXXXXXX1XXX1XXXXX",
               "XXXXXXXXXXXXXXXXXX2X3XXX1X1442X3XX33XX33X2X3X2X3X5",
               "XXX5XXX5XXXXXXXXXXXXXXXXXXXXXXXXXXXXX1111111003000",
               "22222223444413333333")
sml <- strsplit(gsms, split="")[[1]]

sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

gs <- factor(sml)
# groups <- make.names(c("AML","BALL","TALL","GRLC","B","T","MNC","CD34"))
groups <- make.names(c("AML","GRLC","B","T","MNC","CD34"))
levels(gs) <- groups
gset$group <- gs

ex <- exprs(gset)

# 1. do we need to normalize it?
# max(ex) is 13.76154 so it's already normalized
# if it wasn't:
# ex <- log2(ex + 1)
# exprs(gset) <- ex

# 2. quality control
# 2.1 box plot to see if the values are normalized or not
pdf('Results/boxplot.pdf', width = 170)
boxplot(ex) # you can rotate the x-axis labels to visualize it in a prettier way
dev.off()
# as we can see the data is normalized.
# if it wasn't:
# ex <- normalizeQuantiles(ex)
# exprs(gset) <- ex

# 2.2 correlation heatmap between each pair of samples
pdf('Results/CoreHeatmap2.pdf', width = 30, height = 30)
pheatmap(cor(ex),
         # labels_row = gs,
         # labels_col = gs,
         color = greenred(256),
         border_color = NA,)
dev.off()

#2.3 PCA
pc <- prcomp(ex)
pdf('Results/PC.pdf')
plot(pc)
plot(pc$x[, 1:2]) # x column are genes (so each point in this plot represents a gene)
dev.off()

ex.scale <- t(scale(t(ex), scale = F))
pc <- prcomp(ex.scale)
pdf('Results/PC_scaled2.pdf')
plot(pc)
plot(pc$x[, 1:2]) # x column are genes (so each point in this plot represents a gene)
dev.off()

# PCA samples
pcr <- data.frame(pc$rotation[, 1:3], Group=gs)
pdf('Results/PC_samples2.pdf')
ggplot(pcr, aes(x = PC1, y = PC2, color = Group)) + geom_point() + theme_dark()
dev.off()

# #3 find differential expression analysis
# design <- model.matrix(~group + 0, gset)
# colnames(design) <- levels(gs)
# 
# fit <- lmFit(gset, design) # fit linear model
# 
# # set up contrasts of interest and recalculate model coefficients
# cts <- paste(groups, c(tail(groups, -1), head(groups, 1)), sep="-")
# cont.matrix <- makeContrasts(AML-B, levels=design)
# fit2 <- contrasts.fit(fit, cont.matrix)
# 
# # compute statistics and table of top significant genes
# fit2 <- eBayes(fit2, 0.01)
# tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
# 
# tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
# write.table(tT, 'Results/', row.names=F, sep="\t", quote = F)
# 
# aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
# aml.up.genes <- unique(aml.up$ID)

# tSNE dimension reduction
# install.packages("Rtsne")
library(Rtsne)
library(gridExtra)
library(ggpubr)

tsne_results <- list(Rtsne(t(ex), perplexity=5, check_duplicates = FALSE),
                     Rtsne(t(ex), perplexity=10, check_duplicates = FALSE),
                     Rtsne(t(ex), perplexity=15, check_duplicates = FALSE),
                     Rtsne(t(ex), perplexity=20, check_duplicates = FALSE),
                     Rtsne(t(ex), perplexity=30, check_duplicates = FALSE),
                     Rtsne(t(ex), perplexity=40, check_duplicates = FALSE),
                     Rtsne(t(ex), perplexity=50, check_duplicates = FALSE))

plots.list <- list()
plot_data = function(data) {
  ggplot(data, aes(X1, X2, color = Group)) + geom_point() + theme_dark()
}
for(i in seq_along(tsne_results)) {
  tsne <- data.frame(tsne_results[[i]]$Y[, 1:2], Group=gs)
  plot.t <- plot_data(tsne)
  plots.list[[i]] <- plot.t
  # plot(tsne$X1, tsne$X2, col='blue', bg=tsne$Group)
}

plots <- ggarrange(plotlist = plots.list, 
          labels = c(5, 10, 15, 20, 30, 40, 50),
          ncol = 1,
          nrow = 2)
ggexport(plots, filename = 'Results/tSNE_samples.pdf')
# plot(tsne_results$Y, col='blue', bg=gs)
