library("RColorBrewer")
library("gplots")
library("pheatmap")

args<-commandArgs(trailingOnly = TRUE)
input<-args[1]

d<-read.table(input,header=T,stringsAsFactors=FALSE)
data <- d[,seq(9,ncol(d),1)]
rownames(data)<-d[,1]
rownames(data)
anno <- data.frame(d[,1])
rownames(anno) <- rownames(data)
data.m <- data.matrix(data)
data.n <- as.numeric(data.m)
vmax <- max(data.n)
vmin <- min(data.n)

pdf("matrix/figure/Sites.heatmap1.pdf",30,30)

pheatmap(data,color = colorRampPalette(c("#daeffc","#FFFFFF","#FFCCCC","#FF9999","#FF6666","#FF3300"))(50),cluster_cols=T,legend_breaks=vmax:vmin,cluster_rows=T,cellwidth = 30, cellheight = 1,show_rownames = F,show_colnames = T,border_color="grey",fontsize = 25)

dev.off()

pdf("matrix/figure/Sites.heatmap2.pdf",30,30)

pheatmap(data,color = colorRampPalette(c("#daeffc","#FFFFFF","#FFCCCC","#FF9999","#FF6666","#FF3300"))(50),cluster_cols=F,legend_breaks=vmax:vmin,cluster_rows=T,cellwidth = 30, cellheight = 1,show_rownames = F,show_colnames = T,border_color="grey",fontsize = 25)

dev.off()
