

library(limma)
library(pheatmap)
library(gplots)



plotdata=read.table("heatmap.txt",sep="\t",header=T,check.names=F,row.names = 1)
up=read.table("up.xls",header=T,stringsAsFactors = F,sep = "\t")

# rownames(up) <- up[, 1]
# up<-up[,-1]

down=read.table("down.xls",header=T,stringsAsFactors = F,sep = "\t")

annCol <- data.frame(Group = rep(c("Control","Model"),c(3,3)),
                     row.names = colnames(plotdata),
                     stringsAsFactors = F)

annRow <- data.frame(Direct = rep(c("Up","Down"),c(length(up$X),length(down$X))),
                     row.names = c(up$X,down$X),
                     stringsAsFactors = F)

annColors <- list("Group"=c("Control"="blue",
                            "Model"="red"),
                  "Direct"=c("Up"="yellow",
                             "Down"="green"))

plotdata <- plotdata[c(up$X,down$X),]
plotdata <- t(scale(t(plotdata)))
plotdata[plotdata > 5] <- 5
plotdata[plotdata < -5] <- -5

pdf(file="2.Heatmap.pdf",width = 8,height = 8)
pheatmap(plotdata, 
         scale = "none",
         annotation_row=annRow, 
         annotation_col=annCol, 
         annotation_colors = annColors,
         color = colorRampPalette(c("navy", "white", "palevioletred1"))(20),
         #color = greenred(64),
         fontsize_row=7,
         fontsize_col=5,
         fontsize=9,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = F,
         show_rownames = F)
dev.off()
