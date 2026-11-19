

logFoldChange=0.25
adjustP=0.05

library(limma)
rt=read.table("../../data/RNA/symbol.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
rt=normalizeBetweenArrays(as.matrix(rt))

rt=log2(rt+1)
rt <- rt[,c(4:9)]

#differential
modType=c(rep("Model",3),rep("ICA",3))           ####### replicate
design <- model.matrix(~0+factor(modType))
colnames(design) <- c("ICA","Model")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(ICA-Model,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)

#write table
diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange &P.Value < adjustP )), ]
write.table(diffSig,file="diff.xls",sep="\t",quote=F)
diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & P.Value < adjustP )), ]
write.table(diffUp,file="up.xls",sep="\t",quote=F)
diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & P.Value < adjustP )), ]
write.table(diffDown,file="down.xls",sep="\t",quote=F)

#write expression level of diff gene
hmExp=rt[rownames(diffSig),]
diffExp=rbind(id=colnames(hmExp),hmExp)
write.table(diffExp,file="heatmap.txt",sep="\t",quote=F,col.names=F)

# tiff(file="vol.tiff",
#      width = 12,            #ͼƬ?Ŀ???
#      height =12,            #ͼƬ?ĸ߶?
#      units ="cm",
#      compression="lzw",
#      bg="white",
#      res=600)
pdf("1.Volcano.pdf")
xMax=max(abs(allDiff$logFC))
yMax=max(-log10(allDiff$P.Value))
nonZeroPValues <- allDiff$P.Value[allDiff$P.Value > 0]
if (length(nonZeroPValues) > 0) {
  yMax <- max(-log10(nonZeroPValues))
} else {
  yMax <- 0
}
plot(allDiff$logFC, -log10(allDiff$P.Value), xlab="log2FC",ylab="-log10(P.Value)",
     main="Volcano",xlim=c(-xMax,xMax),ylim=c(0,yMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(allDiff, P.Value<adjustP & logFC>logFoldChange)
points(diffSub$logFC, -log10(diffSub$P.Value), pch=20, col="palevioletred1",cex=0.8)
diffSub=subset(allDiff, P.Value<adjustP & logFC<(-logFoldChange))
points(diffSub$logFC, -log10(diffSub$P.Value), pch=20, col="dodgerblue",cex=0.8)
abline(v=0,lty=2,lwd=3)
dev.off()
