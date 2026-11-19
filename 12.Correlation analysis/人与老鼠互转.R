
#人转鼠     注意列名为Gene name
geneCard <- read.table("ID.txt",sep="\t",header=T,check.names=F) ##
colnames(geneCard)[1]="Gene name"
genehuman=read.table("mart_export2.txt",sep="\t",header=T,check.names=F)
genehuman=genehuman[genehuman$`Rat gene name` != "",]
geneCard=dplyr::inner_join(genehuman,geneCard,by="Gene name")

write.table(geneCard,"aa2.txt",sep="\t",row.names = F,quote = F)





# 鼠转人   注意列名为Rat gene name
geneCard <- read.table("aa1.txt",sep="\t",header=T,check.names=F) ##
#colnames(geneCard)[1]="Gene name"
genemouse=read.table("mart_export2.txt",sep="\t",header=T,check.names=F)
genemouse=genemouse[genemouse$`Gene name` != "",]
geneCard=dplyr::inner_join(genemouse,geneCard,by="Rat gene name")

write.table(geneCard,"aa2.txt",sep="\t",row.names = F,quote = F)