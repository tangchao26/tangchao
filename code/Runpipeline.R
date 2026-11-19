source("SeuratV4CellAnnotateFunc.R")
source("KeyGenePipeline.R")
# 细胞注释 --------------------------------------------------------------------
samples = list.files("../data/scRNA/GSE203499_RAW/")
samples
dir <- file.path('../data/scRNA/GSE203499_RAW/', samples)
names(dir) <- samples

Object_list = list()

for (sample in dir){
  filedir = sample
  project <- stringr::str_split(filedir,"/",simplify = T)[, 6] # 设置每个样本的名字
  message(project)
  scrna_data <- Read10X(filedir)
  Seurat_object <- CreateSeuratObject(
    counts = scrna_data,
    project = project,
    min.cells = 3,
    min.features = 500)
  Object_list[[sample]] = Seurat_object
}

scRNA <- merge(Object_list[[1]],Object_list[2:length(Object_list)])
scRNA$group <- ifelse(substr(scRNA$orig.ident,9,10) %in% c(62,64,70,72),"Disease","Control")
table(scRNA$group,scRNA$orig.ident)

# nFeature_RNA >= 500 & percent.mt <= 35.66733 & nFeature_RNA <= 4099.425 & nCount_RNA <= 9924.445


load("refdata/Homo/DatabaseImmuneCell.se.Rdata")
load("refdata/Homo/BlueprintEncode.Rdata")
load("refdata/Homo/MonacoImmune.Rdata")
load("refdata/Homo/HumanPrimaryCellAtlas.Rdata")
load("refdata/Homo/NovershternHematopoietic.Rdata")

# scRNA <- readRDS("../1.QCFilter/doubletFinder_V4.rds")
# scRNA <- merge(scRNA[[1]], scRNA[2:length(scRNA)])
# scRNA <- readRDS("../.Cluster/data.merge.harmony.2000.rds")
scRNA <- scRNAAutoAnno(
  Path = "../",
  SeuratObject = scRNA,
  PC = 13,
  Multi = TRUE,
  ref = list(
    BlueprintEncode,
    DatabaseImmuneCell.se,
    NovershternHematopoietic,
    HumanPrimaryCellAtlas,
    MonacoImmune
  ),
  labels = list(
    BlueprintEncode$label.main,
    DatabaseImmuneCell.se$label.main,
    NovershternHematopoietic$label.main,
    HumanPrimaryCellAtlas$label.main,
    MonacoImmune$label.main
  ),
  Idents = "RNA_snn_res.0.8"
)
# Idents(scRNA) <- "RNA_snn_res.0.2"
# scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 
# write.csv(scRNA.markers,"../.Cluster/Findall.markers.csv")

# 输入细胞
# 0:Müller glia; 1:Bipolar cells; 2:Rod photoreceptors; 3:Bipolar cells; 4:Müller glia; 5:Müller glia; 6:Rod photoreceptors; 7:Bipolar cells; 8:Rod photoreceptors; 9:Rod photoreceptors; 10:Bipolar cells; 11:Bipolar cells; 12:Rod photoreceptors; 13:Microglial cells; 14:T cells; 15:Bipolar cells; 16:Müller glia; 17:Bipolar cells; 18:Endothelial cells; 19:Bipolar cells; 20:Horizontal cell; 21:Amacrine cells; 22:Bipolar cells; 23:Astrocytes; 24:Bipolar cells
# 输入基因
# EGF,RLBP1,PCP2,OTX2,PDE6A,RHO,HLA-DRA,C1QB,CD3D,CD3E,CFH,ONECUT2,ONECUT1,GAD1,GFAP,GLUL



# scPagwas----
library(scPagwas)

gwas_data <- read.table("../data/MR/finngen_R12_H7_GLAUCOMA.txt",header = T)
Single_data <- readRDS("../result.rds")
Pagwas_data<-scPagwas_main(Pagwas = NULL,
                           gwas_data =gwas_data, # The GWAS Summary statistics files 
                           Single_data =Single_data,# scRNA-seq data in seruat format with "RNA" assays and normalized.
                           output.prefix="BLCA", # the prefix name for output files
                           output.dirs="../4.scPagwas/",# the directory file's name for output
                           block_annotation = block_annotation,# gene position in chromosome is provided by package.
                           assay="RNA", # the assays for scRNA-seq data to use.
                           Pathway_list=Genes_by_pathway_kegg,# pathway list is provided by package, including gene symbols.
                           n.cores=1,
                           iters_singlecell = 10,
                           chrom_ld = chrom_ld,# The LD data is provided by package.
                           singlecell=T, # Whether to run the singlecell process.
                           celltype=T# Whether to run the celltype process.
)
save(Pagwas_data,file="../4.scPagwas/BLCA_scPagwas.RData")

require("RColorBrewer")
require("Seurat")
require("SeuratObject")
require("ggsci")
require("dplyr")
require("ggplot2")
require("ggpubr")
#check the objects
color26 <- c("#548CA8","#00A19D","#D54C4C","#8D2828","#FDD2BF","#E98580","#DF5E5E","#492F10","#334257","#476072","#D9DD6B","#ECEFA4","#ECD662","#5D8233","#284E78","#3E215D","#835151","#F08FC0","#C6B4CE","#BB8760","#FFDADA","#3C5186","#558776","#E99497","#FFBD9B","#0A1D37")
# my_pal2 <- c(
#   "#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
#   "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A","#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
#   "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A","#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
#   "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A","#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
#   "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A"
# )

pdf("../4.scPagwas/plot/Pagwas_umap_cellType_1.pdf",width = 10,height = 8)
Seurat::DimPlot(Pagwas_data,pt.size=1,reduction="umap",label = T, repel=TRUE)+
  scale_colour_manual(name = "cellType_1", values = color26)+
  umap_theme()+ggtitle("cellType_1")+labs(x="umap",y="")+theme(aspect.ratio=1)
dev.off()

scPagwas_Visualization(Single_data=Pagwas_data,
                       p_thre = 0.05,
                       FigureType = "umap",
                       width = 7,
                       height = 7,
                       lowColor = "white", 
                       highColor = "red",
                       output.dirs="../4.scPagwas/plot/",
                       size = 0.2,
                       do_plot = T)

pdf("../4.scPagwas/plot/Plot_bar_positie.pdf",width = 10,height = 8)
plot_bar_positie_nagtive(seurat_obj=Pagwas_data,
                         var_ident="cellType_1",#you should change this to you specific columns
                         var_group="positiveCells",
                         vec_group_colors=c("#E8D0B3","#7EB5A6"),
                         do_plot = T)
dev.off()

# 绘制正负单元格中单元类型比例的条形图
pdf("../4.scPagwas/plot/Plot_bar_positie_nagtive.pdf")
plot_bar_positie_nagtive(seurat_obj=Pagwas_data,
                         var_ident="positiveCells",
                         var_group="cellType_1", #you should change this to you specific columns
                         p_thre = 0.01,
                         vec_group_colors=NULL,
                         f_color=colorRampPalette(brewer.pal(n=10, name="RdYlBu")),
                         do_plot = T)
dev.off()

# 绘制遗传相关基因
pdf("../4.scPagwas/plot/Heritability_cor_scatterplot.pdf")
heritability_cor_scatterplot(gene_heri_cor=Pagwas_data@misc$PCC,
                             topn_genes_label=10,
                             color_low="#035397",
                             color_high ="#F32424",
                             color_mid = "white",
                             text_size=2,
                             do_plot=T,
                             max.overlaps =20,
                             width = 7,
                             height = 7)
dev.off()

# 细胞类型的条形图
Bootstrap_P_Barplot(p_results=Pagwas_data@misc$bootstrap_results$bp_value[-1],
                    p_names=rownames(Pagwas_data@misc$bootstrap_results)[-1],
                    figurenames = "../4.scPagwas/plot/Bootstrap_P_Barplot.pdf",
                    width = 5,
                    height = 7,
                    do_plot=T,
                    title = "BLCA")

# 用于估计值的森林图
Bootstrap_estimate_Plot(bootstrap_results=Pagwas_data@misc$bootstrap_results,
                        figurenames = "../4.scPagwas/plot/Estimateplot.pdf",
                        width = 9,
                        height = 7,
                        do_plot=T)


# irGSEA enrichment analysis-------------
dir <- "../5.irGSEA enrichment analysis/"
dir.create(dir)
library(UCell)
library(irGSEA)
scRNA <- readRDS("../result.rds")
scRNA <- irGSEA.score(object = scRNA, assay = "RNA",
                      slot = "data", seeds = 123, ncores = 1,
                      min.cells = 3, min.feature = 200,
                      custom = F, geneset = NULL, msigdb = T,
                      species = "Homo sapiens", category = "H",    #Mus musculus
                      subcategory = NULL, geneid = "symbol",
                      method = c("AUCell", "UCell", "singscore",
                                 "ssgsea"),
                      aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                      kcdf = 'Gaussian')
result.dge <- irGSEA.integrate(object = scRNA,
                               group.by = "cellType_1",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea"))

class(result.dge)
saveRDS(result.dge,"../5.irGSEA enrichment analysis/result.dge.rds")
set.seed(123)
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge,
                                      method = "RRA",
                                      top = 50,
                                      show.geneset = NULL)
irGSEA.heatmap.plot
ggplot2::ggsave(filename = paste0(dir,'1.IrGSEA.heatmap.pdf'),plot = irGSEA.heatmap.plot,height = 12,width = 10)
# Cellchat----
dir.create("../6.Cellchat")
setwd("../6.Cellchat")
library(ComplexHeatmap)
library(SeuratData)
library(ggalluvial)
library(patchwork)
library(tidyverse)
library(CellChat)
library(ggplot2)
library(svglite)
library(future)
library(Seurat)
library(NMF)
options(stringsAsFactors = FALSE)
library(Seurat)
library(magrittr)
library(openxlsx)
library(GSVA)
library(ggplot2)

scRNA<- readRDS("../result.rds")

# a <- subset(pbmc.plot2,cellType_1 == "Treg")
# pbmc.plot2$Scoring <- data.frame(row.names = colnames(pbmc.plot2),ifelse(pbmc.plot2$Scoring > median(a$Scoring) , "Hsco", "Lsco"))

# sub.NK <- subset(pbmc.plot2,subset = cellType_1 == "Treg")
# sub.NK$scissor <- ifelse(sub.NK$scissor != "Scissor+","other","Scissor+")
# 
# identical(colnames(subset(pbmc.plot2,subset = cellType_1 == "Treg")),colnames(sub.NK))
# Idents(pbmc.plot2, cells = colnames(sub.NK)) <- paste0(sub.NK$scissor,"_",sub.NK$cellType_1)
# pbmc.plot2$cellType_1 <- Idents(pbmc.plot2)
# table(pbmc.plot2$cellType_1 )

meta <- scRNA@meta.data[,"cellType_1",drop=F];colnames(meta) <- "cellType"
cellchat <- createCellChat(scRNA@assays$RNA@data, meta = meta, group.by = "cellType")
rm(scRNA)

CellChatDB <- CellChatDB.human
str(CellChatDB)
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor);head(CellChatDB$complex);head(CellChatDB$geneInfo)
pdf(file="CellChatDB.pdf",width = 10,height = 5)
showDatabaseCategory(CellChatDB)
dev.off()

unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 8)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

options(future.globals.maxSize= +Inf)
cellchat <- computeCommunProb(cellchat, seed.use = 123456)

cellchat <- filterCommunication(cellchat)
df.net <- subsetCommunication(cellchat)
write.table(df.net, file="net_lr.txt", quote=F, sep="\t", row.names=F)

cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat)
write.table(df.netp, file="net_pathway.txt", quote=F, sep="\t", row.names=F)

#####-----------show----------###
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat,"cellchat.rds")
cellchat <- readRDS("cellchat.rds")
groupSize <- as.numeric(table(cellchat@idents))

p1 <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                       label.edge = F, title.name = "Number of interactions")

p2 <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                       label.edge = F, title.name = "Interaction weigths/strength")
pdf(file="TIL.net_number_strength.pdf",width = 10,height = 5)
CombinePlots(plots = list(p1, p2))
dev.off()

# mat <- cellchat@net$count
# pdf(file="TIL.net_number_individual.pdf",width=8,height=8)
# par(mfrow=c(3,3), mar=c(1,1,1,1))
# for (i in 1:nrow(mat)) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i,] <- mat[i,]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
#                    arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }
# dev.off()
# 
# mat <- cellchat@net$weight
# pdf(file="TIL.net_strength_individual.pdf",width=8,height=8)
# par(mfrow=c(3,3), mar=c(1,1,1,1))
# for (i in 1:nrow(mat)) {
#   #i = 1
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i,] <- mat[i,]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
#                    arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }
# dev.off()
# 
# cellchat@netP$pathways
# (pathways.show <- cellchat@netP$pathways[1])
# levels(cellchat@idents)
# 
# pathways.show.all <- cellchat@netP$pathways
# levels(cellchat@idents)
# vertex.receiver = c(1,2,4,6)
# dir.create("all_pathways_com_circle")
# setwd("all_pathways_com_circle")
# for (i in 1:length(pathways.show.all)) {
#   netVisual(cellchat, signaling = pathways.show.all[i], out.format = c("pdf"),
#             vertex.receiver = vertex.receiver, layout = "circle")
#   gg <- netAnalysis_contribution(cellchat, pathways.show.all[i])
#   ggsave(filename = paste0(pathways.show.all[i],"_L-R_contribution.pdf"),
#          plot = gg, width = 5, height = 2.5, units = 'in', dpi = 300)
# }
# setwd("../")

levels(cellchat@idents)

pdf(file="Mye_Lymph_buble.pdf",width=10,height=8)
netVisual_bubble(cellchat, sources.use = c(2:9), targets.use = c(1), remove.isolate = TRUE)
dev.off()

pdf(file="Mye_Lymph_buble_2.pdf",width=10,height=8)
netVisual_bubble(cellchat, sources.use = c(7), targets.use = c(1:6,8,9), remove.isolate = TRUE)
dev.off()

scRNA <- subset(scRNA,cellType_1 %in% c("Müller glia","Horizontal cell"))
KeyGene <- c("VEGFA","VEGFR2","KDR")
p <- DotPlot(scRNA, group.by = "cellType_1", features = unique(KeyGene), cols = "RdYlBu", ) + scale_size_continuous(range = c(0, 10)) + theme(
  panel.border = element_rect(colour = "black"),
  axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5
  ),
  legend.position = "right",
  text = element_text(size = 10)
)

ggsave_fun("1.DotPlot.VEGFA", plot = p, width = 5 + ceiling(2 / 3), height = 5)

# pdf(file="Mye_Lymph_buble_KeyP.pdf",width=5,height=4)
# netVisual_bubble(cellchat, sources.use = c(8), targets.use = c(1),
#                  signaling = pathways.show.all[c(1)], remove.isolate = TRUE)
# dev.off()
# 
# pdf(file="Mye_Lymph_buble_KeyPC.pdf",width=6,height=5)
# pairLR.use <- extractEnrichedLR(cellchat, signaling = pathways.show.all)
# netVisual_bubble(cellchat, sources.use = c(8), targets.use = c(1),
#                  pairLR.use = pairLR.use,
#                  remove.isolate = TRUE)
# dev.off()
# 
# p <- plotGeneExpression(cellchat, signaling = cellchat@netP$pathways[c(1)])
# ggsave(p, filename = paste0(cellchat@netP$pathways[c(1)],"_GeneExpression_vln.pdf"), width = 8, height = 6)
# 
# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# 
# pdf(file="TIL.SNA_CXCL_signalingROLE.pdf",width=8,height=6)
# netAnalysis_signalingRole_network(cellchat, signaling = cellchat@netP$pathways[c(1)],
#                                   width = 15, height = 6, font.size = 10)
# dev.off()
# 
# p1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 5)
# p2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 5)
# pdf(file="TIL.SNA_CXCL_signalingPattern.pdf",width = 12,height = 5)
# print(p1+p2)
# dev.off()

# pdf(file="TIL.NMF_outgoing_selectK.pdf",width = 10,height = 5)
# selectK(cellchat, pattern = "outgoing")
# dev.off()
# 
# nPatterns = 3
# pdf(file="TIL.NMF_outgoing_comPattern_heatmap.pdf",width = 8,height = 7)
# cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns,
#                                           width = 5, height = 9, font.size = 6)
# 
# dev.off()
# 
# pdf(file="TIL.NMF_outgoing_comPattern_river.pdf",width = 8,height = 7)
# netAnalysis_river(cellchat, pattern = "outgoing")
# dev.off()
# 
# pdf(file="TIL.NMF_outgoing_comPattern_dotplot.pdf",width = 12,height = 5)
# netAnalysis_dot(cellchat, pattern = "outgoing")
# dev.off()
# 
# pdf(file="TIL.NMF_incoming_selectK.pdf",width = 10,height = 5)
# selectK(cellchat, pattern = "incoming")
# dev.off()
# 
# nPatterns = 4
# pdf(file="TIL.NMF_incoming_comPattern_heatmap.pdf",width = 10,height = 7)
# cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns,
#                                           width = 5, height = 9, font.size = 6)
# 
# dev.off()
# 
# pdf(file="TIL.NMF_incoming_comPattern_river.pdf",width = 8,height = 5)
# netAnalysis_river(cellchat, pattern = "incoming")
# dev.off()
# 
# pdf(file="TIL.NMF_incoming_comPattern_dotplot.pdf",width = 8,height = 5)
# netAnalysis_dot(cellchat, pattern = "incoming")
# dev.off()
# 
# 
# pdf(file="TIL.NMF_incoming_comPattern_river.pdf",width = 8,height = 5)
# netAnalysis_signalingRole_scatter(cellchat)
# dev.off()
# 
# pdf(file="Chord diagram.pdf",width = 8,height = 5)
# netVisual_chord_gene(cellchat, sources.use = c(2:17), targets.use = c(1), signaling = cellchat@netP$pathways[c(1)],legend.pos.x = 8)
# dev.off()
# 
# a=read.table("net_lr.txt",sep="\t",check.names=F,header=T)
# data <- as.data.frame(table(c(a$source,
#                               a$target)))
# colnames(data) <- c("Cell_Type","all_sum")
# data <- data[order(data$all_sum, decreasing = T),];data$Cell_Type <- factor(data$Cell_Type,levels = data$Cell_Type)
# head(data)
# sample_color <- c("#FB040B","#F6A717","#BA06FA","#172D7A")
# sample_color2 <- ifelse(data$all_sum > 10,ifelse(data$all_sum > 15,ifelse(data$all_sum > 20,sample_color[1],sample_color[2]),sample_color[3]),sample_color[4])
# sample_color1 <- ifelse(data$all_sum > 10,ifelse(data$all_sum > 15,ifelse(data$all_sum > 20,alpha(sample_color[1],0.9),alpha(sample_color[2],0.9)),alpha(sample_color[3],0.9)),alpha(sample_color[4],0.9))
# 
# pB2 <- ggplot(data = data, aes(x = Cell_Type, y = all_sum, fill = Cell_Type, colour = Cell_Type)) +
#   geom_bar(stat = "identity", width=0.8) + 
#   scale_fill_manual(values=sample_color1) +
#   scale_colour_manual(values=sample_color1) +
#   theme_bw()+
#   theme(panel.grid =element_blank()) +
#   labs(x="Cell Type",y="Counts")+
#   theme(axis.text.y = element_text(size=12, colour = "black"))+
#   theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1, vjust = 1 , colour = "black"))
# pdf(file="Interaction Count.pdf",width=10,height=8)
# print(pB2)
# dev.off()

# 细胞干性----
set.seed(2025)
library(CytoTRACE)
library(Seurat)
dir.create("../7.CytoTRACE")
scRNA <- readRDS("../result.rds")
Idents(scRNA) <- "cellType_1"

phe <- scRNA$cellType_1
phe = as.character(phe)
names(phe) <- rownames(scRNA@meta.data)
mat_3k <- as.matrix(scRNA@assays$RNA@counts)
mat_3k[1:4,1:4]
results <- CytoTRACE(mat = mat_3k)
emb <- scRNA@reductions[["umap"]]@cell.embeddings
# plotCytoGenes(results, numOfGenes = 10,outputDir = "/home/data/zxg/1/backup/4.Sub.Myeloid/3.Pseudo-time/2.Neu/")
plotCytoTRACE(results, phenotype = phe, emb = emb, outputDir = '../7.CytoTRACE/')

# 拟时序分析----
setwd("../code/")
source("SeuratV4CellAnnotateFunc.R")
# library(ClusterGVis)
# library(org.Hs.eg.db)
scRNA <- readRDS("../result.rds")

# 检查是否有 cellType_1 这一列
if ("cellType_1" %in% colnames(scRNA@meta.data)) {
  # 从每种细胞类型中抽取10000个细胞
  cell_types <- unique(scRNA@meta.data$cellType_1)
  selected_cells <- lapply(cell_types, function(ct) {
    cells_in_type <- WhichCells(scRNA, expression = cellType_1 == ct)
    # 检查该类细胞数是否大于10000
    if (length(cells_in_type) > 10000) {
      sample(cells_in_type, 10000)
    } else {
      cells_in_type
    }
  })
  # 将列表合并成一个向量
  selected_cells <- unlist(selected_cells)
  
  # 创建一个新的 Seurat 对象包含选中的细胞
  sub_pbmcSub <- subset(scRNA, cells = selected_cells)
  
  # 输出信息
  print("新的 Seurat 对象已创建，包含每种细胞类型随机抽取的10000个细胞。")
} else {
  print("没有找到 'cellType_1' 这一列，请检查元数据。")
}
scRNA <- sub_pbmcSub

# 运行6
scRNA <- scRNAAutoAnno(
  Path = "../",
  KeyCell = "Müller glia",
  SeuratObject = scRNA
)

library(monocle)
load("../8.Pseudo-time/cds.rda")

p <- plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 1);p
ggsave("../8.Pseudo-time/1.Trajectory.Pseudotime.pdf", plot = p,width = 8,height = 8)

p <- plot_cell_trajectory(cds, color_by = "group",cell_size = 1);p
ggsave("../8.Pseudo-time/2.Trajectory.group.pdf", plot = p,width = 8,height = 8)

p <- plot_cell_trajectory(cds, color_by = "Cluster",cell_size = 1)
ggsave( "../8.Pseudo-time/3.Trajectory.Cluster.pdf",plot = p,width = 8,height = 8)

# BEAM_res <- BEAM(cds, branch_point = 1, progenitor_method ="duplicate",cores = 36)
# BEAM_res <- BEAM_res[order(BEAM_res$qval),]
# BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
# 
# anno_color<-list(Cluster = c("1"="red","2"="blue","3"="grey","4"="green"))
# 
# m_ps_hm<-plot_genes_branched_heatmap2(cds[row.names(subset(BEAM_res,
#                                                            qval < 0.01)),],
#                                       branch_point = 3, #绘制的是哪个分支
#                                       num_clusters = 4, #分成几个cluster，根据第一次的结果看看是不是cluster太多了
#                                       cores = 36,
#                                       use_gene_short_name = T,
#                                       hmcols = colorRampPalette(c("#343DAE","#0072B5FF", "#FFFFBF", "#F46D43","#D73027"))(62),
#                                       show_rownames = F)
# save(m_ps_hm, file ="../11.Pseudo-time/m_ps_hm.rda")
# 
# ps_kegg_enrich <- enrichCluster(object = m_ps_hm, 
#                                 OrgDb = org.Hs.eg.db, 
#                                 type = "KEGG", 
#                                 organism = 'hsa', 
#                                 pvalueCutoff = 0.01, 
#                                 topn = 5, 
#                                 seed =123)
# 
# pdf('../8.Pseudo-time/4.Branch-enrich.PDF',height = 9,width = 12,onefile = F)
# visCluster(object = m_ps_hm, 
#            plot.type = "both", 
#            column_names_rot = 45, 
#            show_row_dend = FALSE, 
#            markGenes.side = "right",
#            go.col =rep(jjAnno::useMyCol("calm", n = 4), each = 5), #本例子里cluster1仅富集到3个KEGG通路，四个cluster总共18个通路
#            annoTerm.data = ps_kegg_enrich, 
#            go.size = "pval", 
#            line.side = "left")
# dev.off()


# 关键细胞与基因的相关性----
dir.create("../12.Correlation analysis/")
setwd("../12.Correlation analysis/")
library(pheatmap)
library(GSEABase)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(dplyr)
library(GSVA)
library(limma)
library(gplots)
library(ggExtra)

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}

Exper <- read.table("../data/RNA/symbol.txt",header = T,sep = "\t",row.names = 1,check.names = F)
Exper <- as.matrix(Exper)
geneSet=getGmt("Müller_glia_marker_logFC=1.gmt", geneIdType=SymbolIdentifier())
ssgseaScore=gsva(Exper, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
ssgseaOut=normalize(ssgseaScore)
save(ssgseaOut,file = "ssgseaOut.rda")#load("ssgseaOut.rda")
write.table(rbind(id=colnames(ssgseaOut),ssgseaOut),file="ssgseaOutAll.txt",sep="\t",quote=F,col.names=F)

plotdata=read.table("ssgseaOutAll.txt",sep="\t",header=T,check.names=F,row.names = 1)
gene <- read.table("../11.Venny/Intersection genes.txt")
rownames(gene) <- gene$V1
ARGexp <- Exper[rownames(gene),]

identical(colnames(ARGexp),colnames(plotdata))
ARGexp <- rbind(ARGexp,plotdata)

gene="Muller_glia"
corFilter=0
pFilter=0.05

rt=ARGexp

x=as.numeric(rt[gene,])
gene1=unlist(strsplit(gene,"\\|",))[1]
outputFile=paste(gene1,".cor.xls",sep="")
outTab=data.frame()

for(j in rownames(rt)){
  y=as.numeric(rt[j,])
  gene2=unlist(strsplit(j,"\\|",))[1]
  corT=cor.test(x,y)
  
  z=lm(y~x)
  cor=corT$estimate
  cor=round(cor,3)
  pvalue=corT$p.value
  if(pvalue<0.001){
    pval=signif(pvalue,4)
    pval=format(pval, scientific = TRUE)
  }else{
    pval=round(pvalue,3)}
  
  #??????????ͼƬ
  if((abs(cor)>corFilter) & (pvalue<pFilter)){
    pdfFile=paste(gene1,"_",gene2,".cor.pdf",sep="")
    pdf(file=pdfFile,width =6,height = 6)
    plot(x,y, type="p",pch=16,col="black",main=paste("Cor=",cor," (p-value=",pval,")",sep=""),
         cex=1, cex.lab=1, cex.main=1,cex.axis=1,
         xlab=paste(gene1,"expression"),
         ylab=paste(gene2,"expression") )
    lines(x,fitted(z),col=2)
    dev.off()
    outTab=rbind(outTab,cbind(gene1,gene2,cor,pvalue))
  }
}
write.table(file=outputFile,outTab,sep="\t",quote=F,row.names=F)


# GSEA----
setwd("../code/")
source("KeyGenePipeline.R")
mat <- read.table("../data/RNA/symbol.txt",header = T,row.names = 1,check.names = F,sep = "\t") %>% as.matrix()
mat <- log2(mat + 1)
# 运行1
KeyGenePipeline(
  Path = "../",
  species = "mouse",
  mat = mat,
  KeyGene = c("Gnal","Mir384","Rpl21")
)
# 输入通路：
# Gnal：
# rno03015
# rno04310
# rno04390

# Mir384：
# rno04020
# rno04621
# rno04390

# Rpl21
# rno04623
# rno04064
# rno03320



# GSVA----
library(msigdbr)
setwd("../code/")
source("KeyGenePipeline.R")

KEGG_df_all <-  msigdbr(species = "Rattus norvegicus", # Homo sapiens or Mus musculus
                        category = "H") 
KEGG_df <- dplyr::select(KEGG_df_all,gs_name,gs_exact_source,gene_symbol)
kegg_list <- split(KEGG_df$gene_symbol, KEGG_df$gs_name) ##按照gs_name给gene_symbol分组

mat <- read.table("../data/RNA/symbol.txt",header = T,row.names = 1,check.names = F,sep = "\t") %>% as.matrix()
mat <- log2(mat + 1)
# 运行1
KeyGenePipeline(
  Path = "../",
  species = "mouse",
  mat = mat,
  KeyGene = c("Gnal","Mir384","Rpl21")
)



# 疾病相关性 --------------------------------------------------------------------
setwd("../code/")
source("KeyGenePipeline.R")

mat <- read.table("../data/RNA/symbol.txt", header = T, sep = "\t", row.names = 1, check.names = F)
mat <- log2(mat + 1)

KeyGene <- c("Gnal","Mir384","Rpl21")
# select:3
KeyGenePipeline(Path = "../", 
                mat, 
                KeyGene = KeyGene)

# 疾病基因：
# Col4a1
# Col18a1
# Six6
# Ltbp2
# Lmx1b
# Foxc1
# Loxl1
# Tek
# Efemp1
# Asb10
# Atoh7
# Optn
# Foxe3
# Wdr36
# Opa1
# Tbk1
# Fbn1
# Ntf4
# Sbf2
# Pax6
# Myoc
# Best1
# Pxdn
# Cyp1b1
# Optc


# 免疫代谢通路----
source("SeuratV4CellAnnotateFunc.R")
library(clusterProfiler)
scRNA <- readRDS("../result.rds")
# 运行9
scRNA <- scRNAAutoAnno(
  Path = "../",
  SeuratObject = scRNA,
  KeyGene = c("GNAL","MIR384","RPL21")
)

# 疾病基因共表达----
source("SeuratV4CellAnnotateFunc.R")
scRNA <- readRDS("../result.rds")
# 运行9
scRNA <- scRNAAutoAnno(
  Path = "../",
  SeuratObject = scRNA,
  KeyGene = c("GNAL")
)

# 疾病基因
# COL4A1
# FBN1