require(Seurat)
require(ggplot2)
require(harmony)
require(dplyr)
require(reshape2)
require(scRNAtoolVis)
require(ggpubr)
# require(clusterProfiler)
options(dplyr.print_max = 1e9)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(123456)

ggsave_fun <- function(filename, plot = last_plot(), plot_Device = c(".pdf"), width = 7, height = 7, dpi = 600, ...) {
  filenames <- paste0(filename, plot_Device)
  for (filename_tmp in filenames) {
    ggsave(filename = filename_tmp, plot = plot, width = width, height = height, dpi = dpi, ...)
  }
}

paramSweep_V4 <- function(seu, PCs = 1:10, sct = FALSE, num.cores = 1) {
  require(Seurat)
  require(fields)
  require(parallel)
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
  pN <- seq(0.05, 0.3, by = 0.05)

  ## Remove pK values with too few cells
  min.cells <- round(nrow(seu@meta.data) / (1 - 0.05) - nrow(seu@meta.data))
  pK.test <- round(pK * min.cells)
  pK <- pK[which(pK.test >= 1)]

  ## Extract pre-processing parameters from original data analysis workflow
  orig.commands <- seu@commands

  ## Down-sample cells to 10000 (when applicable) for computational effiency
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 10000, replace = FALSE)]
    data <- seu@assays$RNA@counts[, real.cells]
    n.real.cells <- ncol(data)
  }

  if (nrow(seu@meta.data) <= 10000) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
    n.real.cells <- ncol(data)
  }

  ## Iterate through pN, computing pANN vectors at varying pK
  # no_cores <- detectCores()-1
  if (num.cores > 1) {
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)),
      FUN = parallel_paramSweep,
      n.real.cells,
      real.cells,
      pK,
      pN,
      data,
      orig.commands,
      PCs,
      sct, mc.cores = num.cores
    )
    stopCluster(cl)
  } else {
    output2 <- lapply(as.list(1:length(pN)),
      FUN = parallel_paramSweep,
      n.real.cells,
      real.cells,
      pK,
      pN,
      data,
      orig.commands,
      PCs,
      sct
    )
  }

  ## Write parallelized output into list
  sweep.res.list <- list()
  list.ind <- 0
  for (i in 1:length(output2)) {
    for (j in 1:length(output2[[i]])) {
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }

  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_"))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}

doubletFinder_V4 <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE, annotations = NULL) {
  require(Seurat)
  require(fields)
  require(KernSmooth)

  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if (reuse.pANN != FALSE) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, sep = "_")] <- classifications
    return(seu)
  }

  if (reuse.pANN == FALSE) {
    ## Make merged real-artifical data
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells / (1 - pN) - n_real.cells)
    print(paste("Creating", n_doublets, "artificial doublets...", sep = " "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2]) / 2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    # Keep track of the types of the simulated doublets
    if (!is.null(annotations)) {
      stopifnot(typeof(annotations) == "character")
      stopifnot(length(annotations) == length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    ## Store important pre-processing information
    orig.commands <- seu@commands

    ## Pre-process Seurat object
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets,
        normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
        scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
        margin = orig.commands$NormalizeData.RNA@params$margin
      )

      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
        selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
        loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
        clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
        mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
        dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
        num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
        binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
        nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
        mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
        dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff
      )

      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets,
        features = orig.commands$ScaleData.RNA$features,
        model.use = orig.commands$ScaleData.RNA$model.use,
        do.scale = orig.commands$ScaleData.RNA$do.scale,
        do.center = orig.commands$ScaleData.RNA$do.center,
        scale.max = orig.commands$ScaleData.RNA$scale.max,
        block.size = orig.commands$ScaleData.RNA$block.size,
        min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block
      )

      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets,
        features = orig.commands$ScaleData.RNA$features,
        npcs = length(PCs),
        rev.pca = orig.commands$RunPCA.RNA$rev.pca,
        weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
        verbose = FALSE
      )
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc() # Free up memory
    }

    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)

      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }

    ## Compute PC distance matrix
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)

    ## Compute pANN
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    if (!is.null(annotations)) {
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells)) / k
      if (!is.null(annotations)) {
        for (ct in unique(annotations)) {
          neighbors_that_are_doublets <- neighbors[neighbors > n_real.cells]
          if (length(neighbors_that_are_doublets) > 0) {
            neighbor_types[i, ] <-
              table(doublet_types1[neighbors_that_are_doublets - n_real.cells]) +
              table(doublet_types2[neighbors_that_are_doublets - n_real.cells])
            neighbor_types[i, ] <- neighbor_types[i, ] / sum(neighbor_types[i, ])
          } else {
            neighbor_types[i, ] <- NA
          }
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet", n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data), 1]
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, sep = "_")] <- classifications
    if (!is.null(annotations)) {
      colnames(neighbor_types) <- levels(doublet_types1)
      for (ct in levels(doublet_types1)) {
        seu@meta.data[, paste("DF.doublet.contributors", pN, pK, nExp, ct, sep = "_")] <- neighbor_types[, ct]
      }
    }
    return(seu)
  }
}

human2mouse <- function(x) {
  x1 <- substr(x, 1, 1)
  x2 <- substr(x, 2, nchar(x))
  return(paste0(x1, tolower(x2)))
}

getScatterplot <- function(object, gene1, gene2, cor.method = "pearson", jitter.num = 0.15, pos = TRUE) {
  if (!gene1 %in% rownames(object)) {
    print("gene1 was not found")
    if (!gene2 %in% rownames(object)) {
      print("gene2 was not found")
    }
  } else {
    exp.mat <- GetAssayData(object = object, assay = "RNA") %>%
      .[c(gene1, gene2), ] %>%
      as.matrix() %>%
      t() %>%
      as.data.frame()
    if (pos) {
      if (nrow(exp.mat[which(exp.mat[, 1] > 0 & exp.mat[, 2] > 0), ]) > (nrow(exp.mat) * 0.01)) {
        exp.mat <- exp.mat[which(exp.mat[, 1] > 0 & exp.mat[, 2] > 0), ]
      } else {
        exp.mat <- exp.mat[which(exp.mat[, 1] > 0 | exp.mat[, 2] > 0), ]
      }
    }
    colnames(exp.mat) <- c("Var1", "Var2")
    plots <- ggplot(data = exp.mat, mapping = aes_string(x = "Var1", y = "Var2")) +
      geom_smooth(method = "lm", se = T, color = "red", size = 1) +
      stat_cor(method = cor.method) +
      labs(x = gene1, y = gene2) +
      geom_jitter(width = jitter.num, height = jitter.num, color = "black", size = 1, alpha = 1) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        legend.text = element_text(colour = "black", size = 10),
        axis.text = element_text(colour = "black", size = 10),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(size = 1, linetype = "solid", colour = "black"),
        panel.background = element_rect(fill = "white")
      )
    return(plots)
  }
}

capitalize_first <- function(string) {
  words <- strsplit(string, " ")[[1]]
  capitalized_words <- toupper(substring(words, 1, 1))
  rest_of_words <- substring(words, 2)
  capitalized_string <- paste0(capitalized_words, rest_of_words, collapse = " ")
  return(capitalized_string)
}



#' Main Function
#'
#' @param Path 保存结果的路径
#' @param SeuratObject 创建好seuratV4对象且拥有分组的rds，分组名字为group，肿瘤样本为Normal vs Tumor，非肿瘤样本为Control vs Disease
#' @param Multi 是否为多样本
#' @param ref SingleR参考数据集
#' @param labels 参考数据集的分类，label.main为主要细胞亚群，label.fine为细分亚群
#' @param min.features  细胞所含的最低基因数量
#' @param Idents 注释分辨率选择
#' @param DoubletFinder 是否过滤双细胞
#' @param species 物种，人无需填写，鼠填写"mouse"
#' @param KeyCell 关键细胞，仅在拟时序分析中使用
#' @param KeyGene 关键基因
#' @param logFCfilter logFC阈值
#' @param adjPvalFilter p值阈值
#' @param assay 降维聚类使用的assay，默认RNA
#' @param set.resolutions 分群
#' @param PC 选择主成分个数
#' @param nfeatures 高变基因个数
#' @param npcs 碎石图的主成分个数
#'
#' @return
#' @export
#'
#' @examples
scRNAAutoAnno <- function(Path = ".", SeuratObject = NA, Multi = TRUE, ref, labels,min.features = 500, Idents = "RNA_snn_res.0.2", DoubletFinder = TRUE, species = NA, KeyCell = NA, KeyGene = NA, logFCfilter = 0.585, adjPvalFilter = 0.05,assay = "RNA", set.resolutions = seq(0.2, 1.2, by = 0.1), PC = 30, nfeatures = 3000, npcs = 50) {
  cat("1.质量控制\n")
  cat("2.数据标准化\n")
  cat("3.自动化注释\n")
  cat("4.Cellchat\n")
  cat("5.亚群贡献度\n")
  cat("6.拟时序分析\n")
  cat("7.关键基因表达丰度\n")
  cat("8.疾病基因共表达网络\n")
  cat("9.免疫代谢通路\n")
  cat("10.代谢通路热图\n")
  cat("11.量化热点机制\n")
  cat("12.富集分析\n")
  cat("\033[0;34m请选择要分析的内容，按序号填写，英文逗号分隔")
  AnalysisIndex <- scan(sep = ",", quiet = TRUE)

  # AdvanceIndex <- length(which(AnalysisIndex %in% c(1:12) == FALSE)) # 判断是否有其他高级分析
  # if (is.na(KeyCell) & AdvanceIndex > 0) {
  #   message("是否需要在分析过程中挑选关键细胞？yes/no")
  #   isKeyCell <- scan(sep = ",", what = "character", quiet = TRUE)
  #   if (isKeyCell == "yes") {
  #     cat("1.细胞含量变化最大\n")
  #     cat("2.Cellchat\n")
  #     cat("3.亚群贡献度\n")
  #     cat("\033[0;34m从哪一步分析中挑选关键细胞？\n")
  #     KeyCellIndex <- scan(sep = ",", quiet = TRUE)
  #   }
  # }

  scRNA <- SeuratObject

  if (!is.na(species)) {
    mouse <- ifelse(species == "mouse", T, F)
  } else {
    mouse <- F
  }

  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes

  if (mouse) {
    s.genes <- human2mouse(s.genes)
    g2m.genes <- human2mouse(g2m.genes)
  }

  n <- 0
  KeyCells <- KeyCell
  for (step in AnalysisIndex) {
    n <- n + 1

    # 质量控制
    if (step == 1) {
      message("*************************质量控制*************************")

      dir <- paste0(Path, "/.QCFilter")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      grep("^[M,m][T,t]-", rownames(scRNA), value = T)
      grep("^R[P,p][SL,sl]", rownames(scRNA), value = T)
      scRNA[["percent.mt"]] <- PercentageFeatureSet(object = scRNA, pattern = "^[M,m][T,t]-")
      scRNA[["percent.ribo"]] <- PercentageFeatureSet(scRNA, pattern = "^R[P,p][SL,sl]") # 核糖体基因占比

      nomt <- ifelse(max(scRNA[["percent.mt"]]) == 0, T, F)
      if (nomt) {
        MT_genes <- c("ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB")
        if (mouse) {
          MT_genes <- human2mouse(MT_genes)
        }
        MT_genes <- intersect(MT_genes, rownames(scRNA))
        scRNA[["percent.mt"]] <- PercentageFeatureSet(object = scRNA, features = MT_genes)
      }

      mt_value <- stats::mad(scRNA$percent.mt)
      lower_mt <- median(scRNA$percent.mt) - 3 * mt_value
      upper_mt <- median(scRNA$percent.mt) + 3 * mt_value
      
      nF_value <- stats::mad(scRNA$nFeature_RNA)
      lower_nF <- median(scRNA$nFeature_RNA) - 3 * nF_value
      upper_nF <- median(scRNA$nFeature_RNA) + 3 * nF_value
      
      nC_value <- stats::mad(scRNA$nCount_RNA)
      lower_nC <- median(scRNA$nCount_RNA) - 3 * nC_value
      upper_nC <- median(scRNA$nCount_RNA) + 3 * nC_value
      
      # nR_value <- stats::mad(scRNA$percent.ribo)
      # lower_nR <- median(scRNA$percent.ribo) - 3 * nR_value
      # upper_nR <- median(scRNA$percent.ribo) + 3 * nR_value

      scRNA <- subset(x = scRNA, subset = nFeature_RNA >= min.features & percent.mt <= upper_mt & nFeature_RNA <= upper_nF & nCount_RNA <= upper_nC)
      dim(scRNA)

      if (DoubletFinder) {
        require(DoubletFinder)
        scRNA_list <- SplitObject(scRNA, split.by = "orig.ident")
        newscRMA_list <- list()
        for (i in 1:length(scRNA_list)) {
          data <- scRNA_list[[i]]
          data <- NormalizeData(data)
          data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
          data <- ScaleData(data)
          data <- RunPCA(data)
          data <- RunUMAP(data, dims = 1:20)

          sweep.res.list <- paramSweep_V4(data, PCs = 1:20, sct = FALSE) # 若使用SCT方法标准化则'sct=T'
          sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
          bcmvn <- find.pK(sweep.stats)
          p <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC == max(bcmvn$MeanBC), ]$pK))

          # 期望doublet数量
          homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters) # 可使用注释好的细胞类型
          Doubletrate <- ncol(data) * 8 * 1e-6
          nExp_poi <- round(Doubletrate * ncol(data))
          nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

          # 鉴定doublets
          data <- doubletFinder_V4(data, PCs = 1:20, pN = 0.25, pK = p, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
          colnames(data@meta.data)[ncol(data@meta.data)] <- "doublet_info"
          newscRMA_list[[i]] <- data
        }

        saveRDS(newscRMA_list, paste0(dir, "/doubletFinder_V4.rds"))

        plot_list <- list()
        for (i in 1:length(newscRMA_list)) {
          plot_list[[i]] <- DimPlot(newscRMA_list[[i]], group.by = "doublet_info") + ggtitle(newscRMA_list[[i]]$orig.ident[1]) + theme(plot.title = element_text(hjust = 0.5))
        }


        for (i in 1:length(plot_list)) {
          plot_list[[i]] <- plot_list[[i]] + NoAxes() + NoLegend()
        }

        p <- CombinePlots(plot_list, legend = "right")

        ggsave_fun(filename = paste0(dir, "/1.DoubletFinder"), plot = p, width = 12, height = 10)

        if (Multi) {
          scRNA <- merge(newscRMA_list[[1]], newscRMA_list[2:length(newscRMA_list)])
        } else {
          scRNA <- newscRMA_list[[1]]
        }

        scRNA <- subset(scRNA, subset = doublet_info == "Singlet")

        p <- VlnPlot(scRNA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), group.by = "orig.ident", ncol = 4)
        # p
        ggsave_fun(filename = paste0(dir, "/2.VlnPlot"), plot = p, width = 20, height = 4)

        plot1 <-
          FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = F) +
          NoLegend()
        plot2 <-
          FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = F) +
          NoLegend()
        plot3 <-
          FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.ribo", raster = F) +
          NoLegend()

        plot1 + plot2 + plot3
        ggsave_fun(filename = paste0(dir, "/3.FeatureScatter"), plot = plot1 + plot2 + plot3, height = 7, width = 21)
      } else {
        p <- VlnPlot(scRNA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), group.by = "orig.ident", ncol = 4)
        # p
        ggsave_fun(filename = paste0(dir, "/1.VlnPlot"), plot = p, width = 20, height = 4)

        plot1 <-
          FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = F) +
          NoLegend()
        plot2 <-
          FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = F) +
          NoLegend()
        plot3 <-
          FeatureScatter(scRNA, feature1 = "percent.mt", feature2 = "percent.ribo", raster = F) +
          NoLegend()

        plot1 + plot2 + plot3
        ggsave_fun(filename = paste0(dir, "/2.FeatureScatter"), plot = plot1 + plot2 + plot3, height = 7, width = 21)
      }
    }

    # 数据标准化
    if (step == 2) {
      message("*************************数据标准化*************************")

      dir <- paste0(Path, "/.Cluster")

      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      DefaultAssay(scRNA) <- assay

      scRNA <- NormalizeData(object = scRNA)
      scRNA <- CellCycleScoring(object = scRNA, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
      scRNA <- FindVariableFeatures(object = scRNA, selection.method = "vst", nfeatures = nfeatures)
      top10 <- head(VariableFeatures(scRNA), 10)
      plot1 <- VariableFeaturePlot(scRNA)
      plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
      plot1 + plot2
      ggsave_fun(filename = paste0(dir, "/1.VariableFeaturePlot"), plot = plot1 + plot2, width = 14)
      scRNA <- ScaleData(object = scRNA, vars.to.regress = c("percent.mt", "percent.ribo", "S.Score", "G2M.Score"))
      scRNA <- RunPCA(scRNA, verbose = T, npcs = npcs, features = VariableFeatures(object = scRNA))
      p <- ElbowPlot(object = scRNA, ndims = npcs)
      ggsave_fun(filename = paste0(dir, "/2.ElbowPlot"), plot = p)

      if (Multi) {
        scRNA <- RunHarmony(object = scRNA, group.by.vars = "orig.ident", assay.use = assay, verbose = FALSE)
        scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:PC, verbose = T)
        # scRNA <- RunTSNE(scRNA, reduction = "harmony", dims = 1:PC, verbose = T)
        scRNA <- FindNeighbors(scRNA, dims = 1:PC, reduction = "harmony", verbose = T) # 使用harmony替代PCA
        scRNA <- FindClusters(scRNA, resolution = set.resolutions, verbose = T)
        p <- DimPlot(object = scRNA, reduction = "pca", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/3.PcaPlot"), plot = p)
        p <- DimPlot(object = scRNA, reduction = "harmony", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/4.HarmonyPlot"), plot = p)
      } else {
        scRNA <- RunUMAP(scRNA, reduction = "pca", dims = 1:PC, verbose = T)
        # scRNA <- RunTSNE(scRNA, reduction = "pca", dims = 1:PC, verbose = T)
        scRNA <- FindNeighbors(scRNA, dims = 1:PC, reduction = "pca", verbose = T)
        scRNA <- FindClusters(scRNA, resolution = set.resolutions, verbose = T)
        p <- DimPlot(object = scRNA, reduction = "pca", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/3.PcaPlot"), plot = p)
      }

      pdf(paste0(dir, "/7.data.merge.harmony.pdf"))
      # p <- clustree::clustree(scRNA)
      # print(p)
      p <- DimPlot(object = scRNA, reduction = "umap", label = TRUE, group.by = "orig.ident", raster = F) + NoLegend()
      print(p)
      merge.res <- sapply(set.resolutions, function(x) {
        p <- DimPlot(object = scRNA, reduction = "umap", label = TRUE, group.by = paste0(assay, "_snn_res.", x), raster = F) + NoLegend()
        print(p)
      })
      dev.off()
      saveRDS(scRNA, paste0(dir, "/data.merge.harmony.2000.rds"))
    }

    # 自动化注释
    if (step == 3) {
      message("*************************自动化注释*************************")

      dir <- paste0(Path, "/.CellAnnotate")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      # p <- DimPlot(scRNA,group.by = Idents,label = T,reduction = "tsne")+NoLegend();p
      # ggsave_fun(filename = "2.Cluster/5.Tsne.Cluster",plot = p)

      Cluster.dir <- list.files(Path, pattern = "^.Cluster$", all.files = TRUE)

      p <- DimPlot(scRNA, group.by = Idents, label = T, reduction = "umap", raster = F ,pt.size = 0.1) + NoLegend()
      print(p)
      ggsave_fun(filename = paste0(Path, "/", Cluster.dir, "/5.Umap.Cluster"), plot = p)

      CELL <- SingleR::SingleR(
        test = as.matrix(scRNA@assays$RNA@data),
        ref = ref,
        labels = labels,
        clusters = scRNA@meta.data[[Idents]]
      )
      
      cells <- CELL[, ncol(CELL):1]
      cells <- data.frame(cells)
      cell_id <- cells$pruned.labels
      names(cell_id) <- levels(scRNA@meta.data[[Idents]])
      Idents(scRNA) <- Idents
      
      if (anyNA(cell_id)) {
        scRNA <- RenameIdents(scRNA, cell_id)
        scRNA$cellType_1 <- Idents(scRNA)
        print(cell_id)
        p <- DimPlot(scRNA, group.by = "cellType_1", raster = F)
        print(p)
        message("某些细胞簇没注释到，是否强行为NA的细胞命名为unknow？yes/no")
        exit <- scan(what = "character", sep = ",", quiet = TRUE)
        
        if (exit == "yes") {
          cell_id[which(is.na(cell_id))] <- "unknow"
        }
        
        if (exit == "no") {
          message("某些细胞簇没注释到，具体请看右侧umap图，建议选择其他分辨率或参考数据集重新运行此函数，跳过1.质量控制和2.数据标准化，选择3重新注释。")
          return(scRNA)
        }
      }
      
      if (!anyNA(cell_id)) {
        p <- DimPlot(scRNA, group.by = Idents, raster = F,label = T)+NoLegend()
        print(p)
        print(paste0(paste0(names(cell_id),":",as.character(cell_id)), collapse = "; "))
        
        # a <- scan(what = "character")
        while (T) {
          message("请输入基因进行检验，跳过该步骤请按回车键：")
          a <- scan(what = "character")
          gene <- intersect(a,rownames(scRNA))
          if(length(gene) > 0){
            p <- FeaturePlot(scRNA,features = gene,raster = F,label = T,cols = c("lightgrey","red"))
            print(p)
          }
          if(length(a) == 0){break}
        }
        
        cat("\033[0;34m由于使用了多个参考数据，参考数据细胞名字不同，请按顺序对下述细胞进行重命名，英文;分隔。\n例：将\033[31mCD8+ T cells,Neutrophils,CD8+ T-cells\033[0;34m更改为\033[31mCD8+ T cells,Neutrophils,CD8+ T cells\033[0;34m以统一细胞名字\n")
        print(paste0(paste0(names(cell_id),":",as.character(cell_id)), collapse = "; "))
        cell_id <- scan(what = "character", sep = ";", quiet = FALSE)
        cell_id <- gsub(".*:","",cell_id)
        names(cell_id) <- levels(scRNA@meta.data[[Idents]])
        Idents(scRNA) <- Idents
        scRNA <- RenameIdents(scRNA, cell_id)
        scRNA$cellType_1 <- Idents(scRNA)
        print(table(scRNA$cellType_1))
        p <- DimPlot(scRNA, group.by = "cellType_1", raster = F)
        print(p)
        saveRDS(scRNA, paste0(Path, "/result.rds"))

        Idents(scRNA) <- "cellType_1"
        scRNA.markers <- FindAllMarkers(scRNA, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
        scRNA.markers <- scRNA.markers[scRNA.markers$p_val_adj < 0.05,]
        write.csv(scRNA.markers, paste0(dir, "/Findall.markers.cellType_1.csv"))

        # scRNA.markers <- read.csv("../.CellAnnotate/Findall.markers.cellType_1.csv",header = T,row.names = 1)
        scRNA.markers %>%
          group_by(cluster) %>%
          top_n(n = 5, wt = avg_log2FC) -> features

        print((features))

        features <- features$gene
        cat("\033[0;34m是否要改写基因？yes/no\n")
        RenameGene <- scan(what = "character", sep = ",", quiet = TRUE)

        if (RenameGene == "yes") {
          print(paste0(as.character(features), collapse = ","))
          cat("\033[0;34m请改写\n")
          features <- scan(what = "character", sep = ",", quiet = FALSE)
        }

        # plot_list <- FeaturePlot(
        #   scRNA,
        #   features = unique(features),
        #   reduction = "umap",
        #   raster = F,
        #   combine = FALSE, cols = c("lightgrey", "red")
        # )
        #
        # for (i in 1:length(plot_list)) {
        #   plot_list[[i]] <- plot_list[[i]] + NoLegend() + NoAxes()
        # }
        #
        # p <- CombinePlots(plot_list)
        # # p
        # ggsave_fun(filename = paste0(dir, "/1.Featureplot.umap"), plot = p, width = 12, height = 12)

        # p <- DimPlot(scRNA, group.by = "cellType_1", label = T, label.size = 3, reduction = "umap", raster = F) + NoLegend()
        # p
        # ggsave_fun(filename = paste0(dir, "/2.Umap.plot"), plot = p, width = 7, height = 7)
        my_pal2 <- c(
          "#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
          "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A","#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
          "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A","#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
          "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A","#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
          "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A"
        )

        umap <- scRNA@reductions$umap@cell.embeddings %>%
          as.data.frame() %>%
          cbind(cellType = scRNA@meta.data$cellType_1)

        ## 计算标签中心位置区域
        celltypepos <- umap %>%
          group_by(cellType) %>%
          summarise(
            umap_1 = median(UMAP_1),
            umap_2 = median(UMAP_2)
          )

        p <- clusterCornerAxes(object = scRNA, reduction = "umap", pSize = 0.1, clusterCol = "cellType_1", noSplit = T) +
          ggrepel::geom_label_repel(aes(x = umap_1, y = umap_2, label = cellType, color = cellType),
            fontface = "bold",
            data = celltypepos,
            box.padding = 0.5, show.legend = FALSE
          )
        print(p)
        ggsave_fun(filename = paste0(dir, "/1.Umap1.plot"), plot = p, width = 10, height = 7)

        # add circle
        # p <- clusterCornerAxes(object = scRNA, reduction = "umap", pSize = 0.1, clusterCol = "cellType_1", noSplit = T, cornerTextSize = 3.5, addCircle = TRUE)
        # print(p)
        # ggsave_fun(filename = paste0(dir, "/3.Umap2.plot"), plot = p, width = 10, height = 7)

        p <- DotPlot(scRNA, features = unique(features), cols = "RdYlBu", group.by = "cellType_1") + scale_size_continuous(range = c(0, 10)) + theme(
          panel.border = element_rect(colour = "black"),
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
          ),
          legend.position = "top",
          legend.key.height = unit(0.3, "cm"),
          legend.key.width = unit(0.8, "cm"),
        )
        print(p)
        if (length(unique(features)) >= 50) {
          ggsave_fun(filename = paste0(dir, "/2.Dotplot"), plot = p, width = 20, height = 9)
        }
        if (length(unique(features)) < 50) {
          ggsave_fun(filename = paste0(dir, "/2.Dotplot"), plot = p, width = 17, height = 9)
        }

        tryCatch(
          {
            cellRatio <- cellRatioPlot(
              object = scRNA,
              sample.name = "group",
              celltype.name = "cellType_1",
              flow.curve = 0.5,
              fill.col = my_pal2
            ) + theme(
              axis.text.x = element_text(
                angle = 45,
                hjust = 1,
                vjust = 1
              )
            )
            print(cellRatio)

            ggsave_fun(filename = paste0(dir, "/3.CellAnnoted.cellType.group.ratio.plot"), plot = cellRatio, width = 7, height = 7)

            data_p <- scRNA@meta.data[, c("group", "cellType_1")]

            p <- ggstatsplot::ggbarstats(
              data = data_p,
              x = group,
              y = cellType_1,
              palette = "Set3"
            ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
            print(p)
            leg <- length(unique(data_p$cellType_1))
            ggsave_fun(filename = paste0(dir, "/4.Dysfunction"), plot = p, width = 2+1.5*leg, height = 7)

            EachSampleClusterDis <- lapply(unique(scRNA@meta.data$orig.ident), function(Sample) {
              seur <- subset(scRNA, cells = rownames(scRNA@meta.data[scRNA@meta.data$orig.ident %in% Sample, ]))
              seur@meta.data$cellType_1 %>%
                table() %>%
                data.frame() %>%
                magrittr::set_colnames(c("CellTypes", "Number")) %>%
                dplyr::mutate(Per = 100 * Number / sum(Number))
            })
            names(EachSampleClusterDis) <- unique(scRNA@meta.data$orig.ident)
            EachSampleClusterDisA <- dplyr::bind_rows(EachSampleClusterDis) %>% dplyr::mutate(Sample = rep(names(EachSampleClusterDis), times = unlist(lapply(EachSampleClusterDis, nrow))))
            p1 <- ggplot(data = EachSampleClusterDisA, aes(x = reorder(CellTypes, Number), fill = Sample, y = Per)) + # log2(as.numeric(CopyNumber)))
              geom_col(position = "fill", width = 0.8) +
              scale_y_continuous(expand = c(0, 0)) +
              labs(y = "Percentage of cell") +
              scale_fill_manual(values = my_pal2, name = "") +
              theme(
                axis.text = element_text(color = "black"),
                panel.background = element_blank(), # panel.grid=element_blank(),
                legend.title = element_blank(),
                axis.text.y = element_text(color = "black"),
                axis.line = element_line(color = "black"),
                axis.title.y = element_blank(), # 轴标题
                legend.position = "right",
                # legend.key.height = unit(0.3, "cm"),
                # legend.key.width = unit(0.3, "cm"),
                legend.direction = "vertical", # vertical  legend.direction = 'horizontal'
                legend.text = element_text(size = 10)
              ) +
              coord_flip()

            EachSampleClusterDis <- lapply(unique(scRNA@meta.data$group), function(Sample) {
              seur <- subset(scRNA, cells = rownames(scRNA@meta.data[scRNA@meta.data$group %in% Sample, ]))
              seur@meta.data$cellType_1 %>%
                table() %>%
                data.frame() %>%
                magrittr::set_colnames(c("CellTypes", "Number")) %>%
                dplyr::mutate(Per = 100 * Number / sum(Number))
            })

            my_pal1 <- c("#D51F26", "#272E6A", "#208A42", "#89288F", "#6387C5")

            names(EachSampleClusterDis) <- unique(scRNA@meta.data$group)
            EachSampleClusterDisA <- dplyr::bind_rows(EachSampleClusterDis) %>% dplyr::mutate(Sample = rep(names(EachSampleClusterDis), times = unlist(lapply(EachSampleClusterDis, nrow))))
            p2 <- ggplot(data = EachSampleClusterDisA, aes(x = reorder(CellTypes, Number), fill = Sample, y = Per)) + # log2(as.numeric(CopyNumber)))
              geom_col(position = "fill", width = 0.8) +
              scale_y_continuous(expand = c(0, 0)) +
              labs(y = "Percentage of cell") +
              scale_fill_manual(values = my_pal1, name = "") +
              theme(
                axis.text = element_text(color = "black", ),
                panel.background = element_blank(), panel.grid = element_blank(),
                legend.title = element_blank(),
                axis.text.y = element_blank(), axis.title.x = element_blank(), # 删除轴文本
                axis.ticks.y.left = element_blank(),
                axis.title.y = element_blank(), # 轴标题
                axis.line.x = element_line(color = "black"),
                legend.position = "right",
                legend.direction = "vertical", # vertical  legend.direction = 'horizontal'
                # legend.key.height = unit(0.3, "cm"),
                # legend.key.width = unit(0.3, "cm"),
                legend.text = element_text(size = 10)
              ) +
              coord_flip()

            # 计算每种细胞类型的细胞数量
            cell_counts <- table(scRNA@meta.data$cellType_1)

            # 创建数据框
            cell_counts_df <- data.frame(
              CellType = names(cell_counts),
              Count = as.numeric(cell_counts / 1000)
            )

            # 创建柱状图
            cell_counts_df <- cell_counts_df[order(cell_counts_df$Count, decreasing = T), ]

            p3 <- ggplot(cell_counts_df, aes(x = reorder(CellType, Count), y = Count)) +
              scale_y_continuous(expand = c(0, 0)) +
              labs(y = "Number of cell(10^3)") +
              geom_bar(stat = "identity", fill = "#6387C5") +
              theme(
                axis.text = element_text(color = "black"),
                panel.background = element_blank(), panel.grid = element_blank(),
                legend.title = element_blank(),
                axis.text.y = element_blank(),
                # axis.title.x = element_blank(),#删除轴文本
                # plot.tag=element_blank(),
                axis.ticks.y.left = element_blank(),
                # axis.title.x.top = T,
                axis.title.y = element_blank(), # 轴标题
                axis.line.x = element_line(color = "black"),
                legend.text = element_text(size = 10)
              ) +
              coord_flip()
            print(p1 + p2 + p3)
            ggsave_fun(filename = paste0(dir, "/5.MergeRatio"), plot = p1 + p2 + p3, width = 15, height = 6)
            # if (isKeyCell == "yes" && is.na(KeyCell) && 1 %in% KeyCellIndex) {
            pB2_df <- table(scRNA@meta.data$cellType_1, scRNA@meta.data$group) %>% melt()
            colnames(pB2_df) <- c("Cluster", "Sample", "Number")
            pB2_df$Cluster <- factor(pB2_df$Cluster)

            pB2_df <- pB2_df %>%
              group_by(Sample) %>%
              mutate(Percentage = Number / sum(Number))

            pB2_df_tumor <- pB2_df %>%
              filter(Sample %in% c("Disease", "Tumor", "Cancer")) %>%
              select(Cluster, Tumor_Number = Number)

            pB2_df_control <- pB2_df %>%
              filter(Sample == c("Control", "Normal")) %>%
              select(Cluster, Control_Number = Number)

            pB2_df_combined <- merge(pB2_df_tumor, pB2_df_control, by = "Cluster")
            # pB2_df_combined <- pB2_df_combined[pB2_df_combined$Control_Number >= 50 && pB2_df_combined$Tumor_Number >= 50,]
            pB2_df_combined$rate <- (pB2_df_combined$Tumor_Number - pB2_df_combined$Control_Number) / pB2_df_combined$Control_Number
            KeyCell <- pB2_df_combined[abs(pB2_df_combined$rate) == max(abs(pB2_df_combined$rate)), ]$Cluster %>% as.character()
            print(paste0("变化最大的细胞：", KeyCell))
            write.csv(pB2_df_combined, paste0(dir, "/Cell.rate.csv"))
            # }
          },
          error = function(e) {
            message("没有分组信息，跳过此步骤\n", e)

            cellRatio <- cellRatioPlot(
              object = scRNA,
              sample.name = "orig.ident",
              celltype.name = "cellType_1",
              flow.curve = 0.5,
              fill.col = my_pal2
            ) + theme(
              axis.text.x = element_text(
                angle = 45,
                hjust = 1,
                vjust = 1
              )
            )
            print(cellRatio)

            ggsave_fun(filename = paste0(dir, "/3.Sample.ratio.plot"), plot = cellRatio, width = 7 + 0.25 * (length(unique(scRNA$orig.ident))), height = 7)
          }
        )
      }
    }

    # 细胞通讯
    if (step == 4) {
      message("*************************细胞通讯*************************")
      require(CellChat)

      dir <- paste0(Path, "/.Cellchat")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      options(stringsAsFactors = FALSE)
      cellchat <- createCellChat(object = scRNA, group.by = "cellType_1")
      cellchat
      groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

      if (mouse) {
        CellChatDB <- CellChatDB.mouse
      } else {
        CellChatDB <- CellChatDB.human
      }

      CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
      cellchat@DB <- CellChatDB.use

      cellchat <- subsetData(cellchat)
      # future::plan("multisession", workers = 4)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)
      cellchat <- projectData(cellchat, PPI.human)

      cellchat <- computeCommunProb(cellchat,seed.use = 123456)
      cellchat <- filterCommunication(cellchat)
      df.net <- subsetCommunication(cellchat)
      write.table(df.net, file = paste0(dir, "/net_lr.txt"), quote = F, sep = "\t", row.names = F)

      cellchat <- computeCommunProbPathway(cellchat)
      df.netp <- subsetCommunication(cellchat)
      write.table(df.netp, file = paste0(dir, "/net_pathway.txt"), quote = F, sep = "\t", row.names = F)

      ##### -----------show----------#####
      cellchat <- aggregateNet(cellchat)
      saveRDS(cellchat,paste0(dir,"/cellchat.rds"))
      groupSize <- as.numeric(table(cellchat@idents))
      
      # p1 <- netVisual_circle(cellchat@net$count,
      #   vertex.weight = groupSize, weight.scale = T,
      #   label.edge = F, title.name = "Number of interactions"
      # )
      # 
      # p2 <- netVisual_circle(cellchat@net$weight,
      #   vertex.weight = groupSize, weight.scale = T,
      #   label.edge = F, title.name = "Interaction weigths/strength"
      # )
      # # pdf(file="1.Net_number_strength.pdf",width = 10,height = 5)
      # p <- CombinePlots(plots = list(p1, p2), ncol = 2)
      # ggsave_fun(paste0(dir, "/1.Net_number_strength"), plot = p, width = 14)
      pdf(file=paste0(dir, "/1.Net_number_strength.pdf"),width = 10,height = 5)
      par(mfrow = c(1,2), xpd=TRUE)
      netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
      netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
      dev.off()

      df.net <- read.table(paste0(dir, "/net_lr.txt"), sep = "\t", check.names = F, header = T)
      data <- as.data.frame(table(c(
        df.net$source,
        df.net$target
      )))
      colnames(data) <- c("Cell_Type", "all_sum")
      data <- data[order(data$all_sum, decreasing = T), ]
      data$Cell_Type <- factor(data$Cell_Type, levels = data$Cell_Type)
      head(data)
      sample_color <- c("#FB040B", "#F6A717", "#BA06FA", "#172D7A")
      sample_color2 <- ifelse(data$all_sum > summary(data$all_sum)[2], ifelse(data$all_sum > summary(data$all_sum)[3], ifelse(data$all_sum > summary(data$all_sum)[5], sample_color[1], sample_color[2]), sample_color[3]), sample_color[4])
      sample_color1 <- ifelse(data$all_sum > summary(data$all_sum)[2], ifelse(data$all_sum > summary(data$all_sum)[3], ifelse(data$all_sum > summary(data$all_sum)[5], alpha(sample_color[1], 0.9), alpha(sample_color[2], 0.9)), alpha(sample_color[3], 0.9)), alpha(sample_color[4], 0.9))

      pB2 <- ggplot(data = data, aes(x = Cell_Type, y = all_sum, fill = Cell_Type, colour = Cell_Type)) +
        geom_bar(stat = "identity", width = 0.8) +
        scale_fill_manual(values = sample_color1) +
        scale_colour_manual(values = sample_color2) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        labs(x = "Cell Type", y = "Counts") +
        theme(axis.text.y = element_text(size = 12, colour = "black")) +
        theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, colour = "black"))

      ggsave_fun(paste0(dir, "/2.Interaction Count"), plot = pB2, width = 10, height = 8)

      # if (isKeyCell == "yes" && is.na(KeyCell) && 2 %in% KeyCellIndex) {
      #   KeyCell <- as.character(data[1, 1])
      # }
    }

    # 亚群贡献
    if (step == 5) {
      message("*************************亚群贡献*************************")

      dir <- paste0(Path, "/.Subpopulation contribution")
      if (!dir.exists(paste0(dir))) {
        dir.create(paste0(dir))
      } else {
        print("Dir already exists!")
      }

      scRNA$group <- factor(scRNA$group, c(grep("Control|Normal", unique(scRNA$group), value = T), grep("Disease|Tumor|Cancer", unique(scRNA$group), value = T)))
      print(names(table(scRNA$group)))
      if (names(table(scRNA$group))[2] %in% c("Disease","Tumor","Cancer") & names(table(scRNA$group))[1]  %in% c("Control","Normal")){
        Bulk.DEGs <- FindMarkers(scRNA,
                                 group.by = "group", min.pct = 0.25, logfc.threshold = logFCfilter,
                                 ident.1 = names(table(scRNA$group))[2], ident.2 = names(table(scRNA$group))[1]
        )
        Bulk.DEGs <- Bulk.DEGs[Bulk.DEGs$p_val_adj < adjPvalFilter, ]
        Bulk.DEGs <- Bulk.DEGs[Bulk.DEGs$avg_log2FC > 0,]
        print(nrow(Bulk.DEGs))
        Bulk.DEGs$symbol <- rownames(Bulk.DEGs)
        # Bulk.DEGs <- arrange(Bulk.DEGs, Bulk.DEGs$p_val_adj, -Bulk.DEGs$avg_log2FC)
        # Bulk.DEGs <- Bulk.DEGs[1:100, ]
        subset.DEGs <- lapply(SplitObject(scRNA, split.by = "cellType_1"), function(subset) {
          deg <- FindMarkers(subset,
                             features = Bulk.DEGs$symbol, min.pct = 0, logfc.threshold = 0,
                             group.by = "group", ident.1 = names(table(scRNA$group))[2], ident.2 = names(table(scRNA$group))[1]
          )
          deg$symbol <- rownames(deg)
          deg$celltype <- unique(subset$cellType_1)
          return(deg)
        })
        subset.DEGs <- do.call(rbind, subset.DEGs)
        subset.DEGs$FCexp <- 2^subset.DEGs$avg_log2FC
        subset.DEGs$FCprop <- subset.DEGs$pct.1 / subset.DEGs$pct.2
        subset.DEGs$FCscore <- sqrt(subset.DEGs$FCexp * subset.DEGs$FCprop)
        subset.DEGs$FCscore[is.infinite(subset.DEGs$FCscore)] <- NA
        FCscore <- dcast(subset.DEGs, symbol ~ celltype, measure.var = "FCscore")
        write.table(FCscore, paste0(dir, "/output_FCscore.txt"),
                    sep = "\t", row.names = F, col.names = T, quote = F
        )
        
        
        plot.data <- data.frame(
          "subtypes" = colnames(FCscore[, -1]),
          "FCscore" = colMeans(FCscore[, -1], na.rm = T)
        )
        plot.data <- arrange(plot.data, plot.data$FCscore)
        plot.data$subtypes <- factor(plot.data$subtypes, levels = plot.data$subtypes)
        my_pal2 <- c(
          "#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
          "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00"
        )
        p <- ggplot(plot.data, aes(x = subtypes, y = FCscore, fill = subtypes)) +
          geom_bar(stat = "identity") +
          scale_fill_manual(values = my_pal2[1:nrow(plot.data)]) +
          geom_hline(
            yintercept = median(plot.data$FCscore),
            color = "grey", linetype = "dashed"
          ) +
          coord_polar() +
          theme_classic() +
          theme(
            axis.text = element_blank(), axis.title = element_blank(),
            axis.line = element_blank(), axis.ticks = element_blank()
          )
        # p
        # ggsave_fun(filename = paste0(dir, "/1.ContributionScore.pdf"), plot = p, width = 6, height = 4)
        
        plot.data <- plot.data[order(plot.data$FCscore, decreasing = T), ]
        plot.data$subtypes <- factor(plot.data$subtypes, levels = plot.data$subtypes)
        sample_color <- c("#FB040B", "#F6A717", "#BA06FA", "#172D7A")
        sample_color2 <- ifelse(plot.data$FCscore > summary(plot.data$FCscore)[2], ifelse(plot.data$FCscore > summary(plot.data$FCscore)[3], ifelse(plot.data$FCscore > summary(plot.data$FCscore)[5], sample_color[1], sample_color[2]), sample_color[3]), sample_color[4])
        sample_color1 <- ifelse(plot.data$FCscore > summary(plot.data$FCscore)[2], ifelse(plot.data$FCscore > summary(plot.data$FCscore)[3], ifelse(plot.data$FCscore > summary(plot.data$FCscore)[5], alpha(sample_color[1], 0.9), alpha(sample_color[2], 0.9)), alpha(sample_color[3], 0.9)), alpha(sample_color[4], 0.9))
        
        pB2 <- ggplot(data = plot.data, aes(x = subtypes, y = FCscore, fill = subtypes, colour = subtypes)) +
          geom_bar(stat = "identity", width = 0.8) +
          scale_fill_manual(values = sample_color1) +
          scale_colour_manual(values = sample_color2) +
          theme_bw() +
          theme(panel.grid = element_blank()) +
          labs(x = "Cell Type", y = "Counts") +
          theme(axis.text.y = element_text(size = 12, colour = "black")) +
          theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, colour = "black"))
        ggsave_fun(filename = paste0(dir, "/1.CombinedContribution"), plot = p + pB2, width = 12, height = 6)
        
        write.csv(plot.data, paste0(dir, "/plot.data.csv"), row.names = F)
      }
      # if (isKeyCell == "yes" && is.na(KeyCell) && 3 %in% KeyCellIndex) {
      #   KeyCell <- as.character(plot.data[1, 1])
      # }
    }

    # 拟时序分析
    if (step == 6) {
      message("*************************Pseudo-time*************************")
      require(monocle)

      dir <- paste0(Path, "/.Pseudo-time")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      scRNASub <- subset(scRNA, subset = cellType_1 == KeyCell)
      # scRNASub <- scRNA

      monocle.matrix <- GetAssayData(object = scRNASub, slot = "data", assay = "RNA")
      monocle.sample <- scRNASub@meta.data
      monocle.geneAnn <- data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
      
      data <- as(as.matrix(monocle.matrix), "sparseMatrix")
      pd <- new("AnnotatedDataFrame", data = monocle.sample)
      fd <- new("AnnotatedDataFrame", data = monocle.geneAnn)
      cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
      if ("scissor" %in% names(pData(cds))) {
        names(pData(cds))[names(pData(cds)) == "group"] <- "Cluster"
      } else {
        names(pData(cds))[names(pData(cds)) == "orig.ident"] <- "Cluster"
      }
      pData(cds)[, "Cluster"] <- paste0("cluster", pData(cds)[, "Cluster"])
      
      cds <- estimateSizeFactors(cds)
      cds <- estimateDispersions(cds)
      
      disp_table <- dispersionTable(cds)
      disp_table <- arrange(disp_table, -dispersion_empirical)
      track_gene <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
      
      cds <- setOrderingFilter(cds, track_gene)
      cds <- reduceDimension(cds, max_components = 2, reduction_method = "DDRTree")
      cds <- orderCells(cds)
      
      p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
      ggsave(paste0(dir, "/1.Trajectory.Pseudotime.pdf"), plot = p)
      
      p <- plot_cell_trajectory(cds, color_by = "cellType_1")
      ggsave(paste0(dir, "/2.Trajectory.cellType_2.pdf"), plot = p)
      
      p <- plot_cell_trajectory(cds, color_by = "Cluster")
      ggsave(paste0(dir, "/3.Trajectory.Cluster.pdf"), plot = p)
      
      save(cds, file = paste0(dir, "/cds.rda"))
      
      pseudo_time_diff <- differentialGeneTest(cds, cores = 4, fullModelFormulaStr = "~sm.ns(Pseudotime)")
      pseudo_time_diff <- arrange(pseudo_time_diff, qval)
      write.csv(pseudo_time_diff,file = paste0(dir,'/heatmap_pseudo_time_gene.csv'))
      pseudo_time_genes <- rownames(pseudo_time_diff)
      plot_heatmap_gene <- pseudo_time_genes[1:50]
      len <- length(plot_heatmap_gene)
      p <- monocle::plot_pseudotime_heatmap(cds[plot_heatmap_gene, ], num_clusters = 3, show_rownames = T, return_heatmap = T, hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100))
      ggsave(filename = paste0(dir, "/4.Heatmap.pdf"), plot = ggplotify::as.ggplot(p), width = 5, height = len * 0.12)
      if (!is.null(dev.list())) {
        dev.off()
      }
      
      # clusters <- cutree(p$tree_row, k = 3)
      # clustering <- data.frame(clusters)
      # clustering[,1] <- as.character(clustering[,1])
      # colnames(clustering) <- "Gene_Clusters"
      # clustering <- arrange(clustering,Gene_Clusters)
      # write.csv(clustering,file = paste0(dir,'/','heatmap_cluster_gene.csv'))
      
      if (!class(KeyGene) == "logical") {
        cds_subset <- cds[KeyGene, ]
        p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "Cluster")
        p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "cellType_2")
        p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
        p1 | p2 | p3
        ggsave(filename = paste0(dir, "/5.KeyGene.pseudotime.pdf"), plot = p1 | p2 | p3, width = 12, height = length(KeyGene) * 2)
      }
    }

    # 关键基因表达丰度
    if (step == 7) {
      message("*************************关键基因表达丰度*************************")

      dir <- paste0(Path, "/.KeyGene expression abundance")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      KeyGene <- intersect(rownames(scRNA), KeyGene)
      leg <- length(unique(KeyGene))
      p <- FeatureCornerAxes(
        object = scRNA, reduction = "umap",
        groupFacet = NULL,
        relLength = 0.5,
        relDist = 0.2,
        features = unique(KeyGene), nLayout = 3,
        aspect.ratio = 1, pSize = 0.1
      )
      if (leg <= 3) {
        ggsave_fun(filename = paste0(dir, "/1.Keygene.FeaturePlot"), plot = p, width = leg * 5, height = 5)
      }
      if (leg > 3) {
        ggsave_fun(filename = paste0(dir, "/1.Keygene.FeaturePlot"), plot = p, width = 15, height = ceiling(leg / 3) * 5)
      }

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
      if (leg <= 3) {
        ggsave_fun(filename = paste0(dir, "/2.Keygene.DotPlot"), plot = p, width = 5 + ceiling(leg / 3), height = 7)
      }
      if (leg > 3) {
        ggsave_fun(filename = paste0(dir, "/2.Keygene.DotPlot"), plot = p, width = 5 + 0.5 * leg, height = 7)
      }
    }

    # 疾病基因共表达
    if (step == 8) {
      message("*************************疾病基因共表达网络*************************")

      dir <- paste0(Path, "/.Disease gene co-expression")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      showGenes <- intersect(rownames(scRNA), KeyGene)
      cat("\033[0;34m请复制疾病相关基因：\n")
      geneCard <- scan(sep = "\n", what = "character", quiet = FALSE)
      geneCard <- intersect(geneCard,rownames(scRNA))[1:6]

      for (i in showGenes) {
        PlotNum <- 0
        path <- paste0(dir, "/", which(showGenes == i), ".", i)
        if (!dir.exists(path)) {
          dir.create(path)
        } else {
          print("Dir already exists!")
        }
        for (j in geneCard) {
          PlotNum <- PlotNum + 1
          p1 <- FeaturePlot(scRNA,
            features = c(i, j),
            blend = TRUE, cols = c("gray80", "red", "green"),
            pt.size = 0.5, raster = F, reduction = "umap"
          ) +
            theme(aspect.ratio = 1)
          p2 <- getScatterplot(scRNA,
            gene1 = j, gene2 = i,
            jitter.num = 0.15, pos = TRUE
          ) +
            theme(aspect.ratio = 1)
          p <- CombinePlots(plots = list(p1, p2), ncol = 2, rel_widths = c(4, 1))
          ggsave_fun(paste0(path, "/", PlotNum, ".", j, " ~ ", i), plot = p, width = 15, height = 4)
        }
      }
    }

    # 免疫代谢通路
    if (step == 9) {
      message("*************************免疫代谢通路*************************")

      dir <- paste0(Path, "/.Immunometabolic pathways")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      ShowGen <- intersect(rownames(scRNA), KeyGene)
      Idents(scRNA) <- "cellType_1"
      geneset <- clusterProfiler::read.gmt("refdata/Homo/h.all.v7.5.1.symbols.gmt")
      if (mouse) {
        geneset <- clusterProfiler::read.gmt("refdata/Mus/mh.all.v2023.2.Mm.symbols.gmt")
      }
      geneset <- split(geneset$gene, geneset$term)
      genesetInfo <- read.delim("refdata/GenesetInfo.txt", sep = ",")
      genesetInfo <- subset(genesetInfo, Classification != "")
      levels <- c("Immune", "Metabolism", "Signaling", "Proliferation")
      genesetInfo$Classification <- factor(genesetInfo$Classification, levels)
      genesetInfo <- arrange(genesetInfo, genesetInfo$Classification, genesetInfo$geneset)
      genesetInfo$geneset <- factor(genesetInfo$geneset, levels = genesetInfo$geneset)
      mat <- GetAssayData(object = scRNA, assay = "RNA", slot = "data")
      # rm(scRNA.plot2)
      cells_rankings <- AUCell::AUCell_buildRankings(mat, nCores = 1, plotStats = F)
      score <- AUCell::AUCell_calcAUC(geneset, cells_rankings,
        nCores = 1,
        aucMaxRank = nrow(cells_rankings) * 0.05
      )
      score <- AUCell::getAUC(score) # score@assays@data$AUC

      # compare <<- "Hexp-Lexp"
      # adjust.method <<- "bonferroni"
      DP.list <- base::lapply(ShowGen, function(cellType) {
        # cellType <- ShowGen[1]
        # print(cellType)
        # print(1)
        DatGroup <- FetchData(scRNA, vars = c(cellType), slot = "data")
        # print(2)
        group <- setNames(
          object = ifelse(DatGroup[, 1] > median(DatGroup[, 1]), "Hexp", "Lexp"),
          nm = rownames(DatGroup)
        )
        # print(3)
        if (length(unique(group)) > 1) {
          # print(4)
          design <- model.matrix(~ 0 + factor(group))
          # print(5)
          colnames(design) <- levels(factor(group))
          # print(6)
          rownames(design) <- names(group)
          # print(7)
          contrast.matrix <- limma::makeContrasts("Hexp-Lexp", levels = design) # should be Test-Control
          # print(8)
          fit <- limma::lmFit(score[, names(group)], design)
          # print(9)
          fit2 <- limma::contrasts.fit(fit, contrast.matrix)
          # print(10)
          fit2 <- limma::eBayes(fit2)
          # print(11)
          DPs <- limma::topTable(fit2, coef = 1, n = Inf, adjust.method = "bonferroni")
          # print(12)
          DPs$KeyGenes <- cellType
          # print(13)
          DPs$Pathway <- rownames(DPs)
          return(DPs)
        }
      })
      DPs <- do.call(rbind, DP.list)
      write.table(DPs, file = paste0(dir, "/output_hallmark.txt"), sep = "\t", row.names = F, col.names = T, quote = F)

      plot.data <- DPs
      plot.data <- subset(plot.data, Pathway %in% genesetInfo$geneset)
      plot.data$KeyGenes <- factor(plot.data$KeyGenes)
      plot.data$Pathway <- factor(plot.data$Pathway, levels = genesetInfo$geneset)
      plot.data$FDR <- cut(plot.data$adj.P.Val,
        breaks = c(0, 1e-125, 1e-75, 1e-25, 1),
        include.lowest = T
      )
      plot.data$FDR <- factor(as.character(plot.data$FDR),
        levels = rev(levels(plot.data$FDR))
      )
      levels(plot.data$Pathway) <- tolower(gsub("HALLMARK_", "", levels(plot.data$Pathway)))

      color <- c("#4682B4", "#FFFFFF", "#CD2626")
      class.color <- c(
        "Immune" = "#D58986", "Metabolism" = "#80554C",
        "Signaling" = "#71AC7A", "Proliferation" = "#E8D4B4"
      )
      p1 <- ggplot(plot.data, aes(x = Pathway, y = KeyGenes, color = logFC, size = FDR)) +
        geom_point() +
        scale_color_gradient2(low = color[1], mid = color[2], high = color[3]) +
        geom_hline(
          yintercept = seq(
            min(as.numeric(plot.data$KeyGenes)) - 0.5,
            max(as.numeric(plot.data$KeyGenes)) + 0.5
          ),
          color = "grey80"
        ) +
        geom_vline(
          xintercept = seq(
            min(as.numeric(plot.data$Pathway)) - 0.5,
            max(as.numeric(plot.data$Pathway)) + 0.5
          ),
          color = "grey80"
        ) +
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.line = element_blank(),
          legend.position = "top",
          legend.key.height = unit(0.3, "cm"),
          legend.key.width = unit(0.8, "cm"),
        )

      p2 <- ggplot(genesetInfo, aes(x = geneset, y = 1, fill = Classification)) +
        geom_tile() +
        theme_classic() +
        scale_fill_manual(values = class.color) +
        theme(
          axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), axis.line = element_blank(), legend.position = "bottom"
        )
      p <- cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(10, 2))
      ggsave_fun(filename = paste0(dir, "/1.Hallmark"), width = 15, height = 9, plot = p)
    }

    # 代谢通路热图
    if (step == 10) {
      message("*************************代谢通路热图*************************")

      require(GSVA)
      require(GSEABase)
      require(limma)

      dir <- paste0(Path, "/.Metabolic pathway")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      gmtFile <- "refdata/Homo/immune.gmt"
      cat("\033[0;34m请输入分组名字，例：group\n")
      Metagroup <- scan(what = "character")
      rt <- scRNA@assays$RNA@data
      group <- data.frame(ID = rownames(scRNA@meta.data), Group = scRNA@meta.data[[Metagroup]])
      exp <- as.matrix(rt)
      dimnames <- list(rownames(exp), colnames(exp))
      mat <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
      mat <- avereps(mat)
      mat <- mat[rowMeans(mat) > 0, ]
      geneSet <- getGmt(gmtFile,
        geneIdType = SymbolIdentifier()
      )

      ssgseaScore <- gsva(mat, geneSet, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)
      normalize <- function(x) {
        return((x - min(x)) / (max(x) - min(x)))
      }
      ssgseaOut <- normalize(ssgseaScore)
      ssgseaOut <- rbind(id = colnames(ssgseaOut), ssgseaOut)
      write.table(ssgseaOut, file = paste0(dir, "/ssgseaOut.txt"), sep = "\t", quote = F, col.names = F)

      plotdata <- read.table(paste0(dir, "/ssgseaOut.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
      up <- group
      up <- up[order(up$Group), ]
      plotdata <- plotdata[, up$ID]
      up$group <- up$Group
      down <- read.table("refdata/Homo/type.txt", sep = "\t", header = T, check.names = F, stringsAsFactors = F)
      compath <- intersect(rownames(plotdata), down$pathway)
      down <- down[down$pathway %in% compath, ]
      plotdata <- plotdata[down$pathway, ]
      # down[,2][down[,2] == "5C3 specific metabolism signatures"] <- c("4Other metabolism signatures")##
      annCol <- data.frame(
        Group = up$group,
        row.names = up$ID,
        stringsAsFactors = F
      )

      annRow <- data.frame(
        Direct = down$type,
        row.names = down$pathway,
        stringsAsFactors = F
      )

      annGroup <- c("#E7B800", "#2E9FDF")
      names(annGroup) <- unique(group$Group)

      annColors <- list(
        annGroup,
        "Direct" = c(
          "1Amino acid metabolism relevant signatures" = "#247BA0",
          "2lipid metabolism relevant signatures" = "#70C1B3",
          "3Drug metabolism relevant signatures" = "#B2DBBF",
          "4Other metabolism signatures" = "#F3FFBD",
          "5C3 specific metabolism signatures" = "#6CD3A7"
        )
      )

      # plotdata <- plotdata[c(rownames(up),rownames(down)),]
      plotdata <- t(scale(t(plotdata)))
      plotdata[plotdata > 1] <- 1
      plotdata[plotdata < -1] <- -1

      pdf(file = paste0(dir, "/1.Heatmap.pdf"), width = 25, height = 16)
      pheatmap(plotdata,
        scale = "none",
        annotation_row = annRow,
        annotation_col = annCol,
        annotation_colors = annColors,
        color = colorRampPalette(c("#009BC7", "#F3F3F1", "#F15E4C"))(20),
        # color = greenred(64),
        fontsize_row = 12,
        fontsize_col = 8,
        fontsize = 12,
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        gaps_row = c(13, 22, 25, 38),
        #         gaps_row = c(13,22,25),
        show_colnames = F
      )
      dev.off()

      # png(file = paste0(dir, "/1.Heatmap.png"), width = 2500, height = 1600, res = 150)
      # pheatmap(plotdata,
      #   scale = "none",
      #   annotation_row = annRow,
      #   annotation_col = annCol,
      #   annotation_colors = annColors,
      #   color = colorRampPalette(c("#009BC7", "#F3F3F1", "#F15E4C"))(20),
      #   # color = greenred(64),
      #   fontsize_row = 12,
      #   fontsize_col = 8,
      #   fontsize = 12,
      #   cluster_cols = FALSE,
      #   cluster_rows = FALSE,
      #   gaps_row = c(13, 22, 25, 38),
      #   #         gaps_row = c(13,22,25),
      #   show_colnames = F
      # )
      # dev.off()
    }

    # 量化热点机制
    if (step == 11) {
      message("*************************量化热点机制*************************")

      require(GSVA)
      require(GSEABase)
      require(limma)

      dir <- paste0(Path, "/.Hotspot mechanism")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      normalize <- function(x) {
        return((x - min(x)) / (max(x) - min(x)))
      }
      mat <- GetAssayData(object = scRNA, assay = "RNA")
      mat <- as(mat[rowMeans(mat) > 0.3, ], "matrix")
      geneSet <- GSEABase::getGmt("ssGSEA.gmt", geneIdType = GSEABase::SymbolIdentifier())

      ssgseaScore <- gsva(mat, geneSet, method = "ssgsea", parallel.sz = 10, kcdf = "Gaussian", abs.ranking = TRUE)
      ssgseaOut <- normalize(ssgseaScore)
      save(ssgseaOut, file = paste0(dir, "/ssgseaOut.rda")) # load("ssgseaOut.rda")
      ssgseaOut <- rbind(id = colnames(ssgseaOut), ssgseaOut)
      write.table(ssgseaOut, file = paste0(dir, "/ssgseaOutAll.txt"), sep = "\t", quote = F, col.names = F)

      load(paste0(dir, "/ssgseaOut.rda"))
      ssgseaOut <- as.data.frame(t(ssgseaOut))
      ScoGroup <- ifelse(ssgseaOut[rownames(scRNA[["orig.ident"]]), names(geneSet)] > median(ssgseaOut[rownames(scRNA[["orig.ident"]]), names(geneSet)]), "Hsco", "Lsco")
      scRNA$ScoGroup <- data.frame(row.names = rownames(scRNA[["orig.ident"]]), ScoGroup)

      load(paste0(dir, "/ssgseaOut.rda"))
      meta <- as.data.frame(t(ssgseaOut))
      scRNA <- AddMetaData(scRNA, meta)
      saveRDS(scRNA, paste0(dir, "/result.rds"))

      tryCatch(
        {
          for (i in rownames(ssgseaOut)) {
            tmp <- scRNA@meta.data[, c("cellType_1", "group", i)]
            colnames(tmp)[3] <- "val"
            p1 <- ggplot(
              tmp,
              aes(
                x = cellType_1, y = val,
                fill = group,
                color = group
              )
            ) +
              geom_boxplot(
                notch = F, alpha = 0.95,
                outlier.shape = 16,
                outlier.size = 0.65
              ) +
              xlab("") +
              ylab("") +
              scale_fill_manual(values = c("#D5EBFB", "#FBEEB7", "#B4FBCD", "#F5B3FC")) +
              scale_color_manual(values = c("#0073C2", "#EFC000", "#00C244", "#C501D7")) +
              ggtitle("") +
              theme_classic() +
              theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                axis.text.y = element_text(vjust = 0.5, size = 12),
                axis.title.y = element_text(angle = 90, size = 15)
              ) +
              theme(legend.position = "top") +
              stat_compare_means(method = "wilcox.test", hjust = 0.5, vjust = 0, hide.ns = T, label = "p.signif", show.legend = F)
            p2 <- FeaturePlot(object = scRNA, ncol = 1, features = i, cols = c("green", "red"))
            p <- CombinePlots(plots = list(p2, p1), rel_widths = c(1, 1.6))
            ggsave_fun(paste0(dir, "/1.", i), width = 12, height = 4, plot = p)
          }
        },
        error = function(e) {
          for (i in rownames(ssgseaOut)) {
            tmp <- scRNA@meta.data[, c("cellType_1", i)]
            colnames(tmp)[2] <- "val"
            p1 <- ggplot(
              tmp,
              aes(
                x = cellType_1, y = val,
                fill = cellType_1
              )
            ) +
              geom_boxplot(
                notch = F, alpha = 0.95,
                outlier.shape = 16,
                outlier.size = 0.65
              ) +
              xlab("") +
              ylab("") +
              ggtitle("") +
              theme_classic() +
              theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                axis.text.y = element_text(vjust = 0.5, size = 12),
                axis.title.y = element_text(angle = 90, size = 15)
              ) +
              theme(legend.position = "top")
            p2 <- FeaturePlot(object = scRNA, ncol = 1, features = i, cols = c("green", "red"))
            p <- CombinePlots(plots = list(p2, p1), rel_widths = c(1, 1.6))
            ggsave_fun(paste0(dir, "/1.", i), width = 12, height = 4, plot = p)
          }
        }
      )

      message("是否挑选评分最显著的细胞做火山图？yes/no")
      isSlect <- scan(what = "character", quiet = T)
      if (isSlect == "yes") {
        cat("\033[0;34m请输入细胞名字\n")
        sig.Cell <- scan(what = "character", quiet = T)
        pbmc <- subset(scRNA, subset = cellType_1 %in% sig.Cell)
      }
      if (isSlect == "no") {
        pbmc <- scRNA
      }

      PbmcMarkers <- FindMarkers(pbmc,
        ident.1 = "Hsco",
        min.pct = 0.1,
        group.by = "ScoGroup",
        assay = "RNA",
        logfc.threshold = 0.1
      )

      write.csv(PbmcMarkers, file = paste0(dir, "/DiffGene.csv"), row.names = T, quote = F)

      PbmcMarkers <- read.csv(paste0(dir, "/DiffGene.csv"), header = T, row.names = 1, check.names = F)

      require(ggrepel)
      cols <- c("#DC143C", "#00008B", "#808080")
      names(cols) <- c("Up", "Down", "NoSignifi")

      PbmcMarkers <- PbmcMarkers[-log10(PbmcMarkers$p_val_adj) > 0, ]
      PbmcMarkers$gene <- rownames(PbmcMarkers)
      PbmcMarkers$lab <- ""
      PbmcMarkers[order(PbmcMarkers$avg_log2FC), ][c(1:10, (nrow(PbmcMarkers) - 9):nrow(PbmcMarkers)), ]$lab <- PbmcMarkers[order(PbmcMarkers$avg_log2FC), ][c(1:10, (nrow(PbmcMarkers) - 9):nrow(PbmcMarkers)), ]$gene
      PbmcMarkers$logP <- -log10(PbmcMarkers$p_val_adj)

      p <- ggplot(PbmcMarkers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = avg_log2FC)) +
        geom_point(aes(size = logP), alpha = 0.9) +
        scale_color_gradient2(
          low = "#0500FF", high = "#FF0000", mid = "#EFEFEF",
          midpoint = 0, space = "Lab", # , limit = c(-1, 1)
          name = "logFC"
        ) +
        scale_size("-log10(p)") +
        geom_text_repel(
          data = PbmcMarkers,
          aes(x = avg_log2FC, y = -log10(p_val_adj), label = lab),
          size = 4, box.padding = unit(0.5, "lines"),
          point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = F
        ) +
        theme_bw() +
        ylab("-log10 (p_val_adj)") +
        xlab("avg_log2FC") +
        geom_vline(xintercept = c(-logFCfilter, logFCfilter), lty = 3, col = "black", lwd = 0.5) +
        geom_hline(yintercept = -log10(adjPvalFilter), lty = 3, col = "black", lwd = 0.5)
      p
      ggsave_fun(filename = paste0(dir, "/2.VlnPlot"), plot = p, width = 10, height = 8)
    }
    
    # 富集分析
    if (step == 12) {
      message("*************************富集分析*************************")
      
      dir <- paste0(Path, "/.Enrichment analysis")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      require(enrichplot)
      require(clusterProfiler)
      # cat("\033[0;34m请复制基因：\n")
      # genes <- scan(sep = "\n", what = "character", quiet = FALSE)
      genes <- KeyGene
      go <<- kk <<- NULL
      if(mouse){
        require(org.Mm.eg.db)
        entrezIDs <- mget(genes, org.Mm.egSYMBOL2EG, ifnotfound = NA)
        entrezIDs <- as.character(entrezIDs)
        out <- cbind(genes,entrezID=entrezIDs)
        write.table(out, file = paste0(dir,"/entrezIDs.txt"), sep = "\t", quote = F, row.names = F)
        
        go <<- enrichGO(
          gene = entrezIDs,
          OrgDb = org.Mm.eg.db,
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.05,
          ont = "all",
          readable = T
        )
        write.csv(go, file = paste0(dir, "/GO.csv"), quote = F, row.names = F)
        
        kk <<- enrichKEGG(gene = entrezIDs, organism = "mmu", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
        write.csv(kk, file = paste0(dir, "/KEGG.csv"), quote = F, row.names = F)
        
      }else{
        require(org.Hs.eg.db)
        entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
        entrezIDs <- as.character(entrezIDs)
        out <- cbind(genes,entrezID=entrezIDs)
        write.table(out, file = paste0(dir,"/entrezIDs.txt"), sep = "\t", quote = F, row.names = F)
        go <<- enrichGO(gene = entrezIDs,
                        OrgDb = org.Hs.eg.db, 
                        pvalueCutoff =0.05, 
                        qvalueCutoff = 0.05,
                        ont="all",
                        readable =T)
        write.csv(go, file = paste0(dir, "/GO.csv"), quote = F, row.names = F)
        
        kk <<- enrichKEGG(gene = entrezIDs, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
        write.csv(kk, file = paste0(dir, "/KEGG.csv"), quote = F, row.names = F)
      }
      # view GO
      if(nrow(go@result) > 0){
        go.bar <- barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY",label_format=100) + facet_grid(ONTOLOGY~., scale='free')
        ggsave_fun(filename = paste0(dir,"/1.GOBarplot"),plot = go.bar,width = 12,height = 15)
        
        go.dot <- dotplot(go,showCategory = 10,split="ONTOLOGY",label_format=100) + facet_grid(ONTOLOGY~., scale='free')
        ggsave_fun(filename = paste0(dir,"/2.GODotplot"),plot = go.dot,width = 12,height = 15)
      }
      if(nrow(kk@result) > 0){
        if(mouse){
          kk@result$Description <- gsub(" - Mus musculus \\(house mouse\\)","",kk@result$Description)
        }
        kk.bar <- barplot(kk, drop = TRUE, showCategory = 30,label_format=100)
        ggsave_fun(filename = paste0(dir,"/3.KEGGBarplot"),plot = kk.bar,width = 12,height = 15)
        kk.dot <- dotplot(kk, showCategory = 30,label_format=100)
        ggsave_fun(filename = paste0(dir,"/4.KEGGBarplot"),plot = kk.dot,width = 12,height = 15)
      }
    }
    
  }
  return(scRNA)
}
