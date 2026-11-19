Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(123456)

ggsave_fun <- function(filename, plot = last_plot(), plot_Device = c(".pdf"), width = 7, height = 7, dpi = 600, ...) {
  filenames <- paste0(filename, plot_Device)
  for (filename_tmp in filenames) {
    ggsave(filename = filename_tmp, plot = plot, width = width, height = height, dpi = dpi, ...)
  }
}

RePlot <- function(scRNA = scRNA,Ident,Multi = TRUE) {
  require(Seurat)
  require(ggplot2)
  require(dplyr)
  require(reshape2)
  require(scRNAtoolVis)
  require(ggpubr)
  n <- 0
  for (step in 1:3) {
    n <- n + 1

    # 质量控制
    if (step == 1) {
      message("*************************质量控制*************************")

      dir <- paste0(n, ".QCFilter")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      p <- VlnPlot(scRNA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), group.by = "orig.ident", ncol = 4)
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

    # 数据标准化
    if (step == 2) {
      message("*************************数据标准化*************************")

      dir <- paste0(n, ".Cluster")

      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      assay <- "RNA"
      set.resolutions <- seq(0.2, 1.2, by = 0.1)
      PC <- 20
      nfeatures <- 2000
      npcs <- 50
      
      top10 <- head(VariableFeatures(scRNA), 10)
      plot1 <- VariableFeaturePlot(scRNA)
      plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
      plot1 + plot2
      ggsave_fun(filename = paste0(dir, "/1.VariableFeaturePlot"), plot = plot1 + plot2, width = 14)
      p <- ElbowPlot(object = scRNA, ndims = npcs)
      ggsave_fun(filename = paste0(dir, "/2.ElbowPlot"), plot = p)

      if (Multi) {
        p <- DimPlot(object = scRNA, reduction = "pca", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/3.PcaPlot"), plot = p)
        p <- DimPlot(object = scRNA, reduction = "harmony", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/4.HarmonyPlot"), plot = p)
      } else {
        p <- DimPlot(object = scRNA, reduction = "pca", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/3.PcaPlot"), plot = p)
      }
    }

    # 自动化注释
    if (step == 3) {
      message("*************************自动化注释*************************")
      
      dir <- paste0(n, ".CellAnnotate")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }


      Cluster.dir <- list.files(pattern = ".Cluster$", all.files = TRUE)

      p <- DimPlot(scRNA, group.by = Ident, label = T, reduction = "umap") + NoLegend()
      print(p)
      ggsave_fun(filename = paste0(Cluster.dir, "/5.Umap.Cluster"), plot = p)


      # Idents(scRNA) <- "cellType_1"
      # scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      # write.csv(scRNA.markers, paste0(dir, "/Findall.markers.cellType_1.csv"))
      
      scRNA.markers <- read.csv("Findall.markers.cellType_1.csv",header = T,row.names = 1)

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
      # ggsave_fun(filename = paste0(dir, "/1.Featureplot.umap"), plot = p, width = 12, height = 12)

      my_pal2 <- c(
        "#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
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
}


# scRNA <- readRDS("result.rds")
# RePlot(scRNA = scRNA,Ident = "RNA_snn_res.0.2",Multi = TRUE)
