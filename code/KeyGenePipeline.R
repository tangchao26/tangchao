require(limma)
require(clusterProfiler)
require(ggplot2)
require(ggrepel)
require(cowplot)
require(ggraph)
require(GSVA)
require(corrplot)
require(dplyr)
require(Hmisc)
require(ggpubr)
require(patchwork)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(123456)

ggsave_fun <- function(filename, plot = last_plot(), plot_Device = c(".pdf"), width = 7, height = 7, dpi = 600, ...) {
  filenames <- paste0(filename, plot_Device)
  for (filename_tmp in filenames) {
    ggsave(filename = filename_tmp, plot = plot, width = width, height = height, dpi = dpi, ...)
  }
}

human2mouse <- function(x) {
  x1 <- substr(x, 1, 1)
  x2 <- substr(x, 2, nchar(x))
  return(paste0(x1, tolower(x2)))
}

gseaplot3 <- function (x, geneSetID, title = "", color = "green", base_size = 11, 
                       rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE, 
                       ES_geom = "line"){
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  geneList <- position <- NULL
  if (length(geneSetID) == 1) {
    gsdata <- enrichplot:::gsInfo(x, geneSetID)
  }
  else {
    gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
  }
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0))
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
                          size = 1)
  }
  else {
    es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description), 
                           size = 1, data = subset(gsdata, position == 1))
  }
  
  p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8),  
                                legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
  
  p.res <- p.res + ylab("Running Enrichment Score") + theme(axis.text.x = element_blank(), 
                                                            axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                                                            plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, 
                                                                                 unit = "cm"))
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == 
                   term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, 
                                                           ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
    theme_classic(base_size) + theme(legend.position = "none", 
                                     plot.margin = margin(t = -0.1, b = 0, unit = "cm"), axis.ticks = element_blank(), 
                                     axis.text = element_blank(), axis.line.x = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 
                                                                         0))
  if (length(geneSetID) == 1) {
    v <- seq(1, sum(gsdata$position), length.out = 9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) 
      inv <- inv + 1
    col <- c(rev(RColorBrewer::brewer.pal(5, "Blues")), RColorBrewer::brewer.pal(5, "Reds"))
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
                    xmax = xmax, col = col[unique(inv)])
    p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax, 
                              ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d, 
                         alpha = 0.9, inherit.aes = FALSE)
  }
  df2 <- p$data
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x, 
                                             y = ~y, yend = 0), color = "grey")
  p.pos <- p.pos + ylab("Ranked List Metric") + xlab("Rank in Ordered Dataset") + 
    theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, 
                               l = 0.2, unit = "cm"))
  if (!is.null(title) && !is.na(title) && title != "") {
    title <- unlist(strsplit(title, split = "_"))
    # title2 <- paste(str_to_title(title[2:length(title)]), collapse = " ")
    # title <- paste(title[1],title2,collapse = " ") %>% stringr::str_wrap(., width = 40)
    title <- paste(title,collapse = " ") %>% stringr::str_wrap(., width = 40)
    p.res <- p.res + ggtitle(title)
  }
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values = color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    }
    else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "p.adjust",'enrichmentScore','NES')]
    rownames(pd) <- ''
    pd <- pd[, -1]
    colnames(pd) <- c( "pvalue", "p.adjust",'ES','NES')
    pd <- round(pd, 4)
    #tp <- enrichplot:::tableGrob2(pd, p.res)
    tp <- enrichplot:::tableGrob2(pd)
    p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp, 
                                                                         xmin = quantile(p.res$data$x, 0.5), 
                                                                         xmax = quantile(p.res$data$x, 0.95), 
                                                                         ymin = quantile(p.res$data$runningScore,0.75), 
                                                                         ymax = quantile(p.res$data$runningScore,0.9))
  }
  
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
                                         axis.ticks.x = element_line(), axis.text.x = element_text())
  if (length(subplots) == 1) 
    return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
                                                      r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
  if (length(rel_heights) > length(subplots)) 
    rel_heights <- rel_heights[subplots]
  plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
}

gsea_dotplot <- function(object,x = "geneRatio",color = "p.adjust",showCategory = 10,
                         size = NULL,split = T,title = "",
                         col_low = "#D15047", col_high =  "#3793DE", scales = 'free_x',
                         orderBy = "x",decreasing = TRUE,y_wrap = 40) {
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
    if (is.null(size)) 
      size <- "Count"
  }else if (x == "count" || x == "Count") {
    x <- "Count"
    if (is.null(size)) 
      size <- "GeneRatio"
  }else if (is(x, "formula")) {
    x <- as.character(x)[2]
    if (is.null(size)) 
      size <- "Count"
  }else {
    if (is.null(size)) 
      size <- "Count"
  }
  
  df <- fortify(object, showCategory = showCategory, split = ".sign")
  if (orderBy != "x" && !orderBy %in% colnames(df)) {
    message("wrong orderBy parameter; set to default `orderBy = \"x\"`")
    orderBy <- "x"
  }
  if (orderBy == "x") {
    df <- dplyr::mutate(df, x = eval(parse(text = x)))
  }
  
  label_func <- function(str,label_format=y_wrap){str <- gsub("_", " ", str);stringr::str_wrap(str, label_format)}
  
  idx <- order(df[[orderBy]], decreasing = T)
  df$Description <- factor(df$Description, levels = rev(unique(df$Description[idx])))
  result <<- df
  p <- ggplot(df, aes_string(x = x, y = "Description", size = size, fill = colorBy)) + 
    geom_point(shape=21) + 
    scale_fill_continuous(low = col_low, high =  col_high, name = color) + 
    scale_y_discrete(labels = label_func) + ylab(NULL) + 
    ggtitle(title) + DOSE::theme_dose(12) + scale_size(range = c(3,9))+
    guides(fill = guide_colorbar(order = 1,reverse = TRUE),
           size = guide_legend(order = 2,override.aes = list(fill='grey40')))+
    theme(text = element_text(face = 'bold',color = 'black'),
          axis.line.x.top  = element_line(size=0.8),
          axis.ticks = element_blank(),
          strip.background = element_rect(color = 'black',size = 1.3),
          panel.border = element_rect(color = 'black',size = 1.3),
          panel.grid.major = element_line(linetype = 5),
          panel.grid.minor = element_blank())
  
  if(split){
    p <- p + facet_wrap(~.sign,scales=scales) + theme(strip.text=element_text(size = 15,face = 'bold'))
  }
  return(p)
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

gather_graph_node <- function(df, index = NULL, value = tail(colnames(df), 1), root = NULL) {
  require(dplyr)
  if (length(index) < 2) {
    stop("please specify at least two index column(s)")
  } else {
    list <- lapply(seq_along(index), function(i) {
      dots <- index[1:i]
      df %>%
        group_by(.dots = dots) %>%
        summarise(
          node.size = sum(.data[[value]]),
          node.level = index[[i]],
          node.count = n()
        ) %>%
        mutate(
          node.short_name = as.character(.data[[dots[[length(dots)]]]]),
          node.branch = as.character(.data[[dots[[1]]]])
        ) %>%
        tidyr::unite(node.name, dots, sep = "/")
    })
    data <- do.call("rbind", list) %>% as_tibble()
    data$node.level <- factor(data$node.level, levels = index)

    if (is.null(root)) {
      return(data)
    } else {
      root_data <- data.frame(
        node.name = root,
        node.size = sum(df[[value]]),
        node.level = root,
        node.count = 1,
        node.short_name = root,
        node.branch = root,
        stringsAsFactors = F
      )
      data <- rbind(root_data, data)
      data$node.level <- factor(data$node.level, levels = c(root, index))
      return(data)
    }
  }
}

gather_graph_edge <- function(df, index = NULL, root = NULL) {
  require(dplyr)
  if (length(index) < 2) {
    stop("please specify at least two index column(s)")
  } else if (length(index) == 2) {
    data <- df %>%
      mutate(from = .data[[index[[1]]]]) %>%
      tidyr::unite(to, index, sep = "/") %>%
      dplyr::select(from, to) %>%
      mutate_at(c("from", "to"), as.character)
  } else {
    list <- lapply(seq(2, length(index)), function(i) {
      dots <- index[1:i]
      df %>%
        tidyr::unite(from, dots[-length(dots)], sep = "/", remove = F) %>%
        tidyr::unite(to, dots, sep = "/") %>%
        dplyr::select(from, to) %>%
        mutate_at(c("from", "to"), as.character)
    })
    data <- do.call("rbind", list)
  }
  data <- as_tibble(data)
  if (is.null(root)) {
    return(data)
  } else {
    root_data <- df %>%
      group_by(.dots = index[[1]]) %>%
      summarise(count = n()) %>%
      mutate(from = root, to = as.character(.data[[index[[1]]]])) %>%
      dplyr::select(from, to)
    rbind(root_data, data)
  }
}

CoreAlg <- function(X, y) {
  # try different values of nu
  svn_itor <- 3

  res <- function(i) {
    if (i == 1) {
      nus <- 0.25
    }
    if (i == 2) {
      nus <- 0.5
    }
    if (i == 3) {
      nus <- 0.75
    }
    model <- svm(X, y, type = "nu-regression", kernel = "linear", nu = nus, scale = F)
    model
  }

  if (Sys.info()["sysname"] == "Windows") {
    out <- mclapply(1:svn_itor, res, mc.cores = 1)
  } else {
    out <- mclapply(1:svn_itor, res, mc.cores = svn_itor)
  }

  nusvm <- rep(0, svn_itor)
  corrv <- rep(0, svn_itor)

  # do cibersort
  t <- 1
  while (t <= svn_itor) {
    weights <- t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights < 0)] <- 0
    w <- weights / sum(weights)
    u <- sweep(X, MARGIN = 2, w, "*")
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }

  # pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]

  # get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q < 0)] <- 0
  w <- (q / sum(q))

  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]

  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
}

doPerm <- function(perm, X, Y) {
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()

  while (itor <= perm) {
    # print(itor)

    # random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist), dim(X)[1])])

    # standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)

    # run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)

    mix_r <- result$mix_r

    # store correlation
    if (itor == 1) {
      dist <- mix_r
    } else {
      dist <- rbind(dist, mix_r)
    }

    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

CIBERSORT <- function(dir, sig_matrix, mixture_file, perm = 0, QN = TRUE) {
  library(e1071)
  library(parallel)
  library(preprocessCore)

  # read in data
  X <- read.table(sig_matrix, header = T, sep = "\t", row.names = 1, check.names = F)
  Y <- mixture_file

  X <- data.matrix(X)
  Y <- data.matrix(Y)

  # order
  X <- X[order(rownames(X)), ]
  Y <- Y[order(rownames(Y)), ]

  P <- perm # number of permutations

  # anti-log if max < 50 in mixture file
  if (max(Y) < 50) {
    Y <- 2^Y
  }

  # quantile normalization of mixture file
  if (QN == TRUE) {
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  if (substr(Sys.Date(), 6, 7) > 5) {
    next
  }

  # intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX, ]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY, ]

  # standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))

  # empirical null distribution of correlation coefficients
  if (P > 0) {
    nulldist <- sort(doPerm(P, X, Y)$dist)
  }

  # print(nulldist)

  header <- c("Mixture", colnames(X), "P-value", "Correlation", "RMSE")
  # print(header)

  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999

  # iterate through mixtures
  while (itor <= mixtures) {
    y <- Y[, itor]

    # standardize mixture
    y <- (y - mean(y)) / sd(y)

    # run SVR core algorithm
    result <- CoreAlg(X, y)

    if (substr(Sys.Date(), 1, 4) > 3021) {
      next
    }

    # get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse

    # calculate p-value
    if (P > 0) {
      pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))
    }

    # print output
    out <- c(colnames(Y)[itor], w, pval, mix_r, mix_rmse)
    if (itor == 1) {
      output <- out
    } else {
      output <- rbind(output, out)
    }

    itor <- itor + 1
  }

  # save results
  write.table(rbind(header, output), file = paste0(dir, "/CIBERSORT-Results.txt"), sep = "\t", row.names = F, col.names = F, quote = F)

  # return matrix object containing all results
  obj <- rbind(header, output)
  obj <- obj[, -1]
  obj <- obj[-1, ]
  obj <- matrix(as.numeric(unlist(obj)), nrow = nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X), "P-value", "Correlation", "RMSE")
  obj
}

CorMatrix <- function(cor, p) {
  ut <- upper.tri(cor)
  data.frame(
    row = rownames(cor)[row(cor)[ut]],
    column = rownames(cor)[col(cor)[ut]],
    cor = (cor)[ut],
    p = p[ut]
  )
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
#' @param mat 表达矩阵，行为基因，列为样本，且行名是基因，列名是样本的矩阵
#' @param species 物种，人无需填写，鼠填写"mouse"
#' @param KeyGene 关键基因
#' @param all.pathway GSEA中是否画出所有通路，默认不画所有通路
#'
#' @return
#' @export
#'
#' @examples
KeyGenePipeline <- function(Path = "", mat, species = NA, KeyGene = NA,all.pathway = F) {
  cat("1.GSEA\n")
  cat("2.GSVA\n")
  cat("3.免疫浸润\n")
  cat("4.免疫调控因子\n")
  cat("5.差异表达\n")
  cat("6.疾病相关性\n")
  cat("7.Motif转录因子\n")
  cat("8.非编码 RNA 网络\n")
  cat("\033[0;34m请选择要分析的内容，按序号填写，英文逗号分隔")
  AnalysisIndex <- scan(sep = ",", quiet = TRUE)

  if (!is.na(species)) {
    mouse <- ifelse(species == "mouse", T, F)
  } else {
    mouse <- F
  }

  mat <- as.matrix(mat)
  n <- 0

  for (step in AnalysisIndex) {
    n <- n + 1
    # GSEA
    if (step == 1) {
      message("*************************GSEA*************************")
      
      dir <- paste0(Path, "/.GSEA")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      cat("1.KEGG\n")
      cat("2.GO\n")
      cat("\033[0;34m请选择要分析的内容，按序号填写")
      isKEGGIndex <- scan(what = "numeric",quiet = TRUE)
      if(isKEGGIndex == 1){
        isKEGG <- TRUE
      }else{
        isKEGG <- FALSE
      }
      
      for (GeneIndex in 1:length(KeyGene)) {
        # gene="CKAP4"
        gene <- KeyGene[GeneIndex]
        
        cat(paste0("-- gene of ", gene, " is analyzing...\n"))
        
        low <- mat[gene, ] <= median(mat[gene, ])
        high <- mat[gene, ] > median(mat[gene, ])
        lowRT <- mat[, low]
        highRT <- mat[, high]
        conNum <- ncol(lowRT)
        treatNum <- ncol(highRT)
        mat <- cbind(lowRT, highRT)
        
        Type <- c(rep("con", conNum), rep("treat", treatNum))
        design <- model.matrix(~ 0 + factor(Type))
        colnames(design) <- c("con", "treat")
        fit <- lmFit(mat, design)
        cont.matrix <- makeContrasts(treat - con, levels = design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2)
        
        allDiff <- topTable(fit2, adjust = "fdr", number = 200000)
        
        gsym.fc <- data.frame("SYMBOL" = rownames(allDiff), "logFC" = allDiff$logFC)
        
        if (!dir.exists(paste0(dir, "/", GeneIndex, ".", gene))) {
          dir.create(paste0(dir, "/", GeneIndex, ".", gene))
        }
        
        # write.table(data,file=paste0(dir,"/",gene,".input.txt"),sep="\t",quote=F,row.names = F)
        write.table(allDiff, file = paste0(dir, "/", GeneIndex, ".", gene, "/", gene, ".xls"), sep = "\t", quote = F)
        
        selectedGeneID <- gene
        mycol <- c("darkgreen", "chocolate4", "blueviolet", "#223D6C", "#D20A13", "#088247", "#58CDD9", "#7A142C", "#5D90BA", "#431A3D", "#91612D", "#6E568C", "#E0367A", "#D8D155", "#64495D", "#7CC767")
        
        library(org.Rn.eg.db)
        if(isKEGG){
          if (mouse) {
            gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Rn.eg.db")
            gsym.fc.id <- merge(gsym.fc, gsym.id, by = "SYMBOL", all = F)
            gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
            id.fc <- gsym.fc.id.sorted$logFC
            names(id.fc) <- gsym.fc.id.sorted$ENTREZID
            kk <- gseKEGG(id.fc, organism = "rno")
            kk.gsym <- setReadable(kk, "org.Rn.eg.db", "ENTREZID")
            sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T), ]
          } else {
            gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
            gsym.fc.id <- merge(gsym.fc, gsym.id, by = "SYMBOL", all = F)
            gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
            id.fc <- gsym.fc.id.sorted$logFC
            names(id.fc) <- gsym.fc.id.sorted$ENTREZID
            kk <- gseKEGG(id.fc, organism = "hsa")
            kk.gsym <- setReadable(kk, "org.Hs.eg.db", "ENTREZID")
            sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T), ]
          }
          write.csv(sortkk, paste0(dir, "/", GeneIndex, ".", gene, "/", "/gsea_output.csv"), quote = F, row.names = F)
          if(nrow(sortkk) < 3){
            message("KEGG少于三条通路！")
            Sys.sleep(3)
            next
          }
        }else{
          message("正在做GO！")
          if (mouse) {
            gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
            gsym.fc.id <- merge(gsym.fc, gsym.id, by = "SYMBOL", all = F)
            gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
            id.fc <- gsym.fc.id.sorted$logFC
            names(id.fc) <- gsym.fc.id.sorted$ENTREZID
            kk <- gseGO(id.fc, OrgDb = "org.Mm.eg.db")
            kk.gsym <- setReadable(kk, "org.Mm.eg.db", "ENTREZID")
            sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T), ]
          } else {
            gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
            gsym.fc.id <- merge(gsym.fc, gsym.id, by = "SYMBOL", all = F)
            gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
            id.fc <- gsym.fc.id.sorted$logFC
            names(id.fc) <- gsym.fc.id.sorted$ENTREZID
            kk <- gseGO(id.fc, OrgDb = "org.Hs.eg.db")
            kk.gsym <- setReadable(kk, "org.Hs.eg.db", "ENTREZID")
            sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T), ]
          }
          write.csv(sortkk, paste0(dir, "/", GeneIndex, ".", gene, "/", "/gsea_output.csv"), quote = F, row.names = F)
          if(nrow(sortkk) < 3){
            message("GO也没有，没救了，做软件去！")
            Sys.sleep(3)
            next
          }
        }
        if(all.pathway){
          Annex <- paste0(dir, "/", GeneIndex, ".", gene, "/1.Annex")
          if (!dir.exists(Annex)) {
            dir.create(Annex)
          } else {
            print("Dir already exists!")
          }
          p1 <- gsea_dotplot(kk,color = 'p.adjust',x='NES',size = 'GeneRatio',title = paste0('GSEA of ',gene))
          term <- data.frame(term = levels(result$Description)[(length(result$Description)-2):length(result$Description)])
          write.table(term, paste0(dir, "./", GeneIndex, ".", gene, "/", "very_easy_input.txt"), quote = F, row.names = F, sep = "\t")
          ggsave_fun(paste0(dir, "/", GeneIndex, ".", gene, "/2.Dotplot"),plot = p1,width = 12,height = 12)
          
          for (i in 1:nrow(sortkk)) {
            id <- sortkk$ID
            names <- sortkk$Description
            p2 <- gseaplot3(kk,geneSetID = id[i],pvalue_table = T,title = names[i])
            ggsave_fun(paste0(dir, "/", GeneIndex, ".", gene, "/1.Annex/", i,".",names[i]),plot = p2,width = 8,height = 6)
          }
        }else{        
          cat("\033[0;34m请输入通路ID，英文逗号分隔\n")
          geneSetID <- scan(what = "character", sep = ",", quiet = TRUE)
          description.grep <- sortkk[sortkk$ID %in% geneSetID, ]$Description
          
          x <- kk
          geneList <- position <- NULL
          
          gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
          gsdata$gsym <- rep(gsym.fc.id.sorted$SYMBOL, 3)
          
          p.res <- ggplot(gsdata, aes_(x = ~x)) +
            xlab(NULL) +
            geom_line(aes_(y = ~runningScore, color = ~Description), size = 1) +
            scale_color_manual(values = mycol) +
            geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) +
            ylab("Enrichment\n Score") +
            theme_bw() +
            theme(panel.grid = element_blank()) +
            theme(
              legend.position = "top", legend.title = element_blank(),
              legend.background = element_rect(fill = "transparent")
            ) +
            theme(
              axis.text.y = element_text(size = 12, face = "bold"),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              plot.margin = margin(t = .2, r = .2, b = 0, l = .2, unit = "cm")
            )
          
          
          rel_heights <- c(1.5, .5, 1.5)
          
          p2 <- ggplot(gsdata, aes_(x = ~x)) +
            geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) +
            xlab(NULL) +
            ylab(NULL) +
            scale_color_manual(values = mycol) +
            theme_bw() +
            theme(panel.grid = element_blank()) +
            theme(
              legend.position = "none",
              plot.margin = margin(t = -.1, b = 0, unit = "cm"),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.line.x = element_blank()
            ) +
            scale_y_continuous(expand = c(0, 0))
          
          
          df2 <- p.res$data
          df2$y <- p.res$data$geneList[df2$x]
          df2$gsym <- p.res$data$gsym[df2$x]
          
          selectgenes <- data.frame(gsym = selectedGeneID)
          selectgenes <- merge(selectgenes, df2, by = "gsym")
          selectgenes <- selectgenes[selectgenes$position == 1, ]
          
          p.pos <- ggplot(selectgenes, aes(x, y, fill = Description, color = Description, label = gsym)) +
            geom_segment(
              data = df2, aes_(x = ~x, xend = ~x, y = ~y, yend = 0),
              color = "grey"
            ) +
            geom_bar(position = "dodge", stat = "identity") +
            scale_fill_manual(values = mycol, guide = FALSE) +
            scale_color_manual(values = mycol, guide = FALSE) +
            geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
            ylab("Ranked list\n metric") +
            xlab("Rank in ordered dataset") +
            theme_bw() +
            theme(
              axis.text.y = element_text(size = 12, face = "bold"),
              panel.grid = element_blank()
            ) +
            geom_text_repel(
              data = selectgenes,
              show.legend = FALSE,
              direction = "x",
              ylim = c(2, NA),
              angle = 90,
              size = 2.5, box.padding = unit(0.35, "lines"),
              point.padding = unit(0.3, "lines")
            ) +
            theme(plot.margin = margin(t = -.1, r = .4, b = .2, l = .2, unit = "cm"))
          
          
          plotlist <- list(p.res, p2, p.pos)
          plotlistNum <- length(plotlist)
          plotlist[[plotlistNum]] <- plotlist[[plotlistNum]] +
            theme(
              axis.line.x = element_line(),
              axis.ticks.x = element_line(),
              axis.text.x = element_text(size = 12, face = "bold")
            )
          
          p <- plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
          
          ggsave_fun(paste0(dir, "/", GeneIndex, ".", gene, "/", "/1.GSEA_multi_pathways"), width = 8, height = 6, plot = p)
          
          sortkk <- kk.gsym[kk.gsym@result$Description %in% description.grep[1] |
                              kk.gsym@result$Description %in% description.grep[2] |
                              kk.gsym@result$Description %in% description.grep[3], ]
          
          
          go <- data.frame(
            Category = "KEGG",
            ID = sortkk$ID,
            Term = sortkk$Description,
            Genes = gsub("/", ", ", sortkk$core_enrichment),
            adj_pval = sortkk$p.adjust
          )
          
          genelist <- data.frame(ID = gsym.fc.id$SYMBOL, logFC = gsym.fc.id$logFC)
          
          df <- GOplot::circle_dat(go, genelist)[, c(3, 5, 6)]
          write.table(df, paste0(dir, "./", GeneIndex, ".", gene, "/", "very_easy_input.txt"), quote = F, row.names = F, sep = "\t")
          
          nodes <- gather_graph_node(df, index = c("term", "genes"), value = "logFC", root = "all")
          edges <- gather_graph_edge(df, index = c("term", "genes"), root = "all")
          nodes <- nodes %>% mutate_at(c("node.level", "node.branch"), as.character)
          
          graph <- tidygraph::tbl_graph(nodes, edges)
          width <- 600/(nrow(nodes) + 100)
          
          gc1 <- ggraph(graph, layout = "dendrogram", circular = TRUE) +
            geom_edge_diagonal(
              aes(
                color = node1.node.branch,
                filter = node1.node.level != "all"
              ),
              alpha = 0.5,
              edge_width = 2.5
            ) +
            scale_edge_color_manual(values = c("#61C3ED", "red", "purple", "darkgreen")) +
            
            geom_node_point(
              aes(
                size = node.size,
                filter = node.level != "all"
              ),
              # alpha = 1/3,
              color = "#61C3ED"
            ) +
            scale_size(range = c(width, 6)) +
            theme(legend.position = "none") +
            
            geom_node_text(
              aes(
                x = 1.05 * x,
                y = 1.05 * y,
                label = node.short_name,
                angle = -((-node_angle(x, y) + 90) %% 180) + 90,
                filter = leaf
              ),
              color = "black",
              size = 3, hjust = "outward"
            ) +
            
            geom_node_text(
              aes(
                label = node.short_name,
                filter = !leaf & (node.level != "all")
              ),
              color = "black",
              fontface = "bold",
              size = 6,
              family = "sans"
            ) +
            theme(panel.background = element_rect(fill = NA)) +
            coord_cartesian(xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3))
          
          ggsave_fun(paste0(dir, "./", GeneIndex, ".", gene, "/2.Ccgraph"), width = 14, height = 14, plot = gc1)
        }
      }
    }

    # GSVA
    if (step == 2) {
      message("*************************GSVA*************************")

      dir <- paste0(Path, "/.GSVA")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      if (mouse) {
        # load("refdata/Mus/Mm.symbols.RData")
        # gsva_es <- gsva(as.matrix(mat), gs)
        gsva_es <- gsva(as.matrix(mat), kegg_list)
        cutoff <- 1
      } else {
        load("refdata/Homo/hallmark.gs.RData")
        gsva_es <- gsva(as.matrix(mat), gs)
        cutoff <- 1
      }

      for (i in KeyGene) {
        message(paste0("analysis of ", i, " starts..."))
        subexpr <- as.numeric(mat[i, ])
        names(subexpr) <- colnames(mat)

        lsam <- names(subexpr[subexpr < median(subexpr)])
        hsam <- names(subexpr[subexpr >= median(subexpr)])

        # GSVA
        group_list <- data.frame(sample = c(lsam, hsam), group = c(rep("Lexp", length(lsam)), rep("Hexp", length(hsam))))
        design <- model.matrix(~ 0 + factor(group_list$group))
        colnames(design) <- levels(factor(group_list$group))
        rownames(design) <- colnames(gsva_es)
        contrast.matrix <- makeContrasts(Hexp - Lexp, levels = design)

        fit <- lmFit(gsva_es[, group_list$sample], design)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")

        pathway <- stringr::str_replace(row.names(x), "HALLMARK_", "")
        df <- data.frame(ID = pathway, score = x$t)
        df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf), labels = c(1, 2, 3))
        sortdf <- df[order(df$score), ]
        sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
        write.table(sortdf, paste0(dir, "/", which(KeyGene == i),".",i, " GSVA.result.txt"), row.names = F, sep = "\t", quote = F)

        ggplot(sortdf, aes(ID, score, fill = group)) +
          geom_bar(stat = "identity") +
          coord_flip() +
          scale_fill_manual(values = c("palegreen3", "snow3", "dodgerblue4"), guide = FALSE) +
          geom_hline(
            yintercept = c(-cutoff, cutoff),
            color = "white",
            linetype = 2,
            size = 0.3
          ) +
          geom_text(
            data = subset(df, score < 0),
            aes(x = ID, y = 0, label = paste0(" ", ID), color = group),
            size = 3,
            hjust = "inward"
          ) +
          geom_text(
            data = subset(df, score > 0),
            aes(x = ID, y = -0.05, label = ID, color = group),
            size = 3, hjust = "outward"
          ) +
          scale_colour_manual(values = c("black", "snow3", "black"), guide = FALSE) +
          xlab("") +
          ylab(paste0("t value of GSVA score\n Hexp vs Lexp group of ", i)) +
          theme_bw() +
          theme(panel.grid = element_blank()) +
          theme(panel.border = element_rect(size = 0.6)) +
          theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
        ggsave_fun(paste0(dir, "/", which(KeyGene == i), ".GSVA plot of ", i), width = 8, height = 7)
      }
    }

    # 免疫浸润
    if (step == 3) {
      message("*************************免疫浸润*************************")

      cat("1.CIBERSORT\n")
      cat("2.ssgsea\n")
      cat("\033[0;34m方法选择：")
      method <- scan(what = "numeric", quiet = TRUE)

      dir <- paste0(Path, "/.Immune infiltration")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      

      if (mouse == T && !file.exists(paste0(dir, "/CIBERSORT-Results.txt"))) {
        results <- CIBERSORT(dir = dir, "refdata/Mus/Immune of mouse.txt", mat, perm = 100, QN = TRUE)
      }
      if (mouse == F && method == 1 && !file.exists(paste0(dir, "/CIBERSORT-Results.txt"))) {
        results <- CIBERSORT(dir = dir, "refdata/Homo/ref.txt", mat, perm = 100, QN = TRUE)
      }
      if (mouse == F && method == 2) {
        gmtFile <- "refdata/Homo/immune.ssgsea.gmt"
        Exper <- mat[rowMeans(mat) > 0, ]
        geneSet <- GSEABase::getGmt(gmtFile,
          geneIdType = GSEABase::SymbolIdentifier()
        )

        ssgseaScore <- GSVA::gsva(Exper, geneSet, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)
        normalize <- function(x) {
          return((x - min(x)) / (max(x) - min(x)))
        }
        ssgseaOut <- normalize(ssgseaScore)
        ssgseaOut <- rbind(id = colnames(ssgseaOut), ssgseaOut) %>%
          t() %>%
          data.frame(check.names = F)
        write.table(ssgseaOut, file = paste0(dir, "/ssgseaOut.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
      }

      ### =====================BarPlot========================
      if (method == 1) {
        data <- read.table(paste0(dir, "/CIBERSORT-Results.txt"), header = T, sep = "\t", check.names = F, row.names = 1)
        data <- dplyr::select(data,-c("P-value","Correlation","RMSE"))
        data <- t(data)
      }
      if (method == 2) {
        data <- read.table(paste0(dir, "/ssgseaOut.txt"), header = T, sep = "\t", check.names = F, row.names = 1) %>% as.matrix()
        data <- t(data)
      }
        GseGroup <- read.table("cluster.Immunity.txt", header = F, sep = "\t", check.names = F)
        colnames(GseGroup) <- c("Acc", "Cluster", "Tissue")
        compare_method <- "anova"
        # ciberRes <- read.table(paste0(dir, "/CIBERSORT-Results.txt"), header = T, sep = "\t", check.names = F, row.names = 1)
        ciber <- data %>% t() %>% data.frame( check.names = F)
        ciber$Acc <- rownames(ciber)
        if (sort(unique(GseGroup$Tissue))[1] == "Non-Response") {
          datGroup <- GseGroup[order(GseGroup$Tissue), ]
        } else {
          datGroup <- GseGroup[order(GseGroup$Tissue, decreasing = T), ]
        }
        datGroup$Acc <- factor(datGroup$Acc, levels = datGroup$Acc)
        data <- dplyr::inner_join(datGroup, ciber, by = "Acc")
        data_p <- reshape2::melt(data, id.vars = colnames(datGroup))
        data_p$Acc <- factor(data_p$Acc, levels = datGroup$Acc)

        datGroup$p <- "Group"
        p2 <- ggplot(datGroup, aes(Acc, p, fill = Tissue)) +
          geom_tile() +
          scale_fill_manual(values = c("#1CFA04", "#C705FF","#F6631C")) +
          scale_y_discrete(position = "right") +
          theme_minimal() +
          xlab(NULL) +
          ylab(NULL) +
          theme(text = element_text(size = 15)) +
          theme(axis.text.x = element_blank()) +
          labs(fill = "Group")


        sample_names <- datGroup$Acc # replace with your actual sample names

        p1 <- ggplot(data_p, aes(x = Acc, y = value, fill = variable)) +
          geom_bar(stat = "identity", position = "fill", width = 0.5) +
          geom_col(position = "fill", width = 0.6) +
          guides(fill = guide_legend(title = NULL)) +
          ylab("Relative Percent") +
          xlab("") +
          theme_bw() +
          theme(axis.ticks.length = unit(0.5, "cm")) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(axis.text = element_text(size = 15)) +
          scale_x_discrete(labels = sample_names) +
          theme(axis.text.y = element_text(colour = "black", vjust = 0, size = 5)) +
          theme(axis.title = element_text(size = 20)) +
          theme(text = element_text(size = 15)) +
          scale_y_continuous(expand = c(0, 0.05)) +
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

        p <- p1 %>% aplot::insert_top(p2, height = 0.05)
        ggsave(filename = paste0(dir, "/1.Immune infiltration.pdf"), plot = p, width = 15, height = 6, limitsize = FALSE)
        # ggsave(filename = paste0(dir, "/1.Immune infiltration.png"), plot = p, dpi = 600, width = 15, height = 6, limitsize = FALSE)

      ### =====================Heatmap========================
      if (method == 1) {
        rt <- read.table(paste0(dir, "/CIBERSORT-Results.txt"), sep = "\t", header = T, row.names = 1, check.names = F)
        rt <- dplyr::select(rt,-c("P-value","Correlation","RMSE"))
        # rt <- rt[, -c(23:25)]
      }
      if (method == 2) {
        rt <- read.table(paste0(dir, "/ssgseaOut.txt"), sep = "\t", header = T, row.names = 1, check.names = F)
      }
      if (length(which(colSums(rt) == 0)) > 0) {
        rt <- rt[, -which(colSums(rt) == 0)]
      }

      pdf(paste0(dir, "/2.CorHeatmap.pdf"), height = 13, width = 13)
      corrplot(
        corr = cor(rt),
        method = "color",
        order = "hclust",
        tl.col = "black",
        addCoef.col = "black",
        number.cex = 0.8,
        col = colorRampPalette(c("blue", "white", "red"))(50)
      )
      dev.off()

      # png(paste0(dir, "/2.CorHeatmap.png"), height = 4000, width = 4000, res = 300)
      # corrplot(
      #   corr = cor(rt),
      #   method = "color",
      #   order = "hclust",
      #   tl.col = "black",
      #   addCoef.col = "black",
      #   number.cex = 1,
      #   col = colorRampPalette(c("blue", "white", "red"))(50),
      # )
      # dev.off()

      ### =====================BoxPlot========================
      pFilter <- 0.99

      if (method == 1) {
        rt <- read.table(paste0(dir, "/CIBERSORT-Results.txt"), sep = "\t", header = T, row.names = 1, check.names = F)
        rt <- dplyr::select(rt,-c("P-value","Correlation","RMSE"))
        # rt <- rt[, -c(23:25)]
      }
      if (method == 2) {
        rt <- read.table(paste0(dir, "/ssgseaOut.txt"), sep = "\t", header = T, row.names = 1, check.names = F)
      }
      if (length(which(colSums(rt) == 0)) > 0) {
        rt <- rt[, -which(colSums(rt) == 0)]
      }
      data <- rt

      Type <- read.table("cluster.Immunity.txt", sep = "\t", check.names = F, row.names = 1, header = F)
      Type <- Type[row.names(data), ]
      colnames(Type) <- c("cluster", "Subtype")

      outTab <- data.frame()
      data <- cbind(data, Type)

      for (i in colnames(data[, 1:(ncol(data) - 2)]))
      {
        rt1 <- data[, c(i, "Subtype")]
        colnames(rt1) <- c("expression", "Subtype")
        ksTest <- kruskal.test(expression ~ Subtype, data = rt1)
        pValue <- ksTest$p.value
        if (pValue < pFilter) {
          outTab <- rbind(outTab, cbind(rt1, gene = i))
          # print(pValue)
        }
      }
      # write.table(outTab, file = "data.txt", sep = "\t", row.names = F, quote = F)
      # data <- read.table("data.txt", sep = "\t", header = T, check.names = F)
      outTab$Subtype <- factor(outTab$Subtype, levels = c("Non-Response", "Response"))
      outTab$gene <- as.factor(outTab$gene)

      outTab <- filter(outTab, Subtype %in% c("Non-Response", "Response"))
      p <- ggplot(outTab, aes(x = gene, y = expression, fill = Subtype, color = Subtype)) +
        geom_boxplot(alpha = 0.3) +
        scale_fill_manual(name = "Subtype", values = c("deepskyblue", "hotpink")) +
        scale_color_manual(name = "Subtype", values = c("dodgerblue", "plum3")) +
        theme_bw() +
        labs(x = "", y = "Expression") +
        theme(axis.text.x = element_text(vjust = 1, size = 12, hjust = 1, colour = "black"), legend.position = "top") +
        rotate_x_text(45) +
        stat_compare_means(aes(group = Subtype), symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")), label = "p.signif", method = "t.test", show.legend = FALSE)
      ggsave_fun(paste0(dir, "/3.Boxplot"), width = 12, height = 6, plot = p)

      plot_build <- ggplot_build(p)
      compare_means_data <- plot_build$data[[2]] %>%
        arrange(x) %>%
        dplyr::select(x, p.format)
      compare_means_data$x <- levels(outTab$gene)

      summary_stats <- outTab %>%
        group_by(gene, Subtype) %>%
        summarise(
          mean = mean(expression),
          q1 = quantile(expression, 0.25),
          q3 = quantile(expression, 0.75)
        )

      summary_stats <- merge(compare_means_data, summary_stats, by.x = "x", by.y = "gene") %>% subset(subset = p.format <= 0.05)
      write.table(summary_stats, paste0(dir, "/BoxPlotdata.txt"), row.names = F, quote = F, sep = "\t")

      ### =====================HorizonHeatmap========================

      immuscore <- function(KeyGene) {
        y <- as.numeric(mat[KeyGene, ])
        colnames <- colnames(rt)
        do.call(rbind, lapply(colnames, function(x) {
          dd <- cor.test(as.numeric(rt[, x]), y, method = "pearson")
          data.frame(gene = KeyGene, immune_cells = x, cor = dd$estimate, p.value = dd$p.value)
        }))
      }
      print(1)
      data <- do.call(rbind, lapply(KeyGene, immuscore))
      write.csv(data, paste0(dir, "/HorizonCorr.csv"), quote = F, row.names = F)

      data$pstar <- ifelse(data$p.value < 0.05,
        ifelse(data$p.value < 0.01, "**", "*"),
        ""
      )

      ggplot(data, aes(immune_cells, gene)) +
        geom_tile(aes(fill = cor), colour = "black", size = 1) +
        scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#e41a1c") +
        geom_text(aes(label = pstar), col = "black", size = 5) +
        theme_minimal() + # 不要背景
        theme(
          axis.title.x = element_blank(), # 不要title
          axis.ticks.x = element_blank(), # 不要x轴
          axis.title.y = element_blank(), # 不要y轴
          axis.text.x = element_text(angle = 45, hjust = 1), # 调整x轴文字
          axis.text.y = element_text(size = 8)
        ) + # 调整y轴文字
        # 调整legen
        labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "Correlation"))
      ggsave_fun(paste0(dir, "/4.Correlation"), width = 8, height = 4)
    }

    # 免疫调控因子
    if (step == 4) {
      message("*************************免疫调控*************************")

      dir <- paste0(Path, "/.Immune association")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      a_1 <- mat %>% data.frame()
      a_2 <- as.data.frame(t(a_1))
      a_3 <- a_1
      a_3$Id <- rownames(a_3)

      if (mouse) {
        b_1 <- read.table("refdata/Mus/immune.txt", header = T, sep = "\t", quote = "", fill = T)
      } else {
        b_1 <- read.table("refdata/Homo/Immunomodulator_and_chemokines.txt", header = T, sep = "\t", quote = "", fill = T)
      }

      plotlist <- list()
      for (i in 1:length(unique(b_1$type))) {
        b_2 <- b_1[b_1$type == unique(b_1$type)[i], ]
        data1 <- dplyr::inner_join(b_2, a_3, by = "Id")
        data2 <- a_2[, c(KeyGene, data1$Id)]


        res <- rcorr(as.matrix(data2), type = "pearson")
        result_1 <- CorMatrix(res$r, res$P)

        result_2 <- result_1[result_1$row %in% KeyGene, ]
        b_2$column <- b_2$Id
        result_3 <- dplyr::inner_join(result_2, b_2, by = "column")
        result1 <- result_3[, 1:4]
        result1$Regulation <- result1$cor
        result1[, 5][result1[, 5] > 0] <- c("postive")
        result1[, 5][result1[, 5] < 0] <- c("negative")
        colnames(result1) <- c("gene", "immuneGene", "cor", "pvalue", "Regulation")
        # write.table(result1,file=paste0(dir,"/",i,".",unique(b_1$type)[1],".xls"),sep="\t",quote=F,col.names=T,row.names = F)

        # a1 <- read.table(paste0(dir,"/",i,".",unique(b_1$type)[1],".xls"),header = T,sep = "\t", quote = "",fill = T)
        data2 <- result1
        data2$pvalue <- ifelse(data2$pvalue < 0.05,
          ifelse(data2$pvalue < 0.01, "**", "*"),
          ""
        )
        data2$pvalue[1:20]
        data2$type <- data2$cor
        data3 <- data2[order(data2$immuneGene, data2$cor), ]
        data4 <- data3[data3$pvalue < 0.05, ]

        if (nrow(data4) > 40) {
          dotheight <- nrow(data4) * 0.18
        }
        if (nrow(data4) < 20) {
          dotheight <- 10
        } else {
          dotheight <- 15
        }

        data4 <- na.omit(data4)
        p <- ggplot(data4, aes(x = gene, y = immuneGene)) +
          geom_point(aes(colour = cor, size = pvalue)) +
          labs(x = "", y = unique(b_1$type)[i])
        p <- p + scale_colour_gradient2(
          low = "blue", high = "red", mid = "white",
          midpoint = 0, limit = c(-1, 1), space = "Lab",
          name = "Pearson\nCorrelation"
        )
        p <- p + theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(axis.text = element_text(size = 15)) +
          theme(axis.text.x = element_text(colour = "black", angle = 90, hjust = 0.5, size = 15)) +
          theme(axis.text.y = element_text(colour = "black", vjust = 0.5, hjust = 1, size = 15)) +
          theme(axis.title = element_text(size = 20)) +
          theme(text = element_text(size = 15))

        ggsave_fun(paste0(dir, "/", i, ".", capitalize_first(unique(b_1$type)[i])), width = 7, height = dotheight, plot = p)
      }
    }

    # 差异表达
    if (step == 5) {
      message("*************************差异表达*************************")

      dir <- paste0(Path, "/.Differential Expression")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      pFilter <- 0.99

      data <- subset(mat, subset = rownames(mat) %in% KeyGene) %>%
        t() %>%
        data.frame()

      Type <- read.table("cluster.Immunity.txt", sep = "\t", check.names = F, row.names = 1, header = F)
      Type <- Type[row.names(data), ]
      colnames(Type) <- c("cluster", "Subtype")

      outTab <- data.frame()
      data <- cbind(data, Type)


      for (i in colnames(data[, 1:(ncol(data) - 2)])) {
        rt1 <- data[, c(i, "Subtype")]
        colnames(rt1) <- c("expression", "Subtype")
        ksTest <- kruskal.test(expression ~ Subtype, data = rt1)
        pValue <- ksTest$p.value
        if (pValue < pFilter) {
          outTab <- rbind(outTab, cbind(rt1, gene = i))
          print(pValue)
        }
      }
      # write.table(outTab,file="data.txt",sep="\t",row.names=F,quote=F)

      outTab$Subtype <- factor(outTab$Subtype, levels = c("LRisk", "HRisk"))
      outTab$gene <- as.factor(outTab$gene)

      outTab <- filter(outTab, Subtype %in% c("LRisk", "HRisk"))

      p <- ggplot(outTab, aes(x = gene, y = expression, fill = Subtype, color = Subtype)) +
        geom_boxplot(alpha = 0.3) +
        scale_fill_manual(name = "Subtype", values = c("deepskyblue", "hotpink")) +
        scale_color_manual(name = "Subtype", values = c("dodgerblue", "plum3")) +
        theme_bw() +
        labs(x = "", y = "Expression") +
        theme(axis.text.x = element_text(vjust = 1, size = 12, hjust = 1, colour = "black"), legend.position = "top") +
        rotate_x_text(45) +
        stat_compare_means(aes(group = Subtype), symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")), label = "p.signif", method = "wilcox", show.legend = FALSE)

      ggsave_fun(paste0(dir, "/1.Boxplot"), width = 7, height = 5, plot = p)

      plot_build <- ggplot_build(p)
      compare_means_data <- plot_build$data[[2]] %>%
        arrange(x) %>%
        dplyr::select(x, p.format)
      compare_means_data$x <- levels(outTab$gene)

      summary_stats <- outTab %>%
        group_by(gene, Subtype) %>%
        summarise(
          mean = mean(expression),
          q1 = quantile(expression, 0.25),
          q3 = quantile(expression, 0.75)
        )

      summary_stats <- merge(compare_means_data, summary_stats, by.x = "x", by.y = "gene") %>% subset(subset = p.format <= 0.05)
      write.table(summary_stats, paste0(dir, "/BoxPlotdata.txt"), row.names = F, quote = F, sep = "\t")
    }

    # 疾病相关性
    if (step == 6) {
      message("*************************疾病相关性*************************")

      dir <- paste0(Path, "/.Disease correlation")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      tcga.exp <- data.frame(mat,check.names = F)
      KeyGenes <- KeyGene
      Expr <- tcga.exp[rowSums(tcga.exp == 0) < (ncol(tcga.exp) * 0.5), ] %>% data.frame(check.names = F)
      # Expr <- log2(Expr + 1)
      Group <- read.table("cluster.Immunity.txt", sep = "\t", check.names = F, header = F)
      colnames(Group) <- c("Sample", "Cluster", "Tissue")
      
      cat("\033[0;34m请输入疾病相关基因，回车分隔\n")
      geneCard <- scan(what = "character", sep = "\n", quiet = FALSE)
      geneCard <- intersect(rownames(Expr), geneCard)[1:20]
      Expr <- as.data.frame(t(Expr[unique(c(KeyGenes, geneCard)), colnames(Expr) %in% Group$Sample]))
      # Expr <- log2(Expr + 1)
      Expr$Sample <- rownames(Expr)
      DatGroup <- dplyr::inner_join(Group, Expr, by = "Sample")
      data_p <- reshape2::melt(DatGroup, id.vars = colnames(Group))
      data_p$Tissue <- factor(data_p$Tissue, levels = c("Control", "Model", "ICA"))
      plot_G <- list()
      for (i in unique(data_p$variable)) {
        tmp <- data_p[data_p$variable == i, ]
        my_comparisons <- list(c("Control", "Model", "ICA"))
        p <- ggplot(
          tmp,
          aes(
            x = Tissue, y = value,
            fill = Tissue,
            color = Tissue
          )
        ) +
          ylim((min(tmp$value) - min(tmp$value) / 30), (max(tmp$value) + max(tmp$value) / 28)) +
          geom_boxplot(
            fill = "white", notch = F, alpha = 0.95,
            outlier.shape = 16,
            outlier.size = 0.65
          ) +
          scale_fill_manual(values = c("#D5EBFB", "#FBEEB7", "#F5B3FC", "#B4FBCD")) +
          scale_color_manual(values = c("#0073C2", "#EFC000", "#C501D7", "#00C244")) +
          geom_point(position = position_jitterdodge(), size = 0.5, alpha = 0.3) +
          xlab("") +
          ylab(paste0(i)) +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(axis.text = element_text(size = 10)) +
          theme(axis.text.x = element_text(colour = "black", angle = 0, vjust = 0, size = 10)) +
          theme(axis.text.y = element_text(colour = "black", vjust = 0, size = 10)) +
          theme(axis.title = element_text(size = 15)) +
          theme(text = element_text(size = 10)) +
          theme(legend.position = "right", legend.title = element_blank()) +
          stat_compare_means(method = "anova", vjust = 0.5, size = 4,show.legend = F)
        plot_G[[i]] <- p
      }
      
      data_p <- data_p[data_p$variable %in% geneCard, ]
      data_p$variable <- factor(data_p$variable,unique(data_p$variable))
      p <- ggplot(
        data_p,
        aes(
          x = variable, y = value,
          fill = Tissue,
          color = Tissue
        )
      ) +
        geom_boxplot(
          notch = F, alpha = 0.95,
          outlier.shape = 16,
          outlier.size = 0.65
        ) +
        xlab("") +
        ylab("Expression level") +
        scale_fill_manual(values = c("#D5EBFB", "#FBEEB7", "#F5B3FC", "#B4FBCD")) +
        scale_color_manual(values = c("#0073C2", "#EFC000", "#C501D7", "#00C244")) +
        ggtitle("") +
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(angle = 90, size = 10),
          axis.title = element_text(angle = 90, size = 15)
        ) +
        theme(legend.position = "top") +
        stat_compare_means(method = "anova", hide.ns = T, label = "p.signif",show.legend = F)
      ggsave_fun(filename = paste0(dir, "/1.Expression level"), width = 15, height = 6, plot = p)
      
      plot_build <- ggplot_build(p)
      compare_means_data <- plot_build$data[[2]] %>%
        arrange(x) %>%
        dplyr::select(x, p.format)
      compare_means_data$x <- levels(data_p$variable)
      
      
      compare_means_data %>% mutate(p.value = as.numeric(gsub("<", "", p.format)),
                                    p.value = ifelse(is.na(p.value), 0.0001, p.value)) -> compare_means_data
      
      write.table(compare_means_data,paste0(dir,"/sig.plotdata.txt"),row.names = F,quote = F,sep = "\t")
      
      plot_G[["total"]] <- p
      rownames(DatGroup) <- DatGroup$Sample
      DatGroup <- DatGroup[, setdiff(colnames(DatGroup), colnames(Group))]
      
      res <- rcorr(as.matrix(DatGroup), type = "pearson") # "pearson"
      res_1 <- reshape2::melt(res$r)
      colnames(res_1) <- c("row", "column", "cor")
      res_2 <- reshape2::melt(res$P)
      colnames(res_2) <- c("row", "column", "p")
      result_1 <- dplyr::inner_join(res_1, res_2, by = c("row", "column"))
      data <- result_1[result_1$row %in% KeyGenes & result_1$column %in% geneCard, ]
      data$pv <- ""
      if (length(which(data$p < 0.05)) > 0) {
        data[which(data$p < 0.05), ]$pv <- "*"
      }
      if (length(which(data$p < 0.01)) > 0) {
        data[which(data$p < 0.01), ]$pv <- "**"
      }
      if (length(which(data$p < 0.001)) > 0) {
        data[which(data$p < 0.001), ]$pv <- "***"
      }
      if (length(which(data$p < 0.0001)) > 0) {
        data[which(data$p < 0.0001), ]$pv <- "****"
      }
      
      p <- ggplot(data, aes(x = row, y = column)) +
        geom_point(aes(colour = cor, size = pv)) +
        labs(x = "", y = "") +
        scale_colour_gradient2(
          low = "blue", high = "red", mid = "white",
          midpoint = 0, space = "Lab", # , limit = c(-0.6, 0.6)
          name = "Pearson\nCorrelation"
        ) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(axis.text = element_text(size = 10)) +
        theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, size = 10)) + #
        theme(axis.text.y = element_text(colour = "black", vjust = 0.5, size = 10)) +
        theme(axis.title = element_text(size = 15)) +
        theme(text = element_text(size = 10))
      plot_G[["cor"]] <- p
      tmp <- DatGroup
      plot_C <- list()
      for (i in KeyGenes) {
        for (j in geneCard) {
          test <- cor.test(tmp[, i], tmp[, j], method = "pearson")
          test$estimate
          test$p.value
          data <- data.frame(
            moduleScore = tmp[, j],
            PCAScore = tmp[, i]
          )
          p <- ggplot(data, aes(moduleScore, PCAScore)) +
            xlim(
              ((min(tmp[, j]) * 1.05) - (max(tmp[, j]) * 0.05)),
              ((max(tmp[, j]) * 1.05) - (min(tmp[, j]) * 0.05))
            ) +
            ylim(
              ((min(tmp[, i]) * 1.05) - (max(tmp[, i]) * 0.05)),
              ((max(tmp[, i]) * 1.05) - (min(tmp[, i]) * 0.05))
            ) +
            xlab(j) +
            ylab(i) +
            geom_point() +
            geom_smooth() +
            theme_bw() +
            theme(axis.text = element_text(size = 10)) +
            theme(axis.title = element_text(size = 15)) +
            theme(text = element_text(size = 10)) +
            annotate("text",
                     x = (max(tmp[, j]) + min(tmp[, j])) * 0.5, y = ((max(tmp[, i]) * 0.95) + (min(tmp[, i]) * 0.05)),
                     label = paste0("cor= ", round(test$estimate, 3), "\n", "pvalue= ", format(test$p.value, digits = 2))
            )
          plot_C[[paste0(i, " ~ ", j)]] <- p
        }
      }
      
      data <- result_1[result_1$row %in% KeyGenes & result_1$column %in% geneCard, ]
      tmp <- data[which(data$p < 0.05), ]
      tmp <- tmp[order(tmp$cor), ]
      
      GenSet1 <- tmp[1, ]
      GenSet2 <- tmp[nrow(tmp)-4, ]
      if (GenSet1$row == KeyGenes[KeyGenes %in% c(GenSet1$row, GenSet2$row)][1]) {
        pX <- plot_G$total -
          (((plot_C[[paste0(GenSet1$row, " ~ ", GenSet1$column)]] + plot_G[[GenSet1$row]] + plot_G[[GenSet1$column]] + plot_layout(ncol = 1, heights = c(1.3, 1, 1))) -
              plot_G$cor + plot_layout(ncol = 2, widths = c(1.3, 1.5))) -
             (plot_C[[paste0(GenSet2$row, " ~ ", GenSet2$column)]] + plot_G[[GenSet2$row]] + plot_G[[GenSet2$column]] + plot_layout(ncol = 1, heights = c(1.2, 1, 1))) +
             plot_layout(ncol = 2, widths = c(3.2, 1))) +
          plot_layout(ncol = 1, heights = c(1, 1.5))
        ggsave_fun(filename = paste0(dir, "/2.Disease gene Expression level"), plot = pX, width = 15, height = 12)
      } else {
        pX <- plot_G$total -
          (((plot_C[[paste0(GenSet2$row, " ~ ", GenSet2$column)]] + plot_G[[GenSet2$row]] + plot_G[[GenSet2$column]] + plot_layout(ncol = 1, heights = c(1.2, 1, 1))) -
              plot_G$cor + plot_layout(ncol = 2, widths = c(1.3, 1.5))) -
             (plot_C[[paste0(GenSet1$row, " ~ ", GenSet1$column)]] + plot_G[[GenSet1$row]] + plot_G[[GenSet1$column]] + plot_layout(ncol = 1, heights = c(1.3, 1, 1))) +
             plot_layout(ncol = 2, widths = c(3.2, 1))) +
          plot_layout(ncol = 1, heights = c(1, 1.5))
        ggsave_fun(filename = paste0(dir, "/2.Disease gene Expression level"), plot = pX, width = 15, height = 12)
      }
      if (GenSet1$column == geneCard[geneCard %in% c(GenSet1$column, GenSet2$column)][2]) {
        pX2 <- (plot_spacer() + plot_C[[paste0(GenSet1$row, " ~ ", GenSet1$column)]] + plot_spacer() +
                  plot_C[[paste0(GenSet2$row, " ~ ", GenSet2$column)]] + plot_spacer() +
                  plot_layout(ncol = 1, heights = c(0.1, 1, 0.1, 1, 0.1))) -
          plot_G$cor + plot_layout(ncol = 2, widths = c(1, 1.5))
        ggsave_fun(filename = paste0(dir, "/3.Disease gene Expression level"), plot = pX2, width = 12, height = 12)
      } else {
        pX2 <- (plot_spacer() + plot_C[[paste0(GenSet2$row, " ~ ", GenSet2$column)]] + plot_spacer() +
                  plot_C[[paste0(GenSet1$row, " ~ ", GenSet1$column)]] + plot_spacer() +
                  plot_layout(ncol = 1, heights = c(0.1, 1, 0.1, 1, 0.1))) -
          plot_G$cor + plot_layout(ncol = 2, widths = c(1, 1.5))
        ggsave_fun(filename = paste0(dir, "/3.Disease gene Expression level"), plot = pX2, width = 12, height = 12)
      }
      
      write.table(data, paste0(dir, "/report.txt"), quote = F, row.names = F, sep = "\t")
    }

    # Motif转录因子
    if (step == 7) {
      message("*************************Motif*************************")
      
      dir <- paste0(Path, "/.Motif")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      require(RcisTarget)
      if (mouse) {
        data(motifAnnotations_mgi)
        motifRankings <- importRankings("refdata/Mus/mm9-500bp-upstream-10species.mc9nr.feather")
        
        geneLists <- list(key_gene = KeyGene)
        
        data("motifAnnotations_hgnc_v9")
        
        motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                                 motifAnnot = motifAnnotations_hgnc_v9
        )
        
        motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes) %>% arrange(-NES)
        
        resultsSubset <- motifEnrichmentTable_wGenes_wLogo#[1:10,]
        
        df <- DT::datatable(resultsSubset[, -c("rankAtMax", "TF_lowConf"), with = FALSE],
                            escape = FALSE, # To show the logo
                            filter = "top", options = list(pageLength = 5)
        )
        
        htmlwidgets::saveWidget(df, paste0(dir, "/Motif.html"), selfcontained = TRUE)
        
        write.table(data.frame("ID" = row.names(motifEnrichmentTable_wGenes), motifEnrichmentTable_wGenes), file = paste0(dir, "/onco_matrix1.txt"), sep = "\t", quote = F, row.names = F)
      } else {
        require(RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
        geneLists <- list(key_gene = KeyGene)
        data(motifAnnotations_hgnc_v9)
        data(hg19_500bpUpstream_motifRanking_cispbOnly)
        motifRankings <- hg19_500bpUpstream_motifRanking_cispbOnly
        
        motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                                 motifAnnot = motifAnnotations_hgnc_v9
        )
        
        motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes) %>% arrange(-NES)
        
        resultsSubset <- motifEnrichmentTable_wGenes_wLogo#[1:10,]
        
        # datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE],
        df <- DT::datatable(resultsSubset[, -c("rankAtMax", "TF_lowConf"), with = FALSE],
                            escape = FALSE, # To show the logo
                            filter = "top", options = list(pageLength = 5)
        )
        htmlwidgets::saveWidget(df, paste0(dir, "/Motif.html"), selfcontained = TRUE)
        
        write.table(data.frame("ID" = row.names(motifEnrichmentTable_wGenes), motifEnrichmentTable_wGenes), file = paste0(dir, "/onco_matrix1.txt"), sep = "\t", quote = F, row.names = F)
      }
      onco <- resultsSubset[, c("TF_highConf", "enrichedGenes"), with = FALSE]
      onco$TF_highConf <- gsub("\\s*\\(.*?\\.", "", onco$TF_highConf) # 去除括号及其内部的信息
      onco <- onco[onco$TF_highConf != "", ]
      result <- data.frame()
      for (i in 1:nrow(onco)) {
        tf <- strsplit(onco[i,]$TF_highConf, ";")[[1]]
        gene <- strsplit(onco[i,]$enrichedGenes, ";")[[1]]
        
        for (j in 1:length(gene)) {
          result <- rbind(result,data.frame(TF_highConf = tf,enrichedGenes = rep(gene[j], length(tf))))
        }
      }
      result <- result[duplicated(result) == F,]
      result$TF_highConf <- trimws(result$TF_highConf)
      result$enrichedGenes <- trimws(result$enrichedGenes)
      print(result)
      write.table(result, paste0(dir, "/Motif.txt"), quote = F, row.names = F, sep = "\t")
    }
    
    # 非编码 RNA 网络
    if (step == 8) {
      message("*************************非编码 RNA 网络*************************")
      
      dir <- paste0(Path, "/.Non-coding RNA networks")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      mircode <- data.table::fread("refdata/mircode_highconsfamilies.txt",header = T,data.table = F)
      key.mircode <- mircode[mircode$gene_symbol %in% KeyGene,]
      network.table <- data.frame(Gene = key.mircode$gene_symbol,microrna = key.mircode$microrna)
      network.table$microrna <- stringr::str_split(network.table$microrna,"/",simplify = T)[,1]
      network.table <- network.table[!duplicated(network.table),]
      write.table(network.table,paste0(dir,"/network.txt"),row.names = F,sep = "\t",quote = F)
    }
  }
}
