#' save_rds
#'
#' @details function to save data as .rds and create m5sum digest
#' 
#' @ based on https://stackoverflow.com/questions/65403794/compare-md5-of-an-object-in-memory-to-the-md5-of-it-as-an-rds
#' 
#' @param ADD
#' @param ADD
save_rds <- function(x, rdsname = paste0(substitute(x), ".rds"), digestname = paste0(rdsname,".digest")) {
  hash <- digest::digest(x, algo = "md5")
  saveRDS(x, rdsname)
  writeLines(hash, digestname)
}

#' digest_match
#'
#' @details function to check digest md5sum of object saved in .rds file
#' 
#' @ based on https://stackoverflow.com/questions/65403794/compare-md5-of-an-object-in-memory-to-the-md5-of-it-as-an-rds
#' 
#' @param ADD
#' @param ADD
digest_match <- function(x, digestfile) {
  hash <- digest::digest(x, algo = "md5")
  orig_hash <- readLines(digestfile)
  return( hash==orig_hash )
}

#' create_log2FC_rnk_list
#'
#' @details function to create ranked lists from list of raw DESeq result datasets
#'   ranking is done based on log2FC
#'   following columns are required: Gene.Name.HGNC, log2FoldChange, padj
#'
#' @param ADD
#' @param ADD
#' @return ADD
create_log2FC_rnk_list <- function(de_results_df=NULL,
                                   gene_name_column=NULL,
                                   log2FC_column=NULL,
                                   padj_column=NULL){
  import::here("%>%", .from = magrittr)
  import::here(select, filter, rename, arrange, desc, group_by, summarise, ungroup, .from = dplyr)
  
  
  gene_name_column <- as.name(gene_name_column)
  log2FC_column <- as.name(log2FC_column)
  padj_column <- as.name(padj_column)
  
  # filtering of de_results_df
  clean_data <- de_results_df %>%
    rename(gene_symbol = {{ gene_name_column }},
           log2FC = {{ log2FC_column }},
           padj = {{ padj_column }}) %>%
    filter(!is.na(padj)) %>%  # removing NA padj (DESeq filter for lowly/noise expressed genes)
    filter(!is.na(gene_symbol)) %>%  # removing NA gene symbols
    select(gene_symbol, log2FC) %>%
    group_by(gene_symbol) %>%
    summarise(log2FC = log2FC[abs(log2FC) == max(abs(log2FC))]) %>%  # duplicated gene_symbols keep one with max(abs(log2FC))
    #filter(!duplicated(gene_symbol)) %>%  # removing duplicated gene_symbols; randomly keeping one (check log2FC values!)
    ungroup() %>%
    arrange(desc(log2FC)) 
  
  # creating ranked list
  rnk_list <- clean_data$log2FC
  names(rnk_list) <- clean_data$gene_symbol
  
  return(rnk_list)
}


#' gsea_function
#'
#' @details gsea wrapper calling clusterProfiler::GSEA
#'
#' @param ADD
#' @param ADD
#' @return ADD
gsea_function <- function(geneset=NULL, geneset_name=NULL, ranked_list=NULL, ncores=6){
  import::here(SnowParam, .from = BiocParallel)
  import::here(GSEA, .from = clusterProfiler)
  # added geneset_name mainly for debugging and printing in log
  # pvalueCutoff = 1, no filtering at this stage
  # For some pathways, in reality P-values are less than 1e-10. You can set the `eps` argument to zero for better estimation.
  # consider eps = 0
  message("...enriching: ", geneset_name)
  GSEA(geneList = ranked_list,
       TERM2GENE = geneset,
       minGSSize = 10,
       maxGSSize = 500,
       pvalueCutoff = 1,
       seed = 123,
       verbose = TRUE,
       by = "fgsea",
       BPPARAM=SnowParam(workers = ncores),
       eps = 0,
       nPermSimple = 100000)
}

#' gseaPlotESprofile
#'
#' @details function based on enrichplot::gseaplot2 [all credit goes to the original author]
#'
#' @param ADD
#' @param ADD
#' @return ADD
gseaPlotESprofile <- function(x, geneSetID, title = "", color = "green", base_size = 11, 
                               rel_heights = c(1.5, 0.5, 1), subplots = 1:3, ES_geom = "line",
                               pvalue_table = FALSE, 
                               NES_table = FALSE,
                               remove_substring = NULL) {
  # function based on enrichplot::gseaplot2 [all credit goes to the original author]
  #  - adding option for NES_table that contains NES, pval, padj
  #  - added horizontal dashed line at yintercept = 0 for running ES
  #  - remove_substring = option to shorted geneset name(s) by fixed strink e.g. HALLMARK_ 
  
  # require(RColorBrewer) # brewer.pal
  # require(ggplot2)
  import::here(brewer.pal, .from = RColorBrewer)
  import::here(ggplot2)
  import::here(str_wrap, str_replace_all, .from = stringr)
  import::here(plot_grid, .from=cowplot)
  
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  geneList <- position <- NULL
  if (length(geneSetID) == 1) {
    gsdata <- gsInfo(x, geneSetID)
    
    # cleaning and wrapping long strings
    #gsdata$Description <- stringr::str_wrap(stringr::str_replace_all(gsdata$Description, pattern = "_", replacement = " "), width=20)
  }
  else {
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
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
  # adding horizontal line to distinguish between neg and pos running enrichment scores
  p.res <- p.res + geom_hline(yintercept = 0, color="grey", linetype = "dashed")
  
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
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
  if (length(geneSetID) == 1) {
    v <- seq(1, sum(gsdata$position), length.out = 9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) 
      inv <- inv + 1
    col = c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[unique(inv)])
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
  p.pos <- p.pos + ylab("Ranked list metric") + xlab("Rank in Ordered Dataset") + 
    theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, 
                               l = 0.2, unit = "cm"))
  if (!is.null(title) && !is.na(title) && title != "") 
    p.res <- p.res + ggtitle(title)
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
    pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
    names(pd) <- c("Description", "NES", "pval", "padj")
    pd <- pd[order(pd[, 1], decreasing = FALSE), ]
    rownames(pd) <- pd$Description
    pd <- pd[, -1]
    pd <- round(pd, 4)
    tp <- tableGrob2(pd, p.res)
    p.res <- p.res + theme(legend.position = "none") + 
      annotation_custom(tp, xmin = quantile(p.res$data$x, 0.5), 
                        xmax = quantile(p.res$data$x, 0.95), 
                        ymin = quantile(p.res$data$runningScore, 0.75), 
                        ymax = quantile(p.res$data$runningScore, 0.9))
  }
  
  if (NES_table) {
    pd <- x[geneSetID, c("Description", "NES", "pvalue", "p.adjust")]
    names(pd) <- c("Description", "NES", "pval", "padj")
    pd <- pd[order(pd[, 1], decreasing = FALSE), ]
    # adjust length of geneset names
    # make it dynamic
    # hardcoded removal of hallmark in string
    if (!is.null(remove_substring)) {
      #pd$Description <- gsub(pattern = "HALLMARK_", replacement = "", pd$Description) 
      pd$Description <- gsub(pattern = remove_substring, replacement = "", pd$Description) 
    }
    rownames(pd) <- str_wrap(str_replace_all(pd$Description, pattern = "_", replacement = " "), width=30) # original: pd$Description
    
    pd <- pd[, -1]
    #pd <- round(pd, 4)
    pd$NES <- round(pd$NES, 2)
    pd$pval <- round(pd$pval, 4)
    pd$padj <- round(pd$padj, 4)
    tp <- tableGrob2(pd, p.res)
    
    p.res <- p.res + theme(legend.position = "none") +
      annotation_custom(tp,
                        xmin = quantile(p.res$data$x, 0.5),
                        xmax = quantile(p.res$data$x, 0.95),
                        ymin = quantile(p.res$data$runningScore, 0.75),
                        ymax = quantile(p.res$data$runningScore, 0.9))
    
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

#' tableGrob2
#'
#' @details function based on enrichplot:::tableGrob2() [all credit goes to the original author]
#'
#' @param ADD
#' @param ADD
#' @return ADD
tableGrob2 <- function(d, p = NULL) {
  d <- d[order(rownames(d)),]
  tp <- gridExtra::tableGrob(d) # from gridExtra
  if (is.null(p)) {
    return(tp)
  }
  pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
  j <- which(tp$layout$name == "rowhead-fg")
  
  for (i in seq_along(pcol)) {
    tp$grobs[j][[i+1]][["gp"]] = grid::gpar(col = pcol[i])
  }
  return(tp)
}

#' tableGrob2
#'
#' @details function based on enrichplot:::gsInfo [all credit goes to the original author]
#'
#' @param ADD
#' @param ADD
#' @return ADD
gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin=0
  df$ymax=0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}

import:::here(gseaScores, .from = DOSE)
#gseaScores <- DOSE:::gseaScores # getFromNamespace("gseaScores", "DOSE")

#' DOSE:::gseaScores
#'
#' @details function based on DOSE:::gseaScores [all credit goes to the original author]
#'
#' @param ADD
#' @param ADD
#' @return ADD
# gseaScores <- function (geneList, geneSet, exponent = 1, fortify = FALSE) 
# {
#   geneSet <- intersect(geneSet, names(geneList))
#   N <- length(geneList)
#   Nh <- length(geneSet)
#   Phit <- Pmiss <- numeric(N)
#   hits <- names(geneList) %in% geneSet
#   Phit[hits] <- abs(geneList[hits])^exponent
#   NR <- sum(Phit)
#   Phit <- cumsum(Phit/NR)
#   Pmiss[!hits] <- 1/(N - Nh)
#   Pmiss <- cumsum(Pmiss)
#   runningES <- Phit - Pmiss
#   max.ES <- max(runningES)
#   min.ES <- min(runningES)
#   if (abs(max.ES) > abs(min.ES)) {
#     ES <- max.ES
#   }
#   else {
#     ES <- min.ES
#   }
#   df <- data.frame(x = seq_along(runningES), runningScore = runningES, 
#                    position = as.integer(hits))
#   if (fortify == TRUE) {
#     return(df)
#   }
#   df$gene = names(geneList)
#   res <- list(ES = ES, runningES = df)
#   return(res)
# }