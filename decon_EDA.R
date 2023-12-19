# ***********************************************
# Title       : Deconvolution project on bulk RNA-seq for Daria Lazic (Sabine Taschner-Mandl)
# Subtitle    : featureCounts
# Description : 
#  deconvolution using marker genes signatures
#  deconvolution using scRNA-seq
#
# Author      : Peter Repiscak (peter.repiscak@gmail.com)
# Date        : 25/02/19
# Version     : v1.0
# ***********************************************



# TO-DO:
# sequencing depth and other basic QC (also add QC from fastsQC etc.)
# PCA and variance partition of key variables

# Performing analysis for material and DE subset for material
# (1|material) + (1|rnaseq_batch) 
# disease, performed_at - ignore, 
# rnaseq_batch confounding with library_type
experiment_design="1"
cond_interest <- "material" # column of interest for PCA plots etc.
cond_interest_varPart <- c("library_type", "rnaseq_batch", "material", "gender", "timepoint")
fitform_partVariance <- ~ (1|library_type) + (1|rnaseq_batch) + (1|material) + (1|gender) + (1|timepoint) + (1|patient_oid)  
padj_cutoff = 0.05
log2FC_cutoff = log2(2) #0.58 #(FC=1.5); log2FC=1.0 # (FC=2)
abs_filt_samples <- 5          # filtering for expression in at least N samples (defaults to 3)
var_expl_needed <- 0.6         # at least 60% variance explained needed

generateRNAobjects=TRUE           # Generate dds, log2_norm, vsd, rld objects or load them if precalculated
EDA_QC=TRUE                   # Perform exploratory data analysis (EDA) and Quality control

# library(skimr)
# #install.packages("inspectdf")
# my_skim <- skimr::skim_with(factor = skimr::sfl(top_counts = ~top_counts(., max_char = 100, max_levels = 6)))
# my_skim(select(coldata, cond_interest_varPart))
# 
# inspectdf::inspect_imb(coldata) %>% inspectdf::show_plot(col_palette = 1)
# inspectdf::inspect_types(coldata) %>% inspectdf::show_plot(col_palette = 1)
# inspectdf::inspect_cat(coldata) %>% inspectdf::show_plot(col_palette = 1)
# inspectdf::inspect_cat(coldata) %>% inspectdf::show_plot(high_cardinality = 1, col_palette = 1, label_thresh = 0.01)
# inspectdf::inspect_num(coldata) #%>% inspectdf::show_plot()

# again multiple samples from the same patient - check if only timepoint DE is used!

##### END: Experiment description ----

# Do not change following if defaults apply
#GENOME_DIR="/media/prepisca/DATA/Work/Bioinformatics/GenomesAnnotations/Transcriptome/ensembl_v92/" # directory that contains annotation files, transcriptome files etc.


# Run control options ----


# Loading functions ----
# [ ] replace with packages made from this functions

# REMOVE: 
# plot_pca <- function(transf_object=NULL, intgroup=c("condition", "batch")){
#   # Generating PCA plot 
#   message("first entry in intgroup taken for color and second for shape, third (if specified) as labels")
#   
#   pcaData <- plotPCA(transf_object, intgroup=intgroup, returnData=TRUE)
#   percentVar <- round(100 * attr(pcaData, "percentVar"))
#   
#   # taking first entry as color and second as shape
#   pca_plot <- ggplot(pcaData, aes_string("PC1", "PC2", color=intgroup[[1]], shape=intgroup[[2]])) +
#     geom_point(size=3) +
#     xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#     ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#     coord_fixed() + 
#     theme_bw() +
#     guides(colour = guide_legend(order = 1), 
#            shape = guide_legend(order = 2))
#   
#   if(length(intgroup) == 3){
#     pcaData <- plotPCA(transf_object, intgroup=intgroup, returnData=TRUE)
#     percentVar <- round(100 * attr(pcaData, "percentVar"))
#     pca_plot <- ggplot(pcaData, aes_string("PC1", "PC2", color=intgroup[[1]], shape=intgroup[[2]], label=intgroup[[3]])) +
#       geom_point(size=3) +
#       xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#       ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#       coord_fixed() + 
#       ggrepel::geom_label_repel() +
#       theme_bw() +
#       guides(colour = guide_legend(order = 1), 
#              shape = guide_legend(order = 2))
#   }
#   
#   return(pca_plot)
#   
# }


# functions ----
# REMOVE: 
# remove_batch <- function(vst_transf = NULL, batch_variable=NULL, batch_variable2=NULL, cond_interest = NULL){
#   # account for design if unbalanced
#   # For unbalanced batches (e.g. the condition groups are not distributed balanced across batches), the design argument should be used
#   # https://support.bioconductor.org/p/122069/
#   # explicitly add package::function call
#   message("batch_variable: ", batch_variable, " cond_interest: ", cond_interest)
#   
#   vst_transf_batch <- vst_transf
#   vst_transf_batch_count <- assay(vst_transf_batch)
#   if (!is.null(cond_interest)) {
#     message("using design when removing batch effects...")
#     design0 <- model.matrix(~vst_transf_batch[[cond_interest]])
#     
#     if (!is.null(batch_variable2)) {
#       message("batch_variable2: ", batch_variable2)
#       vst_transf_batch_count_NObatch <- limma::removeBatchEffect(x = vst_transf_batch_count, 
#                                                                  batch = vst_transf_batch[[batch_variable2]],
#                                                                  batch2 = vst_transf_batch[[batch_variable]],
#                                                                  design=design0)
#     } else {
#       vst_transf_batch_count_NObatch <- limma::removeBatchEffect(x = vst_transf_batch_count, 
#                                                                  batch = vst_transf_batch[[batch_variable]],
#                                                                  design=design0)
#     }
#     
#   } else {
#     if (!is.null(batch_variable2)) {
#       message("batch_variable2: ", batch_variable2)
#       message("NOT using design when removing batch effects... Are condition groups balanced across batches?")
#       vst_transf_batch_count_NObatch <- limma::removeBatchEffect(x = vst_transf_batch_count, 
#                                                                  batch = vst_transf_batch[[batch_variable]],
#                                                                  batch2 = vst_transf_batch[[batch_variable]])
#     } else {
#       message("NOT using design when removing batch effects... Are condition groups balanced across batches?")
#       vst_transf_batch_count_NObatch <- limma::removeBatchEffect(x = vst_transf_batch_count, 
#                                                                  batch = vst_transf_batch[[batch_variable]])
#       
#     }
#     
#   }
#   
#   vst_transf_NObatch <- vst_transf_batch
#   assay(vst_transf_NObatch) <- vst_transf_batch_count_NObatch
#   
#   return(vst_transf_NObatch)
# }

partition_variance <- function(vst_transf = NULL, fitform_partVariance = NULL, ncores = 10, max_ntop_genes=2000){
  require("dplyr")
  require("variancePartition")
  require('doParallel')
  require("SummarizedExperiment") # for assay
  
  # requires also rnaSelectTopVarGenes to be loaded!!! (fix and include in separate package)
  
  if (exists("rnaSelectTopVarGenes")) {
    
    # max_ntop_genes - this is maximum most variable genes to consider, but also doing 500 genes by default;
    
    # load library
    # optional step to run analysis in parallel on multicore machines
    # Here use 4 threads
    # This is strongly recommended since the analysis
    # can be computationally intensive
    
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    
    # load simulated data:
    # geneExpr: matrix of gene expression values
    # info: information/metadata about each sample
    # Specify variables to consider
    # Age is continuous so model it as a fixed effect
    # Individual and Tissue are both categorical,
    # so model them as random effects
    # Note the syntax used to specify random effects
    #form <- ~ Age + (1|Individual) + (1|Tissue) + (1|Batch)
    
    # Interpretation:
    #These conclusions are based on the genome-wide median across all genes, but
    #the same type of statements can be made at the gene-level. Moreover, care
    #must be taken in the interpretation of nested variables. For example, Age is
    #nested within Individual since the multiple samples from each individual are
    #taken at the same age. Thus the effect of Age removes some variation from
    #being explained by Individual. 
    #This nesting/summing of effects is common for variables that are properties of the
    #individual rather than the sample.
    
    # Fit model and extract results
    # 1) fit linear mixed model on gene expression
    # If categorical variables are specified,
    # a linear mixed model is used
    # If all variables are modeled as fixed effects,
    # a linear model is used
    # each entry in results is a regression model fit on a single gene
    # 2) extract variance fractions from each model fit
    # for each gene, returns fraction of variation attributable
    # to each variable
    # Interpretation: the variance explained by each variables
    # after correcting for all other variables
    # Note that geneExpr can either be a matrix,
    # and EList output by voom() in the limma package,
    # or an ExpressionSet
    
    # Pre-filter genes by removing non-expressing ones
    #geneExpr_partVariance <- rld_counts_filt 
    geneExpr_partVariance <- rnaSelectTopVarGenes(assay(vst_transf), ntop = max_ntop_genes) #using top 500 most variable genes, alternatively using all genes rld_counts
    info_partVariance <- as.data.frame(colData(vst_transf))
    fitform_partVariance <- fitform_partVariance
    varPart_fit <- variancePartition::fitExtractVarPartModel(geneExpr_partVariance, fitform_partVariance, info_partVariance)
    #stopCluster(cl)
    
    # sort variables (i.e. columns) by median fraction
    # of variance explained
    varPart_fit_sorted <- variancePartition::sortCols( varPart_fit ) #showMethods("sortCols")
    
    # Bar plot of variance fractions for the first 10 genes
    #plotPercentBars( vp[1:10,] )
    
    # violin plot of contribution of each variable to total variance
    varPart_plot <- variancePartition::plotVarPart(varPart_fit_sorted)
    
    # Calculate medians
    varPart_stats <- varPart_plot$data %>%
      dplyr::select(variable, value) %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(median_varExplained = median(value),
                       mean_varExplained = mean(value),
                       IQR_varExplained = IQR(value),
                       max_varExplained = max(value),
                       min_varExplained = min(value)) %>%
      dplyr::mutate_if(is.numeric, round, 1)
    
    varPart_plot_annot <- varPart_plot + geom_text(data = varPart_stats, aes(x = variable, y = 100, label = paste0("mean: ", mean_varExplained)), 
                                                   size = 3, vjust = -0.5)
    
    #Assess correlation between all pairs of variables
    #fitform_partVariance_corr <- ~ mouse_gender + date_dissected + rin + condition # formula for the partition variance corr matrix
    #fitform_partVariance_corr <- fitform_partVariance_corr
    # Compute Canonical Correlation Analysis (CCA)
    # between all pairs of variables
    # returns absolute correlation value
    #varPart_corr = canCorPairs( fitform_partVariance_corr, info_partVariance)
    # Plot correlation matrix
    #varPart_corr_plot <- plotCorrMatrix(varPart_corr)
    # plot squared correlations
    #plotCorrMatrix( varPart_corr^2, dendrogram="none" )
    
    # Check if colinearity is a problem
    #Alternatively, the user can use the colinearityScore function to evaluate
    #whether this is an issue for a single model fit:
    
    # fit model
    #varPart <- fitVarPartModel( geneExpr_partVariance[1:500,], fitform_partVariance, info_partVariance )
    
    # Top500 most variable genes
    geneExpr_partVariance500 <- rnaSelectTopVarGenes(assay(vst_transf), ntop = 500) #using top 500 most variable genes, alternatively using all genes rld_counts
    varPart_fit500 <- variancePartition::fitExtractVarPartModel(geneExpr_partVariance500, fitform_partVariance, info_partVariance)
    stopCluster(cl)
    
    # sort variables (i.e. columns) by median fraction
    # of variance explained
    varPart_fit_sorted500 <- variancePartition::sortCols( varPart_fit500 ) #showMethods("sortCols")
    
    # Bar plot of variance fractions for the first 10 genes
    #plotPercentBars( vp[1:10,] )
    
    # violin plot of contribution of each variable to total variance
    varPart_plot500 <- variancePartition::plotVarPart(varPart_fit_sorted500)
    
    # Calculate medians
    varPart_stats500 <- varPart_plot500$data %>%
      dplyr::select(variable, value) %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(median_varExplained = median(value),
                       mean_varExplained = mean(value),
                       IQR_varExplained = IQR(value),
                       max_varExplained = max(value),
                       min_varExplained = min(value)) %>%
      dplyr::mutate_if(is.numeric, round, 1)
    
    varPart_plot_annot500 <- varPart_plot500 + geom_text(data = varPart_stats500, aes(x = variable, y = 100, label = paste0("mean: ", mean_varExplained)), 
                                                         size = 3, vjust = -0.5)
    
    varPart_results <- list(varPart_fit = varPart_fit,
                            varPart_stats = varPart_stats,
                            varPart_plot_annot =varPart_plot_annot,
                            varPart_fit500 = varPart_fit500,
                            varPart_stats500 = varPart_stats500,
                            varPart_plot_annot500 = varPart_plot_annot500)
    
  } else{
    stop("rnaSelectTopVarGenes function needs to be loaded separately!")
  }
  
  return(varPart_results)
}

# Dataset preparation functions ----
# constructing DESeq2 dds object from Salmon count files, condition table and design
# loading support scripts ----
datasetPreparation_scriptsDIR="/core_bioinformatics/scripts/datasetPreparation/"
#datasetPreparation_scriptsDIR="/media/prepisca/DATA/Work/Bioinformatics/Packages/build_functions_R/R/datasetPreparation/"
# source(paste0(datasetPreparation_scriptsDIR, "generateDatasets.R"))

# Exploratory data analysis functions ----
EDA_scriptsDIR="/core_bioinformatics/scripts/EDAfunctions/"


# generating results 
RNAseqPipeline_scriptsDIR="/core_bioinformatics/scripts/RNAseqPipeline/"
# Gene expression plotting
source(paste0(RNAseqPipeline_scriptsDIR, "generateResults.R")) # change contains(description) from contains(description_) -> upload to github
source(paste0(RNAseqPipeline_scriptsDIR, "extractLimmaResults.R"))
source(paste0(RNAseqPipeline_scriptsDIR, "compareConstrasts.R"))
source(paste0(RNAseqPipeline_scriptsDIR, "plotConstrasts.R"))
source(paste0(RNAseqPipeline_scriptsDIR, "plotConstrasts.R"))

# functional enrichment
FEscripts_scriptsDIR="/core_bioinformatics/scripts/funcionalEnrichment/"

# fetch org. database
source(paste0(FEscripts_scriptsDIR, "getOrgDB.R"))

# Run clusterProfiler enrichment
# source(paste0(FEscripts_scriptsDIR, "cluster_profiler/enrichSpecies.R"))

# Gene expression plotting
# source(paste0(FEscripts_scriptsDIR, "cluster_profiler/enrichGOmodule.R"))
# # enrichment plots ----
# source(paste0(FEscripts_scriptsDIR, "cluster_profiler/enrichmentBarplot.R"))
# source(paste0(FEscripts_scriptsDIR, "cluster_profiler/enrichmentDotplot.R"))
# source(paste0(FEscripts_scriptsDIR, "cluster_profiler/enrichmentCnetplot.R"))

source(paste0(FEscripts_scriptsDIR, "generateRankedList.R"))

# GSEA module
#source(paste0(FEscripts_scriptsDIR, "gsea/runGSEA.R"))
#source(paste0(FEscripts_scriptsDIR, "gsea/gseaPlotESprofile.R"))


# creating gsea genesets ----
run_this=FALSE
if (run_this) {
if(!file.exists(file.path(EXPERIMENT_DIR, "msigdbr_gs_entreIDs_collection.rds"))) {
  msigdbr::msigdbr_species()
  packageVersion("msigdbr")
  msigdbr::msigdbr_collections()
  gs_hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = c("H")) %>%
    dplyr::select(gs_name, entrez_gene)
  gs_C1 <- msigdbr::msigdbr(species = "Homo sapiens", category = c("C1")) %>%
    dplyr::select(gs_name, entrez_gene)
  gs_C2_kegg <- msigdbr::msigdbr(species = "Homo sapiens", category = c("C2"), subcategory = "CP:KEGG") %>%
    dplyr::select(gs_name, entrez_gene)
  gs_C2_reactome <- msigdbr::msigdbr(species = "Homo sapiens", category = c("C2"), subcategory = "CP:REACTOME") %>%
    dplyr::select(gs_name, entrez_gene)
  gs_C5_GOBP <- msigdbr::msigdbr(species = "Homo sapiens", category = c("C5"), subcategory = "GO:BP") %>%
    dplyr::select(gs_name, entrez_gene)
  gs_C5_GOCC <- msigdbr::msigdbr(species = "Homo sapiens", category = c("C5"), subcategory = "GO:CC") %>%
    dplyr::select(gs_name, entrez_gene)
  gs_C5_GOMF <- msigdbr::msigdbr(species = "Homo sapiens", category = c("C5"), subcategory = "GO:MF") %>%
    dplyr::select(gs_name, entrez_gene)
  gs_C5_HPO <- msigdbr::msigdbr(species = "Homo sapiens", category = c("C5"), subcategory = "HPO") %>%
    dplyr::select(gs_name, entrez_gene)
  # gs_C7 <- msigdbr::msigdbr(species = "Homo sapiens", category = c("C7")) %>%
  #   dplyr::select(gs_name, entrez_gene)
  # gs_C8 <- msigdbr::msigdbr(species = "Homo sapiens", category = c("C8")) %>%
  #   dplyr::select(gs_name, entrez_gene)
  
  gs_entreIDs_collection <- list(hallmark = gs_hallmark,
                                 C1 = gs_C1,
                                 C2_kegg = gs_C2_kegg,
                                 C2_reactome = gs_C2_reactome,
                                 C5_GOBP = gs_C5_GOBP,
                                 C5_GOCC = gs_C5_GOCC,
                                 C5_GOMF = gs_C5_GOMF,
                                 C5_HPO = gs_C5_HPO)
  
  #, C7 = gs_C7,
  #C8 = gs_C8
  
  saveRDS(gs_entreIDs_collection, file = file.path(EXPERIMENT_DIR, "msigdbr_gs_entreIDs_collection.rds"))
} else {
  gs_entreIDs_collection <- readRDS(file.path(EXPERIMENT_DIR, "msigdbr_gs_entreIDs_collection.rds"))
}
}

# Preparing datasets ----
# saving key data so it does not need to be prepared again
save(cond_data, count_data,
     file = paste0(OUTPUT_DATASET_DIR, experiment_name, "_coldata_count_data.RData"))

load(file = paste0(OUTPUT_DATASET_DIR, experiment_name, "_coldata_count_data.RData"))

count_data_ensID <- count_data %>%
  tibble::rownames_to_column(var = "ensembl_id")

readr::write_tsv(x = cond_data, file = paste0(OUTPUT_DATASET_DIR, experiment_name, "_metadata.tsv"))  # does not save rownames
readr::write_tsv(x = count_data_ensID, file = paste0(OUTPUT_DATASET_DIR, experiment_name, "_countdata.tsv"))

digest::digest(file = paste0(OUTPUT_DATASET_DIR, experiment_name, "_metadata.tsv"), algo = "sha256")
digest::digest(file = paste0(OUTPUT_DATASET_DIR, experiment_name, "_countdata.tsv"), algo = "sha256")

# reading in files
# save also coltypes!! or specify here!
# coldata <- readr::read_tsv(file = paste0(OUTPUT_DATASET_DIR, experiment_name, "_metadata.tsv"))
# count_data <- readr::read_tsv(file = paste0(OUTPUT_DATASET_DIR, experiment_name, "_countdata.tsv"))
# 
# rownames(coldata) <- coldata$omics_id
#all(rownames(coldata) == colnames(count_data))
# 
# cond_data <- coldata

# prepare_dds <- function(){
#   
# }

if (generateRNAobjects) {
  # Aim: prepare initial dds object
  # [] make running of rlog, vsd optional as this may change if EDA_QC discovers that additional variables need to be added to design!
  # Generate dds, log2_norm, vsd, rld objects or load them if precalculated
  # Checking OUTPUT_DATASET_DIR for dds file existence
  
  # Check if dds file exists
  if(!file.exists(paste0(OUTPUT_DATASET_DIR, experiment_name, "_dds_objects.RData"))){
    # Generating dds and filtered dds objects ----
    dds <- DESeqDataSetFromMatrix(countData = count_data,
                                  colData = coldata,
                                  design= as.formula(paste0("~", experiment_design))) 
      
    dds <- estimateSizeFactors(dds)
    #design(dds) <- as.formula(~condition)
    dds <- DESeq(dds, minReplicatesForReplace=Inf) # do not replace outliers based on replicates

    ensemblAnnot <- generateEnsemblAnnotation(ensembl_ids = rownames(dds),
                                              host="http://nov2020.archive.ensembl.org",
                                              version="Ensembl Genes 102",
                                              dataset="hsapiens_gene_ensembl")
    
    # Filtering
    dds_filt <- filterDatasets(dds, abs_filt = TRUE, abs_filt_samples = abs_filt_samples)# at least in N samples, which is a smallest group size
    dds_filt <- estimateSizeFactors(dds_filt)
    
    save(dds, dds_filt, ensemblAnnot, file = paste0(OUTPUT_DATASET_DIR, experiment_name, "_dds_objects.RData"))
    
  } else{
    load(paste0(OUTPUT_DATASET_DIR, experiment_name, "_dds_objects.RData"))
  }
}

#transformation after filtering
if(!file.exists(paste0(OUTPUT_DATASET_DIR, experiment_name, "_log2_vst_filt.RData"))){
  # transformation
  log2_norm_filt <- normTransform(dds_filt)
  vsd_filt <- vst(dds_filt, blind = TRUE) # blind = TRUE for QC
  rld_filt <- vsd_filt # replacing; but rewrite code
  #rld_filt <- rlog(dds_filt, blind = TRUE)
  #rld_filt <- rlog(dds_filt, blind = TRUE) # not blind to batch effects; rlog too computationally expensive
  save(log2_norm_filt, vsd_filt, file = paste0(OUTPUT_DATASET_DIR, experiment_name, "_log2_vst_filt.RData"))
} else{
  load(paste0(OUTPUT_DATASET_DIR, experiment_name, "_log2_vst_filt.RData"))
}

# by default usign rlog transformation as it should be more robust for varying size factors, but it is much slower than vst!
# using vst - remove confusing rld_counts_filt!!! (rename just to transf variable!)
rld_counts_filt <- assay(vsd_filt) # rld_filt # assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2

if (EDA_QC) {
  # some basic assumptions about condition and filenames!!!
  # [ ] check if other assumptions are needed!
  
  SequencingDepth=TRUE
  NormalizationCheck=TRUE
  Transformation_check=FALSE
  Additional_check=TRUE
  
  # Sequencing depth histogram and size factor histogram ----
  if (SequencingDepth) {
    cat("Checking sequencing depth...\n")
    # Plotting sequencing depth from raw counts
    # Plotting size factor distribution - create FUNCTION

    size_factors <- as.data.frame(sizeFactors(dds)) %>% dplyr::rename(size_factors = "sizeFactors(dds)")# extractSizeFactors(dds_object = dds)
    barfill <- "#4271AE"
    barlines <- "#1F3552"
    size_factors_plot <- ggplot(size_factors, aes(x = size_factors)) +
      geom_histogram(aes(y = ..count..), binwidth = 0.05,
                     colour = barlines, fill = barfill) +
      #scale_x_continuous(name = "Size factors",
      #                   breaks = seq(0.7, 1.31, 0.1), 
      #                   limits=c(0.7, 1.31)) +
      scale_y_continuous(name = "Count") +
      ggtitle(paste0("Histogram of size factors for ", experiment_name," ", ncol(dds), " samples")) + 
      geom_vline(color = "red", linetype = 2, xintercept = 1.0) +
      theme_minimal()
    
    #https://cran.r-project.org/web/packages/ggplotify/vignettes/ggplotify.html
    #SizeFactors_histogram <- ggplotify::as.ggplot(ggplotify::as.grob(SizeFactors_histogram))
      # plotSizeFactors(dds_object = dds, size_factors=size_factors, experiment_name = experiment_name)
    sequencingDepthCheck <- function(dds_object=NULL){
      require(DESeq2)
      require(SummarizedExperiment)
      require(dplyr)
      require(ggplot2)
      
      # need to change condition to respective condition: e.g. material
      
      if (!is.null(dds_object)) {
        
        metadata <- as.data.frame(SummarizedExperiment::colData(dds_object))
        
        # returns data.frame and plots for library size depth before and after normalization
        raw_counts_dds <- base::colSums(DESeq2::counts(dds_object, normalized=FALSE))
        raw_counts_df <- data.frame(sample_name=names(raw_counts_dds), 
                                    library_size=raw_counts_dds) %>%
          dplyr::mutate(library_size_M = library_size/(1e+06))
        
        norm_counts_dds <- base::colSums(DESeq2::counts(dds_object, normalized=TRUE))
        norm_counts_df <- data.frame(sample_name=names(norm_counts_dds), 
                                     library_size=norm_counts_dds) %>%
          dplyr::mutate(library_size_M = library_size/(1e+06))
        
        raw_counts_libSize <- raw_counts_df %>%
          dplyr::left_join(., metadata, by=c("sample_name"="omics_id")) %>%
          #dplyr::mutate(sample_name = factor(sample_name, levels=unique(sample_name))) %>%
          dplyr::arrange(material) %>%
          dplyr::mutate(sample_name = factor(sample_name, levels = sample_name))

        norm_counts_libSize <- norm_counts_df %>%
          dplyr::left_join(., metadata, by=c("sample_name"="omics_id")) %>%
          #dplyr::mutate(sample_name = factor(sample_name, levels=unique(sample_name))) %>%
          dplyr::arrange(material) %>%
          dplyr::mutate(sample_name = factor(sample_name, levels = sample_name))
        
        raw_counts_libSize_plot <- ggplot(raw_counts_libSize, aes(x=sample_name, y=(library_size)/10^6)) + 
          geom_bar(stat="identity", aes(fill=material)) +
          theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          geom_hline(yintercept = 25, linetype="dashed", colour="red") +
          labs(title = "Raw counts library depth", x="Sample names", y="Library size [millions reads]")
        
        norm_counts_libSize_plot <- ggplot(norm_counts_libSize, aes(x=sample_name, y=(library_size)/10^6)) + 
          geom_bar(stat="identity", aes(fill=material)) +
          theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          geom_hline(yintercept = 25, linetype="dashed", colour="red") +
          labs(title = "Normalized counts library depth", x="Sample names", y="Library size [millions reads]")
        
        librarySize_merged_plot <- ggpubr::ggarrange(raw_counts_libSize_plot,
                                                     norm_counts_libSize_plot,
                                                     common.legend = TRUE)
        
        library_size_data <- list(raw_counts_libSize=raw_counts_libSize,
                                  norm_counts_libSize=norm_counts_libSize,
                                  raw_counts_libSize_plot=raw_counts_libSize_plot,
                                  norm_counts_libSize_plot=norm_counts_libSize_plot,
                                  librarySize_merged_plot=librarySize_merged_plot)
        
        return(library_size_data)
        
      } else {
        stop("Need to provide dds object!")
      }
    }
    library_size_data <- sequencingDepthCheck(dds_object = dds)
    # go for at least 5M reads!
    # quantile(library_size_data$raw_counts_libSize$library_size_M, probs = seq(0,1,0.1))
    
    # saving files
    ggsave(filename = paste0(FIGURE_DIR, "SequencingDepth_sizeFactors.png"), plot=size_factors_plot, width = 30, height = 25, units = "cm") 
    ggsave(filename = paste0(FIGURE_DIR, "SequencingDepth.png"), plot=library_size_data$librarySize_merged_plot, width = 30, height = 25, units = "cm")
    ggsave(filename = paste0(FIGURE_DIR, "SequencingDepth_raw.png"), plot=library_size_data$raw_counts_libSize_plot, width = 30, height = 25, units = "cm")
    ggsave(filename = paste0(FIGURE_DIR, "SequencingDepth_norm.png"), plot=library_size_data$norm_counts_libSize_plot, width = 30, height = 25, units = "cm")
    save(size_factors, size_factors_plot, library_size_data, file = paste0(OUTPUT_DIR, experiment_name, "_SeqDepth.RData"))
    
    min_depth=5.0 # 5M
    # boxplot(library_size_data$raw_counts_libSize$library_size_M)
    low_depth_samples <- library_size_data$raw_counts_libSize %>%
      dplyr::filter(library_size_M < min_depth)
    openxlsx::write.xlsx(library_size_data$raw_counts_libSize, file = paste0(OUTPUT_DIR, "sequencing_depth.xlsx")) 
    
    cat("Done\n")
  }
  
  # Check normalization ----
  if (NormalizationCheck) {
    cat("Checking normalization...\n")
    # Plotting normalization effect
    normalizationCheck <- function(dds_object=NULL, dds_filt_object=NULL){
      require(DESeq2)
      require(SummarizedExperiment)
      require(dplyr)
      require(ggplot2)
      require(ggpubr) # for arranging multiple plots
      
      # [ ] add check for dds_filt_object!
      
      if (!is.null(dds_object)) {
        
        metadata <- as.data.frame(SummarizedExperiment::colData(dds_object))
        
        filenames_count_annot <- metadata %>%
          dplyr::select(omics_id, material) 

        # Extracting raw and normalised read counts ----
        raw_counts <- as.data.frame(counts(dds, normalized = FALSE)) %>% 
          tidyr::gather(., key = "omics_id", value = "count") %>%
          dplyr::left_join(filenames_count_annot, by = "omics_id") %>%
          dplyr::arrange(material) %>%
          dplyr::mutate(omics_id = factor(omics_id, levels = unique(omics_id)))
          # dplyr::mutate(omics_id = factor(omics_id, levels=unique(omics_id)))
        
        norm_counts <- as.data.frame(counts(dds, normalized = TRUE)) %>% 
          tidyr::gather(., key = "omics_id", value = "count") %>%
          dplyr::left_join(filenames_count_annot, by = "omics_id") %>%
          dplyr::arrange(material) %>%
          dplyr::mutate(omics_id = factor(omics_id, levels = unique(omics_id)))
          #dplyr::mutate(omics_id = factor(omics_id, levels=unique(omics_id)))
        
        norm_counts_filt <- as.data.frame(counts(dds_filt, normalized = TRUE)) %>% 
          tidyr::gather(., key = "omics_id", value = "count") %>%
          left_join(filenames_count_annot, by = "omics_id") %>%
          dplyr::arrange(material) %>%
          dplyr::mutate(omics_id = factor(omics_id, levels = unique(omics_id)))
          #dplyr::mutate(omics_id = factor(omics_id, levels=unique(omics_id)))
        
        .countDensityPlot <- function(count_data=NULL, title="Change title!"){
          require(ggplot2)
          
          ggplot(count_data, aes(x = log2(count + 1), colour = omics_id)) + 
            geom_density(alpha=.3) + ggtitle(title) +  
            #ggsci::scale_color_jco() +
            theme_bw()
        }
        
        .countBoxPlot <- function(count_data=NULL, title="Change title!"){
          require(ggplot2)
          
          ggplot(count_data, aes(x = omics_id, y = log2(count + 1), fill = material)) + 
            geom_boxplot() + ggtitle(title) + 
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
          
        }
        
        raw_counts_boxplot <- .countBoxPlot(count_data=raw_counts, title="Raw counts")
        norm_counts_boxplot <- .countBoxPlot(count_data=norm_counts, title="Normalized counts")
        normFilt_counts_boxplot <- .countBoxPlot(count_data=norm_counts_filt, title="Normalized counts filtered")

        raw_counts_density <- .countDensityPlot(count_data=raw_counts, title="Raw counts")
        norm_counts_density <- .countDensityPlot(count_data=norm_counts, title="Normalized counts")  # remove sample names or otherwise mark! (e.g. material or other group)
        normFilt_counts_density <- .countDensityPlot(count_data=norm_counts_filt, title="Normalized counts filtered")

        norm_raw_counts_boxplot_comparison <- ggpubr::ggarrange(raw_counts_boxplot, 
                                                                norm_counts_boxplot, 
                                                                normFilt_counts_boxplot, 
                                                                common.legend = TRUE, legend="bottom")
        norm_raw_counts_density_comparison <- ggpubr::ggarrange(raw_counts_density, 
                                                                norm_counts_density, 
                                                                normFilt_counts_density, 
                                                                common.legend = TRUE, legend="bottom")
   
        norm_effect_data <- list(norm_raw_counts_boxplot_comparison=norm_raw_counts_boxplot_comparison, 
                                 norm_raw_counts_density_comparison=norm_raw_counts_density_comparison)
        
        return(norm_effect_data)
        
      } else {
        stop("Need to provide dds object!")
      }
    }
    
    norm_effect_data <- normalizationCheck(dds_object=dds, dds_filt_object=dds_filt)
    
    ggsave(filename = paste0(FIGURE_DIR, "NormEffect_boxplot.png"), plot=norm_effect_data$norm_raw_counts_boxplot_comparison, width = 30, height = 30, units = "cm")
    ggsave(filename = paste0(FIGURE_DIR, "NormEffect_density.png"), plot=norm_effect_data$norm_raw_counts_density_comparison, width = 30, height = 30, units = "cm") 
    
    save(norm_effect_data, file = paste0(OUTPUT_DIR, experiment_name, "_NormEffect.RData"))
    cat("Done\n")
  }
  
  if (Transformation_check) {
    cat("Checking transformation...\n")
    # compare effects of transformation 
    transf_effect_data <- transformationCheck(dds_filt_object=dds_filt, vsd_filt_object=vsd_filt) #, rld_filt_object=rld_filt
    
    ggsave(filename = paste0(FIGURE_DIR, "TransfEffect_boxplot.pdf"), plot=transf_effect_data$norm_vsd_rlog_counts_boxplot_comparison)
    ggsave(filename = paste0(FIGURE_DIR, "TransfEffect_msdPlot.pdf"), plot=transf_effect_data$msd_merged_plot) 
    save(transf_effect_data, file = paste0(OUTPUT_DIR, experiment_name, "_TransfEffect.RData"))
    
    cat("Done\n")
  }
  
  if(Additional_check){
    cat("Checking cooks distance and pval distribution...\n")
    additionalCheck <- function(dds_object=NULL, dds_filt_object=NULL){
      require(DESeq2)
      require(dplyr)
      require(ggplot2)
      require(ggpubr) # for arranging multiple plots
      
      # [ ] add check for dds_filt_object!
      
      if (!is.null(dds_object)) {
        
        metadata <- as.data.frame(SummarizedExperiment::colData(dds_object))
        
        filenames_count_annot <- metadata %>%
          dplyr::select(omics_id, material) 
        
        # Cook's distance ----
        cooksDist_dds <- as.data.frame(log10(assays(dds_object)[["cooks"]])) %>% 
          tidyr::gather(., key = "omics_id", value = "cooks_dist") %>%
          dplyr::left_join(filenames_count_annot, by = "omics_id") %>%
          dplyr::mutate(omics_id = factor(omics_id, levels=unique(omics_id)))
        
        cooksDist_dds_filt <- as.data.frame(log10(assays(dds_filt_object)[["cooks"]])) %>% 
          tidyr::gather(., key = "omics_id", value = "cooks_dist") %>%
          dplyr::left_join(filenames_count_annot, by = "omics_id") %>%
          dplyr::mutate(omics_id = factor(omics_id, levels=unique(omics_id)))
        
        # plotting
        cooksDist_dds_boxplot <- ggplot(cooksDist_dds, aes(x = omics_id, y = cooks_dist, fill = material)) + 
          geom_boxplot() + ggtitle("Cook's distance dds") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ ylim(0, 20)
        
        cooksDist_dds_filt_boxplot <- ggplot(cooksDist_dds_filt, aes(x = omics_id, y = cooks_dist, fill = material)) + 
          geom_boxplot() + ggtitle("Cook's distance dds filtered") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ ylim(0, 20)
        
        #boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
        
        cooksDist_comparison <- ggpubr::ggarrange(cooksDist_dds_boxplot, 
                                                  cooksDist_dds_filt_boxplot,
                                                  common.legend = TRUE, legend="bottom")
        
        # p-value checks ----
        res_dds <- as.data.frame(results(dds))
        res_dds_filt <- as.data.frame(results(dds_filt))
        
        pval_hist <- ggplot(res_dds, aes(x = pvalue)) +
          theme_minimal() +
          geom_histogram(binwidth = 0.025, boundary = 0 ) +
          ggtitle("p-value histogram") #+ geom_hline(yintercept = 40, color = "chartreuse3")
        
        pval_hist_filt <- ggplot(res_dds_filt, aes(x = pvalue)) +
          theme_minimal() +
          geom_histogram(binwidth = 0.025, boundary = 0 ) +
          ggtitle("p-value histogram filtered") # + geom_hline(yintercept = 40, color = "chartreuse3")
        
        pval_hist_comparison <- ggpubr::ggarrange(pval_hist, 
                                                  pval_hist_filt,
                                                  common.legend = TRUE, legend="bottom")
        
        # Dispersion plot ----

        finished_ggDispPlot=FALSE
        if(finished_ggDispPlot){
          ggplotDispEsts <- function(object, ymin, CV = FALSE, genecol = "black", fitcol = "red", 
                                     finalcol = "dodgerblue", legend = TRUE, xlab, ylab, log = "xy",
                                     cex = 0.45, ...) {
            require("ggplot2")
            if (missing(xlab)) 
              xlab <- "mean of normalized counts"
            if (missing(ylab)) {
              if (CV) {
                ylab <- "coefficient of variation"
              }
              else {
                ylab <- "dispersion"
              }
            }
            
            px = mcols(object)$baseMean
            sel = (px > 0)
            px = px[sel]
            f <- if (CV){ 
              sqrt} else {I
              }
            
            py = f(mcols(object)$dispGeneEst[sel])
            
            if (missing(ymin)){ 
              ymin = 10^floor(log10(min(py[py > 0], na.rm = TRUE)) - 0.1)
            }
            
            disp_df <- data.frame(px=px, py=py, 
                                  pch_outlier=as.factor(ifelse(mcols(object)$dispOutlier[sel], 1, 16)),
                                  cex_outlier=as.factor(ifelse(mcols(object)$dispOutlier[sel], 2 * cex, cex)))
            ggplot(disp_df, aes(x=px, y=pmax(py, ymin))) + geom_point(aes(color=pch_outlier, shape=cex_outlier)) + labs(x=xlab, y=ylab) + 
              scale_y_continuous(trans = "log10") + scale_x_continuous(trans = "log10") + theme_minimal()
            
            plotDispEsts(object)  
            plot(px, pmax(py, ymin), xlab = xlab, ylab = ylab, log = log, 
                 pch = ifelse(py < ymin, 6, 20), col = genecol, cex = cex)
            pchOutlier <- ifelse(mcols(object)$dispOutlier[sel], 1, 16)
            cexOutlier <- ifelse(mcols(object)$dispOutlier[sel], 2 * 
                                   cex, cex)
            lwdOutlier <- ifelse(mcols(object)$dispOutlier[sel], 2, 1)
            if (!is.null(dispersions(object))) {
              points(px, f(dispersions(object)[sel]), col = finalcol, 
                     cex = cexOutlier, pch = pchOutlier, lwd = lwdOutlier)
            }
            if (!is.null(mcols(object)$dispFit)) {
              points(px, f(mcols(object)$dispFit[sel]), col = fitcol, 
                     cex = cex, pch = 16)
            }
            if (legend) {
              legend("bottomright", c("gene-est", "fitted", "final"), 
                     pch = 16, col = c(genecol, fitcol, finalcol), bg = "white")
            }
          }
        }
        
        additional_check_data <- list(cooksDist_comparison=cooksDist_comparison, 
                                      pval_hist_comparison=pval_hist_comparison)
        
        return(additional_check_data)
        
      } else {
        stop("Need to provide dds object!")
      }
    }
    additional_check_data <- additionalCheck(dds_object=dds, dds_filt_object=dds_filt)

    # [ ] add: check also  DESeq2:::plotSparsity for validation of negative binomial distribution
    #https://support.bioconductor.org/p/98514/
    ggsave(filename = paste0(FIGURE_DIR, "CooksDistance.png"), plot=additional_check_data$cooksDist_comparison, width = 30, height = 30, units = "cm")
    ggsave(filename = paste0(FIGURE_DIR, "PvalHistogram.png"), plot=additional_check_data$pval_hist_comparison, width = 30, height = 30, units = "cm") 
    save(additional_check_data, file = paste0(OUTPUT_DIR, experiment_name, "_AdditionalChecks.RData"))
    
    cat("Done\n")
  }

  #save(pca_scree_cumVarExpl, PCA_matrix_plot, file = paste0(OUTPUT_DATASET_DIR, experiment_name, "_PCAplots.RData"))
  PCA_analysis <- TRUE
  if (PCA_analysis) {
    # PCA ----
    # [x] add scree plot + cummulative variance explained
    # [] add loadings and rotations (PCs)
    # [] plot 1st N PCs
    ntop_genes=nrow(rld_counts_filt) # default 500
    #rld_counts_ntop <- rnaSelectTopVarGenes_var(rld_counts, ntop = ntop_genes)
    rld_counts_filt_ntop <- rnaSelectTopVarGenes(rld_counts_filt, ntop = ntop_genes)
    
    colnames(rld_counts_filt_ntop) <- paste0(rld_filt[[cond_interest]], 
                                             "_",
                                             colnames(rld_counts_filt_ntop))
    
    res_pca <- stats::prcomp(t(rld_counts_filt_ntop), center = TRUE)
    pca_results <- pcaExtractVariance(res_pca)
    pca_scree_cumVarExpl <- pcaPlotVariance(pca_results, var_expl_needed = var_expl_needed)
    # plot(pca_scree_cumVarExpl)
    
    min_components_needed <- which(pca_results$cummulative_variance > var_expl_needed)[1]
    
    # if less than 5 components needed plot otherwise skip; maybe update up to 5 components?
    if (min_components_needed < 5) {
      # [x] add only 1 legend for the overall plot to save space
      
      plot_PCAs_list <- pcaPlotPCs(res_pca, PCs_range=1:min_components_needed, 
                                  groups = rld_filt[[cond_interest]], 
                                  legend_title = cond_interest, common_legend = TRUE,
                                  add_ellipse = FALSE)
      
      if ("common_legend" %in% names(plot_PCAs_list)) {
        #extract legend from PC plots and plot commong legend
        m_size <- (min_components_needed - 1) # (min_components_needed - 1); 3x3 for 4 PCs
        m <- matrix(NA, m_size, m_size)
        # n!/(k!(n-k)!)
        PC_combinations <- factorial(min_components_needed)/(factorial(2)*factorial(min_components_needed-2))
        m[lower.tri(m, diag = T)] <- 1:PC_combinations # 6 for 4 PCs
        # make it more automatic similar to combinations calculation above so it alco places correctly column, no just row!
        if (ncol(m) == 2) {
          m[ceiling(m_size/2), 2] <- (PC_combinations+1) # 7 for 4 PCs legend position; 
        } else {
          m[ceiling(m_size/2), 3] <- (PC_combinations+1) # 7 for 4 PCs legend position; 
        }
        PCA_matrix_plot <- gridExtra::grid.arrange(grobs = plot_PCAs_list, layout_matrix = m)
        
        # another option is to directly use, however it does not produce nice lower triangular matrix shape
        #ggpubr::ggarrange(plotlist = plot_PCAs_list, common.legend = TRUE, legend="right")
        
      } else {
        m_size <- 3
        m <- matrix(NA, m_size, m_size)
        m[lower.tri(m, diag = T)] <- 1:(min_components_needed+2)
        PCA_matrix_plot <- gridExtra::grid.arrange(grobs = plot_PCAs_list, layout_matrix = m)
      }
    } else {
      print(paste0("Only plot for PC2.vs.PC1 generated since more than 4 PCs needed to explain more than: ", var_expl_needed))
      plot_PCAs_list <- pcaPlotPCs(res_pca, PCs_range=1:2, groups = rld_filt[[cond_interest]], legend_title = "Condition")
      PCA_matrix_plot <- plot_PCAs_list[[1]]
      #plot(plot_PCAs_list[[1]])
    }
    
    # Correlation of PC to covariates: compute and plot significance of the correlation of each covariate versus a principal component.
    # based on pcaExplorer
    sample_table_forCorr <- cond_data
    
    PC_corr <- pcaCorrPCs(res_pca, sample_table_forCorr[cond_interest_varPart], pcs = 1:min_components_needed, pval_exact = TRUE) 
    PC_corr_plot <- pcaCorrPCsPlot(PC_corr)
    
    # function - batch removal testing on PCA
    
    PCA_investigation=TRUE
    
    if(PCA_investigation) {
      # Initial plotting of PCA and highlighting different variables
      # assessing PCAs based on PC_corr_plot
      # also check tables
      #ntop_genes=nrow(vsd_filt) # default 500
      table(coldata$library_type)
      table(coldata$gender)
      table(coldata$timepoint)
      table(coldata$material)
      table(coldata$material, coldata$timepoint)
      table(coldata$material, coldata$disease)
      table(coldata$material, coldata$library_type)
      table(coldata$rnaseq_batch)
      table(coldata$library_type)
      table(coldata$library_type, coldata$rnaseq_batch) # cofounded
      chisq.test(coldata$library_type, coldata$rnaseq_batch)
      chisq.test(coldata$library_type, coldata$material)
      chisq.test(coldata$rnaseq_batch, coldata$gender)
      table(coldata$material, coldata$gender)
      
      vst_transf = rld_filt
      
      # check correlation between depth and library_type
      load(paste0(OUTPUT_DIR, experiment_name, "_SeqDepth.RData"))
      raw_counts_libSize_filt <- library_size_data$raw_counts_libSize %>%
        dplyr::filter(!is.na(library_type)) %>%
        dplyr::mutate(library_type = as.factor(library_type))
      # wilcox.test(library_type ~ library_size, data = raw_counts_libSize_filt)
      library_type_depth_effect <- ggplot(raw_counts_libSize_filt, aes(x=library_type, y=library_size, fill=library_type)) +
        geom_boxplot() + theme_bw()
      
      rnaseqBatch_depth_effect <- ggplot(raw_counts_libSize_filt, aes(x=rnaseq_batch , y=library_size, fill=rnaseq_batch)) +
        geom_boxplot() + theme_bw()
      # table(raw_counts_libSize_filt$rnaseq_batch, raw_counts_libSize_filt$library_type)
      
      # pca_all_cond_libType <- plot_pca(transf_object = vst_transf,
      #                                     intgroup=c(cond_interest, "library_type"))
      # pca_all_material_gender <- plot_pca(transf_object = vst_transf,
      #                                     intgroup=c("material", "gender"))
      
      pca_batch_variables <- cond_interest_varPart[!(cond_interest_varPart %in% cond_interest)]
      names(pca_batch_variables) <- pca_batch_variables
      pca_plots_all <- purrr::map(.x = pca_batch_variables, .f = function(batch_variable){
        message("variable:", batch_variable)
        plot_pca(transf_object = vst_transf,
                 intgroup=c(batch_variable, cond_interest, "patient_oid"))
      })
      
      # patchwork::wrap_plots(pca_plots_all, ncol = 2)
      
      # pca_all_cond_libType_label <- plot_pca(transf_object = vst_transf,
      #                                        intgroup=c(cond_interest, "library_type", "patient_oid"))
      # pca_all_cond_rnaseqBatch_label <- plot_pca(transf_object = vst_transf,
      #                                            intgroup=c("rnaseq_batch", cond_interest, "patient_oid"))
      # pca_all_cond_material_label <- plot_pca(transf_object = vst_transf,
      #                                         intgroup=c("material", cond_interest, "patient_oid"))
      # pca_all_cond_gender_label <- plot_pca(transf_object = vst_transf,
      #                                       intgroup=c("gender", cond_interest, "patient_oid"))
      # pca_all_cond_disease_label <- plot_pca(transf_object = vst_transf,
      #                                        intgroup=c("disease", cond_interest, "patient_oid"))
      
      
      # checking QC on subsets ---- 
      check_QC <- function(dds_object=NULL, cond_subset=NULL, cond_interest=NULL, cond_interest_varPart=NULL, fitform_subset=NULL){
        # using default design on dds_object: ~1 design
        # using same number of PCAs as used in dds_object analysis
        # [] remove dependency on other loaded functions
        
        dds_subset <- dds_object[, colnames(dds_object) %in% cond_subset$omics_id]
        dds_subset <- DESeq(dds_subset)
        vst_subset <- vst(dds_subset, blind = TRUE) # blind = TRUE for QC
        ntop_genes <- nrow(vst_subset)
        vst_subset_ntop <- rnaSelectTopVarGenes(assay(vst_subset), ntop = ntop_genes)
        
        # pca
        pca_subset <- stats::prcomp(t(vst_subset_ntop), center = TRUE)
        PC_corr_df <- pcaCorrPCs(pca_subset, cond_subset[cond_interest_varPart], pcs = 1:6, pval_exact = TRUE) 
        PC_corr_plot <- pcaCorrPCsPlot(PC_corr_df)
        
        pca_batch_variables <- cond_interest_varPart[!(cond_interest_varPart %in% cond_interest)]
        names(pca_batch_variables) <- pca_batch_variables
        pca_plots_all <- purrr::map(.x = pca_batch_variables, .f = function(batch_variable){
          message("variable:", batch_variable)
          plot_pca(transf_object = vst_subset,
                   intgroup=c(batch_variable, cond_interest))
        })
        
        varPart_subset <- partition_variance(vst_transf = vst_subset, fitform_partVariance = fitform_subset, ncores = 40, max_ntop_genes = 5000)
        
        check_QC_results <- list(PC_corr_plot = PC_corr_plot,
                                 pca_plots_all = pca_plots_all,
                                 varPart_subset = varPart_subset)
        
        return(check_QC_results)
      }
      
      # library size QC
      #boxplot(library_size_data$raw_counts_libSize$library_size_M)
      
      library_size_data_subset <- library_size_data$raw_counts_libSize %>%
        dplyr::select(omics_id = sample_name, library_size_M, sizeFactor)
      
      cond_data$
      
      # rnaseq_batch - ABC
      cond_data_ABC <- cond_data %>% filter(rnaseq_batch %in% c("A", "B", "C"))
      cond_interest_varPart_ABC <- c("rnaseq_batch", "material", "gender", "timepoint")
      fitform_partVariance_ABC <- ~(1 | rnaseq_batch) + (1 | material) + (1 | gender) + (1 | timepoint) + (1 | patient_oid)
      rnaseq_batchABC_QC <- check_QC(dds_object = dds_filt, 
                                     cond_subset = cond_data_ABC, 
                                     cond_interest = cond_interest, 
                                     cond_interest_varPart = cond_interest_varPart_ABC,
                                     fitform_subset = fitform_partVariance_ABC)
      
      # rnaseq_batch - GE
      cond_data_GE <- cond_data %>% filter(rnaseq_batch %in% c("G", "E"))
      cond_interest_varPart_GE <- c("rnaseq_batch", "material", "gender", "timepoint")
      fitform_partVariance_GE <- ~(1 | rnaseq_batch) + (1 | material) + (1 | gender) + (1 | timepoint) + (1 | patient_oid)
      rnaseq_batchGE_QC <- check_QC(dds_object = dds_filt, 
                                     cond_subset = cond_data_GE, 
                                     cond_interest = cond_interest, 
                                     cond_interest_varPart = cond_interest_varPart_GE,
                                     fitform_subset = fitform_partVariance_GE)
      

      #fitform_partVariance_filt <- ~(1 | gender) + (1 | rnaseq_batch) + (1 | library_type) + (1 | material)
      vst_transf_varPart <- partition_variance(vst_transf = vst_transf, fitform_partVariance = fitform_partVariance, ncores = 40, max_ntop_genes = 5000)
      
      # removing key batch effect(s) 
      # [] fix error when removing batch effect!
      vst_transf_noBatch <- remove_batch(vst_transf = vst_transf, 
                                         batch_variable = "library_type", 
                                         cond_interest = cond_interest) # material, timepoint? or just intercept?
      
      pca_plots_all_noBatch <- purrr::map(.x = pca_batch_variables, .f = function(batch_variable){
        message("variable:", batch_variable)
        plot_pca(transf_object = vst_transf_noBatch,
                 intgroup=c(batch_variable, cond_interest, "patient_oid"))
      })
      
      vst_transf_noBatch_ntop <- rnaSelectTopVarGenes(assay(vst_transf_noBatch), ntop = ntop_genes)

      res_pca_noBatch <- stats::prcomp(t(vst_transf_noBatch_ntop), center = TRUE)
      PC_corr_noBatch <- pcaCorrPCs(res_pca_noBatch, sample_table_forCorr[cond_interest_varPart], pcs = 1:min_components_needed, pval_exact = TRUE) 
      PC_corr_noBatch_plot <- pcaCorrPCsPlot(PC_corr_noBatch)
      
      vst_transf_noBatch_varPart <- partition_variance(vst_transf = vst_transf_noBatch, fitform_partVariance = fitform_partVariance, ncores = 40, max_ntop_genes = 5000)
      
      # extracting colors from varPart_plot
      varPart_colors <- unique(ggplot_build(vst_transf_varPart$varPart_plot_annot)$data[[1]]["fill"])[["fill"]]
      names(varPart_colors) <- c(cond_interest_varPart, "Residuals")
      varPart_plot_annot <- vst_transf_varPart$varPart_plot_annot + 
        labs(title = "original") + 
        scale_fill_manual(name = "variable", values = varPart_colors) # + theme(axis.text.x=element_blank())
      varPart_plot_annot_noBatch <- vst_transf_noBatch_varPart$varPart_plot_annot + 
        labs(title = "batch effect removed") + 
        scale_fill_manual(name = "variable", values = varPart_colors) #+ theme(axis.text.x=element_blank())
      
      varPart_noBatch_compare <- ggpubr::ggarrange(varPart_plot_annot, 
                                                   varPart_plot_annot_noBatch,
                                                   legend = "bottom",
                                                   common.legend = TRUE) 

      # include library_type, gender and rnasaeq_batch in analysis???
      # table(coldata$rnaseq_batch, coldata$library_type)
      # table(coldata$rnaseq_batch, coldata$library_type)
      # control for library_type and gender!
      # library_type effect is probably coming from library_depth differences between U, SR rather than strandness?

    # add tSNE
    
    # save PCA and varPart plots
    save(vst_transf_varPart, vst_transf_noBatch_varPart, file = paste0(OUTPUT_DIR, experiment_name, "_varPart.RData"))
    
    # Plotting figures into file 
    ggsave(filename = paste0(FIGURE_DIR, "pca_scree_cumVarExpl.png"), plot=pca_scree_cumVarExpl,
           width = 20,
           height = 20,
           units = "cm")

    ggsave(filename = paste0(FIGURE_DIR, "PC_corr_plot.png"), plot=PC_corr_plot,
           width = 20,
           height = 20,
           units = "cm")

    ggsave(filename = paste0(FIGURE_DIR, "library_type_depth_effect.png"), plot=library_type_depth_effect,
           width = 20,
           height = 20,
           units = "cm")
    
    ggsave(filename = paste0(FIGURE_DIR, "rnaseqBatch_depth_effect.png"), plot=rnaseqBatch_depth_effect,
           width = 20,
           height = 20,
           units = "cm")
    
    # saving batch related plots
    ggsave(filename = paste0(FIGURE_DIR, "PC_corr_noBatch_plot.png"), plot=PC_corr_noBatch_plot,
           width = 20,
           height = 20,
           units = "cm")
    
    ggsave(filename = paste0(FIGURE_DIR, "varPart_plot_annot.png"), plot=varPart_plot_annot,
           width = 35,
           height = 25,
           units = "cm")
    
    ggsave(filename = paste0(FIGURE_DIR, "varPart_noBatch_compare.png"), plot=varPart_noBatch_compare,
           width = 35,
           height = 25,
           units = "cm")
    
    # saving PCA plots
    purrr::map2(.x = pca_plots_all, .y = names(pca_plots_all), .f = function(pca_batch, batch_variable){
      message("saving PCA plot:", batch_variable)
      ggsave(pca_batch,
             filename = paste0(FIGURE_DIR, "pca_cond_",batch_variable,".png"),
             width = 25,
             height = 25,
             units = "cm")
    })
    
    # saving noBatch PCAs
    purrr::map2(.x = pca_plots_all_noBatch, .y = names(pca_plots_all_noBatch), .f = function(pca_batch, batch_variable){
      message("saving PCA plot:", batch_variable)
      ggsave(pca_batch,
             filename = paste0(FIGURE_DIR, "pca_cond_",batch_variable,"_noBatch.png"),
             width = 25,
             height = 25,
             units = "cm")
    })
    
    # pdfs
    # saving PCA plots
    purrr::map2(.x = pca_plots_all, .y = names(pca_plots_all), .f = function(pca_batch, batch_variable){
      message("saving PCA plot:", batch_variable)
      ggsave(pca_batch,
             filename = paste0(FIGURE_DIR, "pca_cond_",batch_variable,".pdf"),
             width = 25,
             height = 25,
             units = "cm")
    })
    
    # saving noBatch PCAs
    purrr::map2(.x = pca_plots_all_noBatch, .y = names(pca_plots_all_noBatch), .f = function(pca_batch, batch_variable){
      message("saving PCA plot:", batch_variable)
      ggsave(pca_batch,
             filename = paste0(FIGURE_DIR, "pca_cond_",batch_variable,"_noBatch.pdf"),
             width = 25,
             height = 25,
             units = "cm")
    })

  }
  
######### END Data quality assessment by sample clustering and visualization
  
  }
}

celltype_estimation=TRUE  
if (celltype_estimation) {
  tumour_purity_estimate=TRUE
  if (tumour_purity_estimate) {
    # Estimate, ConsensusTME, CIBERSORT
    # https://bioinformatics.mdanderson.org/estimate/rpackage.html
    # https://github.com/cansysbio/ConsensusTME
    # https://github.com/icbi-lab/immunedeconv
    # [ ] also see https://shenorrlab.github.io/bseqsc/index.html
    # [ ] or directly using scRNA-seq data in CIBERSORT!?
    # other resources:
    # https://www.biostars.org/p/428905/
    # https://www.nature.com/articles/s41467-020-19015-1
    # https://academic.oup.com/bib/article/22/1/416/5699815
    # https://www.biostars.org/p/439224/
    # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02126-9
    # https://singulomics.com/bulk-rna-seq-with-cell-type-deconvolution/
    # https://www.biorxiv.org/content/10.1101/2020.01.10.897116v1.full
    # https://advances.sciencemag.org/content/6/30/eaba2619
    # https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007510
    # 
    # cibersort - disable quantile normalization for RNA-seq!; requires log2 scale or just normalized data?; ideally TPM (salmon etc.)
    
    # NOTE: single-cell RNA-seq use for deconvolution of bulk-RNA-seq; train unique profiles and improve with new data!
    # use in diagnosis with cheap quant-seq?!
    # multi-omics: assessing heterogeneity at different levels and picking correct level to discriminate!? and also to assess response to treatment!?
    # for single-cell also use single-cell atlas to align to reference?!
    
    #install.packages("estimate", repos="http://r-forge.r-project.org", dependencies=TRUE)
    #remotes::install_github("cansysbio/ConsensusTME") # update none
    #remotes::install_github("icbi-lab/immunedeconv")
    # library(ggplot2)
    # library(ggpubr)
    # library(ggsci)
    
    mcpCounter_results_plot <- function(mcpCounter_results = NULL, id = "omics_id", metadata_df = cond_data, cond_interest = cond_interest, test_groups = list(c("no", "yes")), plot_title = "Add title"){
      # [ ] violin plots
      # [ ] add heatmaps
      # extract_mcpCounter_results: extract results and plot in heatmap
      require(dplyr)
      require(ggpubr)
      
      mcpCounter_results_df <- mcpCounter_results %>%
        tidyr::pivot_longer(names_to = id, values_to = "score", -cell_type) %>%
        dplyr::left_join(., dplyr::select(metadata_df, c("omics_id", tidyselect::all_of(cond_interest))), by = "omics_id") #%>%  dplyr::select(-omics_id)
      
      mcpCounter_results_plot <- mcpCounter_results_df %>%
        ggpubr::ggviolin(., x = cond_interest, y = "score",
                         fill = cond_interest,
                         palette = "jco",
                         facet.by = "cell_type",
                         #scales = "fixed",
                         xlab = cond_interest,
                         ylab = "score",
                         add = c("boxplot", "jitter"), 
                         add.params = list(fill = "white")) +  # to remove outlier point if plotting with jitter 
        ggplot2::facet_wrap(~cell_type, scales = "free") +  # artibtrary unit using free scale
        #ggpubr::facet(TME_PAAD_celltype_enrichments_plot, facet.by = "cell_type", scales = "free") +
        ggpubr::stat_compare_means(comparisons = test_groups,  
                                   method="wilcox.test", paired=FALSE,  
                                   method.args = list(alternative = "two.sided")) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = plot_title)
      
      mcpCounter_results = list(mcpCounter_results_df = mcpCounter_results_df,
                                mcpCounter_results_plot = mcpCounter_results_plot)
      return(mcpCounter_results)
    }
    
    ssgsea_results_plot <- function(ssgsea_results = NULL, id = "omics_id", metadata_df = cond_data, cond_interest = cond_interest, test_groups = list(c("no", "yes")), plot_title = "Add title"){
      # [ ] violin plots
      # [ ] add heatmaps
      # extract_mcpCounter_results: extract results and plot in heatmap
      require(dplyr)
      require(ggpubr)
      require(tibble)

      ssgsea_results_df <- as.data.frame(ssgsea_results) %>%
        tibble::rownames_to_column(., var = "cell_type") %>%
        tidyr::pivot_longer(names_to = id, values_to = "enrichment_score", -cell_type) %>%
        dplyr::left_join(., dplyr::select(metadata_df, c("omics_id", tidyselect::all_of(cond_interest))), by = "omics_id") #%>%  dplyr::select(-omics_id)
      
      ssgsea_results_plot <- ssgsea_results_df %>%
        ggpubr::ggviolin(., x = cond_interest, y = "enrichment_score",
                         fill = cond_interest,
                         palette = "jco",
                         facet.by = "cell_type",
                         #scales = "fixed",
                         xlab = cond_interest,
                         ylab = "Enrichment score",
                         add = c("boxplot", "jitter"), 
                         add.params = list(fill = "white")) +  # to remove outlier point if plotting with jitter 
        ggplot2::facet_wrap(~cell_type, scales = "free") +  # artibtrary unit using free scale
        #ggpubr::facet(TME_PAAD_celltype_enrichments_plot, facet.by = "cell_type", scales = "free") +
        ggpubr::stat_compare_means(comparisons = test_groups,  
                                   method="wilcox.test", paired=FALSE,  
                                   method.args = list(alternative = "two.sided")) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = plot_title)
      
      ssgsea_results = list(ssgsea_results_df = ssgsea_results_df,
                                ssgsea_results_plot = ssgsea_results_plot)
      return(ssgsea_results)
    }
    
    # analyzing immune cells
    # analyzing custom genesets!
    
    # mdsc signatures ----
    mdsc_sign_DIR = "/home/peter_r/reseq_database/datasets/grant_application/MDSC_signatures/"
    G_MDSC_raw <- readxl::read_excel(file.path(mdsc_sign_DIR, "aay6017_Table_S3.xlsx"))
    G_MDSC <- G_MDSC_raw %>%
      dplyr::select(gene, avg_logFC, p_val_adj) %>%
      dplyr::filter(p_val_adj < 0.05) #& abs(avg_logFC) > 2
    # above gives 2956 genes, but should be 642
    
    M_MDSC_raw <- readxl::read_excel(file.path(mdsc_sign_DIR, "aay6017_Table_S4.xlsx"))
    M_MDSC <- M_MDSC_raw %>%
      dplyr::select(gene, avg_logFC, p_val_adj) %>%
      dplyr::filter(p_val_adj < 0.05) #& abs(avg_logFC) > 2
    # above gives 829 genes, but should be 223
    
    combined_MDSC_raw <- readxl::read_excel(file.path(mdsc_sign_DIR, "aay6017_Table_S5.xlsx"))
    combined_MDSC <- combined_MDSC_raw %>%
      dplyr::rename(mouse_gene_symbol=x) %>%
      dplyr::pull(mouse_gene_symbol)
    
    # convert mouse to human
    #ensembl_version="Ensembl Genes 92", 
    #biomart_host="http://apr2018.archive.ensembl.org"
    human_mart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl", 
                              host="http://nov2020.archive.ensembl.org", 
                              version="Ensembl Genes 102")
    mouse_mart <- biomaRt::useMart("ensembl", dataset="mmusculus_gene_ensembl", 
                              host="http://nov2020.archive.ensembl.org", 
                              version="Ensembl Genes 102")
    combined_MDSC_human_all <- biomaRt::getLDS(attributes = c("ensembl_gene_id"), filters = "mgi_symbol", values = combined_MDSC, mart = mouse_mart, 
                           attributesL = c("ensembl_gene_id","hgnc_symbol", "description", "chromosome_name", "start_position", "end_position", "strand"), 
                           martL = human_mart)
    # removing duplicates
    combined_MDSC_human_signature <- combined_MDSC_human_all %>%
      dplyr::rename(gene_symbol = HGNC.symbol) %>%
      dplyr::distinct(gene_symbol) %>%
      dplyr::pull(gene_symbol) 
    combined_MDSC_human_signature = list(combined_MDSC_human_signature = combined_MDSC_human_signature)
    
    # creating input expression matrices ----
    # human  
    vst_transf_counts <- assay(vst_transf)
    vst_transf_noBatch_counts <- assay(vst_transf_noBatch)

    vst_matrix_GeneSymbols <- convertIDsMatrix(mat_obj=vst_transf_counts, annot_df=ensemblAnnot,
                                                       id_from = "ensembl_id", id_to = "hgnc_symbol",
                                                       resolve_duplicates="var") 
    vst_matrix_noBatch_GeneSymbols <- convertIDsMatrix(mat_obj=vst_transf_noBatch_counts, annot_df=ensemblAnnot,
                                                      id_from = "ensembl_id", id_to = "hgnc_symbol",
                                                      resolve_duplicates="mad") 
    
    # vst_transf_GeneSymbols_counts_noDupl <- as.data.frame(vst_transf_counts) %>%
    #   tibble::rownames_to_column(., var = "ensembl_id") %>%
    #   dplyr::left_join(., ensemblAnnot_ensembl_symbol, by = "ensembl_id") %>%
    #   dplyr::select(-ensembl_id) %>%
    #   dplyr::distinct(hgnc_symbol, .keep_all = TRUE) %>%
    #   tibble::column_to_rownames(., var = "hgnc_symbol") %>%
    #   as.matrix(.)
    
    # identical(counts(dds_filt,normalized=TRUE), counts(dds_filt_extractResults,normalized=TRUE)) # TRUE
    # for most of the methods in immunedeconv; check in depth which methods expect what and read original immunedeconv paper (and methods papers)
    # However, log2 for MCP-counter based on this https://github.com/ebecht/MCPcounter/issues/4
    # should be log2 transformed! -> but when compared norm and vst it gives different results
    # norm_counts <- counts(dds_filt_extractResults,normalized=TRUE)
    # norm_matrix_GeneSymbols <- convertIDsMatrix(mat_obj=norm_counts, annot_df=ensemblAnnot,
    #                                            id_from = "ensembl_id", id_to = "hgnc_symbol",
    #                                            resolve_duplicates="var") 
    
    # TME calculations and plot ----
    # https://icbi-lab.github.io/immunedeconv/articles/immunedeconv.html
    
    # try MCP-counter with
    # based on this https://github.com/ebecht/MCPcounter/issues/4
    # should be log2 transformed! -> but when compared norm and vst it gives different results
    # [ ] vst_transf
    # [] normalized values
    # with and without batch effect (for norm. need to substract)
    #cor(norm_matrix_GeneSymbols, vst_matrix_GeneSymbols, method = "pearson") # is this due to duplicates???
    # quantiseq
    # timer
    # cibersort
    # cibersort_abs
    # mcp_counter
    # xcell
    # epic
    
    # should use non-log transformed values; normalized?!
    # use Salmon results to do celltype estimation (import TPM!)
    # TIMER - requires tumor time to be specified
    # cibersort_abs - need to switch off normalization - see consensusTME methods and FAQ on cibersort
    # Xcell - https://xcell.ucsf.edu/
    # xCell uses TPM; uses the expression levels ranking and not the actual values, thus normalization does not have an effect, however normalizing to gene length (RPKM/FPKM/TPM/RSEM) is required.
    # quantiseq uses non-log; https://icbi.i-med.ac.at/software/quantiseq/doc/; non-log!
    # EPIC uses TPM, non-log; https://github.com/GfellerLab/EPIC
    
    # quantiseq
    # ! Expression data must be on non-log scale
    # vst_quantiseq <- immunedeconv::deconvolute(vst_matrix_GeneSymbols, method = "quantiseq", tumor = FALSE)
    # vst_quantiseq_noBatch <- immunedeconv::deconvolute(vst_matrix_noBatch_GeneSymbols, method = "quantiseq", tumor = FALSE)
    # 
    # vst_quantiseq_plot <- mcpCounter_results_plot(mcpCounter_results = vst_quantiseq,
    #                                                id = "omics_id",
    #                                                metadata_df = cond_data,
    #                                                cond_interest = cond_interest,
    #                                                test_groups = list(c("no", "yes")),
    #                                                plot_title = "MCP-Counter")
    # 
    # vst_quantiseq_noBatch_plot <- mcpCounter_results_plot(mcpCounter_results = vst_quantiseq_noBatch,
    #                                                        id = "omics_id",
    #                                                        metadata_df = cond_data,
    #                                                        cond_interest = cond_interest,
    #                                                        test_groups = list(c("no", "yes")),
    #                                                        plot_title = "MCP-Counter no batch")
    # 
    # MCP-counter ----
    # it may actually not be vst - misleading label!!!!
    vst_mcpCounter <- immunedeconv::deconvolute(vst_matrix_GeneSymbols, method = "mcp_counter")
    vst_mcpCounter_noBatch <- immunedeconv::deconvolute(vst_matrix_noBatch_GeneSymbols, method = "mcp_counter")
    
    vst_mcpCounter_plot <- mcpCounter_results_plot(mcpCounter_results = vst_mcpCounter,
                                                   id = "omics_id",
                                                   metadata_df = cond_data,
                                                   cond_interest = cond_interest,
                                                   test_groups = list(c("no", "yes")),
                                                   plot_title = "MCP-Counter")
    
    vst_mcpCounter_noBatch_plot <- mcpCounter_results_plot(mcpCounter_results = vst_mcpCounter_noBatch,
                                                   id = "omics_id",
                                                   metadata_df = cond_data,
                                                   cond_interest = cond_interest,
                                                   test_groups = list(c("no", "yes")),
                                                   plot_title = "MCP-Counter no batch")
    
    # runnig ssgssea on mesenchymal vs adrenergic ----
    # signatures_list
    ssgsea_mes_adr_ncc_noradr <- GSVA::gsva(vst_matrix_GeneSymbols,
                                            signatures_list,
                                            method=c("ssgsea"),
                                            min.sz=1, max.sz=Inf, 
                                            ssgsea.norm=TRUE, verbose=TRUE, parallel.sz=10)
    
    ssgsea_mes_adr_ncc_noradr_gsva <- GSVA::gsva(vst_matrix_GeneSymbols,
                                            signatures_list,
                                            method=c("gsva"),
                                            min.sz=1, max.sz=Inf, 
                                            verbose=TRUE, parallel.sz=10)
    
    # add stats to above and also kruskal for global difference
    # scaling? A <- t(scale(t(ssgsea_mes_adr_ncc_noradr)))
    ssgsea_mes_adr_ncc_noradr_plot <- ssgsea_results_plot(ssgsea_results = ssgsea_mes_adr_ncc_noradr,
                                                               id = "omics_id",
                                                               metadata_df = cond_data,
                                                               cond_interest = cond_interest,
                                                               #test_groups = list(c("no", "yes")),
                                                               plot_title = "mes_adr_ncc_noradr")
    
    ggsave(filename = paste0(FIGURE_DIR, "ssgsea_mes_adr_ncc_noradr_plot.png"), plot=ssgsea_mes_adr_ncc_noradr_plot$ssgsea_results_plot,
           width = 25,
           height = 25,
           units = "cm")
    
    
    heatmap_col_annot_mat_tp <- cond_data %>%
      dplyr::mutate(mat_tp = paste(material, timepoint, sep="_"))
    
    ssgsea_mes_adr_ncc_noradr_mat_tp_plot <- ssgsea_results_plot(ssgsea_results = ssgsea_mes_adr_ncc_noradr,
                                                          id = "omics_id",
                                                          metadata_df = heatmap_col_annot_mat_tp,
                                                          cond_interest = "mat_tp",
                                                          #test_groups = list(c("no", "yes")),
                                                          plot_title = "mes_adr_ncc_noradr")
    
    ggsave(filename = paste0(FIGURE_DIR, "ssgsea_mes_adr_ncc_noradr_mat_tp_plot.png"), plot=ssgsea_mes_adr_ncc_noradr_mat_tp_plot$ssgsea_results_plot,
           width = 25,
           height = 25,
           units = "cm")
    
    # test also gsva and scale or not scale scores?!!
    # also check clustering samples or not!
    # more samples from same patient
    heatmap_col_annot_all <- cond_data %>%
      dplyr::mutate(mat_tp = paste(material, timepoint, sep="_")) %>%
      dplyr::mutate(mat_tp = factor(mat_tp, levels = c("DTC_DX", "DTC_REL", "MNC_DX", "MNC_REL", "TUM_DX", "BMn_DX"))) %>%
      dplyr::select(omics_id, library_type, rnaseq_batch, gender, mat_tp) %>%
      tibble::column_to_rownames(var = "omics_id") %>%
      dplyr::arrange(mat_tp)
    
    
    # order by mat_tp
    identical(colnames(ssgsea_mes_adr_ncc_noradr), rownames(heatmap_col_annot_all))
    ssgsea_mes_adr_ncc_noradr_order <- dplyr::select(as.data.frame(ssgsea_mes_adr_ncc_noradr), rownames(heatmap_col_annot_all)) %>% as.matrix(.)
    identical(colnames(ssgsea_mes_adr_ncc_noradr_order), rownames(heatmap_col_annot_all))
    
    # , filename = paste0(FIGURE_DIR, "ssgsea_mes_adr_ncc_noradr_heatmap_all.png")
    ssgsea_mes_adr_ncc_noradr_heatmap_clust_all <- pheatmap::pheatmap(ssgsea_mes_adr_ncc_noradr_order,
                                                            scale = "row",
                                                            annotation_col = heatmap_col_annot_all,
                                                            cluster_rows = TRUE,
                                                            cluster_cols = TRUE,
                                                            color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                            show_colnames = FALSE, 
                                                            filename = paste0(FIGURE_DIR, "ssgsea_mes_adr_ncc_noradr_heatmap_clust_all.png"))
    
    heatmap_col_annot <- heatmap_col_annot_all %>%
      dplyr::select(mat_tp) 
    ssgsea_mes_adr_ncc_noradr_clust_heatmap_clust <- pheatmap::pheatmap(ssgsea_mes_adr_ncc_noradr_order,
                                                                scale = "row",
                                                                annotation_col = heatmap_col_annot,
                                                                cluster_rows = TRUE,
                                                                cluster_cols = TRUE,
                                                                color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                                show_colnames = FALSE,
                                                            filename = paste0(FIGURE_DIR, "ssgsea_mes_adr_ncc_noradr_clust_heatmap_clust.png"))
    
    # not clustering 
    ssgsea_mes_adr_ncc_noradr_heatmap_all <- pheatmap::pheatmap(ssgsea_mes_adr_ncc_noradr_order,
                                                                      scale = "row",
                                                                      annotation_col = heatmap_col_annot_all,
                                                                      cluster_rows = TRUE,
                                                                      cluster_cols = FALSE,
                                                                      color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                                show_colnames = FALSE,
                                                                filename = paste0(FIGURE_DIR, "ssgsea_mes_adr_ncc_noradr_heatmap_all.png"))
    
    ssgsea_mes_adr_ncc_noradr_heatmap <- pheatmap::pheatmap(ssgsea_mes_adr_ncc_noradr_order,
                                                                  scale = "row",
                                                                  annotation_col = heatmap_col_annot,
                                                                  cluster_rows = TRUE,
                                                                  cluster_cols = FALSE,
                                                                  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                                  show_colnames = FALSE,
                                                                  filename = paste0(FIGURE_DIR, "ssgsea_mes_adr_ncc_noradr_heatmap.png"))
    
    
    # in the paper from Daria's summary do they scale gsva scores or not?
    ssgsea_mes_adr_ncc_noradr_gsva_heatmap <- pheatmap::pheatmap(ssgsea_mes_adr_ncc_noradr_gsva,
                                                            scale = "none",
                                                            annotation_col = heatmap_col_annot,
                                                            cluster_rows = TRUE,
                                                            cluster_cols = FALSE,
                                                            color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                            show_colnames = FALSE)
    
    
    ggsave(filename = paste0(FIGURE_DIR, "varPart_noBatch_compare.png"), plot=varPart_noBatch_compare,
           width = 35,
           height = 25,
           units = "cm")
    
    # ConsensusTME ----
    # running unfiltered
    matchedSigs_preprocessed_genesets <- ConsensusTME::consensusGeneSets

    ssgsea_consensusTME_unfiltered <- GSVA::gsva(vst_matrix_GeneSymbols,
                                          matchedSigs_preprocessed_genesets$Unfiltered,
                                          method=c("ssgsea"),
                                          min.sz=1, max.sz=Inf, 
                                          ssgsea.norm=TRUE, verbose=TRUE, parallel.sz=10)
    
    ssgsea_consensusTME_unfiltered_noBatch <- GSVA::gsva(vst_matrix_noBatch_GeneSymbols,
                                          matchedSigs_preprocessed_genesets$Unfiltered,
                                          method=c("ssgsea"),
                                          min.sz=1, max.sz=Inf, 
                                          ssgsea.norm=TRUE, verbose=TRUE, parallel.sz=10)
    
    ssgsea_consensusTME_unfiltered_plot <- ssgsea_results_plot(ssgsea_results = ssgsea_consensusTME_unfiltered,
                                                           id = "omics_id",
                                                           metadata_df = cond_data,
                                                           cond_interest = cond_interest,
                                                           test_groups = list(c("no", "yes")),
                                                           plot_title = "consensusTME (unfiltered)")
    
    ssgsea_consensusTME_unfiltered_noBatch_plot <- ssgsea_results_plot(ssgsea_results = ssgsea_consensusTME_unfiltered_noBatch,
                                                               id = "omics_id",
                                                               metadata_df = cond_data,
                                                               cond_interest = cond_interest,
                                                               test_groups = list(c("no", "yes")),
                                                               plot_title = "consensusTME (unfiltered) - no batch")
    
    # combined MDSC signature ----
    ssgsea_combined_MDSC <- GSVA::gsva(vst_matrix_GeneSymbols,
                                          combined_MDSC_human_signature,
                                          method=c("ssgsea"),
                                          min.sz=1, max.sz=Inf, 
                                          ssgsea.norm=TRUE, verbose=TRUE, parallel.sz=10)
    
    ssgsea_combined_MDSC_noBatch <- GSVA::gsva(vst_matrix_noBatch_GeneSymbols,
                                                  combined_MDSC_human_signature,
                                                  method=c("ssgsea"),
                                                  min.sz=1, max.sz=Inf, 
                                                  ssgsea.norm=TRUE, verbose=TRUE, parallel.sz=10)
    
    ssgsea_combined_MDSC_plot <- ssgsea_results_plot(ssgsea_results = ssgsea_combined_MDSC,
                                                               id = "omics_id",
                                                               metadata_df = cond_data,
                                                               cond_interest = cond_interest,
                                                               test_groups = list(c("no", "yes")),
                                                               plot_title = "combined MDSC")
    
    ssgsea_combined_MDSC_noBatch_plot <- ssgsea_results_plot(ssgsea_results = ssgsea_combined_MDSC_noBatch,
                                                                       id = "omics_id",
                                                                       metadata_df = cond_data,
                                                                       cond_interest = cond_interest,
                                                                       test_groups = list(c("no", "yes")),
                                                                       plot_title = "combined MDSC - no batch")
    
    # MsigDB: C8:"HAY_BONE_MARROW" ----
    # msigdb: C8 cell_types
    gs_C8_celltypes <- msigdbr::msigdbr(species = "Homo sapiens", category = c("C8")) %>%
      dplyr::select(gs_name, gene_symbol)
    
    gs_C8_celltypes_boneMarrow <- gs_C8_celltypes %>%
      dplyr::filter(., grepl("BONE_MARROW", gs_name))
    table(gs_C8_celltypes_boneMarrow$gs_name)
    gs_C8_celltypes_boneMarrow_list_dfs <- split(x = gs_C8_celltypes_boneMarrow, f = gs_C8_celltypes_boneMarrow$gs_name)
    gs_C8_celltypes_boneMarrow_list <- purrr::map(.x = gs_C8_celltypes_boneMarrow_list_dfs, .f = function(x){
      x$gene_symbol
    })
    names(gs_C8_celltypes_boneMarrow_list) <- gsub(pattern = "HAY_BONE_MARROW_", replacement = "", names(gs_C8_celltypes_boneMarrow_list))
    
    ssgsea_gs_C8_celltypes_boneMarrow <- GSVA::gsva(vst_matrix_GeneSymbols,
                                                    gs_C8_celltypes_boneMarrow_list,
                                                    method=c("ssgsea"),
                                                    min.sz=1, max.sz=Inf, 
                                                    ssgsea.norm=TRUE, verbose=TRUE, parallel.sz=10)
    
    ssgsea_gs_C8_celltypes_boneMarrow_noBatch <- GSVA::gsva(vst_matrix_noBatch_GeneSymbols,
                                                            gs_C8_celltypes_boneMarrow_list,
                                                            method=c("ssgsea"),
                                                            min.sz=1, max.sz=Inf, 
                                                            ssgsea.norm=TRUE, verbose=TRUE, parallel.sz=10)
    
    ssgsea_gs_C8_celltypes_boneMarrow_plot <- ssgsea_results_plot(ssgsea_results = ssgsea_gs_C8_celltypes_boneMarrow,
                                                     id = "omics_id",
                                                     metadata_df = cond_data,
                                                     cond_interest = cond_interest,
                                                     test_groups = list(c("no", "yes")),
                                                     plot_title = "C8: HAY_BONE_MARROW_")
    
    ssgsea_gs_C8_celltypes_boneMarrow_noBatch_plot <- ssgsea_results_plot(ssgsea_results = ssgsea_gs_C8_celltypes_boneMarrow_noBatch,
                                                             id = "omics_id",
                                                             metadata_df = cond_data,
                                                             cond_interest = cond_interest,
                                                             test_groups = list(c("no", "yes")),
                                                             plot_title = "C8: HAY_BONE_MARROW_ - no batch")
    
    
    
    # saving results ----
    save(vst_matrix_GeneSymbols, vst_matrix_noBatch_GeneSymbols,
         vst_mcpCounter, vst_mcpCounter_noBatch,
         ssgsea_consensusTME_unfiltered, ssgsea_consensusTME_unfiltered_noBatch,
         ssgsea_combined_MDSC, ssgsea_combined_MDSC_noBatch,
         ssgsea_gs_C8_celltypes_boneMarrow, ssgsea_gs_C8_celltypes_boneMarrow_noBatch,
         file = paste0(OUTPUT_DIR, experiment_name, "_TME_results.RData"))
    
    save(vst_mcpCounter_plot, vst_mcpCounter_noBatch_plot,
         ssgsea_consensusTME_unfiltered_plot, ssgsea_consensusTME_unfiltered_noBatch_plot,
         ssgsea_combined_MDSC_plot, ssgsea_combined_MDSC_noBatch_plot,
         ssgsea_gs_C8_celltypes_boneMarrow_plot, ssgsea_gs_C8_celltypes_boneMarrow_noBatch_plot,
         file = paste0(OUTPUT_DIR, experiment_name, "_TME_resultsPlots.RData"))
    
    # testing wilcoxon ----
    # ssgsea_gs_C8_celltypes_boneMarrow_noBatch_df_CD34_POS_LMPP <- as.data.frame(ssgsea_gs_C8_celltypes_boneMarrow_noBatch) %>%
    #   tibble::rownames_to_column(., var = "cell_type") %>%
    #   tidyr::pivot_longer(names_to = "omics_id", values_to = "enrichment_score", -cell_type) %>%
    #   dplyr::left_join(., dplyr::select(cond_data, c("omics_id", tidyselect::all_of(cond_interest))), by = "omics_id") %>% #%>%  dplyr::select(-omics_id)
    #   dplyr::filter(cell_type == "CD34_POS_LMPP") #%>%
    #   wilcox.test(enrichment_score ~ material, data = .)
    # 
    # ssgsea_gs_C8_celltypes_boneMarrow_noBatch_df_CD34_POS_LMPP_split <- split(x = ssgsea_gs_C8_celltypes_boneMarrow_noBatch_df_CD34_POS_LMPP, f = ssgsea_gs_C8_celltypes_boneMarrow_noBatch_df_CD34_POS_LMPP$material)
    # wilcox.test(x = ssgsea_gs_C8_celltypes_boneMarrow_noBatch_df_CD34_POS_LMPP_split$yes$enrichment_score, 
    #             y = ssgsea_gs_C8_celltypes_boneMarrow_noBatch_df_CD34_POS_LMPP_split$no$enrichment_score)
    
    # Saving TME plots ----
    # replotting with pdf and svg
    # if exists only load!
    load(file = paste0(OUTPUT_DIR, experiment_name, "_TME_results.RData"))
    load(file = paste0(OUTPUT_DIR, experiment_name, "_TME_resultsPlots.RData"))
    
    celltype_deconvolution_DIR = paste0(OUTPUT_DIR, "celltype_deconvolution/")
    dir.create(celltype_deconvolution_DIR, showWarnings = TRUE)
    
    # MCP-Counter
    ggsave(filename = paste0(celltype_deconvolution_DIR, "mcpCounter_plot.png"), 
           plot=vst_mcpCounter_plot$mcpCounter_results_plot, 
           width = 25, height = 25, units = "cm") 
    ggsave(filename = paste0(celltype_deconvolution_DIR, "mcpCounter_noBatch_plot.png"), 
           plot=vst_mcpCounter_noBatch_plot$mcpCounter_results_plot, 
           width = 25, height = 25, units = "cm") 
    # consensusTME
    ggsave(filename = paste0(celltype_deconvolution_DIR, "ssgsea_consensusTME_unfiltered_plot.png"), 
           plot=ssgsea_consensusTME_unfiltered_plot$ssgsea_results_plot, 
           width = 25, height = 25, units = "cm") 
    ggsave(filename = paste0(celltype_deconvolution_DIR, "ssgsea_consensusTME_unfiltered_noBatch_plot.png"), 
           plot=ssgsea_consensusTME_unfiltered_noBatch_plot$ssgsea_results_plot, 
           width = 25, height = 25, units = "cm") 
    # combined MDSC signature
    ggsave(filename = paste0(celltype_deconvolution_DIR, "ssgsea_combined_MDSC_plot.png"), 
           plot=ssgsea_combined_MDSC_plot$ssgsea_results_plot, 
           width = 10, height = 10, units = "cm") 
    ggsave(filename = paste0(celltype_deconvolution_DIR, "ssgsea_combined_MDSC_noBatch_plot.png"), 
           plot=ssgsea_combined_MDSC_noBatch_plot$ssgsea_results_plot, 
           width = 10, height = 10, units = "cm") 
    # combined MDSC signature
    ggsave(filename = paste0(celltype_deconvolution_DIR, "ssgsea_gs_C8_celltypes_boneMarrow_plot.png"), 
           plot=ssgsea_gs_C8_celltypes_boneMarrow_plot$ssgsea_results_plot, 
           width = 33, height = 33, units = "cm") 
    ggsave(filename = paste0(celltype_deconvolution_DIR, "ssgsea_gs_C8_celltypes_boneMarrow_noBatch_plot.png"), 
           plot=ssgsea_gs_C8_celltypes_boneMarrow_noBatch_plot$ssgsea_results_plot, 
           width = 33, height = 33, units = "cm") 
    
    # add cibersort results!

    # MCP-Counter
    ggsave(filename = paste0(celltype_deconvolution_DIR, "mcpCounter_plot.pdf"), 
           plot=vst_mcpCounter_plot$mcpCounter_results_plot, 
           width = 25, height = 25, units = "cm") 
    ggsave(filename = paste0(celltype_deconvolution_DIR, "mcpCounter_noBatch_plot.pdf"), 
           plot=vst_mcpCounter_noBatch_plot$mcpCounter_results_plot, 
           width = 25, height = 25, units = "cm") 
    # consensusTME
    ggsave(filename = paste0(celltype_deconvolution_DIR, "ssgsea_consensusTME_unfiltered_plot.pdf"), 
           plot=ssgsea_consensusTME_unfiltered_plot$ssgsea_results_plot, 
           width = 25, height = 25, units = "cm") 
    ggsave(filename = paste0(celltype_deconvolution_DIR, "ssgsea_consensusTME_unfiltered_noBatch_plot.pdf"), 
           plot=ssgsea_consensusTME_unfiltered_noBatch_plot$ssgsea_results_plot, 
           width = 25, height = 25, units = "cm") 
    # combined MDSC signature
    ggsave(filename = paste0(celltype_deconvolution_DIR, "ssgsea_combined_MDSC_plot.pdf"), 
           plot=ssgsea_combined_MDSC_plot$ssgsea_results_plot, 
           width = 10, height = 10, units = "cm") 
    ggsave(filename = paste0(celltype_deconvolution_DIR, "ssgsea_combined_MDSC_noBatch_plot.pdf"), 
           plot=ssgsea_combined_MDSC_noBatch_plot$ssgsea_results_plot, 
           width = 10, height = 10, units = "cm") 
    # combined MDSC signature
    ggsave(filename = paste0(celltype_deconvolution_DIR, "ssgsea_gs_C8_celltypes_boneMarrow_plot.pdf"), 
           plot=ssgsea_gs_C8_celltypes_boneMarrow_plot$ssgsea_results_plot, 
           width = 33, height = 33, units = "cm") 
    ggsave(filename = paste0(celltype_deconvolution_DIR, "ssgsea_gs_C8_celltypes_boneMarrow_noBatch_plot.pdf"), 
           plot=ssgsea_gs_C8_celltypes_boneMarrow_noBatch_plot$ssgsea_results_plot, 
           width = 33, height = 33, units = "cm") 
    
  }
  
    
}

exportResults=TRUE
if (exportResults) {
  # using filtering to remove lowly expressed genes! and avoiding 0 count genes?, but this genes may be important in terms of outliers in expression
  # if doing changes in design or using filtered this ensures the final dds to extrac results has this changes without overriding original files
  dds_filt_extractResults <- dds_filt
  
  # adjusting for batch effects/ individual id
  #  the condition effect represents the overall effect controlling for differences due to individual mouse_ids
  # cond_interest_varPart
  # potentially accounting for gender as well; examine gender effects
  #table(coldata$material, coldata$material)
  #table(coldata$library_type, coldata$material)
  # + material - material removed; add gender?!!
  new_design_variables <- c("library_type", "material")
  new_design <- paste0("~", paste(new_design_variables, collapse = " + "))
  if (!(design(dds_filt_extractResults) == as.formula(new_design))) {
    message("changing design to: ", new_design)
    design(dds_filt_extractResults) <- as.formula(new_design)
    
    #design(dds_filt_extractResults) <- ~ library_type + material 
    
    # removing one sample where material == NA
    # not removing sample without material if not considering material!!!
    # material_NA_sample <- coldata %>%
    #   dplyr::filter(is.na(material))
    #dds_filt_extractResults <- dds_filt_extractResults[, !(colnames(dds_filt_extractResults) %in% material_NA_sample$omics_id)]
    
    dds_filt_extractResults <- DESeq(dds_filt_extractResults) # replacing outliers and refitting
  }

  experiment_comparisons <- resultsNames(dds_filt_extractResults) 
  experiment_comparisons <- experiment_comparisons[experiment_comparisons != "Intercept"] # removing intercept
  cat("There are", length(experiment_comparisons),"comparisons in the result table.\n")
  
  
  # [] add function to export results, annotate, write to excel
  
  results_coeff <- experiment_comparisons
  results_summary_df <- data.frame(test = as.character(),
                                   design = as.character(),
                                   signif_genes = as.numeric(),
                                   signif_genes_UP = as.numeric(),
                                   signif_genes_DOWN = as.numeric())
  
  
  cat("Generating results...\n")
  # run also other comparisons!
  # 1:length(results_coeff)
  
  for (res_test in 1:length(results_coeff)) {
   
    condition_test <- results_coeff[[res_test]] # "material_high_vs_no_or_low" # 
    condition_variable <- gsub(pattern = paste0("(",paste(new_design_variables, collapse = "|"),")(_.+_vs_.+)"), 
                               replacement = "\\1",
                               x = condition_test)
    
    cond_numerator <- gsub(pattern = paste0("(",condition_variable,"_)(.+)(_vs_.+)"), replacement = "\\2", condition_test)
    cond_denominator <- gsub(pattern = paste0("(",condition_variable,"_)(.+_vs_)(.+)"), replacement = "\\3", condition_test)
    message("condition_variable: ", condition_variable, 
            " cond_numerator: ", cond_numerator,
            " cond_denominator: ", cond_denominator)
    #message(cond_numerator)
    #message(cond_denominator)
  
    RESULTS_DIR = paste0(OUTPUT_DIR, condition_test, "/")
    cat("Saving results in RESULTS_DIR:", RESULTS_DIR, "\n")
    dir.create(RESULTS_DIR, showWarnings = TRUE)
    
    # Generating results for each of the comparisons
    results_all <- generateResults(dds_object=dds_filt_extractResults, res_extract=condition_test, annot_df=ensemblAnnot, 
                                   variable=condition_variable, padj_cutoff=padj_cutoff, shrinkage_type="apeglm") # "normal"
    
    results_data <- results_all$results_data
    results_data_annot <- results_all$results_data_annot
    
    # extracting only significant results
    results_data_annot_signif <- results_data_annot %>%
      dplyr::filter((!is.na(padj) & (padj < padj_cutoff)) & abs(log2FoldChange) > log2FC_cutoff)
    
    results_data_annot_signif_geneGO <- results_data_annot_signif %>%
      dplyr::select(ensembl_id, log2FoldChange, padj)
    
    temp_results_summary_df <- data.frame(test = condition_test,
                                          design = paste(as.character(design(dds_filt_extractResults)), collapse=""),
                                          signif_genes = nrow(results_data_annot_signif),
                                          signif_genes_UP = sum(results_data_annot_signif$log2FoldChange > log2FC_cutoff),
                                          signif_genes_DOWN = sum(results_data_annot_signif$log2FoldChange < log2FC_cutoff))
    results_summary_df <- rbind(results_summary_df, temp_results_summary_df)
    
    if(nrow(results_data_annot_signif) != 0){
      # Capturing gene_biotype, chromosome, strand distribution in significantly DEGs
      gene_biotype_dist <- as.data.frame(table(results_data_annot_signif$gene_biotype)) %>%
        rename(gene_biotype=Var1) %>%
        arrange(desc(Freq))
      
      gene_biotype_dist_plot <- 
        ggplot(results_data_annot_signif, aes(forcats::fct_rev(forcats::fct_infreq(gene_biotype)))) + #Reorder factors levels by frequency, reverse order - fct_rev
        geom_bar() + 
        labs(x="gene biotype") +
        coord_flip() +
        theme_bw()
      
      chromosomes_dist <- as.data.frame(table(results_data_annot_signif$chromosome_name)) %>%
        rename(chromosome=Var1) %>%
        arrange(desc(Freq))
      
      chromosome_dist_plot <- 
        ggplot(results_data_annot_signif, aes(chromosome_name)) + 
        geom_bar() + 
        scale_x_discrete(limits=gtools::mixedsort(unique(results_data_annot_signif$chromosome_name))) + # natural sort for chromosome names
        labs(x="chromosome") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      theme_bw()
      
      strand_dist <- as.data.frame(table(results_data_annot_signif$strand)) %>%
        rename(strand=Var1) %>%
        arrange(desc(Freq))
      save(results_data_annot, results_data_annot_signif, gene_biotype_dist_plot, chromosome_dist_plot, 
           file = paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_resultsDE.RData"))
      ggsave(filename = paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_signif_GeneBiotypeDist.png"), 
             plot = gene_biotype_dist_plot,
             width = 20, height = 20, units = "cm")
      ggsave(filename = paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_signif_ChromosomeDist.png"),
             plot = chromosome_dist_plot,
             width = 20, height = 20, units = "cm")
    } else {
      save(results_data_annot, results_data_annot_signif,
           file = paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_resultsDE.RData"))
    }
    
    # saving results for use in notebook
    save(dds_filt_extractResults, file = paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_dds_filt_extractResults.RData"))
    
    writeLines(capture.output(summary(results_data, alpha = padj_cutoff)), paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_summary.csv"))
    #write_csv(results_data_annot, file=paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_DE_results.csv"))
    #write_csv(results_data_annot_signif, file=paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_DE_results_signif.csv"))
    write_tsv(results_data_annot_signif_geneGO, file=paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_signifGeneGO.tsv"))
    write_csv(results_summary_df, file=paste0(OUTPUT_DIR, "results_summary.csv"))

    #       dplyr::filter((!is.na(padj) & (padj < padj_cutoff)) & abs(log2FoldChange) > log2FC_cutoff)
    analysis_notes = data.frame(design_formula = new_design,
                                padj_cutoff = padj_cutoff,
                                log2FC_cutoff = log2FC_cutoff,
                                signif_filter = paste0("!is.na(padj) & (padj < padj_cutoff)) & abs(log2FoldChange) > log2FC_cutoff"))
    de_results_data_list <- list(DE_results = results_data_annot,
                                 DE_results_signif = results_data_annot_signif,
                                 analysis_notes = analysis_notes)
    openxlsx::write.xlsx(de_results_data_list, file = paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_DE_results.xlsx")) 
        
    gsea_enrichment = TRUE
    if (gsea_enrichment) {
  
      # log2FC ranked list
      rnk_log2FC_geneList <- generateRankedList(results_data=results_data_annot, ranking="log2FC", convert2entrezIDs = TRUE)
      
      #quantile(rnk_log2FC_geneList)
      rnk_log2FC_geneList_df <- data.frame(log2FoldChange=rnk_log2FC_geneList,
                                           index = seq(1:length(rnk_log2FC_geneList))) 
      
      #plot(rnk_log2FC_geneList)
      rnk_log2FC_geneList_plot <- ggplot(rnk_log2FC_geneList_df, aes(x=index, y=log2FoldChange)) + geom_bar(stat="identity") + theme_bw()
      
      gsea_function <- function(geneset=NULL, geneset_name=NULL, ranked_list=NULL, ncores=6){
        # added geneset_name mainly for debugging and printing in log
        # pvalueCutoff = 1, no filtering at this stage
        # For some pathways, in reality P-values are less than 1e-10. You can set the `eps` argument to zero for better estimation.
        # consider eps = 0
        message("...enriching: ", geneset_name)
        clusterProfiler::GSEA(geneList = ranked_list,
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
      
      furrr_map_options <- furrr_options(seed = 123)
      future::plan(multisession, workers = 6)
      #future::plan(sequential)
      
      gsea_results <- furrr::future_map2(.x = gs_entreIDs_collection,
                                         .y = names(gs_entreIDs_collection),
                                         .f = function(geneset, geneset_name)
                                           gsea_function(geneset=geneset,
                                                         geneset_name=geneset_name,
                                                         ranked_list=rnk_log2FC_geneList),
                                         .options = furrr_map_options,
                                         .progress = FALSE) 
      
      future::plan(sequential)
      
      # converting gsea to data.frames and adding gene symbol core_enrichment
      entrez2symbol_df <- ensemblAnnot %>%
        dplyr::select(entrezgene, hgnc_symbol) %>%
        dplyr::filter(hgnc_symbol != "")
      
      gsea_results_df <- purrr::map(.x = gsea_results, .f = function(x) {
        temp_df <- as.data.frame(x) 
        temp_df_coreSymbols <- temp_df %>%
          tidyr::separate_rows(., core_enrichment, sep = "/", convert = FALSE) %>%
          dplyr::left_join(., entrez2symbol_df, by = c("core_enrichment" = "entrezgene")) %>%
          dplyr::rename(core_enrichment_symbols = hgnc_symbol) %>%
          dplyr::distinct(.) %>%
          dplyr::group_by(ID) %>%
          dplyr::summarise(across(.cols = everything(), .fns = ~ paste(unique(.x), collapse = "/")), .groups = "drop") %>%
          dplyr::mutate(across(c("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalues", "rank"), as.double)) %>%
          dplyr::mutate(across(c("setSize"), as.integer)) %>%     
          dplyr::ungroup(.) %>%
          dplyr::arrange(match(ID, temp_df$ID)) 
        
        return(temp_df_coreSymbols)
        })
      
      # saving data ----
      ggsave(rnk_log2FC_geneList_plot, filename = paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_rnk_log2FC_geneList_plot.png")) 
      openxlsx::write.xlsx(gsea_results_df, file = paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_gsea_results.xlsx")) 
      save(gsea_results, gsea_results_df, file = paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_resultsGSEA.RData"))                     
    }
    
  }
}

# saving data matrix and cond_data
# also plot dispersion to see outliers?
# what is replaceable? where is it coming from
summary(results_all$results_data)
cond_data_extractResults <- as.data.frame(colData(dds_filt_extractResults))

readr::write_tsv(x = cond_data_extractResults,
                 file = paste0(OUTPUT_DIR, experiment_name, "_annotData.tsv"))

vst_filt_extractResults <- rlog(dds_filt_extractResults, blind = FALSE) #vst(dds_filt_extractResults, blind=FALSE)
vst_filt_extractResults_counts <- assay(vst_filt_extractResults) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., var = "ensembl_id")

norm_extractResults_counts <- counts(dds_filt_extractResults, normalized=TRUE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., var = "ensembl_id")

readr::write_tsv(x = vst_filt_extractResults_counts,
                 file = paste0(OUTPUT_DIR, experiment_name, "_vst_counts.tsv"))

readr::write_tsv(x = norm_extractResults_counts,
                 file = paste0(OUTPUT_DIR, experiment_name, "_norm_counts.tsv"))


# generate plots
generate_pub_plots <- TRUE
if (generate_pub_plots) {
# generating key plots for publication
# DE genes plots
#  1a. heatmap of top 30-50 genes - add information about infiltration and library type
#  1b. volcano plot and highlight 30-50 genes
#  2. For these 30-50 genes and boxplot/violin, cloud_plot or half plot betweem yes, no group and separation based on library type
# GSEA:
# top 20 upregulated: KEGG, GOBP, Hallmarks
# ssGSEA plot as pdf (svg, eps)
# save(results_data_annot, results_data_annot_signif, gene_biotype_dist_plot, chromosome_dist_plot, 
#        file = paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_results.RData"))  
# save(dds_filt_extractResults, file = paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_dds_filt_extractResults.RData"))  
# save(gsea_results, gsea_results_df, file = paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_resultsGSEA.RData"))
# TO-DO:
# [x] send plots
# [x] add violin with library_type effect and calculate per gene varPartition!
# [x] use rlog transformed log2 counts! for violin plots; use batch corrected and uncorrected!
# [] plot both splitted and unsplitted
# [] summarise discussion: we want gsea_results_dotPlots_v2; waitint for summary from Daria; not about non-significant genes and 
# [] plot QC for selected samples
# [] try all unstranded
# [] explain main problem is coming from only 2 samples rather than batch for which you can control 
# [] plot violin plot for subset of genes and GSEA plot for subset of genesets
  
# re-check results from previous and 
PUBLICATION_FIGURES_DIR="/home/peter_r/results/Results/material_timepointDE/material_yes_vs_no/pub_figures/"
dir.create(PUBLICATION_FIGURES_DIR)
load("/home/peter_r/results/Results/material_timepointDE/material_yes_vs_no/yes_vs_no_resultsDE.RData")  # results_data_annot_signif
load("/home/peter_r/results/Results/material_timepointDE/material_yes_vs_no/yes_vs_no_dds_filt_extractResults.RData")  # dds_filt_extractResults
load("/home/peter_r/results/Results/material_timepointDE/material_yes_vs_no/yes_vs_no_resultsGSEA.RData")  # GSEA results

# filtering top30 based on padj
# identical(head(results_data_annot_signif$ensembl_id, 30), results_data_annot_signif_top30$ensembl_id)
# generating rlog matrix
# [] plot 
signif_rld <- rlog(dds_filt_extractResults, blind = FALSE) # not blind to batch effects
signif_rld_noBatch <- remove_batch(vst_transf = signif_rld,
                                   batch_variable = "library_type",
                                   cond_interest = cond_interest)
signif_rld_counts <- assay(signif_rld)
signif_rld_counts_noBatch <- assay(signif_rld_noBatch)
signif_normCounts <- DESeq2::counts(object = dds_filt_extractResults, normalized=TRUE)
cond_data <- as.data.frame(colData(dds_filt_extractResults))

# for fun checking expressio nmatrix ----
# especially for project for Johannes test GC and length correction! using:
# http://bioconductor.org/packages/release/bioc/vignettes/EDASeq/inst/doc/EDASeq.html#retrieving-gene-length-and-gc-content
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#sample-gene-dependent-normalization-factors
# see how to obtain GC, length from biomart 
# how to ofset DESeq2
# compare to salmon results and uncorrected results; compare on filtered, unfiltered counts,...

# https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
# Since the majority of genes are not differentially expressed, the majority of genes in each sample should have similar ratios within the sample.

# see DESeq2 size factor calculation: ?estimateSizeFactors 
#  and corresponding equation in description
# default_size_factors <- DESeq2::sizeFactors(dds_filt_extractResults)
# raw_count_matrix <- DESeq2::counts(dds_filt_extractResults, normalized = FALSE)
# 
# # using function from DESeq2 - see 
# #  estimateSizeFactors.DESeqDataSet https://github.com/mikelove/DESeq2/blob/master/R/methods.R
# # poscounts
# geoMeanNZ <- function(x) {
#   if (all(x == 0)) { 0 } else { exp( sum(log(x[x > 0])) / length(x) ) }
# }
# 
# # gm_mean = function(x, na.rm=TRUE){
# #   # https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
# #   exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
# #   #exp(mean(log(x)))
# # }
# 
# # following method on https://github.com/mikelove/DESeq2/blob/5a32b40b563e7db35fc3dafed705269695ae2dca/R/core.R
# # poscounts
# lc <- log(raw_count_matrix)
# lc[!is.finite(lc)] <- 0
# loggeomeans <- rowMeans(lc)
# allZero <- rowSums(raw_count_matrix) == 0
# loggeomeans[allZero] <- -Inf
# 
# # just ratio
# loggeomeans <- rowMeans(log(raw_count_matrix))
# sf_default <- apply(raw_count_matrix, 2, function(cnts) {
#   locfunc <- stats::median
#   exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))
# })
# 
# # my approach...
# raw_count_gene_geomMeans <- apply(raw_count_matrix, 1, geoMeanNZ)
# raw_count_gene_ratio <- sweep(raw_count_matrix, 1, raw_count_gene_geomMeans, FUN = "/")
# raw_count_gene_MedianRatio <- apply(raw_count_gene_ratio, 2, median)


ntop_signif = 50  # nrow(results_data_annot_signif) # 50 
# dim(results_data_annot_signif)

# NA in pvalues????
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
# Cook's distance
# removing NA from annotated results - all
results_data_annot_pvalNA <- results_data_annot %>%
  dplyr::filter(is.na(pvalue))

results_data_annot_noNA <- results_data_annot %>%
  dplyr::filter(!is.na(pvalue)) %>%
  dplyr::mutate(gene_name = if_else(hgnc_symbol == "", ensembl_id, hgnc_symbol)) %>% # if hgnc_symbol not present use ensembl_id
  dplyr::select(gene_name, ensembl_id, hgnc_symbol, description, log2FoldChange, padj, pvalue)

results_data_annot_signif_Ntop <- results_data_annot_signif %>%
  dplyr::filter(hgnc_symbol != "") %>%  # removing genes without symbol
  dplyr::arrange(desc(log2FoldChange)) %>%  # ordering based on log2FC check if only UPREGULATED or limit; abs(log2FoldChange)
  dplyr::filter(log2FoldChange > 0) %>% # only positive
  #dplyr::arrange(padj) %>%  # ordering based on padj
  dplyr::slice_head(n = ntop_signif) %>%
  dplyr::mutate(gene_name = if_else(hgnc_symbol == "", ensembl_id, hgnc_symbol)) %>% # if hgnc_symbol not present use ensembl_id
  dplyr::select(gene_name, ensembl_id, hgnc_symbol, description, log2FoldChange, padj, pvalue)

sum(results_data_annot_signif_Ntop$log2FoldChange < 0)
sum(results_data_annot_signif_Ntop$hgnc_symbol == "")

# cond_data_signif_Ntop <- cond_data %>%
#   dplyr::

# extracting rlog normalized matrix for plotting heatmap ----
signif_Ntop_rld_counts_mat <- as.data.frame(signif_rld_counts) %>%
  tibble::rownames_to_column(var="ensembl_id") %>%
  dplyr::filter(ensembl_id %in% results_data_annot_signif_Ntop$ensembl_id) %>%
  dplyr::left_join(., dplyr::select(results_data_annot_signif_Ntop, gene_name, ensembl_id), by = "ensembl_id") %>%
  dplyr::select(-ensembl_id) %>%
  tibble::column_to_rownames(var = "gene_name") %>%
  as.matrix(.)

# batch effect removed
signif_Ntop_rld_counts_mat_noBatch <- as.data.frame(signif_rld_counts_noBatch) %>%
  tibble::rownames_to_column(var="ensembl_id") %>%
  dplyr::filter(ensembl_id %in% results_data_annot_signif_Ntop$ensembl_id) %>%
  dplyr::left_join(., dplyr::select(results_data_annot_signif_Ntop, gene_name, ensembl_id), by = "ensembl_id") %>%
  dplyr::select(-ensembl_id) %>%
  tibble::column_to_rownames(var = "gene_name") %>%
  as.matrix(.)

# scaling - calculating z-scores ----
signif_Ntop_rld_counts_mat_scaled <- t(scale(t(signif_Ntop_rld_counts_mat), center = TRUE, scale = TRUE))
# same scaling as in ComplexHeatmap/Pheatmap scale = "row"; identical(signif_Ntop_heatmap@matrix, signif_Ntop_rld_counts_mat_scaled)
# restrict ranges 
signif_Ntop_rld_counts_mat_scaled_quantiles <- quantile(signif_Ntop_rld_counts_mat_scaled)
max_scale_limit <- floor(min(abs(signif_Ntop_rld_counts_mat_scaled_quantiles[1]), abs(signif_Ntop_rld_counts_mat_scaled_quantiles[5])))
signif_Ntop_rld_counts_mat_scaled[signif_Ntop_rld_counts_mat_scaled > max_scale_limit] <- max_scale_limit  # clip values > 3.0
signif_Ntop_rld_counts_mat_scaled[signif_Ntop_rld_counts_mat_scaled < -1*max_scale_limit] <- -1*max_scale_limit  # clip values < -3.0
quantile(signif_Ntop_rld_counts_mat_scaled)

# batch effect removed 
signif_Ntop_rld_counts_mat_scaled_noBatch <- t(scale(t(signif_Ntop_rld_counts_mat_noBatch), center = TRUE, scale = TRUE))
signif_Ntop_rld_counts_mat_scaled_noBatch_quantiles <- quantile(signif_Ntop_rld_counts_mat_scaled_noBatch)
max_scale_limit_noBatch <- floor(min(abs(signif_Ntop_rld_counts_mat_scaled_noBatch_quantiles[1]), abs(signif_Ntop_rld_counts_mat_scaled_noBatch_quantiles[5])))
signif_Ntop_rld_counts_mat_scaled_noBatch[signif_Ntop_rld_counts_mat_scaled_noBatch > max_scale_limit_noBatch] <- max_scale_limit_noBatch  # clip values > 3.0
signif_Ntop_rld_counts_mat_scaled_noBatch[signif_Ntop_rld_counts_mat_scaled_noBatch < -1*max_scale_limit_noBatch] <- -1*max_scale_limit_noBatch  # clip values < -3.0
quantile(signif_Ntop_rld_counts_mat_scaled_noBatch)

# plotting heatmap ----
# identical(rownames(cond_data), cond_data$omics_id)
# results_data_annot_signif_Ntop_annot <- cond_data %>%
#   dplyr::select(material, infiltration_percentage, library_type)
# interesting color palletes
#devtools::install_github("karthik/wesanderson")

# generating heatmap annotation ----
# results_data_annot_signif_Ntop_annot <- as_tibble(cond_data) %>%
#   dplyr::select(material, infiltration_percentage, library_type) %>%
#   dplyr::mutate(infiltration_percentage = ifelse(infiltration_percentage == "<1", "0.9", infiltration_percentage)) %>%
#   dplyr::mutate(infiltration_percentage = as.numeric(infiltration_percentage)) %>%
#   #dplyr::ungroup() %>%
#   dplyr::group_by(material) %>%
#   dplyr::arrange(infiltration_percentage, .by_group = TRUE) %>%
#   as.data.frame(.)

results_data_annot_signif_Ntop_annot <- as_tibble(cond_data) %>%
  dplyr::select(material, infiltration_percentage, library_type) %>%
  dplyr::mutate(infiltration_percentage = case_when(
    infiltration_percentage == "0" ~ "0",
    infiltration_percentage %in% c("<1", "1") ~ ">1",
    infiltration_percentage %in% c("10", "15") ~ ">=10",
    infiltration_percentage %in% c("20", "30") ~ ">=20",
    infiltration_percentage %in% c("50", "60","70") ~ ">=50",
    TRUE ~ as.character(infiltration_percentage)
  )) %>%
  dplyr::mutate(infiltration_percentage = factor(infiltration_percentage, levels = c("0", ">1", ">=10", ">=20", ">=50"))) %>%
  #dplyr::ungroup() %>%
  dplyr::group_by(material) %>%
  dplyr::arrange(infiltration_percentage, .by_group = TRUE) %>%
  as.data.frame(.)

# ann_colors = list(
#   material = c(no = "#7570B3", yes = "#E7298A"), #, Path3 = "#66A61E"
#   library_type = c(U = "#1B9E77", SR = "#D95F02")
# )
# https://github.com/karthik/wesanderson
#infiltration_perc_pal <- wesanderson::wes_palette("Zissou1", 10, type = "continuous")
infiltration_perc_pal <- ggsci::pal_material("purple")(4)
#scales::show_col(ggsci::pal_material("purple")(4))

# for consistency using same color for TME, heatmap, gene comparison 
# scales::show_col(ggsci::pal_jco("default")(10))
ann_colors = list(
  material = c(no = "#0073C2FF", yes = "#EFC000FF"), #, Path3 = "#66A61E"
  library_type = c(U = "#1B9E77", SR = "#D95F02"),
  infiltration_percentage = c("0" = "lightgray", 
                              ">1" = infiltration_perc_pal[1], 
                              ">=10" = infiltration_perc_pal[2], 
                              ">=20" = infiltration_perc_pal[3], 
                              ">=50" = infiltration_perc_pal[4])
)

library(ComplexHeatmap)
# improvements:
# [ ] annotation move on top of each other
# [ ] order by percentage

fontsize_setting=14
signif_Ntop_heatmap <- ComplexHeatmap::pheatmap(signif_Ntop_rld_counts_mat_scaled,
                                                #annotation_col = results_data_annot_signif_Ntop_annot,
                                                #annotation_colors = ann_colors,
                                                col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                show_rownames = TRUE,
                                                show_colnames = FALSE,
                                                cluster_cols = FALSE,
                                                cluster_rows = TRUE,
                                                legend_breaks = seq(from = -1*max_scale_limit, to = max_scale_limit, by = 1),
                                                border_color = NA, #  if too many genes; https://github.com/jokergoo/ComplexHeatmap/issues/651
                                                #column_split = results_data_annot_signif_Ntop_annot$material,
                                                scale = "none",
                                                fontsize=fontsize_setting,
                                                top_annotation = HeatmapAnnotation(df = results_data_annot_signif_Ntop_annot,
                                                                                   col = ann_colors,
                                                                                   annotation_name_side = "right",
                                                                                   annotation_name_gp = gpar(fontsize = fontsize_setting),
                                                                                   annotation_legend_param = list(title_gp = gpar(fontsize = fontsize_setting),
                                                                                                                  labels_gp = gpar(fontsize = fontsize_setting))),
                                                heatmap_legend_param = list(title = "Z-score", # "z-score",
                                                                            #legend_direction="horizontal",
                                                                            #title_position="topcenter",
                                                                            #at = seq(from = -1*max_scale_limit, to = max_scale_limit, by = 1),
                                                                            title_gp = gpar(fontsize = fontsize_setting),
                                                                            labels_gp = gpar(fontsize = fontsize_setting),
                                                                            legend_height = unit(6, "cm")))

# batch effect removed 
signif_Ntop_heatmap_noBatch <- ComplexHeatmap::pheatmap(signif_Ntop_rld_counts_mat_scaled_noBatch,
                                                #annotation_col = results_data_annot_signif_Ntop_annot,
                                                #annotation_colors = ann_colors,
                                                col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                show_rownames = TRUE,
                                                show_colnames = FALSE,
                                                cluster_cols = FALSE,
                                                cluster_rows = TRUE,
                                                legend_breaks = seq(from = -1*max_scale_limit, to = max_scale_limit, by = 1),
                                                border_color = NA, #  if too many genes; https://github.com/jokergoo/ComplexHeatmap/issues/651
                                                #column_split = results_data_annot_signif_Ntop_annot$material,
                                                scale = "none",
                                                fontsize=fontsize_setting,
                                                top_annotation = HeatmapAnnotation(df = results_data_annot_signif_Ntop_annot,
                                                                                   col = ann_colors,
                                                                                   annotation_name_side = "right",
                                                                                   annotation_name_gp = gpar(fontsize = fontsize_setting),
                                                                                   annotation_legend_param = list(title_gp = gpar(fontsize = fontsize_setting),
                                                                                                                  labels_gp = gpar(fontsize = fontsize_setting))),
                                                heatmap_legend_param = list(title = "Z-score", # "z-score",
                                                                            #legend_direction="horizontal",
                                                                            #title_position="topcenter",
                                                                            #at = seq(from = -1*max_scale_limit, to = max_scale_limit, by = 1),
                                                                            title_gp = gpar(fontsize = fontsize_setting),
                                                                            labels_gp = gpar(fontsize = fontsize_setting),
                                                                            legend_height = unit(6, "cm")))

# compare batch and noBatch ----
# NOTE to generate this cluster_rows = FALSE for both heatmaps!!! - rerun for final plot with cluster_rows = TRUE
ht_list = signif_Ntop_heatmap + signif_Ntop_heatmap_noBatch

pdf(paste0(PUBLICATION_FIGURES_DIR, "heatmap_signifDE_ntop", ntop_signif, "_compareBatchRemoval.pdf"),  width = (35/2.5), height = (30/2.5)) # ,  width = 25, height = 30
ComplexHeatmap::draw(ht_list, column_title = "Comparison of batch uncorrected (left) and corrected (right) heatmaps for top50 signif. upregulated genes", 
                     column_title_gp = gpar(fontsize = 14),
                     ht_gap = unit(2, "cm"),
                     auto_adjust = FALSE)  # auto_adjust = FALSE to keep rownames
#print(Heatmap(dat))
dev.off()

# saving heatmap
pdf(paste0(PUBLICATION_FIGURES_DIR, "heatmap_signifDE_ntop", ntop_signif, ".pdf"),  width = (25/2.5), height = (30/2.5)) # ,  width = 25, height = 30
ComplexHeatmap::draw(signif_Ntop_heatmap, merge_legend = FALSE)
#print(Heatmap(dat))
dev.off()

#signif_Ntop_heatmap2 <- ComplexHeatmap::rowAnnotation(expression = ComplexHeatmap::anno_barplot(1:nrow(signif_Ntop_rld_counts_mat)))
#signif_Ntop_heatmap + signif_Ntop_heatmap2

# Volcano plot ----
# add a column of NAs
results_data_annot_noNA$signif_DE <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
results_data_annot_noNA$signif_DE[results_data_annot_noNA$log2FoldChange > log2FC_cutoff & results_data_annot_noNA$padj < padj_cutoff] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
results_data_annot_noNA$signif_DE[results_data_annot_noNA$log2FoldChange < -log2FC_cutoff & results_data_annot_noNA$padj < padj_cutoff] <- "DOWN"

#install.packages("gghighlight")
# 
# BiocManager::install('EnhancedVolcano')
# EnhancedVolcano::EnhancedVolcano(results_data,
#                 lab = rownames(results_data),
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 pCutoff = padj_cutoff,
#                 FCcutoff = log2FC_cutoff)

# saving volcano plot
results_data_annot_noNA$signif_DE <- factor(results_data_annot_noNA$signif_DE,
                                            levels = c("NO", "DOWN", "UP"))
signif_Ntop_volcanoPlot <- ggplot(data = results_data_annot_noNA, aes(x = log2FoldChange, y = -log10(pvalue), col=signif_DE)) +
  geom_point() +
  #gghighlight::gghighlight(signif_DE %in% c("DOWN", "UP")) +
  geom_label_repel(data = . %>% filter(gene_name %in% results_data_annot_signif_Ntop$gene_name), aes(label = gene_name),
                   show.legend = FALSE) +
  #geom_vline(xintercept=c(-log2FC_cutoff, log2FC_cutoff), col="red", linetype="dashed") +
  #geom_hline(yintercept=-log10(padj_cutoff), col="red", linetype="dashed") + # need to adjust to match padj_cutoff
  #scale_color_manual(values=c(DOWN="navy", UP="firebrick3")) +
  scale_color_manual(values=c(DOWN="navy", UP="firebrick3", NO = "grey")) +
  theme_minimal(base_size = 14) +
  labs(x = "log2FC")
       #y = "-log10( p-value )",color = "signif. DE") 

ggsave(plot = signif_Ntop_volcanoPlot,
       filename = paste0(PUBLICATION_FIGURES_DIR, "volcanoPlot_signifDE_ntop", ntop_signif, ".pdf"),
       width = 25, height = 25,
       units = "cm")

# top 50 based on significance?
# [] plot aso library effect!
# [-] rather use rlog norm log2 counts - use log2 tranformed normalized counts;
# highlihgt batch effects - https://support.bioconductor.org/p/68859/
# https://support.bioconductor.org/p/76099/#101057
# https://support.bioconductor.org/p/116821/
sum(results_data_annot_signif_Ntop$log2FoldChange > 0)

# norm counts
signif_normCounts_annot <- as.data.frame(signif_normCounts) %>%
  tibble::rownames_to_column(var="ensembl_id") %>%
  dplyr::filter(ensembl_id %in% results_data_annot_signif_Ntop$ensembl_id) %>%
  dplyr::left_join(., dplyr::select(results_data_annot_signif_Ntop, gene_name, ensembl_id), by = "ensembl_id") %>%
  dplyr::select(-ensembl_id) %>%
  tidyr::pivot_longer(-gene_name, names_to="omics_id", values_to="norm_counts") %>%
  dplyr::left_join(., cond_data, by = "omics_id") %>%
  dplyr::mutate(counts_pc = (norm_counts + 1.0)) %>%
  dplyr::mutate(log10_normCounts = log10(counts_pc),
                log2_normCounts = log2(counts_pc))  # or using pc = 0.5; similar to DESeq2::plotCounts()

# using rlog instead!
signif_rlogCounts_annot <- as.data.frame(signif_rld_counts) %>%
  tibble::rownames_to_column(var="ensembl_id") %>%
  dplyr::filter(ensembl_id %in% results_data_annot_signif_Ntop$ensembl_id) %>%
  dplyr::left_join(., dplyr::select(results_data_annot_signif_Ntop, gene_name, ensembl_id), by = "ensembl_id") %>%
  dplyr::select(-ensembl_id) %>%
  tidyr::pivot_longer(-gene_name, names_to="omics_id", values_to="rlog_counts") %>%
  dplyr::left_join(., cond_data, by = "omics_id") 

signif_rlogCounts_noBatch_annot <- as.data.frame(signif_rld_counts_noBatch) %>%
  tibble::rownames_to_column(var="ensembl_id") %>%
  dplyr::filter(ensembl_id %in% results_data_annot_signif_Ntop$ensembl_id) %>%
  dplyr::left_join(., dplyr::select(results_data_annot_signif_Ntop, gene_name, ensembl_id), by = "ensembl_id") %>%
  dplyr::select(-ensembl_id) %>%
  tidyr::pivot_longer(-gene_name, names_to="omics_id", values_to="rlog_counts") %>%
  dplyr::left_join(., cond_data, by = "omics_id") 


# The y-axis is log scale, but the tick marks are on the count scale
#https://support.bioconductor.org/p/105938/

# ggplot(data = signif_normCounts_annot, aes(x = material, y = counts_pc, colour=material)) + 
#   geom_boxplot() +
#   ggpubr::yscale("log10", .format = FALSE) +
#   facet_wrap(~gene_name, scales = "free_y") 

intgroup_plot <- "material"
# Gene_countsPlot <- signif_normCounts_annot %>%
#   #mutate(norm_counts = (norm_counts + 1.0)) %>%
#   ggpubr::ggviolin(., x = intgroup_plot, y = "counts_pc",
#                     color = intgroup_plot, palette = "jco",
#                     #facet.by = "gene_name",
#                     #title=paste0("Expression of ", analysis_gene),
#                     xlab = intgroup_plot,
#                     ylab = "Normalized counts",
#                     add = "jitter",
#                     outlier.shape = NA,) +  # to remove outlier point if plotting with jitter 
#   ggpubr::yscale("log10", .format = FALSE) + 
#   ggplot2::facet_wrap(~gene_name, scales = "free") #+  # artibtrary unit using free scale
#   #theme(axis.text.x = element_text(angle = 45, hjust = 1))

# the following with free_y creates negative y-axis ticks!?
Gene_countsPlot <- signif_normCounts_annot %>%
  ggpubr::ggviolin(., x = intgroup_plot, y = "log2_normCounts",
                   color = intgroup_plot, palette = "jco",
                   #facet.by = "gene_name",
                   #title=paste0("Expression of ", analysis_gene),
                   xlab = intgroup_plot,
                   ylab = "log2(Normalized counts)",
                   add = "jitter",
                   outlier.shape = NA,) +  # to remove outlier point if plotting with jitter
  ggplot2::facet_wrap(~gene_name, scales = "free_y") #+  # artibtrary unit using free scale theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Gene_countsPlot_libSplit <- signif_normCounts_annot %>%
#   ggpubr::ggviolin(., x = "library_type", y = "log2_normCounts",
#                    color = intgroup_plot, palette = "jco",
#                    #facet.by = "gene_name",
#                    #title=paste0("Expression of ", analysis_gene),
#                    xlab = intgroup_plot,
#                    ylab = "log2(Normalized counts)",
#                    add = "jitter",
#                    outlier.shape = NA,) +  # to remove outlier point if plotting with jitter
#   ggplot2::facet_wrap(~gene_name, scales = "free_y") #+  # artibtrary unit using free scale theme(axis.text.x = element_text(angle = 45, hjust = 1))

Gene_countsPlot_rlog <- signif_rlogCounts_annot %>%
  ggpubr::ggviolin(., x = intgroup_plot, y = "rlog_counts",
                   color = intgroup_plot, palette = "jco",
                   #facet.by = "gene_name",
                   #title=paste0("Expression of ", analysis_gene),
                   xlab = intgroup_plot,
                   ylab = "regularized log2(Normalized counts)",
                   add = "jitter",
                   outlier.shape = NA,) +  # to remove outlier point if plotting with jitter
  ggplot2::facet_wrap(~gene_name, scales = "free_y") #+  # artibtrary unit using free scale theme(axis.text.x = element_text(angle = 45, hjust = 1))

Gene_countsPlot_rlog_noBatch <- signif_rlogCounts_noBatch_annot %>%
  ggpubr::ggviolin(., x = intgroup_plot, y = "rlog_counts",
                   color = intgroup_plot, palette = "jco",
                   #facet.by = "gene_name",
                   #title=paste0("Expression of ", analysis_gene),
                   xlab = intgroup_plot,
                   ylab = "regularized log2(Normalized counts)",
                   add = "jitter",
                   outlier.shape = NA,) +  # to remove outlier point if plotting with jitter
  ggplot2::facet_wrap(~gene_name, scales = "free_y") #+  # artibtrary unit using free scale theme(axis.text.x = element_text(angle = 45, hjust = 1))

# example library_split
signif_rlogCounts_annot_SLPI <- signif_rlogCounts_annot %>%
  dplyr::filter(gene_name == "SLPI") #%>%  dplyr::mutate(library_type = paste0("library type - ", library_type))

Gene_countsPlot_rlog_SLPI <- signif_rlogCounts_annot_SLPI %>%
  ggpubr::ggviolin(., x = intgroup_plot, y = "rlog_counts",
                   color = intgroup_plot, palette = "jco",
                   facet.by = "library_type",
                   panel.labs = list(library_type = c("library type - U", "library type - SR")),
                   add = c("boxplot", "jitter"),
                   title=paste0("Expression of SLPI gene"),
                   xlab = intgroup_plot,
                   ylab = "regularized log2(Normalized counts)",
                   outlier.shape = NA,) 


# plotting for selected genes ----
pub_selected_genes <- c("CCL20","IL1A","IL1B","CCL8","CXCL3","CCL3","IL11","CCL4","IL2RA","CD274","CXCL2","CXCL8","ARG2")
results_data_annot_signif_pub_selected <- results_data_annot_signif %>%
  dplyr::filter(hgnc_symbol != "") %>%  # removing genes without symbol
  dplyr::mutate(gene_name = if_else(hgnc_symbol == "", ensembl_id, hgnc_symbol)) %>% # if hgnc_symbol not present use ensembl_id
  dplyr::filter(gene_name %in% pub_selected_genes) %>%
  dplyr::select(gene_name, ensembl_id, hgnc_symbol, description, log2FoldChange, padj, pvalue)

signif_rlogCounts_annot_pub_selected <- as.data.frame(signif_rld_counts) %>%
  tibble::rownames_to_column(var="ensembl_id") %>%
  dplyr::filter(ensembl_id %in% results_data_annot_signif_pub_selected$ensembl_id) %>%
  dplyr::left_join(., dplyr::select(results_data_annot_signif_pub_selected, gene_name, ensembl_id), by = "ensembl_id") %>%
  dplyr::select(-ensembl_id) %>%
  tidyr::pivot_longer(-gene_name, names_to="omics_id", values_to="rlog_counts") %>%
  dplyr::left_join(., cond_data, by = "omics_id") 

sum(!(pub_selected_genes %in% signif_rlogCounts_annot_pub_selected$gene_name))
pub_selected_genes[!(pub_selected_genes %in% signif_rlogCounts_annot_pub_selected$gene_name)]

# highlight library_type
Gene_countsPlot_rlog_pub_selected <- signif_rlogCounts_annot_pub_selected %>%
  ggpubr::ggviolin(., x = intgroup_plot, y = "rlog_counts",
                   fill = intgroup_plot, 
                   #color = "library_type",
                   palette = "jco",
                   #facet.by = "gene_name",
                   #title=paste0("Expression of ", analysis_gene),
                   xlab = intgroup_plot,
                   ylab = "regularized log2(Normalized counts)",
                   # add = c("jitter"),
                   # add.params = list(color = "library_type"),
                   add = c("boxplot", "jitter"),
                   add.params = list(fill = "white"),
                   #facet.by = "gene_name",
                   outlier.shape=NA) +  # to remove outlier point if plotting with jitter
  ggplot2::facet_wrap(vars(gene_name), scales = "free_y") + #+  # artibtrary unit using free scale theme(axis.text.x = element_text(angle = 45, hjust = 1))
  scale_colour_manual(name = "library_type", values = c("red", "black"))

Gene_countsPlot_rlog_pub_selected <- signif_rlogCounts_annot_pub_selected %>%
  ggpubr::ggviolin(., x = intgroup_plot, y = "rlog_counts",
                   fill = intgroup_plot, 
                   #color = "library_type",
                   palette = "jco",
                   #facet.by = "gene_name",
                   #title=paste0("Expression of ", analysis_gene),
                   xlab = intgroup_plot,
                   ylab = "regularized log2(Normalized counts)",
                   #add = c("jitter"),
                   #add.params = list(color = "library_type"),
                   #add = c("boxplot", "dotplot"),
                   #add.params = list(fill = "white"),
                   #facet.by = "gene_name",
                   outlier.shape=NA) +  # to remove outlier point if plotting with jitter
  ggplot2::facet_wrap(vars(gene_name), scales = "free_y") + #+  # artibtrary unit using free scale theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggplot2::geom_boxplot(width=0.2, fill="white", outlier.shape = NA) +
  ggplot2::geom_jitter(aes(color=library_type), size=1, width=0.1, alpha=0.8) +
  scale_colour_manual(name = "library_type",values = c("red", "black"))

# updating deconvolution plots with the above approach and coloring? 




# plot batch corrected!? - to sort of show values used in DE analysis!?  
# Gene_countsPlot_rlog_pub_selected <- signif_rlogCounts_annot_pub_selected %>%
#   ggpubr::ggviolin(., x = intgroup_plot, y = "rlog_counts",
#                    fill = intgroup_plot, 
#                    #color = "library_type",
#                    palette = "jco",
#                    #facet.by = "gene_name",
#                    #title=paste0("Expression of ", analysis_gene),
#                    xlab = intgroup_plot,
#                    ylab = "regularized log2(Normalized counts)",
#                    add = c("jitter"),
#                    add.params = list(color = "library_type"),
#                    #add = c("boxplot", "dotplot"),
#                    #add.params = list(fill = "white"),
#                    #facet.by = "gene_name",
#                    outlier.shape=NA) +  # to remove outlier point if plotting with jitter
#   ggplot2::facet_wrap(vars(gene_name), scales = "free_y") + #+  # artibtrary unit using free scale theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   scale_colour_manual(name = "library_type",values = c("red", "black"))

# generate one plot per gene
signif_rlogCounts_annot_pub_selected_perGene <- base::split(signif_rlogCounts_annot_pub_selected, f = signif_rlogCounts_annot_pub_selected$gene_name)

signif_rlogCounts_annot_pub_selected_perGene_vlnPlot <- purrr::map2(.x = signif_rlogCounts_annot_pub_selected_perGene,
                                                                    .y = names(signif_rlogCounts_annot_pub_selected_perGene),
                                                                    .f = function(gene_exp, gene_name){
                                                                       ggpubr::ggviolin(gene_exp, x = intgroup_plot, y = "rlog_counts",
                                                                                        fill = intgroup_plot, 
                                                                                        palette = "jco",
                                                                                        facet.by = "library_type",
                                                                                        panel.labs = list(library_type = c(U="library type - U", SR="library type - SR")),
                                                                                        add = c("jitter"),
                                                                                        #add = c("boxplot", "jitter"),
                                                                                        #add.params = list(fill = "white"),
                                                                                        title=paste0("Expression of ", gene_name),
                                                                                        xlab = intgroup_plot,
                                                                                        ylab = "regularized log2(Normalized counts)",
                                                                                        outlier.shape = NA) #+ ggplot2::facet_wrap(vars(library_type), scales = "fixed")
                                                                   })

signif_rlogCounts_annot_pub_selected_perGene_vlnPlot2 <- purrr::map2(.x = signif_rlogCounts_annot_pub_selected_perGene,
                                                                    .y = names(signif_rlogCounts_annot_pub_selected_perGene),
                                                                    .f = function(gene_exp, gene_name){
                                                                      ggpubr::ggviolin(gene_exp, x = "library_type", y = "rlog_counts",
                                                                                       fill = "library_type", 
                                                                                       #palette = "jco",
                                                                                       facet.by = intgroup_plot,
                                                                                       panel.labs = list(material = c(no="infiltration - no",yes="infiltration - yes")),
                                                                                       add = c("jitter"),
                                                                                       #add = c("boxplot", "jitter"),
                                                                                       #add.params = list(fill = "white"),
                                                                                       title=paste0("Expression of ", gene_name),
                                                                                       xlab = "library_type",
                                                                                       ylab = "regularized log2(Normalized counts)",
                                                                                       outlier.shape = NA) #+ ggplot2::facet_wrap(vars(library_type), scales = "fixed")
                                                                    })

devtools::source_gist("2a1bb0133ff568cbe28d", filename = "geom_flat_violin.R")
signif_rlogCounts_annot_pub_selected_perGene_cloudPlot <- purrr::map2(.x = signif_rlogCounts_annot_pub_selected_perGene,
                                                                    .y = names(signif_rlogCounts_annot_pub_selected_perGene),
                                                                    .f = function(gene_exp, gene_name){
                                                                      #devtools::source_gist("2a1bb0133ff568cbe28d", filename = "geom_flat_violin.R")
                                                                      
                                                                      pos <- position_jitter(width = 0.15, seed = 1)
                                                                      
                                                                      p0 <- ggplot(data = gene_exp, 
                                                                                   aes(x = material, 
                                                                                       y = rlog_counts, 
                                                                                       fill = material)) +
                                                                        geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
                                                                        labs(title=paste0("Expression of ", gene_name),
                                                                             xlab = intgroup_plot,
                                                                             ylab = "regularized log2(Normalized counts)") +
                                                                        guides(fill = FALSE) +
                                                                        guides(color = FALSE) +
                                                                        scale_color_brewer(palette = "Dark2") +
                                                                        scale_fill_brewer(palette = "Dark2") +
                                                                        theme_classic() 
                                                                      
                                                                      # raincloud plot
                                                                      p1 <- p0 + 
                                                                        geom_point(aes(color = material), 
                                                                                   position = pos, size = 3, alpha = 0.8) +
                                                                        geom_boxplot(width = .1, show.legend = FALSE, outlier.shape = NA, alpha = 0.5) +
                                                                        ggplot2::facet_wrap(vars(library_type), scales = "fixed")
                                                                    })


Gene_countsPlot_rlog_pub_selected_ARG2 <- signif_rlogCounts_annot_pub_selected_perGene$ARG2 %>%
  ggpubr::ggviolin(., x = intgroup_plot, y = "rlog_counts",
                   fill = intgroup_plot, 
                   #color = "library_type",
                   palette = "jco",
                   #facet.by = "gene_name",
                   #title=paste0("Expression of ", analysis_gene),
                   xlab = intgroup_plot,
                   ylab = "regularized log2(Normalized counts)",
                   add = c("jitter"),
                   add.params = list(color = "library_type"),
                   #add = c("boxplot", "dotplot"),
                   #add.params = list(fill = "white"),
                   #facet.by = "gene_name",
                   outlier.shape=NA) +  # to remove outlier point if plotting with jitter
  ggplot2::facet_wrap(vars(gene_name), scales = "free_y") + #+  # artibtrary unit using free scale theme(axis.text.x = element_text(angle = 45, hjust = 1))
  scale_colour_manual(name = "library_type",values = c("red", "black"))

ggarrange(Gene_countsPlot_rlog_pub_selected_ARG2, signif_rlogCounts_annot_pub_selected_perGene_vlnPlot$ARG2)

# raincloud plot
devtools::source_gist("2a1bb0133ff568cbe28d", filename = "geom_flat_violin.R")

pos <- position_jitter(width = 0.15, seed = 1)

p0 <- ggplot(data = signif_rlogCounts_annot_pub_selected, aes(x = material, y = rlog_counts, fill = material)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  ggplot2::facet_grid(~gene_name, scales = "free_y") #+
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() 

# raincloud plot
p1 <- p0 + 
  geom_point(aes(color = material), 
             position = pos, size = 3, alpha = 0.8) +
  geom_boxplot(width = .1, show.legend = FALSE, outlier.shape = NA, alpha = 0.5) 
p1


# Gene_countsPlot_rlog_pub_selected2 <- ggpubr::ggviolin(signif_rlogCounts_annot_pub_selected, x = intgroup_plot, y = "rlog_counts",
#                                                        fill = intgroup_plot,
#                                                        ylab = "regularized log2(Normalized counts)",
#                                                        font.label = list(size = 12, color = "black", face = "bold"),
#                                                        panel.labs.font = list(size = 12, color = "black", face = "bold"),
#                                                        add = c("boxplot", "jitter"),
#                                                        add.params = list(size = 0.2, fill = "white"),
#                                                        trim = TRUE,  # If TRUE (default), trim the tails of the violins to the range of the data. If FALSE, don't trim the tails.
#                                                        facet.by = "gene_name")


# intGene_vln1 <- ggpubr::ggviolin(norm_log2pcounts_interestingGene, x = "orig.ident", y = "log2_normCounts",
#                                  fill = "orig.ident",
#                                  title = selected_genes,
#                                  ylab = "log2(Norm. counts)",
#                                  font.label = list(size = 12, color = "black", face = "bold"),
#                                  panel.labs.font = list(size = 12, color = "black", face = "bold"),
#                                  add = c("boxplot", "jitter"),
#                                  add.params = list(size = 0.2, fill = "white"),
#                                  #legend title = "sample_id",
#                                  trim = TRUE,  # If TRUE (default), trim the tails of the violins to the range of the data. If FALSE, don't trim the tails.
#                                  facet.by = "clust_celltype") #+ ylim(-0.1,8) +  # to start from 0; not needed if using trim=TRUE
#guides(fill=guide_legend(title="sample_id", title.position = "right", label.position = "right")) +
# theme(#text = element_text(size=8),
#   plot.title = element_text(size=16, color = "black", face = "bold", hjust = 0.5),   #hjust = 0.5
#   #axis.title.y =  element_blank(), 
#   axis.title.x = element_blank(),
#   axis.text.x = element_blank(), 
#   axis.text.y = element_text(face = "bold")) #axis.text.y = element_text(face = "bold", size = 6)

# intGene_vln1 <- ggpar(intGene_vln1,
#                       legend = "right", legend.title = "sample_id") + 
#   theme(#text = element_text(size=8),
#     plot.title = element_text(size=16, color = "black", face = "bold", hjust = 0.5),   #hjust = 0.5
#     #axis.title.y =  element_blank(), 
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(), 
#     axis.text.y = element_text(face = "bold")) #axis.text.y = element_text(face = "bold", size = 6)


# ggsave(plot = Gene_countsPlot,
#        filename = paste0(PUBLICATION_FIGURES_DIR, "gene_normCounts_ntop", ntop_signif, ".pdf"),
#        width = 30, height = 30,
#        units = "cm")

ggsave(plot = Gene_countsPlot_rlog,
       filename = paste0(PUBLICATION_FIGURES_DIR, "gene_rlogCounts_ntop", ntop_signif, ".pdf"),
       width = 30, height = 30,
       units = "cm")

ggsave(plot = Gene_countsPlot_rlog_noBatch,
       filename = paste0(PUBLICATION_FIGURES_DIR, "gene_rlogCounts_noBatch_ntop", ntop_signif, ".pdf"),
       width = 30, height = 30,
       units = "cm")

ggsave(plot = Gene_countsPlot_rlog_SLPI,
       filename = paste0(PUBLICATION_FIGURES_DIR, "gene_rlogCounts_SLPI_gene.pdf"),
       width = 30, height = 30,
       units = "cm")

# saving selected genes for publication ----
ggsave(plot = Gene_countsPlot_rlog_pub_selected,
       filename = paste0(PUBLICATION_FIGURES_DIR, "gene_rlogCounts_pubGenes.pdf"),
       width = 20, height = 20,
       units = "cm")


for (gene_name in names(signif_rlogCounts_annot_pub_selected_perGene_vlnPlot)) {
  ggsave(plot = signif_rlogCounts_annot_pub_selected_perGene_vlnPlot[[gene_name]],
         filename = paste0(PUBLICATION_FIGURES_DIR, "gene_rlogCounts_pubGenes_",gene_name,".pdf"),
         width = 15, height = 15,
         units = "cm")
}

for (gene_name in names(signif_rlogCounts_annot_pub_selected_perGene_vlnPlot2)) {
  ggsave(plot = signif_rlogCounts_annot_pub_selected_perGene_vlnPlot2[[gene_name]],
         filename = paste0(PUBLICATION_FIGURES_DIR, "gene_rlogCounts_pubGenes_",gene_name,"_libType.pdf"),
         width = 15, height = 15,
         units = "cm")
}

# variance partition for top N genes ----
# Access first entries
#varPart_fit <- variancePartition::fitExtractVarPartModel(geneExpr_partVariance, fitform_partVariance, info_partVariance)
#stopCluster(cl)

# can be computationally intensive
library('doParallel')
cl <- makeCluster(40)
registerDoParallel(cl)

partVariance_design <- ~ (1|library_type) + (1|material)

signif_rlogCounts_geneName_matrix <- as.data.frame(signif_rld_counts) %>%
  tibble::rownames_to_column(var="ensembl_id") %>%
  dplyr::filter(ensembl_id %in% results_data_annot_signif_Ntop$ensembl_id) %>%
  dplyr::left_join(., dplyr::select(results_data_annot_signif_Ntop, gene_name, ensembl_id), by = "ensembl_id") %>%
  dplyr::select(-ensembl_id) %>%
  tibble::column_to_rownames(var = "gene_name") %>%
  as.matrix(.)

info_partVariance <- as.data.frame(colData(rld_filt))

varPart_signif_rlog <- variancePartition::fitExtractVarPartModel(signif_rlogCounts_geneName_matrix, 
                                                                 partVariance_design, 
                                                                 info_partVariance)
stopCluster(cl)

# sort variables (i.e. columns) by median fraction
# of variance explained
varPart_signif_rlog_sorted <- variancePartition::sortCols( varPart_signif_rlog ) #showMethods("sortCols")

head(varPart_signif_rlog_sorted)
varPart_signif_rlog_sorted[1:10,]
variancePartition::plotPercentBars( varPart_signif_rlog_sorted[1:10,] )
varPart_signif_rlog_barPlot <- variancePartition::plotPercentBars( varPart_signif_rlog_sorted)

ggsave(plot = varPart_signif_rlog_barPlot,
       filename = paste0(PUBLICATION_FIGURES_DIR, "gene_rlogCounts_ntop", ntop_signif, "_varPartition.pdf"),
       width = 30, height = 30,
       units = "cm")


# GSEA plots ----
# top 20 upregulated: KEGG, GOBP, Hallmarks
# preselect signif. and up for each and plot with selected
# plot 'classical' GSEA plots
ntop_signif_gsea <- 20
gsea_df_h <- gsea_results_df$hallmark
gsea_df_kegg <- gsea_results_df$C2_kegg
gsea_df_gobp <- gsea_results_df$C5_GOBP

gsea_results_list <- list(gsea_hallmark = gsea_df_h,
                          gsea_kegg = gsea_df_kegg,
                          gsea_gobp = gsea_df_gobp)

gsea_results_signif_UP_ntop <- purrr::map(.x = gsea_results_list, .f = function(gsea_df) {
  gsea_df %>%
    dplyr::filter(p.adjust < padj_cutoff & NES > 0) %>%  # significant and upregulated
    dplyr::arrange(desc(NES)) %>%
    dplyr::slice_head(n = ntop_signif_gsea) 
})

gsea_results_plotReady <- purrr::map(.x = gsea_results_signif_UP_ntop, .f = function(gsea_df) {
  # see https://github.com/YuLab-SMU/DOSE/issues/20
  # see https://github.com/YuLab-SMU/enrichplot/blob/497b8ba286b9ab8a7664eb81ae467376d0166d23/R/method-fortify.R for gene counting
  # res$Count <- str_count(res$core_enrichment, "/")
  # res$GeneRatio <- res$Count / res$setSize
  # test_A <- gsea_df %>%
  #   tidyr::separate_rows(., core_enrichment, sep = "/") %>%
  #   dplyr::group_by(ID) %>%
  #   dplyr::summarise(Count = n()) %>%
  #   dplyr::arrange(desc(Count))
  gsea_df_GeneCount <- gsea_df %>% group_by(ID) %>% 
    summarise(count = sum(stringr::str_count(core_enrichment, "/")) + 1) # + 1 because counting / not genes!
  dplyr::left_join(gsea_df, gsea_df_GeneCount, by = "ID") %>% 
    dplyr::mutate(GeneRatio = count/setSize) %>%
    dplyr::mutate(log10_padj = -log10(p.adjust)) %>%
    dplyr::mutate(Status=ifelse(NES > 0, "Upregulated", "Downregulated")) %>% # adding status
    dplyr::mutate(Status=factor(Status, levels = c("Upregulated", "Downregulated"))) # converting to factor and changing levels
})

# test_plot ----
# For the three other plots, it would for us make more sense to have:
# *	NES on the y-axis - DONE; but bettern on x-axis!
# *	Circles instead of triangles - DONE
# *	Circle size representing the gene ratio - DONE
# *	Circle color representing the adjusted p-value

# N_variables = 5 ; NES, padj, Count, GeneRatio, Description
# plot_options = 4-6; x-axis, y-axis, colour, size, (shape, label)

# gsea_df <- gsea_results_plotReady$gsea_hallmark
# padj_limit <- 0.1
# min_log_padj = -log10(padj_limit) #-log10(padj_cutoff); padb = 0.1
# max_ceiling_log_padj <- ceiling(max(-log10(gsea_df$p.adjust)))
# max_round_log_padj <- round(max(-log10(gsea_df$p.adjust)))

# ggplot(data = gsea_df, 
#        mapping=aes(x = forcats::fct_reorder(Description, NES), http://localhost:8888/graphics/plot_zoom_png?width=1200&height=900
#                    y = NES,
#                    #shape = Status,
#                    color = -log10(p.adjust),
#                    #fill = NES,
#                    size = GeneRatio)) +
#  geom_point()
# rm(gsea_df)
# END test plot ----

# -log10(padj)
gsea_results_dotPlots_v1 <- purrr::map(.x = gsea_results_plotReady, .f = function(gsea_df) {
  padj_limit <- 0.1
  min_log_padj = -log10(padj_limit) #-log10(padj_cutoff); padb = 0.1
  max_ceiling_log_padj <- ceiling(max(-log10(gsea_df$p.adjust)))
  max_round_log_padj <- round(max(-log10(gsea_df$p.adjust)))
  
  # plotting NES, -log10(padj), GeneRatio
  ggplot(data = gsea_df, 
         mapping=aes(x = NES, 
                     y = forcats::fct_reorder(Description, NES),
                     color = -log10(p.adjust),
                     size = GeneRatio)) +
    geom_point() + 
    scale_colour_gradient(low = "yellow", high = "red", na.value = NA) +
    theme_bw() +
    ggtitle(paste0("GSEA (p.adj < ", padj_cutoff,")")) +
    theme(text = element_text(size=12),
          plot.title = element_text(hjust = 0.5), 
          axis.title.y =  element_blank(), 
          axis.text.y = element_text(face = "bold")) 
})

# count - discrete, geneRatio
#gsea_results_dotPlots_v2 <- purrr::map(.x = gsea_results_plotReady, .f = function(gsea_df) {
#   padj_limit <- 0.1
#   min_log_padj = -log10(padj_limit) #-log10(padj_cutoff); padb = 0.1
#   max_ceiling_log_padj <- ceiling(max(-log10(gsea_df$p.adjust)))
#   max_round_log_padj <- round(max(-log10(gsea_df$p.adjust)))
#   
#   # plotting NES, Count, GeneRatio
#   ggplot(data = gsea_df, 
#          mapping=aes(x = NES, 
#                      y = forcats::fct_reorder(Description, NES),
#                      color = GeneRatio,
#                      size = count)) +
#     geom_point() + 
#     #scale_shape_manual(values = c(24, 25), guide=guide_legend(order=1)) + #24 - up triangle, 25 - down triangle
#     scale_colour_gradient(low = "yellow", high = "red", na.value = NA) +
#     #scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
#     # scale_size(limits = c(min_log_padj, max_ceiling_log_padj), 
#     #            breaks = c(1.3, seq(2.0, 5, 1)),
#     #            labels = c("1.3", "2.0", "3.0", "4.0", ">= 5.0"),
#     #            guide=guide_legend(override.aes = list(shape=17), order=3)) +
#     theme_bw() +
#     ggtitle(paste0("GSEA (p.adj < ", padj_cutoff,")")) +
#     theme(text = element_text(size=12),
#           plot.title = element_text(hjust = 0.5), 
#           axis.title.y =  element_blank(), 
#           axis.text.y = element_text(face = "bold")) 
# })

# count, geneRatio - discrete
gsea_results_dotPlots_v3 <- purrr::map(.x = gsea_results_plotReady, .f = function(gsea_df) {
  padj_limit <- 0.1
  min_log_padj = -log10(padj_limit) #-log10(padj_cutoff); padb = 0.1
  max_ceiling_log_padj <- ceiling(max(-log10(gsea_df$p.adjust)))
  max_round_log_padj <- round(max(-log10(gsea_df$p.adjust)))
  
  # plotting NES, Count, GeneRatio
  ggplot(data = gsea_df, 
         mapping=aes(x = NES, 
                     y = forcats::fct_reorder(Description, NES),
                     color = count,
                     size = GeneRatio)) +
    geom_point() + 
    #scale_shape_manual(values = c(24, 25), guide=guide_legend(order=1)) + #24 - up triangle, 25 - down triangle
    scale_colour_gradient(low = "yellow", high = "red", na.value = NA) +
    #scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
    # scale_size(limits = c(min_log_padj, max_ceiling_log_padj), 
    #            breaks = c(1.3, seq(2.0, 5, 1)),
    #            labels = c("1.3", "2.0", "3.0", "4.0", ">= 5.0"),
    #            guide=guide_legend(override.aes = list(shape=17), order=3)) +
    theme_bw() +
    ggtitle(paste0("GSEA (p.adj < ", padj_cutoff,")")) +
    theme(text = element_text(size=12),
          plot.title = element_text(hjust = 0.5), 
          axis.title.y =  element_blank(), 
          axis.text.y = element_text(face = "bold")) 
})

# count, GeneRatio - discrete
#gsea_results_dotPlots_v4 <- purrr::map(.x = gsea_results_plotReady, .f = function(gsea_df) {
#   padj_limit <- 0.1
#   min_log_padj = -log10(padj_limit) #-log10(padj_cutoff); padb = 0.1
#   max_ceiling_log_padj <- ceiling(max(-log10(gsea_df$p.adjust)))
#   max_round_log_padj <- round(max(-log10(gsea_df$p.adjust)))
#   
#   # plotting NES, -log10(padj), Count, GeneRatio
#   ggplot(data = gsea_df, 
#          mapping=aes(x = NES, 
#                      y = -log10(p.adjust),
#                      color = count,
#                      size = GeneRatio,
#                      label=Description)) +
#     geom_point() + 
#     ggrepel::geom_label_repel(show.legend = FALSE) +
#     scale_colour_gradient(low = "yellow", high = "red", na.value = NA) +
#     theme_bw() +
#     ggtitle(paste0("GSEA (p.adj < ", padj_cutoff,")")) +
#     theme(plot.title = element_text(hjust = 0.5))
# })

# count - discrete, GeneRatio
# gsea_results_dotPlots_v5 <- purrr::map(.x = gsea_results_plotReady, .f = function(gsea_df) {
#   padj_limit <- 0.1
#   min_log_padj = -log10(padj_limit) #-log10(padj_cutoff); padb = 0.1
#   max_ceiling_log_padj <- ceiling(max(-log10(gsea_df$p.adjust)))
#   max_round_log_padj <- round(max(-log10(gsea_df$p.adjust)))
#   
#   # plotting NES, -log10(padj), Count, GeneRatio
#   ggplot(data = gsea_df, 
#          mapping=aes(x = NES, 
#                      y = -log10(p.adjust),
#                      color = GeneRatio,
#                      size = count,
#                      label=Description)) +
#     geom_point() + 
#     ggrepel::geom_label_repel(show.legend = FALSE) +
#     scale_colour_gradient(low = "yellow", high = "red", na.value = NA) +
#     theme_bw() +
#     ggtitle(paste0("GSEA (p.adj < ", padj_cutoff,")")) +
#     theme(plot.title = element_text(hjust = 0.5))
# })

# triangle plots; GeneRation on x-axis
# gsea_results_trianglePlots <- purrr::map(.x = gsea_results_plotReady, .f = function(gsea_df) {
#   padj_limit <- 0.1
#   min_log_padj = -log10(padj_limit) #-log10(padj_cutoff); padb = 0.1
#   max_ceiling_log_padj <- ceiling(max(-log10(gsea_df$p.adjust)))
#   max_round_log_padj <- round(max(-log10(gsea_df$p.adjust)))
# 
#   ggplot(data = gsea_df,
#          mapping=aes(x = GeneRatio,
#                      y = forcats::fct_reorder(Description, GeneRatio),
#                      shape = Status,
#                      fill = NES,
#                      size = -log10(p.adjust))) +
#     geom_point() +
#     scale_shape_manual(values = c(24, 25), guide=guide_legend(order=1)) + #24 - up triangle, 25 - down triangle; without this there is no fill!; use also color=NES!?
#     scale_colour_gradient(low = "yellow", high = "red", na.value = NA) +
#     scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
#     scale_size(limits = c(min_log_padj, max_ceiling_log_padj),
#                breaks = c(1.3, seq(2.0, 5, 1)),
#                labels = c("1.3", "2.0", "3.0", "4.0", ">= 5.0"),
#                guide=guide_legend(override.aes = list(shape=17), order=3)) +
#     theme_bw() +
#     ggtitle(paste0("GSEA (p.adj < ", padj_cutoff,")")) +
#     theme(text = element_text(size=12),
#           plot.title = element_text(hjust = 0.5),
#           axis.title.y =  element_blank(),
#           axis.text.y = element_text(face = "bold"))
# })

for (gsea_results_name in names(gsea_results_dotPlots_v1)) {
  # triangle plots
  # ggsave(gsea_results_trianglePlots[[gsea_results_name]], 
  #        filename = paste0(PUBLICATION_FIGURES_DIR, gsea_results_name, "_upregulated_ntop", ntop_signif, "_trianglePlot.pdf"),
  #        height = 20, 
  #        width = 25, 
  #        units = "cm")
  
  ggsave(gsea_results_dotPlots_v1[[gsea_results_name]], 
         filename = paste0(PUBLICATION_FIGURES_DIR, gsea_results_name, "_upregulated_ntop", ntop_signif, "_dotPlot_v1.pdf"),
         height = 20, 
         width = 25, 
         units = "cm")
  
  # ggsave(gsea_results_dotPlots_v2[[gsea_results_name]], 
  #        filename = paste0(PUBLICATION_FIGURES_DIR, gsea_results_name, "_upregulated_ntop", ntop_signif, "_dotPlot_v2.pdf"),
  #        height = 20, 
  #        width = 25, 
  #        units = "cm")
  
  ggsave(gsea_results_dotPlots_v3[[gsea_results_name]], 
         filename = paste0(PUBLICATION_FIGURES_DIR, gsea_results_name, "_upregulated_ntop", ntop_signif, "_dotPlot_v3.pdf"),
         height = 20, 
         width = 25, 
         units = "cm")
}

}


# additional tests ----

# # testing means
# d_test <- plotCounts(dds_filt_extractResults, gene="ENSG00000148826", intgroup=c("library_type", "material"), 
#                 returnData=TRUE, transform = FALSE) %>%
#   dplyr::mutate(log2_norm_counts = log2(count + 0.5))
# 
# d_test_summary <- d_test %>%
#   dplyr::group_by(material, library_type) %>%
#   dplyr::summarise(mean_exp = mean(count))
# 
# ggpubr::ggviolin(d_test, x = "material", y = "log2_norm_counts", facet.by = "library_type", add = c("boxplot", "jitter"))
# ggplot(d_test, aes(x=material, y=count)) + 
#   geom_point(position=position_jitter(w=0.1,h=0)) + 
#   scale_y_log10(breaks=c(25,100,400))
# 
# # add material::library_type interaction term
# dds_filt_extractResults_interaction <- dds_filt_extractResults
# dds_filt_extractResults_interaction$group <- factor(paste0(dds_filt_extractResults_interaction$material, 
#                                                            dds_filt_extractResults_interaction$library_type))
# design(dds_filt_extractResults_interaction) <- ~ group
# dds_filt_extractResults_interaction <- DESeq(dds_filt_extractResults_interaction)
# resultsNames(dds_filt_extractResults_interaction)
# results(dds_filt_extractResults_interaction, contrast=c("group", "yesSR", "noSR"))
# results(dds_filt_extractResults_interaction, contrast=c("group", "yesU", "noSR"))
# design(dds_filt_extractResults_interaction) <- ~library_type + material + library_type:material
# dds_filt_extractResults_interaction <- DESeq(dds_filt_extractResults_interaction)
# resultsNames(dds_filt_extractResults_interaction)
# results(dds_filt_extractResults_interaction, contrast=c("group", "yesSR", "noSR"))
# results(dds_filt_extractResults_interaction, contrast=c("group", "yesU", "noSR"))


# test shrinkage
# res_table_unshrunken_test <- DESeq2::results(dds_filt_extractResults, name="material_yes_vs_no", alpha = padj_cutoff)
# res_table_test1 <- DESeq2::lfcShrink(dds_filt_extractResults, coef="material_yes_vs_no", res=res_table_unshrunken_test, type = "normal")
# res_table_test2 <- DESeq2::lfcShrink(dds_filt_extractResults, coef="material_yes_vs_no", res=res_table_unshrunken_test, type = "apeglm")
# summary(res_table_unshrunken_test)
# summary(res_table_test1)
# summary(res_table_test2)
# boxplot(log10(assays(dds_filt_extractResults)[["cooks"]]), range=0, las=2)  # two samples (no; lib SR are outliers on Cook's plot!)
# plotDispEsts(dds_filt_extractResults)
# metadata(res_table_unshrunken_test)$alpha
# metadata(res_table_unshrunken_test)$filterThreshold
# plot(metadata(res_table_unshrunken_test)$filterNumRej, 
#      type="b", ylab="number of rejections",
#      xlab="quantiles of filter")
# lines(metadata(res_table_unshrunken_test)$lo.fit, col="red")
# abline(v=metadata(res_table_unshrunken_test)$filterTheta)
# 
# resNoFilt <- results(dds_filt_extractResults, independentFiltering=FALSE)
# addmargins(table(filtering=(res_table_unshrunken_test$padj < .1),
#                  noFiltering=(resNoFilt$padj < .1)))
# 
# DESeq2::plotMA(res_table_unshrunken_test, ylim=c(-3,3), cex=.8)
# abline(h=c(-1,1), col="dodgerblue", lwd=2)
# 
# DESeq2::plotMA(res_table_test1, ylim=c(-3,3), cex=.8)
# abline(h=c(-1,1), col="dodgerblue", lwd=2)
# 
# DESeq2::plotMA(res_table_test2, ylim=c(-3,3), cex=.8)
# abline(h=c(-1,1), col="dodgerblue", lwd=2)
# 
# # https://support.bioconductor.org/p/77461/#107693
# # https://support.bioconductor.org/p/107680/
# # https://www.biostars.org/p/370268/
# # https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html
# # https://support.bioconductor.org/p/98833/
# # res_table_unshrunken_test3 <- results(dds_filt_extractResults, addMLE=TRUE)
# # DESeq2::plotMA(res_table_unshrunken_test3)
# # DESeq2::plotMA(res_table_unshrunken_test3, MLE=TRUE)
# 
# # pvalue and padj are identical!
# # differences only in
# test0 <- as.data.frame(res_table_unshrunken_test) %>% tibble::rownames_to_column(var="ensembl_id")
# test1 <- as.data.frame(res_table_test1) %>% tibble::rownames_to_column(var="ensembl_id")
# test2 <- as.data.frame(res_table_test2) %>% tibble::rownames_to_column(var="ensembl_id")
# 
# # because we are interested in treated vs untreated, we set 'coef=2'
# resNorm <- lfcShrink(dds, coef=2, type="normal")
# resAsh <- lfcShrink(dds, coef=2, type="ashr")
# 
# par(mfrow=c(1,3), mar=c(4,4,2,1))
# xlim <- c(1,1e5); ylim <- c(-3,3)
# plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
# plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
# plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")



# For selected genes plot batch effects for counts
# plot heatmaps!
# do VarPart to check batch effect

# check modelling of batch effect ----
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# plotDispEsts(dds_filt_extractResults)
# par(mar=c(8,5,2,2))
# boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
# 
# par(mfrow=c(2,2),mar=c(2,2,1,1))
# ylim <- c(-2.5,2.5)
# resGA <- results(dds_filt_extractResults, lfcThreshold=.5, altHypothesis="greaterAbs")
# resLA <- results(dds_filt_extractResults, lfcThreshold=.5, altHypothesis="lessAbs")
# resG <- results(dds_filt_extractResults, lfcThreshold=.5, altHypothesis="greater")
# resL <- results(dds_filt_extractResults, lfcThreshold=.5, altHypothesis="less")
# drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
# DESeq2::plotMA(resGA, ylim=ylim); drawLines()
# DESeq2::plotMA(resLA, ylim=ylim); drawLines()
# DESeq2::plotMA(resG, ylim=ylim); drawLines()
# DESeq2::plotMA(resL, ylim=ylim); drawLines()
# 
# gene_interest <- "ENSG00000166923"
# check_gene <- plotCounts(dds_filt_extractResults, gene = gene_interest, intgroup=c("library_type","material"))
# check_gene <- plotCounts(dds_filt_extractResults, gene = gene_interest, intgroup=c("gender","material"), 
#                 returnData=TRUE)
# library("ggplot2")
# ggplot(check_gene, aes(x=material, y=count)) + 
#   geom_point(position=position_jitter(w=0.1,h=0)) + 
#   scale_y_log10(breaks=c(25,100,400)) +
#   facet_wrap(~library_type)
# 
# ggplot(check_gene, aes(x=material, y=count, fill = material)) + 
#   geom_violin() + 
#   geom_jitter() + 
#   scale_y_log10() +
#   facet_wrap(~gender) +
#   theme_bw() +
#   labs(title = gene_interest)
# 
# # Access first entries
# #varPart_fit <- variancePartition::fitExtractVarPartModel(geneExpr_partVariance, fitform_partVariance, info_partVariance)
# #stopCluster(cl)
# 
# # can be computationally intensive
# library('doParallel')
# cl <- makeCluster(40)
# registerDoParallel(cl)
# 
# fitform_partVariance <- ~ (1|library_type) + (1|material) + (1|material)
# 
# geneExpr_partVarianceAll <- rnaSelectTopVarGenes(rld_counts_filt, ntop = nrow(rld_counts_filt)) #using top 500 most variable genes, alternatively using all genes rld_counts
# varPart_fitAll <- fitExtractVarPartModel(geneExpr_partVarianceAll, fitform_partVariance, info_partVariance)
# stopCluster(cl)
# 
# # sort variables (i.e. columns) by median fraction
# # of variance explained
# varPart_fitAll_sorted <- sortCols( varPart_fitAll ) #showMethods("sortCols")
# 
# head(varPart_fitAll_sorted)
# varPart_fitAll_sorted[1:10,]
# plotPercentBars( varPart_fitAll_sorted[1:10,] )
# gene_interest <- "ENSG00000109265"
# 
# # top5 from
# top5_signif <- list(SPTA1 = "ENSG00000163554",
#                     JCHAIN = "ENSG00000132465",
#                     CXCL11 = "ENSG00000169248",
#                     CENPF = "ENSG00000117724",
#                     CPXM1 = "ENSG00000088882")
# 
# # check my function as there are few things not actually done by this function!
# extract_dds_genes <- function(gene){
#   # extract genes from dds object
#   DESeq2::plotCounts(dds_filt_extractResults, gene = gene, intgroup=c("library_type","material"), 
#                      normalized = TRUE,
#                            returnData=TRUE) %>%
#     dplyr::mutate(ensembl_id = gene)
# }
# genes_interest_df <- purrr::map_dfr(.x = top5_signif, .f = extract_dds_genes) # return data.frame by rowbinding
# 
# gene_interest_plot <- ggplot(genes_interest_df, aes(x=material, y=count, fill = library_type)) + 
#   geom_violin(position = "dodge") + 
#   #geom_jitter() + 
#   scale_y_log10() +
#   facet_wrap(~ensembl_id, scales = "free") +
#   theme_bw() 
# 
# gene_interest_VarPart_plot <- plotPercentBars( varPart_fitAll_sorted[unlist(top5_signif),] )
# 
# ggsave(plot = gene_interest_plot, filename = paste0(OUTPUT_DIR, "example_gene_interest_plot.png"))
# ggsave(plot = gene_interest_VarPart_plot, filename = paste0(OUTPUT_DIR, "example_gene_interest_VarPart_plot.png"))
# 
# head(varPart_fitAll[order(varPart_fitAll$material, decreasing=TRUE),])
# 
# 

# comparing vsd_filt, rld_filt
# and blind=TRUE, blind=FALSE
# vsd_filt_blindTRUE <- vst(dds_filt, blind = TRUE)    # blind = TRUE for QC
# rld_filt_blindTRUE <- rlog(dds_filt, blind = TRUE)   # blind = TRUE for QC
# vsd_filt_blindFALSE <- vst(dds_filt, blind = FALSE)  # not blind to batch effects
# rld_filt_blindFALSE <- rlog(dds_filt, blind = FALSE) # not blind to batch effects
# 
# test_pca_vsd_blindTRUE <- plot_pca(transf_object = vsd_filt_blindTRUE,
#                          intgroup=c(cond_interest, "library_type", "patient_oid"))
# test_pca_rld_blindTRUE <- plot_pca(transf_object = rld_filt_blindTRUE,
#                          intgroup=c(cond_interest, "library_type", "patient_oid"))
# test_pca_vsd_blindFALSE <- plot_pca(transf_object = vsd_filt_blindFALSE,
#                                     intgroup=c(cond_interest, "library_type", "patient_oid"))
# test_pca_rld_blindFALSE <- plot_pca(transf_object = rld_filt_blindFALSE,
#                                     intgroup=c(cond_interest, "library_type", "patient_oid"))
# 
# test_transf_blind_pca_list <- list(test_pca_vsd_blindTRUE = test_pca_vsd_blindTRUE, 
#                                    test_pca_rld_blindTRUE = test_pca_rld_blindTRUE, 
#                                    test_pca_vsd_blindFALSE = test_pca_vsd_blindFALSE, 
#                                    test_pca_rld_blindFALSE = test_pca_rld_blindFALSE)
# test_transf_blind_pca_plot <- ggpubr::ggarrange(plotlist = test_transf_blind_pca_list, 
#                                                 labels = c("vst: blind=TRUE", "rlog: blind=TRUE", "vst: blind=FALSE", "rlog: blind=FALSE"),
#                                                 ncol = 2, nrow = 2, 
#                                                 common.legend = TRUE)
