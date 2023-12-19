# loading libraries ----
library("tidyverse")

# cutNrun metadata

bam_files_md5sum <- "/home/peter_r/cutnrun_data/tala_CutnRun_bamfiles.md5sum"
bam_md5sums <- readr::read_delim(file = bam_files_md5sum, delim=" ", col_names=FALSE) %>%
  dplyr::rename(bam_md5sum = X1, 
                bam_filename = X2) %>%
  dplyr::mutate(bam_filename = stringr::str_trim(bam_filename, side = "both"))

metadata <- bam_md5sums %>%
  dplyr::rename(bsf_bam_filename = bam_filename) %>%
  dplyr::mutate(temp_sample_name = gsub(pattern = "(BSF_0993_HGC2CDRXY_1#S_.+)((Patient|ND).+)(_S.+bam)", replacement = "\\2", bsf_bam_filename)) %>%
  dplyr::mutate(replicate = gsub(pattern = "(Patient_|ND_JBe_)(1|2)(.+)", replacement = "\\2", temp_sample_name),
                marker = gsub(pattern = "(Patient_|ND_JBe_)(1|2)(.+)", replacement = "\\3", temp_sample_name),
                sample_group = gsub(pattern = "(Patient|ND_JBe)(_1|_2)(.+)", replacement = "\\1", temp_sample_name)) %>%
  dplyr::mutate(sample_name = paste0(sample_group, "_", marker, "_", "rep", replicate)) %>%
  dplyr::mutate(bamfilename = paste0(sample_name, ".bam")) %>%
  dplyr::arrange(sample_name)

bam_rename_cmd <- metadata %>%
  dplyr::mutate(rename_cmd = paste("cp", bsf_bam_filename, bamfilename),
                renamed_md5sum = paste(bam_md5sum, bamfilename, sep = "  ")) %>%  # check first then select
  dplyr::select(rename_cmd, renamed_md5sum)

write_lines(x = bam_rename_cmd$rename_cmd, file = "/home/peter_r/cutnrun_data/rename_bams_cutnrun.sh")
write_lines(x = bam_rename_cmd$renamed_md5sum, file = "/home/peter_r/cutnrun_data/renamed_bams_cutnrun.md5sums")


# TO-DO
# 1. rename
# 2. check md5sums
# 3. run find *.bam and convert to fastq files
# 4. run fastq screen - adjust config to check only for e.coli (is it the same as in spike-ins?)
# 5. prepare metadata with 'pseudo-bio' and technical replicates and run pipeline with specified references!
# 0. ! calculate md5sums on fastq files!

# nf-core samplesheet ----
# https://nf-co.re/cutandrun/dev/usage
# group,replicate,fastq_1,fastq_2
# target,1,H3K27me3_S1_L001_R1.fastq.gz,H3K27me3_S1_L001_R2.fastq.gz
# target,2,H3K27me3_S2_L001_R1.fastq.gz,H3K27me3_S2_L001_R2.fastq.gz
# target,3,H3K27me3_S3_L001_R1.fastq.gz,H3K27me3_S3_L001_R2.fastq.gz
# igg,1,IGG_S1_L001_R1.fastq.gz,IGG_S1_L001_R2.fastq.gz
# igg,2,IGG_S2_L001_R1.fastq.gz,IGG_S2_L001_R2.fastq.gz

# initial samplesheet to assess technical replicates
# for the production samplesheet need to change replicates to the same number! (since these are technical!)
# Patient_NFAT_rep2_R1.fastq.gz; Patient_NFAT_rep2_R2.fastq.gz

# fastq_1, fastq_2 - full path! /data_synology_rg2/core_bioinformatics/tala_atac_cut/tala_cutnrun/datasets/
samplesheet_qc <- metadata %>%
  dplyr::mutate(group = paste0(sample_group, "_", marker),
                fastq_1 = gsub(pattern = ".bam", replacement = "_R1.fastq.gz", x = bamfilename),
                fastq_2 = gsub(pattern = ".bam", replacement = "_R2.fastq.gz", x = bamfilename)) %>%
  dplyr::mutate(fastq_1 = paste0("/data_synology_rg2/core_bioinformatics/tala_atac_cut/tala_cutnrun/datasets/", fastq_1),
                fastq_2 = paste0("/data_synology_rg2/core_bioinformatics/tala_atac_cut/tala_cutnrun/datasets/", fastq_2)) %>%  # adding full path to fastq files 
  #dplyr::mutate() %>%  # changing igg to lower case!
  dplyr::select(group, replicate, fastq_1, fastq_2)

write_csv(x = samplesheet_qc, file = "/home/peter_r/cutnrun_data/cutnrun_samplesheet_qc.csv")
  
# metadata for collapsing technical replicates 
# production - collapsing replicates
samplesheet_quick_production <- metadata %>%
  dplyr::mutate(artificial_replicate = if_else(sample_group == "Patient", 1, 2),
                marker = if_else(marker == "IgG", "igg", marker)) %>%
  dplyr::mutate(group = marker,
                fastq_1 = gsub(pattern = ".bam", replacement = "_R1.fastq.gz", x = bamfilename),
                fastq_2 = gsub(pattern = ".bam", replacement = "_R2.fastq.gz", x = bamfilename)) %>%
  dplyr::mutate(fastq_1 = paste0("/data_synology_rg2/core_bioinformatics/tala_atac_cut/tala_cutnrun/datasets/", fastq_1),
                fastq_2 = paste0("/data_synology_rg2/core_bioinformatics/tala_atac_cut/tala_cutnrun/datasets/", fastq_2)) %>%  # adding full path to fastq files
  dplyr::mutate(technical_rep = replicate,
                replicate = artificial_replicate) %>%
  dplyr::select(group, replicate, fastq_1, fastq_2)

write_csv(x = samplesheet_quick_production, file = "/home/peter_r/cutnrun_data/cutnrun_samplesheet_prodHack.csv")

# collapsing technical replicates ----
metadata_collapse_technical <- metadata %>%
  dplyr::mutate(group = paste0(sample_group, "_", marker),
                fastq_1 = gsub(pattern = ".bam", replacement = "_R1.fastq.gz", x = bamfilename),
                fastq_2 = gsub(pattern = ".bam", replacement = "_R2.fastq.gz", x = bamfilename),
                merged_fq1 = gsub(pattern = "(_rep1|_rep2).bam", replacement = "_merged_R1.fastq.gz", x = bamfilename),
                merged_fq2 = gsub(pattern = "(_rep1|_rep2).bam", replacement = "_merged_R2.fastq.gz", x = bamfilename)) %>%
  dplyr::select(group, fastq_1, fastq_2, merged_fq1, merged_fq2) %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(cat_fq1 = paste(fastq_1, collapse = " "),
                   cat_fq2 = paste(fastq_2, collapse = " "),
                   merged_fq1 = paste0(unique(merged_fq1)),
                   merged_fq2 = paste0(unique(merged_fq2))) %>% # cmd to collapse fastq
  dplyr::mutate(cmd_merge_fq1 = paste("cat", cat_fq1, ">",merged_fq1),
                cmd_merge_fq2 = paste("cat", cat_fq2, ">",merged_fq2)) %>%
  dplyr::select(cmd_merge_fq1, cmd_merge_fq2) %>%
  tidyr::pivot_longer(., cols = c("cmd_merge_fq1", "cmd_merge_fq2"),names_to = "merge_name", values_to = "merge_cmds")

write_lines(x = metadata_collapse_technical$merge_cmds, file = "/home/peter_r/cutnrun_data/cutnrun_merge_fastq.sh")  
  
# run only important
# IKAROS and HELIOS with respective IgG
metadata_collapsed <- metadata %>%
  dplyr::mutate(group = paste0(sample_group, "_", marker),
                fastq_1 = gsub(pattern = ".bam", replacement = "_R1.fastq.gz", x = bamfilename),
                fastq_2 = gsub(pattern = ".bam", replacement = "_R2.fastq.gz", x = bamfilename),
                merged_fq1 = gsub(pattern = "(_rep1|_rep2).bam", replacement = "_merged_R1.fastq.gz", x = bamfilename),
                merged_fq2 = gsub(pattern = "(_rep1|_rep2).bam", replacement = "_merged_R2.fastq.gz", x = bamfilename)) %>% 
  dplyr::mutate(artificial_replicate = if_else(sample_group == "Patient", 1, 2),
                marker = if_else(marker == "IgG", "igg", marker)) %>%
  dplyr::mutate(group = marker,
                fastq_1 = merged_fq1,
                fastq_2 = merged_fq2) %>%
  dplyr::mutate(fastq_1 = paste0("/data_synology_rg2/core_bioinformatics/tala_atac_cut/tala_cutnrun/datasets/merged/", fastq_1),
                fastq_2 = paste0("/data_synology_rg2/core_bioinformatics/tala_atac_cut/tala_cutnrun/datasets/merged/", fastq_2)) %>%  # adding full path to fastq files
  dplyr::mutate(technical_rep = replicate,
                replicate = artificial_replicate) %>%
  dplyr::select(group, replicate, fastq_1, fastq_2) %>%
  dplyr::distinct(.)

write_csv(x = metadata_collapsed, file = "/home/peter_r/cutnrun_data/cutnrun_samplesheet_merged_all.csv")

important_markers <- c("Ikaros", "Helios", "igg")
metadata_collapsed_important <- metadata_collapsed %>%
  dplyr::filter(group %in% important_markers)

write_csv(x = metadata_collapsed_important, file = "/home/peter_r/cutnrun_data/cutnrun_samplesheet_important.csv")

# save bed file for plotting heatmaps
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene_cols <- c("gene_id", "tx_id", "tx_name", "tx_chrom", "tx_strand")
txdb_genes <- GenomicFeatures::genes(txdb, columns=gene_cols, single.strand.genes.only=TRUE)
rtracklayer::export.bed(txdb_genes, "/home/peter_r/cutnrun_data/Hsapiens_UCSC_hg38_knownGene.bed")

# testing
txdb_genes2 <- GenomicFeatures::genes(txdb)
rtracklayer::export.bed(txdb_genes2, "/home/peter_r/cutnrun_data/txdb_genes.bed")

txdb_transcripts <- GenomicFeatures::transcripts(txdb)
rtracklayer::export.bed(txdb_transcripts, "/home/peter_r/cutnrun_data/txdb_transcripts.bed")

txdb_promoters <- GenomicFeatures::promoters(txdb)
rtracklayer::export.bed(txdb_transcripts, "/home/peter_r/cutnrun_data/txdb_promoters.bed")

# for postprocessing
exclude_markers <- c("IgG", "Ikaros")
metadata_postprocessing <-  metadata %>%
  dplyr::mutate(group = paste0(sample_group, "_", marker),
                fastq_1 = gsub(pattern = ".bam", replacement = "_R1.fastq.gz", x = bamfilename),
                fastq_2 = gsub(pattern = ".bam", replacement = "_R2.fastq.gz", x = bamfilename)) %>%
  dplyr::mutate(fastq_1 = paste0("/data_synology_rg2/core_bioinformatics/tala_atac_cut/tala_cutnrun/datasets/", fastq_1),
                fastq_2 = paste0("/data_synology_rg2/core_bioinformatics/tala_atac_cut/tala_cutnrun/datasets/", fastq_2)) %>%  # adding full path to fastq files 
  #dplyr::mutate() %>%  # changing igg to lower case!
  dplyr::filter(!marker %in% exclude_markers) %>%
  dplyr::mutate(sample_name = paste0(group, "_R", replicate)) %>%
  dplyr::select(sample_name, replicate, marker)

write_csv(metadata_postprocessing, file = "/home/peter_r/cutnrun_data/sample_table_postprocessing.csv")
