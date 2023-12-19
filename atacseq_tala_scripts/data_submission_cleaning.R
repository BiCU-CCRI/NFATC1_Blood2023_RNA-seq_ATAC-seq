# Tala data submission - filenames cleaning!

# load tala_data_submission_PR.xlsx
# load tala_paper_processed_GEO.md5sum
# load tala_paper_raw_EGA.md5sum

# generate tables of changes and match md5sums

# copy files to new names and verify md5sums
# modify column names in the featureCounts!!!

# Encrypt for EGA
# add metadata and upload metadata to EGA
# prepare access control documents for Irinka and Kaan!

# loading libraries 
library(tidyverse)

# Prepare submission for GEO

input_dir <- "/home/peter_r/datasets/data_submission/"

# loading omics metadata/file overview

metadata_rawdata_raw <- readxl::read_excel(file.path(input_dir, "tala_data_submission_PR.xlsx"),
                                   sheet = "EGA_raw_data")

metadata_processed_raw <- readxl::read_excel(file.path(input_dir, "tala_data_submission_PR.xlsx"),
                                   sheet = "GEO_processed_data")

# md5sums
md5sums_rawdata <- readr::read_table(file.path(input_dir, "tala_paper_raw_EGA.md5sum"), col_names = c("md5sum", "file")) %>%
  dplyr::mutate(filename = gsub(pattern = "(\\./.+/)(.+|.bam|.narrowPeak|.fastq.gz)", replacement = "\\2", x = file),
                filepath = gsub(pattern = "(\\./.+/)(.+|.bam|.narrowPeak|.fastq.gz)", replacement = "\\1", x = file))

md5sums_processed <- readr::read_table(file.path(input_dir, "tala_paper_processed_GEO.md5sum"), col_names = c("md5sum", "file")) %>%
  dplyr::mutate(filename = dplyr::case_when(grepl(pattern = "ND_|P3_", x = file) ~ gsub(pattern = "./processed_data/scrna_seq/", replacement = "", x = file),
                                            TRUE ~ gsub(pattern = "(\\./.+/)(.+|.tsv|.tsv.gz|.mtx.gz|.csv)", replacement = "\\2", x = file)),
                filepath = gsub(pattern = "(\\./.+/)(.+|.tsv|.tsv.gz|.mtx.gz|.csv)", replacement = "\\1", x = file)) %>%
  # adjusting filenames for scrnaseq to contain also folder and to be able to distinguish between patients
  dplyr::mutate(filename = gsub(pattern = "/", replacement = "-", filename))




# add extra metadata about illumina etc. after new sample names are generated and md5sums properly logged in!
metadata_rawdata <- metadata_rawdata_raw %>%
  # expanding concatenated rows (;) file_extensions, filename
  tidyr::separate_rows(., file_extensions, filename, sep = ";") %>%
  dplyr::left_join(., md5sums_rawdata, by = "filename") %>%
  # populating empty md5sums 
  dplyr::mutate(md5sum = dplyr::coalesce(md5sum.x, md5sum.y)) %>%
  dplyr::select(-md5sum.x, -md5sum.y) %>%
  # constructing new file names and paths
  dplyr::mutate(omics_method = gsub(pattern = "-", replacement = "_", x = omics)) %>%
  dplyr::mutate(new_filename = paste(new_sample_name, file_extensions, sep = "_"),
                new_file_path = paste0(input_dir, "data_submission_renamed/raw_data/", omics_method, "/")) %>%
  dplyr::mutate(cp_from = paste0(filepath, filename),
                cp_from_absolutePath = gsub(pattern = "./raw_data", replacement = paste0(input_dir, "raw_data"), x=cp_from),
                cp_to = paste0(gsub(pattern = input_dir, replacement = "./", x=new_file_path), new_filename),
                cp_to_abolutePath = paste0(new_file_path, new_filename),
                cp_cmd = paste0("cp ", cp_from, " ", cp_to),
                new_filename_md5sum = paste0(md5sum, "  ", gsub(pattern = "\\./data_submission_renamed/", replacement = "\\./", cp_to))) %>%
  dplyr::mutate(new_fastq_path = gsub(pattern = "_unaligned.bam", replacement = ".fastq.gz", x = cp_to),
                sed_cmd = dplyr::case_when(file_extensions == "unaligned.bam" ~ paste0("samtools fastq --threads 2 ", cp_to, " | sed 's/",sample_name_BSF,"/",sample_name_article_new,"/g' | pigz -p6 > ", new_fastq_path),
                                           file_extensions == "peaks.narrowPeak" ~ paste0("sed -i 's/",sample_name_BSF,"/",sample_name_article_new,"/g' ", new_fastq_path),
                                           TRUE ~ NA_character_)
                )

# sed command
# samtools fastq --threads 2 atacseq_C6_IKZF2-WT_unaligned.bam | sed 's/ND_02/C6/g' | pigz -p6 > atacseq_C6_IKZF2-WT.fastq.gz
# sed -i 's/ND_01/C5/g' ./data_submission_renamed/raw_data/atac_seq/atacseq_C5_IKZF2-WT_peaks.narrowPeak  # to rename sample_peak names

identical(metadata_rawdata$cp_from, metadata_rawdata$file)

# creating new directories ----
# Note: directories are created here, but the actual copy command and md5sum -c is run on the command line using the following scripts
purrr::map(unique(metadata_rawdata$new_file_path), function(x) {dir.create(path = x, showWarnings = TRUE, recursive=TRUE)})

# copy and rename files ----
# file.copy(from = metadata_rawdata$cp_from,
#          to = metadata_rawdata$cp_to,
#          recursive = FALSE)

# saving copy commands for reference
writeLines(metadata_rawdata$cp_cmd, con = file.path(input_dir, "cp_rename_raw.sh"))

# saving md5sums to check copy process
writeLines(metadata_rawdata$new_filename_md5sum, con = file.path(input_dir, "tala_paper_raw_EGA_renamed.md5sums"))

# sed command write
writeLines(metadata_rawdata$sed_cmd[!is.na(metadata_rawdata$sed_cmd)], con = file.path(input_dir, "ubam2fastq_sed_samplenames.sh"))

# Preparing processed data ----
# for EGA and GEO
# note: featureCounts - need to rename columns! so they match filenames!

metadata_processed <- metadata_processed_raw %>%
  # expanding concatenated rows (;) file_extensions, filename
  tidyr::separate_rows(., file_extensions, filename, sep = ";") %>%
  # adjusting filenames for scrnaseq to contain also folder and to be able to distinguish between patients
  dplyr::mutate(filename = gsub(pattern = "/", replacement = "-", filename)) %>%
  dplyr::left_join(., md5sums_processed, by = "filename") %>%
  # populating empty md5sums 
  dplyr::mutate(md5sum = dplyr::coalesce(md5sum.x, md5sum.y)) %>%
  dplyr::select(-md5sum.x, -md5sum.y) %>%
  # constructing new file names and paths
  dplyr::mutate(omics_method = gsub(pattern = "-", replacement = "_", x = omics)) %>%
  dplyr::mutate(new_filename = paste(new_sample_name, file_extensions, sep = "_"),
                new_file_path = paste0(input_dir, "data_submission_renamed/processed_data/", omics_method, "/")) %>%
  dplyr::mutate(cp_from = file,
                cp_from_absolutePath = gsub(pattern = "./processed_data", replacement = paste0(input_dir, "processed_data"), x=cp_from),
                cp_to = paste0(gsub(pattern = input_dir, replacement = "./", x=new_file_path), new_filename),
                cp_to_abolutePath = paste0(new_file_path, new_filename),
                cp_cmd = paste0("cp ", cp_from, " ", cp_to),
                new_filename_md5sum = paste0(md5sum, "  ", gsub(pattern = "\\./data_submission_renamed/", replacement = "\\./", cp_to))) 


purrr::map(unique(metadata_processed$new_file_path), function(x) {dir.create(path = x, showWarnings = TRUE, recursive=TRUE)})

# saving copy commands for reference
writeLines(metadata_processed$cp_cmd, con = file.path(input_dir, "cp_rename_processed.sh"))

# saving md5sums to check copy process
writeLines(metadata_processed$new_filename_md5sum, con = file.path(input_dir, "tala_paper_processed_GEO_renamed.md5sums"))

# modifying featureCounts columns
smartseq_column_samples <- metadata_rawdata %>%
  dplyr::filter(omics == "smart-seq") %>%
  dplyr::select(sample_name_analysis, new_sample_name)


smartseq_columnames <- smartseq_column_samples$new_sample_name
names(smartseq_columnames) <- smartseq_column_samples$sample_name_analysis

feature_counts_old <- readr::read_csv(file = "/home/peter_r/datasets/data_submission/data_submission_renamed/processed_data/smart_seq/smartseq_P_C1_C5_C6_featureCounts.csv")

smartseq_columnames_ord <- smartseq_columnames[colnames(feature_counts_old)]

feature_counts_new <- feature_counts_old %>% setNames(smartseq_columnames_ord) # TRUE
  
# testing writing if it produces same file before column renaming
#readr::write_csv(feature_counts_old, file = "/home/peter_r/datasets/data_submission/data_submission_renamed/processed_data/smart_seq/smartseq_P_C1_C5_C6_featureCounts_test.csv")

#readr::write_csv(feature_counts_new, file = "/home/peter_r/datasets/data_submission/data_submission_renamed/processed_data/smart_seq/smartseq_P_C1_C5_C6_featureCounts.csv")

# EGA submission ----
# sample registration
sample_registration_columns <- c("title", "alias", "description", "subjectId", "bioSampleId", "caseOrControl", "gender", "organismPart", "cellLine", "region", "phenotype")
ega_samples <- metadata_rawdata %>%
  dplyr::select(omics, sample_name_article_new, new_sample_name, patient_sex, patient_age, tissue, cell_type, sample_description, genotype) %>%
  dplyr::mutate(title = paste0(sample_name_article_new, "-", omics),
                alias = new_sample_name,
    description = sample_description,
    subjectId = sample_name_article_new,
    caseOrControl = if_else(genotype == "I325V-Hom", "case", "control"),
    gender = patient_sex,
    organismPart = tissue,
    cellLine = cell_type,
    phenotype = if_else(genotype == "I325V-Hom", "Patient with defects of immunity", "Healthy donor")) %>%
  dplyr::mutate(bioSampleId = NA,
                region = NA) %>%
  dplyr::select(title, alias, description, subjectId, bioSampleId, caseOrControl, gender, organismPart, cellLine, region, phenotype) %>%
  # reducing to unique samples - 48 total
  dplyr::distinct(.)

readr::write_csv(x = ega_samples, file = file.path(input_dir, "ega_samples.csv"))

# extracting samples per experiment
# should have 62 files
# atac-seq 7 (6 raw 3 fastq and 3 narrowPeak + 1 processed)
# smart-seq 25 (24 raw fastq 2xcell types B-cells/T-cells x 4 patients x 3 replicates + processed featureCounts file)
# scrna-seq 30 (15 raw fastq I1, R1, R2 + 15 processed barcodes, features, matrix)

# combine raw and processed into one data.frame
# add md5sums 
# extract respective parts for e.g. fastq and analysis (phenotype)
rnaseq_columns <- c("Sample alias", "Fastq File", "Checksum", "Unencrypted checksum")
ega_omics_files <- metadata_rawdata %>%
  dplyr::filter(omics == "scrna-seq" | file_extensions == "unaligned.bam") %>%
  dplyr::select(omics, sample_name_article_new, new_fastq_path, new_sample_name, patient_sex, patient_age, tissue, cell_type, sample_description, genotype) %>%
  dplyr::mutate(fastq_filename = gsub(pattern = "(./data_submission_renamed/raw_data/.+/)(atacseq|scrnaseq|smartseq.+.fastq.gz)", replacement = "\\2", x = new_fastq_path))

# md5sums
encrypted_md5sums_raw <- list.files(path = "/home/peter_r/datasets/data_submission/EGA_submission/", pattern = ".md5$", full.names = TRUE)
names(encrypted_md5sums_raw) <- gsub(pattern = "/home/peter_r/datasets/data_submission/EGA_submission//", replacement = "", x= encrypted_md5sums_raw)

encrypted_md5sums_df <- purrr::map_dfr(names(encrypted_md5sums_raw), function(x) {
  data.frame(x, md5sum = readr::read_table(file = encrypted_md5sums_raw[[x]], col_names = "md5sum")$md5sum)
}) 

encrypted_md5sums <- encrypted_md5sums_df %>% 
  dplyr::rename(file = x) %>%
  dplyr::mutate(filename = gsub(pattern = ".md5|.gpg.md5", replacement = "", x = file))

encrypted_md5sums_unencrypt <- encrypted_md5sums %>%
  dplyr::filter(!grepl(pattern = "gpg.md5", x = file)) %>%
  dplyr::rename(file_unencrypted = file,
                md5sum_unencrypted = md5sum)

encrypted_md5sums_encrypt <- encrypted_md5sums %>%
  dplyr::filter(grepl(pattern = "gpg.md5", x = file)) %>%
  dplyr::rename(file_encrypted = file,
                md5sum_encrypted = md5sum)

# [ ] - document submission and generate interactive report for the future!
encrypted_md5sums_final <- encrypted_md5sums_unencrypt %>%
  left_join(., encrypted_md5sums_encrypt, by = "filename") %>%
  dplyr::mutate(alias = gsub(pattern = "(.+)(.fastq.gz|_peaks.narrowPeak|_consensus_peaks_counts.tsv|_featureCounts.csv|_matrix.mtx.gz|_features.tsv.gz|_barcodes.tsv.gz)", replacement = "\\1", x = filename)) %>%
  dplyr::mutate(file_extension = gsub(pattern = "(.+)(.fastq.gz|_peaks.narrowPeak|_consensus_peaks_counts.tsv|_featureCounts.csv|_matrix.mtx.gz|_features.tsv.gz|_barcodes.tsv.gz)", replacement = "\\2", x = filename)) %>%
  dplyr::mutate(alias = gsub(pattern = "(.+)(_S.+_L00.+_.+_001)", replacement = "\\1", x = alias))

sum(!(encrypted_md5sums_final$alias %in% ega_samples$alias))

table(encrypted_md5sums_final$file_extension)
encrypted_md5sums_final[!(encrypted_md5sums_final$alias %in% ega_samples$alias),]

EGA_fastq_files <- encrypted_md5sums_final %>%
  dplyr::filter(file_extension == ".fastq.gz") %>%
  dplyr::select(alias, filename, md5sum_encrypted, md5sum_unencrypted) %>%
  dplyr::rename(`Sample alias` = alias,
                `Fastq File` = filename,
                `Checksum` = md5sum_encrypted,
                `Unencrypted checksum` = md5sum_unencrypted) 

EGA_fastq_files_scRNAseq <- EGA_fastq_files %>%
  dplyr::filter(grepl(pattern = "scrnaseq_.+", x = `Sample alias`))

EGA_fastq_files_RNAseq <- EGA_fastq_files %>%
  dplyr::filter(grepl(pattern = "smartseq_.+", x = `Sample alias`))

EGA_fastq_files_ATACseq <- EGA_fastq_files %>%
  dplyr::filter(grepl(pattern = "atacseq_.+", x = `Sample alias`))

EGA_analysis_files <- encrypted_md5sums_final %>%
  dplyr::filter(file_extension != ".fastq.gz")

readr::write_csv(x = EGA_fastq_files_scRNAseq, file = file.path(input_dir, "EGA_fastq_files_scRNAseq.csv"))
readr::write_csv(x = EGA_fastq_files_RNAseq, file = file.path(input_dir, "EGA_fastq_files_RNAseq.csv"))
readr::write_csv(x = EGA_fastq_files_ATACseq, file = file.path(input_dir, "EGA_fastq_files_ATACseq.csv"))
readr::write_csv(x = EGA_analysis_files, file = file.path(input_dir, "EGA_analysis_files.csv"))

# Experiments descriptions
# for atac-seq and smart-seq it is Illumina HiSeq 3000/4000, however had to pick one
#EGA_experiments <- data.frame(design_name = c("scRNA-seq experiment", "RNA-seq experiment", "ATAC-seq experiment"),
#                               instrument_model = c("Illumina HiSeq 4000", "Illumina HiSeq 4000", "Illumina HiSeq 4000"))

