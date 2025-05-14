library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(microbiome)
library(dada2)

library(ggplot2)
library(gridExtra)
library(readxl)
library(viridis)
library(hrbrthemes)
library(ggthemes)
library(scales)
library(dplyr)
library(vegan)
library(ape)
library(ggnewscale)
library(phyloseq)
library(ggpubr)
library(ggsignif)

# The idea is to use all the 16S full length sequences of individual strains
    # first we filter and trim the sequences from each strain
    # then, infer ASVs - there should be 1, 2 or 4 depending on the species and #unique 16S sequences
    # compare these to make sure they agree with sequences extracted from sequenced genomes
    # found in 99_Analysis/CustomDada2DB_creation/ESL_strains_extracted16S
    # finally, we create a custom database with the ASVs and the full length 16S sequences
    # that can be used by dada2
getwd()
setwd("/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20240399_aprasad_PriorityEffects/")
dir.create("99_Analysis/CustomDada2DB_creation/RDS")
dir.create("99_Analysis/CustomDada2DB_creation/Figures")
dir.create("99_Analysis/CustomDada2DB_creation/Data")
# first, get the list of strains that need to go in the custom DB
strains_files_raw <- Sys.glob(file.path("99_Analysis/CustomDada2DB_creation/00-ESL_strains_raw_pacbio16S/*"))
strains <- as.character( lapply( strains_files_raw, function(x) strsplit(basename(x), ".fastq")[[1]] ))

primer_remove_in_files = strains_files_raw
primer_remove_out_files = as.character(lapply(strains, function(x) paste0("99_Analysis/CustomDada2DB_creation/01-Dada2Preprocessed/", x, "_primers_removed.fastq.gz")))
f_primer <- "AGRGTTYGATYMTGGCTCAG"
r_primer <- "RGYTACCTTGTTACGACTT"
r_primer_comp <- rc(r_primer)
# primer_removal_summary <- removePrimers(primer_remove_in_files, primer_remove_out_files,
#               primer.fwd = f_primer,
#               primer.rev = r_primer_comp,
#               verbose = TRUE)
# saveRDS(primer_removal_summary, "99_Analysis/CustomDada2DB_creation/RDS/primer_removal_summary.rds")
primer_removal_summary <- readRDS("99_Analysis/CustomDada2DB_creation/RDS/primer_removal_summary.rds")

filt_read_files <- as.character(lapply(strains, function(x) paste0("99_Analysis/CustomDada2DB_creation/01-Dada2Preprocessed/", x, "_filt_reads.fastq.gz")))
# trimming_summary <- filterAndTrim(primer_remove_out_files, filt_read_files, minLen=1000, maxLen=1600, rm.phix=FALSE, multithread=T)
# saveRDS(trimming_summary, "99_Analysis/CustomDada2DB_creation/RDS/trimming_summary.rds")
trimming_summary <- readRDS("99_Analysis/CustomDada2DB_creation/RDS/trimming_summary.rds")

# derepd <- derepFastq(filt_read_files, verbose = T)
# errors <- learnErrors(derepd, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)
# saveRDS(errors, "99_Analysis/CustomDada2DB_creation/RDS/errors.rds")
errors <- readRDS("99_Analysis/CustomDada2DB_creation/RDS/errors.rds")
plotErrors(errors)
# dada2 <- dada(derepd, err=errors, BAND_SIZE=32, multithread=T)
# saveRDS(dada2, "99_Analysis/CustomDada2DB_creation/RDS/dada2.rds")
dada2 <- readRDS("99_Analysis/CustomDada2DB_creation/RDS/dada2.rds")
# seqtable <- makeSequenceTable(dada2); dim(seqtable)
# seqtable_nochim <- removeBimeraDenovo(seqtable, method = "consensus", multithread = TRUE, verbose = TRUE); dim(seqtable)
# sum(seqtable_nochim)/sum(seqtable)
# write.csv(file="99_Analysis/CustomDada2DB_creation/seqtable_nochim.csv", seqtable_nochim, quote = F)
seqtable_nochim <- read.csv("99_Analysis/CustomDada2DB_creation/seqtable_nochim.csv", row.names = 1)
# silva_DB_assignTax <- "99_Analysis/Databases/silva_nr99_v138_wSpecies_train_set.fa.gz"
# taxonomySilva <- assignTaxonomy(seqtable_nochim, silva_DB_assignTax, minBoot = 80, multithread=TRUE, verbose = TRUE)
# saveRDS(taxonomySilva, "99_Analysis/CustomDada2DB_creation/RDS/taxonomySilva.rds")
# taxonomySilva <- readRDS("99_Analysis/CustomDada2DB_creation/RDS/taxonomySilva.rds")
# get all numbers from seqtable_nochim and make a histogram
# remove 0s
# as.list(seqtable_nochim)[as.list(seqtable_nochim) > 10 & as.list(seqtable_nochim) < 1000000] %>% unlist() %>% hist(breaks = 1000)
get_samples_of_occurance <- function(x) {
    sub_table <- seqtable_nochim %>% as.data.frame() %>%
        select(x)
    colnames(sub_table) <- c("count")
    names <- sub_table %>% filter(count > 0) %>% rownames_to_column("genome") %>% pull(genome)
    names <- as.character(lapply(names, function(name) strsplit(name, "_filt_reads")[[1]][[1]]))
    if (length(names) == 0) {
        names <- c("None")
    }
    return(names)
}
# make ids for ASVs from taxonomySilva
ASV_nums <- 1:length(colnames(seqtable_nochim))
# ASV_nums <- 1:nrow(taxonomySilva)
ASV_ids <- paste0("ASV", ASV_nums)
asv_num_df <- colnames(seqtable_nochim) %>% as.data.frame()
rownames(asv_num_df) <- ASV_ids
colnames(asv_num_df) <- c("ASV_seq")
asv_num_df <- asv_num_df %>% rownames_to_column("ASV_id")
write.csv(file="99_Analysis/CustomDada2DB_creation/ASV_ids.csv", asv_num_df, quote = F, row.names = F)
# write asv sequences to fasta with seq as asv_num_df$ASV_seq and header as asv_num_df$ASV_id
asv_fasta <- asv_num_df %>% mutate(ASV_seq = paste0(">", ASV_id, "\n", ASV_seq)) %>% pull(ASV_seq)
write(asv_fasta, "99_Analysis/CustomDada2DB_creation/ASV_seqs.fasta")
strain_metadata <- read.csv("Strain_metadata.tsv", sep = "\t") %>%
                    rbind(c("None", "Unknown", "unknown"))
taxonomy_df_info <- data.frame()
for (i in 1:nrow(taxonomySilva)) {
  print(i)
    ASV <- rownames(taxonomySilva)[i]
    strains <- get_samples_of_occurance(ASV)
    print(paste(as.character(sapply(strains, function(x) strain_metadata %>% filter(Strain == x) %>% pull(Species)))))
    for (strain in strains) {
        taxonomy_df_info <- rbind(taxonomy_df_info, data.frame(ASV = ASV, strain = strain, count = seqtable_nochim[paste0(strain, "_filt_reads.fastq.gz"), ASV]))
    }
}
low_coverage_data_strains <- c("ESL0199",
                               "ESL0186",
                               "ESL0185",
                               "ESL0399",
                               "ESL0820",
                               "ESL0197",
                               "ESL0200"
)
taxonomy_df_info <- taxonomy_df_info %>%
                      left_join(strain_metadata, by = c("strain" = "Strain")) %>%
                        mutate(low_coverage = ifelse(strain %in% low_coverage_data_strains, "yes", "no")) %>%
                          group_by(ASV) %>%
                            mutate(is_unique = ifelse(n_distinct(strain) == 1, "yes", "no")) %>%
                              ungroup()
write.csv(file="99_Analysis/CustomDada2DB_creation/ASV_per_strain_info.csv", taxonomy_df_info, quote = F, row.names = F)
get_asv_num <- function(x) {
    return(asv_num_df %>% filter(ASV_seq == x) %>% pull(ASV_id))
}
heatmap_data <- taxonomy_df_info %>%
                  select(ASV, strain, count) %>%
                    group_by(ASV, strain) %>%
                      unique() %>%
                      arrange(strain) %>%
                        pivot_wider(names_from = strain, values_from = count, values_fill = 0) %>%
                          as.data.frame() %>%
                          mutate(ASV = Vectorize(get_asv_num)(ASV)) %>%
                            column_to_rownames("ASV") %>%
                              as.matrix()
species_names <- heatmap_data %>% colnames() %>% as.data.frame() %>% rownames_to_column("strain") %>% left_join(strain_metadata, by = c("." = "Strain")) %>% pull(Species)
col_split <- species_names %>% as.factor()
# highlight names of low coverage strains in red
bottom_anno <- HeatmapAnnotation(df = data.frame(low_coverage = ifelse(colnames(heatmap_data) %in% low_coverage_data_strains, "yes", "no")),
                                 col = list(low_coverage = c("yes" = "red", "no" = "white")),
                                 border = T,
                                 show_legend = T)
genusColors <- c("Bombilactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[1],
                 "Lactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[4],
                 "Bifidobacterium" = brewer.pal(11, "Spectral")[3],
                 "Apibacter" = brewer.pal(11, "Spectral")[4],
                 "Bombella" = brewer.pal(11, "Spectral")[5],
                 "Commensalibacter" = brewer.pal(11, "Spectral")[6],
                 "Bartonella" = brewer.pal(11, "Spectral")[7],
                 "Frischella" = brewer.pal(11, "Spectral")[8],
                 "Apilactobacillus" = brewer.pal(11, "Spectral")[9],
                 "Snodgrassella" = brewer.pal(11, "Spectral")[10],
                 "Gilliamella" = brewer.pal(11, "Spectral")[11]
)
top_anno <- HeatmapAnnotation(Species = anno_text(species_names %>% as.factor() %>% as.character(), location = 0,
                                 rot = 45, just = c("left", "bottom")),
                              Genus = as.character(sapply(species_names, function(x) strsplit(x, " ")[[1]][1])),
                              col = list(Genus = genusColors)
                              )
# annotate unique or not
left_anno <- rowAnnotation(is_unique = ifelse(heatmap_data %>% rowstrainSums() > 1, "no", "yes"),
                            col = list(is_unique = c("yes" = "black", "no" = "white")),
                            show_legend = T
)
heatmap_obj = Heatmap(heatmap_data,
                      name = "Present",
                      col = colorRamp2(c(0, 1, 10, 100, 1000, 10000, 100000), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[2], brewer.pal(9, "PuBu")[3], brewer.pal(9, "PuBu")[4], brewer.pal(9, "PuBu")[7], brewer.pal(9, "PuBu")[8], brewer.pal(9, "PuBu")[9])),
                      heatmap_legend_param = list(at = c(0, 1, 10, 100, 1000, 10000, 100000),
                                                  legend_height = unit(5, "cm"),
                                                  labels = c(0, 1, 10, 100, 1000, 10000, 100000)),
                      column_split = col_split,
                      bottom_annotation = bottom_anno,
                      top_annotation = top_anno,
                      # left_annotation = left_anno,
                      cluster_rows = F,
                      cluster_columns = F,
                      column_title_rot = 40,
                      column_title_gp = gpar(fontsize = 0),
                      cluster_row_slices = T,
                      border = T,
                      show_row_names = T,
                      show_heatmap_legend = T)
pdf("CustomDada2DB_creation/Figures/ASV_per_strain_heatmap_by_genus.pdf", width = 10, height = 15)
draw(heatmap_obj, heatmap_legend_side = "right")
dev.off()
heatmap_obj = Heatmap(heatmap_data,
                      name = "Present",
                      col = colorRamp2(c(0, 1, 10, 100, 1000, 10000, 100000), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[2], brewer.pal(9, "PuBu")[3], brewer.pal(9, "PuBu")[4], brewer.pal(9, "PuBu")[7], brewer.pal(9, "PuBu")[8], brewer.pal(9, "PuBu")[9])),
                      heatmap_legend_param = list(at = c(0, 1, 10, 100, 1000, 10000, 100000),
                                                  legend_height = unit(5, "cm"),
                                                  labels = c(0, 1, 10, 100, 1000, 10000, 100000)),
                      # column_split = col_split,
                      bottom_annotation = bottom_anno,
                      top_annotation = top_anno,
                      # left_annotation = left_anno,
                      cluster_rows = F,
                      cluster_columns = F,
                      column_title_rot = 40,
                      column_title_gp = gpar(fontsize = 0),
                      cluster_row_slices = T,
                      border = T,
                      show_row_names = T,
                      show_heatmap_legend = T)
pdf("CustomDada2DB_creation/Figures/ASV_per_strain_heatmap.pdf", width = 10, height = 15)
draw(heatmap_obj, heatmap_legend_side = "right")
dev.off()
# divide each count by the sum of counts in the column
heatmap_data_norm <- t(t(heatmap_data)/rowSums(t(heatmap_data)))
species_names <- heatmap_data_norm %>% colnames() %>% as.data.frame() %>% rownames_to_column("strain") %>% left_join(strain_metadata, by = c("." = "Strain")) %>% pull(Species)
col_split <- species_names %>% as.factor()
# highlight names of low coverage strains in red
bottom_anno <- HeatmapAnnotation(df = data.frame(low_coverage = ifelse(colnames(heatmap_data_norm) %in% low_coverage_data_strains, "yes", "no")),
                                 col = list(low_coverage = c("yes" = "red", "no" = "white")),
                                 border = T,
                                 show_legend = T)
genusColors <- c("Bombilactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[1],
                 "Lactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[4],
                 "Bifidobacterium" = brewer.pal(11, "Spectral")[3],
                 "Apibacter" = brewer.pal(11, "Spectral")[4],
                 "Bombella" = brewer.pal(11, "Spectral")[5],
                 "Commensalibacter" = brewer.pal(11, "Spectral")[6],
                 "Bartonella" = brewer.pal(11, "Spectral")[7],
                 "Frischella" = brewer.pal(11, "Spectral")[8],
                 "Apilactobacillus" = brewer.pal(11, "Spectral")[9],
                 "Snodgrassella" = brewer.pal(11, "Spectral")[10],
                 "Gilliamella" = brewer.pal(11, "Spectral")[11]
)
top_anno <- HeatmapAnnotation(Species = anno_text(species_names %>% as.factor() %>% as.character(), location = 0,
                                 rot = 45, just = c("left", "bottom")),
                              Genus = as.character(sapply(species_names, function(x) strsplit(x, " ")[[1]][1])),
                              col = list(Genus = genusColors)
                              )
# annotate unique or not
left_anno <- rowAnnotation(is_unique = ifelse(heatmap_data_norm %>% rowstrainSums() > 1, "no", "yes"),
                            col = list(is_unique = c("yes" = "black", "no" = "white")),
                            show_legend = T
)
heatmap_obj = Heatmap(heatmap_data_norm,
                      name = "Present",
                      # col = colorRamp2(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[2], brewer.pal(9, "PuBu")[3], brewer.pal(9, "PuBu")[4], brewer.pal(9, "PuBu")[7], brewer.pal(9, "PuBu")[8], brewer.pal(9, "PuBu")[9])),
                      col = colorRamp2(c(0, 1), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[9])),
                      # column_split = col_split,
                      bottom_annotation = bottom_anno,
                      top_annotation = top_anno,
                      # left_annotation = left_anno,
                      cluster_rows = F,
                      cluster_columns = F,
                      column_title_rot = 40,
                      cluster_row_slices = T,
                      border = T,
                      show_row_names = T,
                      show_heatmap_legend = T)
pdf("CustomDada2DB_creation/Figures/ASV_per_strain_heatmap_norm.pdf", width = 10, height = 15)
draw(heatmap_obj, heatmap_legend_side = "right")
dev.off()
heatmap_obj = Heatmap(heatmap_data_norm,
                      name = "Present",
                      # col = colorRamp2(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[2], brewer.pal(9, "PuBu")[3], brewer.pal(9, "PuBu")[4], brewer.pal(9, "PuBu")[7], brewer.pal(9, "PuBu")[8], brewer.pal(9, "PuBu")[9])),
                      col = colorRamp2(c(0, 1), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[9])),
                      column_split = col_split,
                      bottom_annotation = bottom_anno,
                      top_annotation = top_anno,
                      # left_annotation = left_anno,
                      cluster_rows = F,
                      cluster_columns = F,
                      column_title_gp = gpar(fontsize = 0),
                      column_title_rot = 40,
                      cluster_row_slices = T,
                      border = T,
                      show_row_names = T,
                      show_heatmap_legend = T)
pdf("CustomDada2DB_creation/Figures/ASV_per_strain_heatmap_norm_by_genus.pdf", width = 10, height = 15)
draw(heatmap_obj, heatmap_legend_side = "right")
dev.off()

# divide each count by the sum of counts in the column
# plot AQ and AP separately
strains_AQ <- strain_metadata %>% filter(Source == "AQ" | Source == "Both") %>% pull(Strain)
strains_AP <- strain_metadata %>% filter(Source == "AP" | Source == "Both") %>% pull(Strain)
heatmap_data_norm <- t(t(heatmap_data)/rowSums(t(heatmap_data))) %>% as.data.frame() %>%
                      select(all_of(strains_AQ)) %>%
                        as.matrix()
species_names <- heatmap_data_norm %>% colnames() %>% as.data.frame() %>% rownames_to_column("strain") %>% left_join(strain_metadata, by = c("." = "Strain")) %>% pull(Species)
col_split <- species_names %>% as.factor()
# highlight names of low coverage strains in red
bottom_anno <- HeatmapAnnotation(df = data.frame(low_coverage = ifelse(colnames(heatmap_data_norm) %in% low_coverage_data_strains, "yes", "no")),
                                 col = list(low_coverage = c("yes" = "red", "no" = "white")),
                                 border = T,
                                 show_legend = T)
genusColors <- c("Bombilactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[1],
                 "Lactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[4],
                 "Bifidobacterium" = brewer.pal(11, "Spectral")[3],
                 "Apibacter" = brewer.pal(11, "Spectral")[4],
                 "Bombella" = brewer.pal(11, "Spectral")[5],
                 "Commensalibacter" = brewer.pal(11, "Spectral")[6],
                 "Bartonella" = brewer.pal(11, "Spectral")[7],
                 "Frischella" = brewer.pal(11, "Spectral")[8],
                 "Apilactobacillus" = brewer.pal(11, "Spectral")[9],
                 "Snodgrassella" = brewer.pal(11, "Spectral")[10],
                 "Gilliamella" = brewer.pal(11, "Spectral")[11]
)
top_anno <- HeatmapAnnotation(Species = anno_text(species_names %>% as.factor() %>% as.character(), location = 0,
                                 rot = 45, just = c("left", "bottom")),
                              Genus = as.character(sapply(species_names, function(x) strsplit(x, " ")[[1]][1])),
                              col = list(Genus = genusColors)
                              )
# annotate unique or not
left_anno <- rowAnnotation(is_unique = ifelse(heatmap_data_norm %>% rowstrainSums() > 1, "no", "yes"),
                            col = list(is_unique = c("yes" = "black", "no" = "white")),
                            show_legend = T
)
heatmap_obj = Heatmap(heatmap_data_norm,
                      name = "Present",
                      col = colorRamp2(c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 1), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[2], brewer.pal(9, "PuBu")[3], brewer.pal(9, "PuBu")[4], brewer.pal(9, "PuBu")[7], brewer.pal(9, "PuBu")[8], brewer.pal(9, "PuBu")[9])),
                      # col = colorRamp2(c(0, 1), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[9])),
                      heatmap_legend_param = list(at = c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 1),
                                                  legend_height = unit(5, "cm"),
                                                  labels = c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 1)),
                      bottom_annotation = bottom_anno,
                      top_annotation = top_anno,
                      # left_annotation = left_anno,
                      cluster_rows = T,
                      cluster_columns = T,
                      column_title_rot = 40,
                      cluster_row_slices = T,
                      border = T,
                      show_row_names = T,
                      show_heatmap_legend = T)
pdf("CustomDada2DB_creation/Figures/ASV_per_strain_heatmap_norm_AQ.pdf", width = 10, height = 15)
draw(heatmap_obj, heatmap_legend_side = "right")
dev.off()
heatmap_obj = Heatmap(heatmap_data_norm,
                      name = "Present",
                      col = colorRamp2(c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 1), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[2], brewer.pal(9, "PuBu")[3], brewer.pal(9, "PuBu")[4], brewer.pal(9, "PuBu")[7], brewer.pal(9, "PuBu")[8], brewer.pal(9, "PuBu")[9])),
                      # col = colorRamp2(c(0, 1), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[9])),
                      column_split = col_split,
                      bottom_annotation = bottom_anno,
                      top_annotation = top_anno,
                      heatmap_legend_param = list(at = c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 1),
                                                  legend_height = unit(5, "cm"),
                                                  labels = c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 1)),
                      # left_annotation = left_anno,
                      cluster_rows = F,
                      cluster_columns = F,
                      column_title_gp = gpar(fontsize = 0),
                      column_title_rot = 40,
                      cluster_row_slices = T,
                      border = T,
                      show_row_names = T,
                      show_heatmap_legend = T)
pdf("CustomDada2DB_creation/Figures/ASV_per_strain_heatmap_norm_by_genus_AQ.pdf", width = 10, height = 15)
draw(heatmap_obj, heatmap_legend_side = "right")
dev.off()
# plot AQ and AP separately
heatmap_data_norm <- t(t(heatmap_data)/rowSums(t(heatmap_data))) %>% as.data.frame() %>%
                      select(all_of(strains_AP)) %>%
                        as.matrix()
species_names <- heatmap_data_norm %>% colnames() %>% as.data.frame() %>% rownames_to_column("strain") %>% left_join(strain_metadata, by = c("." = "Strain")) %>% pull(Species)
col_split <- species_names %>% as.factor()
# highlight names of low coverage strains in red
bottom_anno <- HeatmapAnnotation(df = data.frame(low_coverage = ifelse(colnames(heatmap_data_norm) %in% low_coverage_data_strains, "yes", "no")),
                                 col = list(low_coverage = c("yes" = "red", "no" = "white")),
                                 border = T,
                                 show_legend = T)
genusColors <- c("Bombilactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[1],
                 "Lactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[4],
                 "Bifidobacterium" = brewer.pal(11, "Spectral")[3],
                 "Apibacter" = brewer.pal(11, "Spectral")[4],
                 "Bombella" = brewer.pal(11, "Spectral")[5],
                 "Commensalibacter" = brewer.pal(11, "Spectral")[6],
                 "Bartonella" = brewer.pal(11, "Spectral")[7],
                 "Frischella" = brewer.pal(11, "Spectral")[8],
                 "Apilactobacillus" = brewer.pal(11, "Spectral")[9],
                 "Snodgrassella" = brewer.pal(11, "Spectral")[10],
                 "Gilliamella" = brewer.pal(11, "Spectral")[11]
)
top_anno <- HeatmapAnnotation(Species = anno_text(species_names %>% as.factor() %>% as.character(), location = 0,
                                 rot = 45, just = c("left", "bottom")),
                              Genus = as.character(sapply(species_names, function(x) strsplit(x, " ")[[1]][1])),
                              col = list(Genus = genusColors)
                              )
# annotate unique or not
left_anno <- rowAnnotation(is_unique = ifelse(heatmap_data_norm %>% rowSums() > 1, "no", "yes"),
                            col = list(is_unique = c("yes" = "black", "no" = "white")),
                            show_legend = T
)
heatmap_obj = Heatmap(heatmap_data_norm,
                      name = "Present",
                      col = colorRamp2(c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 1), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[2], brewer.pal(9, "PuBu")[3], brewer.pal(9, "PuBu")[4], brewer.pal(9, "PuBu")[7], brewer.pal(9, "PuBu")[8], brewer.pal(9, "PuBu")[9])),
                      # col = colorRamp2(c(0, 1), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[9])),
                      heatmap_legend_param = list(at = c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 1),
                                                  legend_height = unit(5, "cm"),
                                                  labels = c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 1)),
                      bottom_annotation = bottom_anno,
                      top_annotation = top_anno,
                      # left_annotation = left_anno,
                      cluster_rows = T,
                      cluster_columns = T,
                      column_title_rot = 40,
                      cluster_row_slices = T,
                      border = T,
                      show_row_names = T,
                      show_heatmap_legend = T)
pdf("CustomDada2DB_creation/Figures/ASV_per_strain_heatmap_norm_AP.pdf", width = 10, height = 15)
draw(heatmap_obj, heatmap_legend_side = "right")
dev.off()
heatmap_obj = Heatmap(heatmap_data_norm,
                      name = "Present",
                      col = colorRamp2(c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 1), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[2], brewer.pal(9, "PuBu")[3], brewer.pal(9, "PuBu")[4], brewer.pal(9, "PuBu")[7], brewer.pal(9, "PuBu")[8], brewer.pal(9, "PuBu")[9])),
                      # col = colorRamp2(c(0, 1), c(brewer.pal(9, "PuBu")[1], brewer.pal(9, "PuBu")[9])),
                      column_split = col_split,
                      bottom_annotation = bottom_anno,
                      top_annotation = top_anno,
                      heatmap_legend_param = list(at = c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 1),
                                                  legend_height = unit(5, "cm"),
                                                  labels = c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 1)),
                      # left_annotation = left_anno,
                      cluster_rows = F,
                      cluster_columns = F,
                      column_title_gp = gpar(fontsize = 0),
                      column_title_rot = 40,
                      cluster_row_slices = T,
                      border = T,
                      show_row_names = T,
                      show_heatmap_legend = T)
pdf("CustomDada2DB_creation/Figures/ASV_per_strain_heatmap_norm_by_genus_AP.pdf", width = 10, height = 15)
draw(heatmap_obj, heatmap_legend_side = "right")
dev.off()
# formatted_DB_assignSpec <- "/nas/FAC/FBM/DMF/pengel/spirit/D2c/aprasad/20220921_aprasad_PriorityEffectsExperimentPilot/03_PilotExperiment/07_PacBIO_sequencing/database/pilot_experiment_sequences_for_species_assignment.fasta"
# taxonomyCommunity <- assignSpecies(taxonomySilva, formatted_DB_assignSpec, allowMultiple = T, verbose = T)
# saveRDS(taxonomyCommunity, "99_Analysis/CustomDada2DB_creation/RDS/taxonomyCommunity.rds")
# taxonomyCommunity <- readRDS("99_Analysis/CustomDada2DB_creation/RDS/taxonomyCommunity.rds")
# taxonomy_df <- taxonomySilva %>% 
#                 as.data.frame %>%
#                   rownames_to_column("ASV") %>%
#                   left_join(taxonomyCommunity %>%
#                               as.data.frame() %>%
#                                 mutate(SDP = Genus) %>%
#                                 mutate(Strain = Species) %>%
#                                   select(!c(Genus, Species)) %>%
#                                     rownames_to_column("ASV")
#                                 ) %>% 
#                                 mutate(sample = Vectorize(get_sample_of_occurance)(ASV))
# write.csv(file="99_Analysis/CustomDada2DB_creation/Taxonomy_info.csv", taxonomy_df, quote = F, row.names = F)
dir.create("99_Analysis/CustomDada2DB_creation/Indiv_sample")
dir.create("99_Analysis/CustomDada2DB_creation/Indiv_sample/RDS")
# Run dada2 one sample at a time!
for (strain in strain_metadata$Strain) {
  if (strain == "None") {
    next
  }
  if (strain == "ESL0186") {
    next
  }
  print(strain)
  if (file.exists(paste0("99_Analysis/CustomDada2DB_creation/Indiv_sample/RDS/", strain,"_seqtable_nochim.csv"))) {
    next
  }
  filt_file = paste0("99_Analysis/CustomDada2DB_creation/01-Dada2Preprocessed/", strain, "_filt_reads.fastq.gz")
  derepd <- derepFastq(filt_file, verbose = T)
  errors <- learnErrors(derepd, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)
  saveRDS(errors, paste0("99_Analysis/CustomDada2DB_creation/Indiv_sample/RDS/", strain,"_errors.rds"))
  # errors <- readRDS(paste0("99_Analysis/CustomDada2DB_creation/Indiv_sample/RDS/", strain,"_errors.rds"))
  # plotErrors(errors)
  dada2 <- dada(derepd, err=errors, BAND_SIZE=32, multithread=T)
  saveRDS(dada2, paste0("99_Analysis/CustomDada2DB_creation/Indiv_sample/RDS/", strain,"_dada2.rds"))
  # dada2 <- readRDS(paste0("99_Analysis/CustomDada2DB_creation/Indiv_sample/RDS/", strain,"_dada2.rds"))
  seqtable <- makeSequenceTable(dada2); dim(seqtable)
  seqtable_nochim <- removeBimeraDenovo(seqtable, method = "consensus", multithread = TRUE, verbose = TRUE); dim(seqtable)
  sum(seqtable_nochim)/sum(seqtable)
  write.csv(file=paste0("99_Analysis/CustomDada2DB_creation/Indiv_sample/RDS/", strain,"_seqtable_nochim.csv"), seqtable_nochim, quote = F)
  # seqtable_nochim <- read.csv(paste0("99_Analysis/CustomDada2DB_creation/Indiv_sample/RDS/", strain,"_seqtable_nochim.csv"), row.names = 1)
}
df_info = read.csv("/users/aprasad/nas_recherche/20240399_aprasad_PriorityEffects/00_StrainSelection/LongReadGenomes/summaries/strain_spec_info.csv")
heatmap_data <- read.csv("/users/aprasad/nas_recherche/20240399_aprasad_PriorityEffects/00_StrainSelection/LongReadGenomes/summaries/uid_strain_copy_mat.csv") %>%
                  select(!X) %>%
                  column_to_rownames("uid") %>%
                    as.matrix()
species_names <- heatmap_data %>% colnames() %>% as.data.frame() %>% rownames_to_column("strain") %>% left_join(df_info, by = c("." = "strain")) %>% pull(spec_name)
genus_names <- unlist(lapply(species_names, function(x) strsplit(x, " ")[[1]][[2]]))
col_split <- species_names %>% as.factor()
# highlight names of low coverage strains in red
genusColors <- c("Bombilactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[1],
                 "Lactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[4],
                 "Bifidobacterium" = brewer.pal(11, "Spectral")[3],
                 "Apibacter" = brewer.pal(11, "Spectral")[4],
                 "Bombella" = brewer.pal(11, "Spectral")[5],
                 "Commensalibacter" = brewer.pal(11, "Spectral")[6],
                 "Bartonella" = brewer.pal(11, "Spectral")[7],
                 "Frischella" = brewer.pal(11, "Spectral")[8],
                 "Apilactobacillus" = brewer.pal(11, "Spectral")[9],
                 "Snodgrassella" = brewer.pal(11, "Spectral")[10],
                 "Gilliamella" = brewer.pal(11, "Spectral")[11]
)
top_anno <- HeatmapAnnotation(Genus = as.character(genus_names),
                                col = list(Genus = genusColors)
                                )
# top_anno <- HeatmapAnnotation(Species = anno_text(species_names %>% as.factor() %>% as.character(), location = 0,
#                                   rot = 45, just = c("left", "bottom")),
#                                 Genus = as.character(genus_names),
#                                 col = list(Genus = genusColors)
#                                 )
heatmap_obj = Heatmap(heatmap_data,
                      name = "Present",
                      col = colorRamp2(c(0, 1, 2, 3, 4),
                                       c("#ffffff",
                                        "#4393c3",
                                        "#2166ac",
                                        "#fdb863",
                                        "#b2182b"
                                    )
                                ),
                      # heatmap_legend_param = list(at = c(0, 1, 2, 3, 4),
                      #                             legend_height = unit(5, "cm"),
                      #                             labels = c(0, 1, 2, 3, 4)),
                      column_split = col_split,
                      top_annotation = top_anno,
                      # left_annotation = left_anno,
                      cluster_rows = F,
                      cluster_columns = F,
                      column_title_rot = 40,
                      column_title_gp = gpar(fontsize = 12),
                      cluster_row_slices = T,
                      border = T,
                      show_row_names = T,
                      show_heatmap_legend = T)
pdf("/users/aprasad/nas_recherche/20240399_aprasad_PriorityEffects/00_StrainSelection/LongReadGenomes/summaries/strain_ASV_heatmpa.pdf", width = 10, height = 15)
draw(heatmap_obj, heatmap_legend_side = "right")
dev.off()