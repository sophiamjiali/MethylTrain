# ==============================================================================
# Script:           parse_annotations.R
# Purpose:          Consolidates methylation annotations from Bioconductor
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             12/23/2025
#
# Notes:            Fetches cross-reactive, SNP, and sex chromosome annotations
#                   from Bioconductor for TCGA-based workflows, consolidating 
#                   them into a table and saving them as a CSV for compatibility
#                   with python-based workflows; serves as general probe anno.
#
# Platform:         Illumina Methylation Annotation
# Genome Build:     hg19
# Arrays:           27K, 450K, EPIC (850K)
# ==============================================================================

library(minfi)
library(AnnotationDbi)
library(dplyr)
library(tibble)
library(readr)
library(maxprobes)
library(rtracklayer)
library(GenomicRanges)
library(rtracklayer)
library(R.utils)

# =====| Illumina 27K hg19|=====================================================

library(IlluminaHumanMethylation27kanno.ilmn12.hg19)

anno_27k <- getAnnotation(IlluminaHumanMethylation27kanno.ilmn12.hg19)

# SNP, cross-reactive, and multi-map not curated for the 27K array
manifest_27k <- tibble(
  probe_id          = rownames(anno_27k),
  chr               = as.character(anno_27k$chr),
  pos               = anno_27k$pos,
  strand            = anno_27k$strand,
  gene_symbol       = anno_27k$Symbol,
  is_sex_chr        = chr %in% c("chrX", "chrY"),
  has_probe_snp     = FALSE,
  has_cpg_snp       = FALSE,
  has_sbe_snp       = FALSE,
  is_cross_reactive = FALSE,
  is_multi_mapped   = FALSE,
  Distance_to_TSS = as.numeric(anno_27k$Distance_to_TSS),
)

# Gene annotations as functional flags based on Distance_to_TSS
manifest_27k <- manifest_27k %>%
  mutate(
    TSS200 = Distance_to_TSS >= -200 & Distance_to_TSS <= 200,
    TSS1500 = Distance_to_TSS >= -1500 & Distance_to_TSS < -200,
    gene_body = Distance_to_TSS > 200
  )

write_csv(manifest_27k, "illumina27k_annotation_hg19.csv")

# =====| Illumina 450K hg19 |===================================================

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

anno_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

manifest_450k <- tibble(
  probe_id      = rownames(anno_450k),
  chr           = as.character(anno_450k$chr),
  pos           = anno_450k$pos,
  strand        = anno_450k$strand,
  gene_symbol   = anno_450k$UCSC_RefGene_Name,
  is_sex_chr    = chr %in% c("chrX", "chrY"),
  has_cpg_snp   = !is.na(anno_450k$CpG_rs),
  has_sbe_snp   = !is.na(anno_450k$SBE_rs),
  has_probe_snp = !is.na(anno_450k$Probe_rs)
)

# Extract cross-reactive probes via maxprobes
xreactive_450k <- maxprobes::xreactive_probes(array_type = "450K")

manifest_450k <- manifest_450k %>%
  mutate(
    is_cross_reactive = probe_id %in% xreactive_450k
  )

# Extract multi-mapped probes from the bowtie annotation
multi_450k <- read_tsv(
  "HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt",
  col_names = "probe_id",
  show_col_types = FALSE
)

manifest_450k <- manifest_450k %>%
  mutate(
    is_multi_mapped = probe_id %in% multi_450k$probe_id
  )

# Gene annotations as functional flags based on Distance_to_TSS
manifest_450k <- manifest_450k %>%
  mutate(
    TSS200 = grepl("TSS200", anno_450k$UCSC_RefGene_Group, ignore.case = TRUE),
    TSS1500 = grepl("TSS1500", anno_450k$UCSC_RefGene_Group, ignore.case = TRUE),
    gene_body = grepl("Body", anno_450k$UCSC_RefGene_Group, ignore.case = TRUE),
  )

write_csv(
  manifest_450k,
  "illumina450k_annotation_hg19.csv"
)

# =====| Illumina EPIC (850K) hg19 |============================================

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

anno_epic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

manifest_epic <- tibble(
  probe_id      = rownames(anno_epic),
  chr           = as.character(anno_epic$chr),
  pos           = anno_epic$pos,
  strand        = anno_epic$strand,
  gene_symbol   = anno_epic$UCSC_RefGene_Name,
  is_sex_chr    = chr %in% c("chrX", "chrY"),
  has_cpg_snp   = !is.na(anno_epic$CpG_rs),
  has_sbe_snp   = !is.na(anno_epic$SBE_rs),
  has_probe_snp = !is.na(anno_epic$Probe_rs)
)

# Cross‑reactive probes via maxprobes for EPIC
xreactive_epic <- maxprobes::xreactive_probes(array_type = "EPIC")

manifest_epic <- manifest_epic %>%
  mutate(
    is_cross_reactive = probe_id %in% xreactive_epic,
    is_multi_mapped   = FALSE
  )

# Gene annotations as functional flags based on Distance_to_TSS
manifest_epic <- manifest_epic %>%
  mutate(
    TSS200 = grepl("TSS200", anno_epic$UCSC_RefGene_Group, ignore.case = TRUE),
    TSS1500 = grepl("TSS1500", anno_epic$UCSC_RefGene_Group, ignore.case = TRUE),
    gene_body = grepl("Body", anno_epic$UCSC_RefGene_Group, ignore.case = TRUE),
  )

write_csv(manifest_epic, "illuminaEPIC_annotation_hg19.csv")

# =====| Liftover hg38 |========================================================

# Fetch the liftover chain; unzip if necessary
# gunzip("raw/hg19ToHg38.over.chain.gz", 
#        destname = "raw/hg19ToHg38.over.chain", 
#        overwrite = TRUE)
chain <- import.chain("hg19ToHg38.over.chain")

# Initialize each hg19 annotation as a GRanges object
gr_hg19_27k <- GRanges(
  seqnames = manifest_27k$chr,
  ranges   = IRanges(start = manifest_27k$pos, end = manifest_27k$pos),
  strand   = manifest_27k$strand,
  probe_id = manifest_27k$probe_id
)

gr_hg19_450k <- GRanges(
  seqnames = manifest_450k$chr,
  ranges   = IRanges(start = manifest_450k$pos, end = manifest_450k$pos),
  strand   = manifest_450k$strand,
  probe_id = manifest_450k$probe_id
)

gr_hg19_epic <- GRanges(
  seqnames = manifest_epic$chr,
  ranges   = IRanges(start = manifest_epic$pos, end = manifest_epic$pos),
  strand   = manifest_epic$strand,
  probe_id = manifest_epic$probe_id
)

# Liftover each hg19 annotation to hg39
gr_hg38_27k <- liftOver(gr_hg19_27k, chain)
gr_hg38_450k <- liftOver(gr_hg19_450k, chain)
gr_hg38_epic <- liftOver(gr_hg19_epic, chain)

# Flatten each GRanges object
gr_hg38_27k <- unlist(gr_hg38_27k, use.names = TRUE)
gr_hg38_450k <- unlist(gr_hg38_450k, use.names = TRUE)
gr_hg38_epic <- unlist(gr_hg38_epic, use.names = TRUE)

# Finalize the hg38 manifests
manifest_hg38_27k <- manifest_27k %>%
  filter(probe_id %in% mcols(gr_hg38_27k)$probe_id) %>%
  mutate(
    chr    = as.character(GenomicRanges::seqnames(gr_hg38_27k)),
    pos    = GenomicRanges::start(gr_hg38_27k),
    strand = as.character(GenomicRanges::strand(gr_hg38_27k))
  )

manifest_hg38_450k <- manifest_450k %>%
  filter(probe_id %in% mcols(gr_hg38_450k)$probe_id) %>%
  mutate(
    chr    = as.character(GenomicRanges::seqnames(gr_hg38_450k)),
    pos    = GenomicRanges::start(gr_hg38_450k),
    strand = as.character(GenomicRanges::strand(gr_hg38_450k))
  )

manifest_hg38_epic <- manifest_epic %>%
  filter(probe_id %in% mcols(gr_hg38_epic)$probe_id) %>%
  mutate(
    chr    = as.character(GenomicRanges::seqnames(gr_hg38_epic)),
    pos    = GenomicRanges::start(gr_hg38_epic),
    strand = as.character(GenomicRanges::strand(gr_hg38_epic))
  )

# Save each liftover
write_csv(manifest_hg38_27k, "illumina27k_annotation_hg38.csv")
write_csv(manifest_hg38_450k, "illumina450k_annotation_hg38.csv")
write_csv(manifest_hg38_epic, "illuminaEPIC_annotation_hg38.csv")

rm(list = ls())
gc()

# END