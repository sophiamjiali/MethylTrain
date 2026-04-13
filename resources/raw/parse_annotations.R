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
  probe_type        = anno_27k$Type,
  gene_symbol       = anno_27k$Symbol,
  is_sex_chr        = chr %in% c("chrX", "chrY"),
  has_probe_snp     = FALSE,
  has_cpg_snp       = FALSE,
  has_sbe_snp       = FALSE,
  is_cross_reactive = FALSE,
  is_multi_mapped   = FALSE
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
  probe_type    = anno_450k$Type,
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
  "/Volumes/FBI_Drive/MethylCDM-project/resources/raw/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt",
  col_names = "probe_id",
  show_col_types = FALSE
)

manifest_450k <- manifest_450k %>%
  mutate(
    is_multi_mapped = probe_id %in% multi_450k$probe_id
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
  probe_type    = anno_epic$Type,
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

write_csv(manifest_epic, "illuminaEPIC_annotation_hg19.csv")

# =====| Liftover hg38 |========================================================

# Fetch the liftover chain
# chain <- import.chain(gunzip("hg19ToHg38.over.chain.gz"))
chain <- import.chain("hg19ToHg38.over.chain")

# --- Helper: create GRanges ---
make_gr <- function(df) {
  
  df <- df %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      chr      = as.character(chr),
      pos      = as.integer(pos),
      strand   = as.character(strand),
      probe_id = as.character(probe_id)
    )
  
  chr <- df$chr
  pos <- df$pos
  strand <- df$strand
  probe_id <- df$probe_id
  
  stopifnot(
    length(chr) == length(pos),
    length(pos) == length(strand),
    length(strand) == length(probe_id)
  )
  
  GRanges(
    seqnames = chr,
    ranges   = IRanges(start = pos, end = pos),
    strand   = strand,
    probe_id = probe_id
  )
}

# --- Helper: liftover + clean + convert to tibble ---
liftover_to_df <- function(gr, chain) {
  gr_hg38 <- liftOver(gr, chain)
  gr_hg38 <- unlist(gr_hg38, use.names = FALSE)
  
  # Standardize chromosome style (important)
  GenomeInfoDb::seqlevelsStyle(gr_hg38) <- "UCSC"
  
  df <- tibble::tibble(
    probe_id = mcols(gr_hg38)$probe_id,
    chr      = as.character(GenomicRanges::seqnames(gr_hg38)),
    pos      = GenomicRanges::start(gr_hg38),
    strand   = as.character(GenomicRanges::strand(gr_hg38))
  )
  
  # Keep only uniquely mapped probes
  df <- df %>%
    dplyr::group_by(probe_id) %>%
    dplyr::filter(dplyr::n() == 1) %>%
    dplyr::ungroup()
  
  return(df)
}

# --- Helper: apply liftover + join back ---
liftover_manifest <- function(manifest_df) {
  gr <- make_gr(manifest_df)
  df_hg38 <- liftover_to_df(gr, chain)
  
  manifest_hg38 <- manifest_df %>%
    dplyr::inner_join(df_hg38, by = "probe_id")
  
  return(manifest_hg38)
}

# --- Apply to each platform ---
manifest_hg38_27k  <- liftover_manifest(manifest_27k)
manifest_hg38_450k <- liftover_manifest(manifest_450k)
manifest_hg38_epic <- liftover_manifest(manifest_epic)

# Save each liftover
write_csv(manifest_hg38_27k, "illumina27k_annotation_hg38.csv")
write_csv(manifest_hg38_450k, "illumina450k_annotation_hg38.csv")
write_csv(manifest_hg38_epic, "illuminaEPIC_annotation_hg38.csv")

# =====| Liftover full annotation |=============================================

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(GenomicRanges)

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

gr_hg19 <- GRanges(
  seqnames = anno$chr,
  ranges = IRanges(anno$pos, anno$pos),
  strand = "*",
  probe_id = rownames(anno)
)

library(rtracklayer)

chain <- import.chain("hg19ToHg38.over.chain")
gr_list <- liftOver(gr_hg19, chain)

mapped <- elementNROWS(gr_list) == 1
gr_hg38 <- unlist(gr_list[mapped])
names(gr_hg38) <- mcols(gr_hg19)$probe_id[mapped]

library(AnnotationHub)

ah <- AnnotationHub()
query(ah, c("CpG", "hg38"))

cpg_islands <- ah[["AH5086"]]   # CpG islands GRanges

library(GenomicRanges)

islands <- cpg_islands

shores <- flank(islands, 2000, both=TRUE)
shores <- setdiff(shores, islands)

shelves <- flank(islands, 4000, both=TRUE)
shelves <- setdiff(shelves, c(islands, shores))

context <- rep("OpenSea", length(gr_hg38))

context[queryHits(findOverlaps(gr_hg38, islands))] <- "Island"
context[queryHits(findOverlaps(gr_hg38, shores))]  <- "Shore"
context[queryHits(findOverlaps(gr_hg38, shelves))] <- "Shelf"

probe_context <- context
names(probe_context) <- names(gr_hg38)

probe_context

# ==============================================================================
rm(list = ls())
gc()

# END