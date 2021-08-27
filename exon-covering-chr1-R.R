# This script computes the proportion for a single chromosome from 
# Mus musculus (version 10) covered by exons

# Load libraries
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(rtracklayer)
library(GenomicFeatures)
library(ggplot2)

# Load data
txdb <- TxDb.Mmusculus.UCSC.mm39.refGene

# Get the exons
exons <- exons(txdb)

# Select only ranges from chromosome 1
chr1_exons <- exons[seqnames(exons) == "chr1"]

# Merge exons ranges
chr1_exons_collapsed <- reduce(chr1_exons, ignore.strand=TRUE)

# Total exons width
exons_total_width <- sum(width(chr1_exons_collapsed))

# Total chromosome width
chr1_length <- seqlengths(txdb)["chr1"]

# Compute proportion of chromosome 1 covered by exons
exons_total_width/chr1_length
exons_total_width
chr1_length