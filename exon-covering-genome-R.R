# This script computes the proportion for each chromosome from 
# Mus musculus (version 10) from EnsGenecovered by exons and make a bar plot

# Load libraries
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(rtracklayer)
library(GenomicFeatures)
library(ggplot2)

# Load data
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene

# Get the exons
exons <- exons(txdb)

# Split our exons by chromosome
exons_split_chr <- split(exons, seqnames(exons))

# Extract the names of the sequences
chr_names <- names(exons_split_chr)

# Find those sequences names not belonging to a full chromosome
chr_index <- grep(pattern="_",
                  x=chr_names)

# Select chromosome names
chr_names <- chr_names[-chr_index]

# Subset exons by chromosome names
exons_split_chr <- exons_split_chr[names(exons_split_chr) == chr_names] 

# Merge exons ranges
exons_split_chr_collapsed <- lapply(exons_split_chr, reduce, ignore.strand=TRUE)

# Compute exons width for each chromosome
exons_split_chr_length <- sapply(exons_split_chr_collapsed, 
                                 function(x) sum(width(x)))

# Create ranges for chromosomes
chr_ranges <- GRanges(seqnames=chr_names,
                      IRanges(start=1, end=seqlengths(exons)[names(seqlengths(exons))==chr_names]))

# Split these ranges by chromosome
chr_ranges_split <- split(chr_ranges, seqnames(chr_ranges))

# Calculate lengths for each chromosome
chr_length <- sapply(chr_ranges_split, function(x) sum(width(x)))

# Calculate exon proportions for each chromosome
exon_proportions <- exons_split_chr_length/chr_length
exon_proportions

# VISUALIZATION OF RESULTS
# We change the names of our chromosomes to make the visualization simpler
names(exon_proportions) <- sub(pattern="chr([1-9]*|X|Y)", 
                               replacement="\\1", 
                               x=names(exon_proportions))

# Create the dataframe for using ggplot
chr_exon_prop_df <- data.frame(chr=names(exon_proportions), 
                               prop=exon_proportions)

# Remove row names
rownames(chr_exon_prop_df) = NULL

# Create the correct order of the chromosomes
chr_order<-c(paste(1:19,sep=""),"X","Y","M")
chr_exon_prop_df$chr <- factor(chr_exon_prop_df$chr, levels=chr_order)

# Create the bar plot
ggplot(data=chr_exon_prop_df, aes(x=chr, y=prop)) +
  geom_bar(stat="identity", color="black", fill="dodgerblue") +
  labs(x="Chromosomes", y="Proportion of exonic regions") +
  ggtitle("Exonic regions mouse version 10 EnsGene") +
  ylim(c(0, 1))

# Remove mitochondrial chromosome for a new visualization
# Remove mitochondrial chromosome while creating the dataframe for ggplot
chr_exon_prop_df_noM <- data.frame(chr=names(exon_proportions)[names(exon_proportions) !="M"], 
                                   prop=exon_proportions[names(exon_proportions) != "M"])

# Remove row names
rownames(chr_exon_prop_df) = NULL

# Create the correct order of the chromosomes
chr_order<-c(paste(1:19,sep=""),"X","Y")
chr_exon_prop_df_noM$chr <- factor(chr_exon_prop_df_noM$chr, levels=chr_order)

# Create the bar plot
ggplot(data=chr_exon_prop_df_noM, aes(x=chr, y=prop)) +
  geom_bar(stat="identity", color="black", fill="dodgerblue") +
  labs(x="Chromosomes", y="Proportion of exonic regions") +
  ggtitle("Exonic regions mouse version 10 EnsGene (-MT)") +
  ylim(c(0, 0.06))