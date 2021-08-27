#!/bin/bash

# This script computes the proportion of a single chromosome from 
# Mus musculus (version 39) covered by exons

# Extract all records from chromosome 1
grep -P "^1\t|#" ../data/Mus_musculus.GRCm39.104.gtf > ../data/Mus_musculus.GRCm39.104_chr1.gtf

# Extract exons
bioawk -c gff '$feature ~ /exon/ {print $0}' ../data/Mus_musculus.GRCm39.104_chr1.gtf > ../data/exons.gtf

# Convert from gff to bed
sortBed -i ../data/exons.gtf | gff2bed > ../data/exons.bed

# Select chromosome, start and end
cut -f1,2,3 ../data/exons.bed > ../data/exons-ranges.bed

# Sort and merge exons ranges
sort -k1,1 -k2,2n ../data/exons-ranges.bed > ../data/exons-ranges-sorted.bed
bedtools merge -i ../data/exons-ranges-sorted.bed > ../data/exons-merged.bed

# Compute total width of exons and save in a variable
exons_width_chr1=$(bioawk -c bed 'BEGIN{ s = 0 }; { s += $end-$start }; END{ print s }' ../data/exons-merged.bed)
echo "$exons_width_chr1"


# Compute total width for chromosome 1
chr1_width=$(cut -f2 ../data/Mus_musculus.GRCm39_genome.txt | head -n1)
echo "$chr1_width"

# Compute exon proportion for chromosome 1
echo "$exons_width_chr1/$chr1_width" | bc -l

# Remove all intermediate files
rm ../data/Mus_musculus.GRCm39.104_chr1.gtf ../data/exons.gtf ../data/exons.bed
rm ../data/exons-ranges.bed ../data/exons-ranges-sorted.bed ../data/exons-merged.bed