#!/bin/bash

# This script computes the proportion of the genome from 
# Mus musculus (version 39) covered by exons

# Extract exons
bioawk -c gff '$feature ~ /exon/ {print $0}' ../data/Mus_musculus.GRCm39.104.gtf > ../data/exons.gtf

# Convert from gff to bed
sortBed -i ../data/exons.gtf | gff2bed > ../data/exons.bed

# Select chromosome, start and end
cut -f1,2,3 ../data/exons.bed > ../data/exons-ranges.bed

# Sort and merge exons ranges
sort -k1,1 -k2,2n ../data/exons-ranges.bed > ../data/exons-ranges-sorted.bed
bedtools merge -i ../data/exons-ranges-sorted.bed > ../data/exons-merged.bed

# Filter chromosomes from exons and genome
bioawk '$1 ~ "^[1-9XYM]" {print $0}' ../data/exons-merged.bed > ../data/exons-merged-chr.bed
grep -E "^[1-9XYM]T*[^A-Za-z]" ../data/Mus_musculus.GRCm39_genome.txt > ../data/mm_genome-chr.txt

# Compute total width of exons and save in a variable
exons_width=$(bioawk -c bed 'BEGIN{ s = 0 }; { s += ($end-$start) }; END{ print s }' ../data/exons-merged-chr.bed)
echo "$exons_width"
# Compute total width for the whole genome
genome_width=$(bioawk 'BEGIN{ s = 0 }; { s += $2 }; END{ print s }' ../data/mm_genome-chr.txt)

# Compute exon proportion of the genome
echo "$exons_width/$genome_width" | bc -l

# Remove all intermediate files
rm ../data/exons.gtf ../data/exons.bed ../data/exons-ranges.bed ../data/exons-ranges-sorted.bed 
rm ../data/exons-merged.bed ../data/exons-merged-chr.bed ../data/mm_genome-chr.txt