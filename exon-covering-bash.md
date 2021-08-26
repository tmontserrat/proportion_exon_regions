# Percentage of genome covered by exons not in pseudogenes

## Preparing files

We will need the GFF file for mm10/GR38.75 and the txt file with the sequence information.  

- *Mus_musculus.GRCm38.75_chr1.gtf.gz*  
- *Mus_musculus.GRCm38_genome.txt*

Both files comes from <https://github.com/vsbuffalo/bds-files/tree/master/chapter-09-working-with-range-data> and, as the README file in this repository says, you can get the data following this instructions:  

- The first file comes from <http://ftp.ensembl.org/pub/release-75/gtf/mus_musculus/> and then we have selected those registers of the chromosome 1:  

        # Ref.: https://github.com/vsbuffalo/bds-files/tree/master/chapter-09-working-with-range-data
        gzcat Mus_musculus.GRCm38.75.gtf.gz | egrep "^(1\t|#)" | gzip > Mus_musculus.GRCm38.75_chr1.gtf.gf

- The second file comes from:  

        # Ref.: https://github.com/vsbuffalo/bds-files/tree/master/chapter-09-working-with-range-data
        curl ftp://ftp.ensembl.org/pub/release-75/fasta/mus_musculus/dna/Mus_musculus.GRCm38.75.dna_rm.toplevel.fa.gz \
        bioawk -c fastx '{print $name"\t"length($seq)}' > Mus_musculus.GRCm38_genome.txt

## Unzip GFF files

First of all, we have to unzip our GFF file with the features of mm10:  

    $ gunzip Mus_musculus.GRCm38.75_chr1.gtf.gz

## Extract exons ranges

Now, we will extract the exon ranges (we are not incorporating pseudogenes):  

    $ bioawk -c gff '$feature ~ /exon/ {print $seqname,$start,$end,$score,$filter,$9}' Mus_musculus.GRCm38.75_chr1.gtf | grep -v "pseudogene" | cut -f1,2,3,4,5 | head
    1	3102016	3102125	.	+
    1	3213609	3216344	.	-
    1	3205901	3207317	.	-
    1	3213439	3215632	.	-
    1	3206523	3207317	.	-
    1	3670552	3671498	.	-
    1	3421702	3421901	.	-
    1	3214482	3216968	.	-
    1	3466587	3466687	.	+
    1	3513405	3513553	.	+

    $ bioawk -c gff '$feature ~ /exon/ {print $seqname,$start,$end,$score,$filter,$9}' Mus_musculus.GRCm38.75_chr1.gtf | grep -v "pseudogene" | cut -f1,2,3,4,5 > exons.bed

## Merge exons that overlap

We first have to sort the file:  

    $ sort -k1,1 -k2,2n exons.bed > exons-sorted.bed

    $ bedtools merge -i exons-sorted.bed | head
    1	3102016	3102125
    1	3205901	3207317
    1	3213439	3216968
    1	3421702	3421901
    1	3466587	3466687
    1	3513405	3513553
    1	3670552	3671498
    1	3783876	3783933
    1	4343507	4350091
    1	4351910	4352081

    $ bedtools merge -i exons-sorted.bed > exons-merged.bed

## Compute the width of all the merged exons

We can do that using `bioawk`:  

    $ bioawk -c bed 'BEGIN{ s = 0}; {s += ($3-$2)}; END{ print "width: " s }' exons-merged.bed
    width: 5410233

We can save it in a variable to use later:  

    $ exons_width=$(bioawk -c bed 'BEGIN{ s = 0 }; {s += ($3-$2)}; END{ print s }' exons-merged.bed)

    $ echo $exons_width
    5410233


## Width of chromosome 1 genome

We have the information in the *Mus_musculus.GRCm38_genome.txt* file:  

    $ head Mus_musculus.GRCm38_genome.txt
    1	195471971
    10	130694993
    11	122082543
    12	120129022
    13	120421639
    14	124902244
    15	104043685
    16	98207768
    17	94987271
    18	90702639

    $ cut -f2 Mus_musculus.GRCm38_genome.txt | head -n1
    195471971

Again, we can save it as a variable:  

    $ chr1_width=$(cut -f2 Mus_musculus.GRCm38_genome.txt | head -n1)

    $ echo $chr1_width
    195471971

## Computing the proportion of the genome covered by exons not in pseudogenes

    $ echo "$exons_width/$chr1_width" | bc -l
    .02767779427568160143

# Percentage of the whole genome covered by exons not in pseudogenes

# Preparing files

We will need the GFF file for mm10/GR38.75 and the txt file with the sequence information.  

- *Mus_musculus.GRCm39.104.gtf.gz*  
- *Mus_musculus.GRCm39_genome.txt*

The first file comes from 
    $ curl http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus/Mus_musculus.GRCm39.104.gtf.gz \
    > Mus_musculus.GRCm39.104.gtf.gz

The second file comes from:  

    $ curl http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz | bioawk -c fastx '{print $name"\t"length($seq)}' > Mus_musculus.GRCm39_genome.txt  

## Extracting exons ranges

We select all exons not belonging to pseudogenes:  

    $ bioawk -c gff '$feature ~ /exon/ {print $seqname,$start,$end,$score,$filter,$9}' Mus_musculus.GRCm39.104.gtf | grep -v "pseudogene" | cut -f1,2,3,4,5 | head -n2
    1	150956201	150958296	.	+
    1	108344807	108347562	.	+

    $ bioawk -c gff '$feature ~ /exon/ {print $seqname,$start,$end,$score,$filter,$9}' Mus_musculus.GRCm39.104.gtf | grep -v "pseudogene" | cut -f1,2,3,4,5 > exons-genome.bed

Note that now we are working in the **mm_39** genome version. 

## Sorting and mergin exons ranges

First, we sort by chromosome and by start position:   

    $ sort -k1,1 -k2,2n exons-genome.bed > exons-genome-sorted.bed

Then, we can merge the exons ranges:  

    $ bedtools merge -i exons-genome-sorted.bed > exons-genome-merged.bed

## Filtering for chromosomes

We select only those registers about chromosome, not scaffolds or other sequences:  

    $ bioawk '$1 ~ "^[1-9XYM]" {print $0}' exons-genome-merged.bed > exons-genome-merged-chr.bed

Now we'll do the same but for sequences lengths:  

    $ grep -E "^[1-9XYM]T*[^A-Za-z]" Mus_musculus.GRCm39_genome.txt > mm_genome-chr.txt

## Computing the total width of exons not in a pseudogene for the whole genome  

    $ bioawk -c bed 'BEGIN{ s = 0 }; {s += ($end-$start)}; END{ print "width: " s}' exons-genome-merged-chr.bed
    width: 117833254

We can save it in a variable to use later:  

    $ exons_width=$(bioawk -c bed 'BEGIN{ s = 0 }; {s += ($end-$start)}; END{ print s}' exons-genome-merged-chr.bed)

## Computing the total width for the whole genome

    $ awk 'BEGIN{s = 0 }; {s += $2}; END{print "width: " s}' mm_genome-chr.txt
    width: 2.72343e+09

Again, we can save it in a variable for later use:  

    $ genome_width=$(awk 'BEGIN{s = 0 }; {s += $2}; END{print s}' mm_genome-chr.txt)

And we change the scientific notation format for another that `bc` can understant:  

    $ echo "$genome_width" | sed -e 's/e+/\*10\^/'
    2.72343*10^09

    $ genome_width=$(echo "$genome_width" | sed -e 's/e+/\*10\^/')

    $ echo "$genome_width" | bc
    2723430000.00000

    $ genome_width=$(echo "$genome_width" | bc)

## Compute the result

    $ echo "$exons_width/$genome_width" | bc -l
    .04326648894959664834