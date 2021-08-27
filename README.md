# Proportion of exon regions
A little study using R and Bioconductor of the proportion of exon regions among chromosomes in mouse:  

- exon-covering-chr1-R.R  
- exon-covering-genome-R.R
- exons-covering-genome-nucleotide-R.R

Also, you can find two shell scripts for computing the proportion of exon regions among chromosome 1 and among the whole genome in mouse from a .gtf file:  

- exon-covering-chr1.sh
- exon-covering-genome.sh

For the shell scripts the data was got with the following command:  
  
    # Download the file
    $ curl http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus/Mus_musculus.GRCm39.104.gtf.gz\
    > Mus_musculus.GRCm39.104.gtf.gz
    
    # Unzip the file
    gunzip Mus_musculus.GRCm39.104.gtf.gz
    
    # Get the length of each chromosome
    curl http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz | \
    bioawk -c fastx '{print $name"\t"length($seq)}' > Mus_musculus.GRCm39_genome.txt 
    

  
