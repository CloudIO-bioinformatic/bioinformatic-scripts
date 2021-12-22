#!/bin/bash

## Download from: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP230751&o=acc_s%3Aa


wget https://sra-download.ncbi.nlm.nih.gov/traces/sra57/SRR/010323/SRR10571719
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra59/SRR/010323/SRR10571721
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra25/SRR/010323/SRR10571737
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra17/SRR/010323/SRR10571755
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra25/SRR/010323/SRR10571756





for run in SRR10571719 SRR10571721 SRR10571737 SRR10571755 SRR10571756
	do echo $run
	#fastq-dump -X 5 -I --split-files $run 
	fastq-dump -I --split-files $run
done




## Download from: https://www.ensembl.org/Homo_sapiens/Info/Index
## wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/

## Download from: https://www.ncbi.nlm.nih.gov/genome/?term=human
## gtf
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz
## fna
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
gunzip *.gz



## Extract exons and splicesites

hisat2_extract_exons.py GCF_000001405.39_GRCh38.p13_genomic.gtf > exons
hisat2_extract_splice_sites.py GCF_000001405.39_GRCh38.p13_genomic.gtf > splicesites


## Creating index
hisat2-build -p 4 --ss splicesites --exon exons GCF_000001405.39_GRCh38.p13_genomic.fna index_hisat2



## Alignment

for run in SRR10571719 SRR10571721 SRR10571737 SRR10571755 SRR10571756
	do hisat2 -x index_hisat2 -1 $run'_1.fastq' -2 $run'_2.fastq' | samtools sort -O BAM -o $run'_align.bam'
	samtools index $run'_align.bam'
	## Rsubread arreglado
	r -e 'library(Rsubread);featurecounts_output<-featureCounts("'$run'_align.bam",annot.ext="GCF_000001405.39_GRCh38.p13_genomic.gtf", isGTFAnnotationFile = TRUE);write.table(featurecounts_output$counts,file=paste("'$run'_align.bam","counts",sep="."),quote=FALSE,sep="\t");'
done



rm *.fastq
rm *.bam
rm *.ht2
rm *sites
rm exons
rm *rf
mkdir counts
mv *counts counts	
