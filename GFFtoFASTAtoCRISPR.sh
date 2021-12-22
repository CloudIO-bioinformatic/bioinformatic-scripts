#!/bin/bash

########################################################################
###   This is a simple tutorial of what commands I use to execute,   ###
###   filter better results and generate .fasta files from .gff,     ###
###   using CRISPRDetect software, awk and bedtools. I               ###
###   hope it helps you                                              ###
########################################################################
### This run with perl 5, version 26, subversion 1 (v5.26.1) built for x86_64-linux-gnu-thread-multi
### GNU Awk 4.1.4, API: 1.1 (GNU MPFR 4.0.1, GNU MP 6.1.2)
### bedtools v2.30.0

### this takes about 4 hours with all 126 genomes (Intel 5 7th Gen, 8GB RAM, SSD NvME M2)###

for genome in $(ls PATH/genomes/); # here you need to put your path where the downloaded genomes are
do echo perl PATH/CRISPRDetect.pl -f PATH/data/$genome/$genome* -o PATH/CRISPRDetectResult/$genome;done #here you need to create the CRISPRDetectResult folder in the same directory where you are executing the command

### for example in my case I ran it this way ###
#for genome in $(ls ../../Lihui/genomes/); I use this
#perl CRISPRDetect.pl -f /home/cquevedo/tesis/sara_cuadros/Lihui/genomes/$genome/$genome* -o /home/cquevedo/tesis/sara_cuadros/Lihui/CRISPRDetectResult/$genome;done I use this



### commands to separate and sort the resulting files from CRISPRDetect ###
mkdir gff  # create directory
mkdir fp  # create directory
mkdir 1 # create directory
mv *.gff gff/
mv *.fp fp/
mv G* 1/


### commands to filter the best scores (> = 4) in the CRISPRDetect GFF results, inside the CRISPRDetectResult directory ###
mkdir best_gff_repeatregion # create directory
mkdir best_gff_directrepeat # create directory
mkdir best_gff_bindingsite  # create directory
for genome in $(ls gff/*.gff);
do awk 'BEGIN{OFS=FS="="}($NF>=4){print $1"="$2"="$3"="$4"="$5"="$6$7"\t"$8}' $genome|awk 'BEGIN{OFS=FS="\t"}($3=="repeat_region"){print $0}' > best_gff_repeatregion/filtered_$genome;
awk 'BEGIN{OFS=FS="="}($NF>=4){print $1"="$2"="$3"="$4"="$5"="$6$7"\t"$8}' $genome|awk 'BEGIN{OFS=FS="\t"}($3=="direct_repeat"){print $0}' > best_gff_directrepeat/filtered_$genome;
awk 'BEGIN{OFS=FS="="}($NF>=4){print $1"="$2"="$3"="$4"="$5"="$6$7"\t"$8}' $genome|awk 'BEGIN{OFS=FS="\t"}($3=="binding_site"){print $0}' > best_gff_bindingsite/filtered_$genome;done

# $NF = last column
# $n = column
# BEGIN{OFS=FS="="} = indicates how the columns are separated

### commands to generate .fasta from filtered .gff, inside the CRISPRDetectResult directory ###
mkdir gfftofasta_repeatregion # create directory
mkdir gfftofasta_directrepeat # create directory
mkdir gfftofasta_bindingsite # create directory
for genome in $(ls PATH/data/);
do bedtools getfasta -fi PATH/data/$genome/*.fna -bed gff/best_gff_repeatregion/*$genome* > gfftofasta_repeatregion/$genome.fasta;
bedtools getfasta -fi PATH/data/$genome/*.fna -bed gff/best_gff_directrepeat/*$genome* > gfftofasta_directrepeat/$genome.fasta;
bedtools getfasta -fi PATH/data/$genome/*.fna -bed gff/best_gff_bindingsite/*$genome* > gfftofasta_bindingsite/$genome.fasta;done

### for example in my case I ran it this way ###
#for genome in $(ls /home/cquevedo/tesis/sara_cuadros/Lihui/data/); I use this
#do bedtools getfasta -fi /home/cquevedo/tesis/sara_cuadros/Lihui/data/$genome/*.fna -bed gff/best_gff_repeatregion/*$genome* > gfftofasta_repeatregion/$genome.fasta; I use this
#bedtools getfasta -fi /home/cquevedo/tesis/sara_cuadros/Lihui/data/$genome/*.fna -bed gff/best_gff_directrepeat/*$genome* > gfftofasta_directrepeat/$genome.fasta; I use this
#bedtools getfasta -fi /home/cquevedo/tesis/sara_cuadros/Lihui/data/$genome/*.fna -bed gff/best_gff_bindingsite/*$genome* > gfftofasta_bindingsite/$genome.fasta;done I use this
