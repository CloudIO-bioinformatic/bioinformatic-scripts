#!/bin/bash

##########################################################################################
###   This is a simple tutorial of what commands I use to install and execute blast.   ###
###   I hope it helps you.                                                             ###
##########################################################################################
### This run with BLAST 2.11.0+ built for x86_64-linux-gnu-thread-multi
### GNU Awk 4.1.4, API: 1.1 (GNU MPFR 4.0.1, GNU MP 6.1.2)
### GNU Wget 1.19.4 built on linux-gnu.

# Step 0, download all 126 genomes.
./datasets download genome taxon 28889 --filename PATH-genome/


# The first step is to download it, you can download it from here (if I remember correctly you use MAC)
# >  https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
# Inside the "/bins" folder you will find the executables

# The second step is to download the CRISPR amino acid sequence database (in PADS, https://bigd.big.ac.cn/padsarsenal/download.php)
wget https://bigd.big.ac.cn/padsarsenal/download/sequence/PADS_Arsenal_Archaea_CRISPR-CAS_v1_2019.09.09.faa.tar.gz
wget https://bigd.big.ac.cn/padsarsenal/download/sequence/PADS_Arsenal_Bacteria_CRISPR-CAS_v1_2019.09.09.faa.tar.gz


# Then use "makeblastb" to create a database. If you're wondering why we should pass a fasta to the database, it's because alignment makes it so much faster.
PATH-ncbi-blast/makeblastdb -in PATH-pads/PADS_Arsenal_Archaea_CRISPR-CAS_v1_2019.09.09.faa \
-dbtype prot -max_file_sz 4080218932 \
-out PATH-pads/pads.db
# PATH-ncbi-blast: path where you have downloaded BLAST.
# PATH-pads: path where you have the PADS files.

# Finally, perform out the alignments and save the information in tabular format (without header).
for genome in $(ls PATH-genome/genomes/);
do blastx -query PATH-genome/genomes/$genome/*.fna \
 -outfmt 6 -num_threads 4 -db PATH-pads/pads.db > PATH-result/$genome.blastx;
 echo $genome;echo "Done!";done
# PATH-genome: path where you have the 126 genome files.
# PATH-result:path where you want to save the result files.
# blastx: search protein databases using a translated nucleotide query. (Genomes = nucleotide query and databases(PADS) = protein database).
# -outfmt 6: tabular without header.
# -num_threads 4: you can change it!.
