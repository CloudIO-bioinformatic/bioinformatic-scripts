from Bio.Seq import Seq
from Bio import SeqIO
import warnings
import re
import sys
import os
from os import remove
from os import path

def parse_fasta(file):

    return SeqIO.parse(file, "fasta")

def calculate_orf_without_partial_codon(seq,id):

    warnings.filterwarnings("ignore")

    table = 1

    min_pro_len = 100

    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:

        file = open("TSS_pos.txt", "a")

        for frame in range(3):

            for pro in nuc[frame:].translate(table).split("*"):

                if len(pro) >= min_pro_len:

                    first_Met = pro.find('M')

                    if first_Met>-1:

                        print ("Gene: %s - %s...[ %s ]...%s - Length: %i, Strand: %i, Frame: %i, M-pos: %i" % (id,pro[:10],pro[first_Met:first_Met+1], pro[-10:], len(pro), strand, frame,first_Met))

                        file.write(id+"\t"+str(first_Met*3)+"\n")

                        return
        file.close()


def main():

    if(len(sys.argv)==1):
        print("\n***\n")
        print("ERROR: missing .fasta file")

        return

    print ("Looking for TSS's in fasta: ",sys.argv[1])
    print("\n***\n")
    fasta_file=sys.argv[1]

    if path.exists("TSS_pos.txt"):
        remove("TSS_pos.txt")

    data = []

    for record in parse_fasta(fasta_file):

        calculate_orf_without_partial_codon(record.seq,record.id)

    print("\n***\n")
    print("An output has been generated with the TSS positions of each sequence: TSS_pos.txt")

if __name__ == "__main__":

    main()
