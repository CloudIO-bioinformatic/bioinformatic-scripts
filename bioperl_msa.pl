#!usr/bin/perl

use strict;
use warnings;

use Bio::AlignIO;
my $file = shift;
my $format = shift;

my $alignment_object = new Bio::AlignIO(-file=>"$file",-format=>"$format");

my $alignment = $alignment_object->next_aln;

#print $alignment->num_sequences."\n";
#print $alignment->length."\n";
#si tiene o o no gaps
#print $alignment->is_flush."\n";
my $all="all <- c(";
my $count = 0;

# imprime las secuencias de un MSA, para procesarlals en un ggseqlogo

#para no mostrar el alineamiento completo, es recomendable cortar
#my $column1 = $alignment->column_from_residue_number("sp|P23246|SFPQ_HUMAN",442);
#my $column2 = $alignment->column_from_residue_number("sp|P23246|SFPQ_HUMAN",538);
#my $subaligment = $alignment->slice($column1,$column2,1);
##y se cambia el for por $subaligment


foreach my $sequence ($alignment->each_seq){
  #print $sequence->seq."\n";
  my $seq = $sequence->seq;
  print "seq$count <- '".$seq."'\n";
  $all= $all."seq$count,";
  $count = $count +1;
}
$all=$all.")";

print $all."\n";
