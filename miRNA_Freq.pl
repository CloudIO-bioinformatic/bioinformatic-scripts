#/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
my $file = shift;
my $file2 = shift;

my $genoma = `cat $file2`;

my $object_seqIO = new Bio::SeqIO(-file=>"$file",-format=>"fasta",-alphabet=>"dna");
my $cont = 1;
while ( my $object_secuencia = $object_seqIO->next_seq()){
    my $miRNA = $object_secuencia->seq;
    if (my $freq = $genoma=~ m/$miRNA/){
      print "miRNA ".$cont." fue encontrado : ".$freq." veces en el genoma\n";
    }else{
      print "miRNA ".$cont." no fue encontrado con RegExp en el genoma\n";
    }
    $cont++;
}
