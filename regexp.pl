#/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
my $file = shift;
my $object_seqIO = new Bio::SeqIO(-file=>"$file",-format=>"fasta",-alphabet=>"protein");
while ( my $object_secuencia = $object_seqIO->next_seq()){
  my $secuencia = $object_secuencia->seq;
  #if ($secuencia =~ m/(D.{86,107}D.{35,43}E)/){
  #comentar el if de arriba y descomentar el if de abajo para efectuar la busqueda de dominio DDD
  if ($secuencia =~ m/(D.{86,107}D.{35,43}D)/){
    print $object_secuencia->id."\t";
    print "p\n";
  }else{
    print $object_secuencia->id."\t";
    print "n\n";
  }
}
#my $object_seqIO2 = new Bio::SeqIO(-file=>"$file",-format=>"fasta",-alphabet=>"protein");
#while ( my $object_secuencia2 = $object_seqIO2->next_seq()){
#  my $secuencia = $object_secuencia2->seq;
#  if ($secuencia =~ m/(D.{86,107}D.{35,43}D)/){
#    print $object_secuencia2->id."\t";
#    print "p\tDDD\n";
#  }else{
#    print $object_secuencia2->id."\t";
#    print "n\tDDD\n";
#  }
#}
