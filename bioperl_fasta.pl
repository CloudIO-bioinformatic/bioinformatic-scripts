#/usr/bin/perl
use strict;
use warnings;
#revolver un arreglo
use List::Util qw(shuffle);
use Bio::Seq;
use Bio::SeqIO;
my $file = shift;
my $object_seqIO = new Bio::SeqIO(-file=>"$file",-format=>"fasta",-alphabet=>"protein");
my @secuencia;
while ( my $object_secuencia = $object_seqIO->next_seq()){
  my $secuencia = $object_secuencia->seq;
  push @secuencia,$secuencia;
}
my $num = @secuencia;
for (my $i=0;$i<$num;$i++){
  print $secuencia[$i]."\n";
}

#hacer un msa

system("./clustalo -i $file -o $file".".msa");
