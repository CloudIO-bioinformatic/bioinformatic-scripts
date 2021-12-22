#!/usr/bin/perl
use strict;
use warnings;
use Bio::SearchIO;
use Bio::SeqIO;
use Bio::Seq;
#Escriba un script en BioPerl, que reciba como entrada un largo L, y genere un genoma “al azar” de
#ese largo, en formato FASTA. Utilice su script para generar 4 genomas de largos 50.000, 100.000,
#1.000.000 y 3.000.000.
my $length = shift;
my $i = 1;
my $genoma = "";
if ($length == 50000){}elsif ($length == 100000){}elsif ($length == 1000000){}elsif ($length == 3000000){}else{print "error!\n usa 50000 o 100000 o 1000000 o 3000000\n";exit;}
while($i<=$length){
  #print $i."\n";
  my $random_number = int(rand(4))+1;
  #print $random_number."\n";
  if ( $random_number == 1 ){
    $genoma = $genoma."A";
  }
  elsif ( $random_number == 2 ){
    $genoma = $genoma."C";
  }
  elsif ( $random_number == 3 ){
    $genoma = $genoma."T";
  }
  elsif ( $random_number == 4 ){
    $genoma = $genoma."G";
  }

  $i++;
}
#print "el largo es de : ".$length."\n";

if ($length == 50000){
  my $seq_obj = Bio::Seq->new(-seq=>"$genoma",-display_id => "genoma50k",-desc => "Genoma al azar con 50000 nucleótidos",-alphabet => "dna" );
  my $seqio_obj = Bio::SeqIO->new(-file => '>genoma50k.fasta',-format => 'fasta' );
  $seqio_obj->write_seq($seq_obj);
}
elsif ($length == 100000){
  my $seq_obj = Bio::Seq->new(-seq=>"$genoma",-display_id => "genoma100k",-desc => "Genoma al azar con 100000 nucleótidos",-alphabet => "dna" );
  my $seqio_obj = Bio::SeqIO->new(-file => '>genoma100k.fasta',-format => 'fasta' );
  $seqio_obj->write_seq($seq_obj);
}
elsif ($length == 1000000){
  my $seq_obj = Bio::Seq->new(-seq=>"$genoma",-display_id => "genoma1m",-desc => "Genoma al azar con 1000000 nucleótidos",-alphabet => "dna" );
  my $seqio_obj = Bio::SeqIO->new(-file => '>genoma1m.fasta',-format => 'fasta' );
  $seqio_obj->write_seq($seq_obj);
}
elsif ($length == 3000000){
  my $seq_obj = Bio::Seq->new(-seq=>"$genoma",-display_id => "genoma3m",-desc => "Genoma al azar con 3000000 nucleótidos",-alphabet => "dna" );
  my $seqio_obj = Bio::SeqIO->new(-file => '>genoma3m.fasta',-format => 'fasta' );
  $seqio_obj->write_seq($seq_obj);
}else{
  print "Largos incorrectos!, deben ser 50000, 100000, 1000000 o 3000000";
}
