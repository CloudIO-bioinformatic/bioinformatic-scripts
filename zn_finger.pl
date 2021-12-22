#!/usr/bin/perl
		
use	strict;
use	warnings;
use	Bio::SeqIO;
my $file =	shift or die "$!\n";
my $seqIOobj =	new	Bio::SeqIO(-file=>"$file",-format=>"fasta");
while (my $seqObj =	$seqIOobj->next_seq){
								if ($seqObj->seq =~	m/C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H/){
																print ">".$seqObj->id."	".$seqObj->desc."\n".$seqObj->seq."\n";
								}
}
