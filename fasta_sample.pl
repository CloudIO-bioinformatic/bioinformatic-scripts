#!/usr/bin/perl
		
use	strict;
use	warnings;
use	Bio::SeqIO;
use	Bio::DB::Fasta;

my $file =	shift or die "$!\n";
my $n = shift or die "$!\n";
my $BioDBobj =	new	Bio::DB::Fasta("$file");
my @idArray =	$BioDBobj->get_all_primary_ids;
while(	$n>0 ){
								my $selected_seq_index =	int(rand(scalar(@idArray)));
								my $seqObj =	$BioDBobj->get_Seq_by_id($idArray[$selected_seq_index]);
								print ">".$seqObj->id."	".$seqObj->desc."\n".$seqObj->seq."\n";
								splice(@idArray,$selected_seq_index,1);
								$n--;
}
