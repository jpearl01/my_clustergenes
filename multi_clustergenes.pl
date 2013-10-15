#!/usr/bin/perl

use strict;
use warnings;
#Currently have to change "genes.fna" to your current nucleotide gene annotations file (and location if you are running this anywhere but 
#the directory which contains that file) "genes_vs_peptides.tfasty" to your tfasty alignment outputs, "dnaGenes_vs_allContigs.fasta" to 
#your fasta contig/na alignment output.


my $fastaSeqs	= 'genes.fna';
my $tfasty	= 'genes_vs_peptides.tfasty';
my $allFasta	= 'dnaGenes_vs_allContigs.fasta';
my $aa_iden 	= .5;
my $nuc_iden	= .5;
my $ratio	= .5;
my $cg 		= 'clusterGenes';
my $cr		= 'report_output';


for my $i (0..5){
	for my $j (0..5){
		my $aa_arg = $aa_iden+($i/10);
		my $nuc_arg= $nuc_iden+($i/10);
		my $ratio_arg =$ratio+($j/10);
		my $out = "clusterGenes_".($aa_arg*100)."_".($nuc_arg*100)."_".($ratio_arg*100);
		print "$out\n";
		`perl clusterGenes3.pl $fastaSeqs $tfasty $allFasta $aa_arg $nuc_arg $ratio_arg >$out`;
		my $cr_out = $cr."_".($aa_arg*100)."_".($nuc_arg*100)."_".($ratio_arg*100);
		`perl clusterReport.pl <$out >$cr_out`;
	}
}
