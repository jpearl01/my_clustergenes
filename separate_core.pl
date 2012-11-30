#!/usr/bin/perl

use strict;
use warnings;

##This program is to be run after extract_fasta_clusterGenes.pl
##The idea is after you get your clusters into individual fasta files
##You then want to separate them into core and distributed folders
##This should do something like that.

open CR, $ARGV[0] or die "Can't open the $ARGV[0] file: $!\n";

my %core_clusters;
my $bool = 0;

while(<CR>){

	if (/<report stop>/){
		last;
	}
	if ($bool==1){
	    chomp;
	    $core_clusters{$_} = 1;
	}
	if (/<report start>/){
		$bool = 1;
	}
}

`mkdir core` unless (-d "core");

for my $i (keys %core_clusters){
    `mv $i.fasta core/`;
}
