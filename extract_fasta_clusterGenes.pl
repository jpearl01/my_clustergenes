#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

##Creates fasta files of all the genes from clusterGenes output
##
##Usage: perl extract_fasta_clusterGenes.pl seqs.fna clusterGenes 
##Gene sequences file = seqs.fna
##clusterGenes output = clusterGenes

open SEQS, $ARGV[0] or die "Can't open the sequence file: $!\n";
open C, $ARGV[1] or die "Can't open the clustergenes file: $!\n";

my %seqs;
my %clusts;

while (<SEQS>){
	chomp;
	s/ //;
	if (/^>/){$seqs{$_}=1};	
}

my $bool =0;
my $curr_clust;
while (<C>){
	if (/<cluster start>\s(\S+)/){
		$curr_clust = ">".$1;
		chomp $curr_clust;
	}
	if (/<genes stop>/){
		$bool = 0;
	}
	if ($bool==1){
		$clusts{$curr_clust} .= $_;
	}
	if (/<genes start>/){
		$bool = 1;
	}
}
`mkdir fastas` unless (-d "fastas");

my $CG_size = scalar(keys %clusts);
my $seqs_size = scalar(keys %seqs);
#print Dumper(%seqs);
print "The number of Clusters in $ARGV[1] is $CG_size, the number of seqs in $ARGV[0] is $seqs_size\n";

for my $k (keys %clusts){
	if (exists $seqs{$k}){
		my ($out) = $k =~/\>(.+)/;
		open O, ">","fastas/$out.fasta" or die "Can't open output file $k : $!\n";
		print O $clusts{$k};
		close O; 
	}
	else{
	    #print "$k does not exist in the fasta file\n";
	}
}
