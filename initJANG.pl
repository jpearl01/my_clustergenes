#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
######################################
#This formating program was designed specifically for RAST style clusterGenes
#Although, so far it hasn't broken on PGAAP style that I can tell
#
#Usage:
#perl initJang.pl clusterGenes report_output
#
#Where clusterGenes is the output of the clustergenes.pl program and report_output is the output
#from the clusterReport.pl program (must have the cluster names for core/distributed!)
#Output:
#CoreGeneList
#This file is used in the JANG.pl program to produce the distance matrix used for tree building
######################################

die "We need two arguments.  Usage: perl iniJang.pl clusterGenes report_output\n" unless ($ARGV[0] && $ARGV[1]);

open CR, $ARGV[0] or die "clusterGenes did not successfully open: $!";
open RO, $ARGV[1] or die "report_output did not successfully open: $!";
open MCR, '>CoreGeneList' or die "Couldn't make CoreGeneList: $!";

my %hashOfGenes;
my %hashOfGenomes;
my $bool = 0;
my $currCluster;

while (<RO>){
	chomp;
	if (($bool==1) && ($_=~ /report stop/)){
	    last;
	}
	if ($bool==1){
	    $hashOfGenes{$_}=1;
	}
	if($_=~ /<report start> total core clusters:/){
	    $bool=1;
	}
}

#print Dumper(%hashOfGenes);
print "The number of Core Genes is: ".(keys %hashOfGenes)."\n";
my @array1;
my @array2;
$bool=0;
while (<CR>){
	chomp;
	s/\s+$//;
	if ($_ =~ /(<cluster start> )([A-Za-z0-9]+_\d+)/){
		if (exists $hashOfGenes{$2}){
				$bool = 1;
		}
		else {
			$bool = 0;
		}
	}

	if ($bool){
		if ($_ =~ /(<cluster start> )([A-Za-z0-9]+_\d+)/){
			print MCR ">cluster_$2"."\n";
			push (@array1, $_)
		}
		#Match reps we actually have sequence data for.
		elsif($_ =~ />([A-Za-z0-9]+_\d+)/){
			push(@array2, $1);
			print MCR $1."\n";
		}
	}
	elsif($_ eq ""){
		print "There is an empty string in your data at line $. Continuing normally.\n"
	}
}
if (@array1==0){
	print "Warning!  Your new file contains no data!  Probably a matching problem\n";
}

print "Size of incoming cluster list (should be the same as above): ".(scalar (@array1))."\n";
print "Size of incoming core genes ".(scalar @array2)."\n";
print "Program completed successfully\n"
