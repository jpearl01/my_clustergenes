#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

#This formating program is specifically for RAST clusterGenes.  It outputs
#MyClusterReport, which is the input for the JNG.pl neighbor grouping program.
#This file is poorly named.
#
#Usage:
#
# perl formatCR.pl clusterGenes report_output

open CR, $ARGV[0] or die "clusterGenes did not successfully open: $!";
open RO, $ARGV[1] or die "report_output did not successfully open: $!";
open MCR, '>MyClusterReport' or die "Couldn't make MyClusterReport: $!";

my %hashOfGenes; #This holds all the distributed (contingency) clusters
my $bool = 0;

while (<RO>){
	chomp;
	if ($_=~ /report stop/){
	    $bool = 0;
	}
	if ($bool==1){
	    $hashOfGenes{$_}=1;
	}
	if($_=~ /<report start> total contingency clusters:/){
	    $bool=1;
	}
}

#print $hashOfGenes{'HPAG1_0028'}.": So there!\n";
#print Dumper(%hashOfGenes);
for my $g(keys %hashOfGenes){
    if ($g =~ /HPAG1_0028/){
	print $hashOfGenes{$g}."\n";
    }
}
print "The number of Distributed Gene Clusters is: ".(keys %hashOfGenes)."\n";
my @array1; #The first array holds all the names of the clusters (from "<cluster start> ...")
my @array2; #The second array holds the total number of genes in all distributed clusters
while (<CR>){
	chomp;
	if ($_ =~ /(<cluster start> )([A-Za-z0-9]+_[0-9]+)/){
	    my $gene = $2;
#	    print $gene."\n";
	    if (exists $hashOfGenes{$gene}){
	    	$bool = 1;
	    }
	    else {
	    	$bool = 0;
	    }
	}

	if ($bool){
		if ($_ =~ /(<cluster start> )([A-Za-z0-9]+_[0-9]+)/){
		    my $gene2 = $2;
		    print MCR ">cluster_$gene2"."\n";
		    push (@array1, $_)
		}
		elsif($_ =~ /(\s+)([A-Za-z0-9]+)(\s+)(1)/){
		    my $strain = $2;
		    print MCR $strain."_1"."\n"; #the one at the end is to make it compatible with previous version
		    push @array2, $strain;
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
print "Number of distributed genes ".(scalar @array2)."\n";
print "Program completed successfully\n"
