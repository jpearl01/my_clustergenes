#!/usr/bin/perl

use strict;
use warnings;

my @files = `ls report_output*`;
my @hists; #Array to store the histograms
my $numFiles = 0;

for my $f (@files){
	open F, $f or die "Can't open the report_output file: $!\n";
	my ($col_name) = $f =~ /(\d+_\d+)$/;
	my $bool =0;
	my $counter=0;
	$hists[$numFiles][$counter] =$col_name;
	$counter++;
	while (<F>){
		chomp;
		if (/\<report stop\>/){
			$bool=0;
		}
		if ($bool==1){
			$counter++;
			$hists[$numFiles][$counter] = (split ("\t",$_))[1];
		}
		if (/\<report start\> gene copy histogram/){
			$bool = 1;
		}
	}
	$numFiles++;
}
foreach my $i (0..@files){
	foreach my $j (0..17){
		print $hists[$i][$j]."\t";		
		#open O, ">","$k.hist" or die "Can't open output file $k : $!\n";
		#print O $clusts{$k};
		#close O; 
	}
	print "\n";
}
