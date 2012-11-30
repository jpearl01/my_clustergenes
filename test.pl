#!/usr/bin/perl

use strict;
use warnings;

my $cr = "report_output";

open C, $cr or die "Can't open the clustergenes file: $!\n";


my @hist; #Array to store the histogram
my $bool =0;
my $counter=0;
while (<C>){
	chomp;
	if (/\<report stop\>/){
		$bool=0;
	}
	if ($bool==1){
		$counter++;
		$hist[$counter] = (split ("\t",$_))[1];
	}
	if (/\<report start\> gene copy histogram/){
		$hist[0]= $cr;
		$bool = 1;
	}
}
for my $k (@hist){
	print $k."\n";		
	#open O, ">","$k.hist" or die "Can't open output file $k : $!\n";
	#print O $clusts{$k};
	#close O; 
}
