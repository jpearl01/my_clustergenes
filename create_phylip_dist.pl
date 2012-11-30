#!/usr/bin/perl -w
use strict;

#I need a program that will take as input the distance matrices generated from the NG
#algorithm(s) and make the changes required to be compatable w/ phylip. Note however that
#the naming scheme you use will probably have to be altered by hand, phylip requires 
#that the name be "exactly" 10 characters, and if it is more it must be trimmed, if it
#is less, then x "padding" spaces must be added on the end so that name+x = 10
#It also doesn't seem to like more than 5 decimal places. This program takes care of
#everything but the names that are too long

#Be aware the input file should have a beginning line with just a tab and the number
#of genomes in the tree, and everything else should jut be tab delimited.

open M, $ARGV[0] or die "Can't open the matrix: \n$!\n";
open O, ">", "out_matrix" or die "Can't open the output file for the new matrix: \n$!\n";

my @inFile=<M>;
my $size = @inFile;
my $is_first = 1;

print O "\t$size\n";
foreach my $var (@inFile){
	chomp $var;
	my @a = split("\t",$var);
	if (@a<2){next}
	if (length($a[0])<10){$a[0] = $a[0]." " until(length($a[0])==10)}
	print O $a[0]."\t";
	foreach my $v (@a){
		if ($is_first == 0){
			if (length($v)>5){
				$v =~ s/^\s//;
				print O substr $v, 0, 7;
				print O "\t";
			}
		}
		else{
			$is_first=0;
		}
	}
	print O "\n"; 
	$is_first = 1;

}