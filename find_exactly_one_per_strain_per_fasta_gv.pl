#!/usr/bin/perl
use warnings;
use strict;

my @strains = qw(B476
B477
B478
B473
AMD
B483
B475
B472
14019
B479
B512
14018
Gv5
B482
B474
B513
409);
my %strainHash = map { $_ => 0 } @strains;

`mkdir exactly_one` unless (-d exactly_one);

FILES: while(<*>){
	my $fname = $_;
	open I, $fname or die "Can't open $fname: \n$!\n";
	my @fasta=<I>;
	chomp(@fasta);
	if (@fasta > 20){next FILES;}
	else{
	    for my $s(@strains){
		my @count=grep {$_=~/$s/}@fasta;    
		if(@count!=1)(next FILES}
	    }
	}
#	`cp $fname exactly_one/`;
	print "File $fname has exactly one per strain\n";
}
