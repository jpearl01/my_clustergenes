#!/usr/bin/perl

use warnings;
use strict;
use Statistics::Descriptive;
use Data::Dumper;
use Graph::Undirected;
use Algorithm::EquivalenceSets; 
########################################
#
#Version 1.1
#
#This program calculates a distance matrix based on allelic differences in core genes from a supragenome project
#
#Syntax:
#perl JANG_fasta.pl CoreGenesList /mnt/p/o_drive/Homes/jearl/mcat/all.fasta 12
#
#Where CoreGenesList is created from initJANG.pl and all.fasta is the conglomerated fasta matches from a normal
#supragenome project run.  I've added the requirement to have the "strains_list" in the directory, so we can 
#look at subsets of the full run.
#
#Added in a file so that we can examine the gene names of the identities.
#############################################
my $cl;					#File of the list of core clusters
my $tf;					#File of the tfasty alignments
my $supraSize;			#number of genomes in supragenome
our %nearestNeighbor;	#Hash that holds nearest neighbor, and the score
my $validNeighbor;
my %singletsHash;		#probably more accurate to call this NOT singlets, since its a list of all the strains in groups
my %groupHash;			#hash that holds all the groups (except singlets)
my %updatedSingletsHash;
my @equiv_sets;

#stat object for mean, sdev and stand-err
my $stat = Statistics::Descriptive::Full->new();

#Create a graph to hold the matrix
my $graph = Graph::Undirected->new;  

#missing sequence, and warnings log
open MS, '>', 'MissingSequence.log' or die "Couldn't open the missing sequence logfile :$!";
open W, '>', 'warnings.log' or die "Couldn't open the warnings.log file: $!"; 

die "We need three inputs, USAGE: perl JANG_fasta.pl CoreGeneList all-against-all.fasta36 integer_value_of_strains\n" unless ($ARGV[0] && $ARGV[1] && $ARGV[2]);

#Check if Args are initialized (command line)
if (defined $ARGV[0]){
	$cl = $ARGV[0];
}
else {
	die "Please enter a core clusters file first, followed by the FASTA alignments file, followed by the total number of strains in the project (core clusters error)"
}

if (defined $ARGV[1]){
	$tf = $ARGV[1];
}
else {
	die "Please enter a core clusters file first, followed by the FASTA alignments file, followed by the total number of strains in the project (fasta error)"
}

if (defined $ARGV[2]){
	$supraSize = $ARGV[2];
}
else {
	die "Please enter a valid integer value of the number of strains";
}

#Added by Josh, the functionality to read in a strain list from a file (woot)
open (STRAIN, "strain_list") or die "can't open strain_list file:$!\n";
our @STRAIN_LIST ="";
our %strain_hash;

while (<STRAIN>)
{
	chomp $_;
	push(@STRAIN_LIST, $_);
	$strain_hash{$_}=1;
}

my $cgHash = load_coreGenes('CoreGeneList');

#Some core clusters do not have sequence data available, i.e. they matched
#against a contig instead of another protein, we can't use those, so get rid of them.
my $tot_core_clusters = scalar(keys %$cgHash);
my $coreCounter = 0;
foreach my $k (keys %$cgHash){
	if (scalar(@{${$cgHash}{$k}})<$supraSize){
		print MS "Not enough sequence data for cluster $k\n";
		print MS "Size of cluster is ".scalar (@{${$cgHash}{$k}})."\n";
		$coreCounter++;
		delete ${$cgHash}{$k};
	}
#The other problem is some clusters have too many genes. Exclude these too
	elsif(scalar(@{${$cgHash}{$k}})>2*$supraSize){
		print MS "Too many sequence data for cluster $k (more than twice as many as expected)\n";
		print MS "Size of cluster is ".scalar (@{${$cgHash}{$k}})."\n";
		$coreCounter++;
		delete ${$cgHash}{$k};
	}
}
print "The total number of clusters excluded due to lack of sequence info is $coreCounter\n";
print "The total number of core clusters from Justin's analysis is: $tot_core_clusters\n";
my $used_clusters = $tot_core_clusters-$coreCounter;
print "The total clusters used for this analysis is $used_clusters\n";

my ($fastaIdenRef, $strainRef) = load_fasta($tf);

my $matrixRef = build_matrix($cgHash, $fastaIdenRef, $strainRef);

build_groups($strainRef);

#print Dumper($matrixRef);

print_matrix($matrixRef, $strainRef);

print_jgraph($matrixRef);

print "Finished successfully\n";

###/=============\###
###| SUBROUTINES |###
###\=============/###

sub load_coreGenes{
	my $core_file = shift;
	my %coreGenes;
	my $cluster;
	my @genes;
	my $bool;
	open CORE, $core_file or die "ERROR: Couldn't open the file $core_file: $!";
	while (my $line_input = <CORE>){
		chomp $line_input;
		if ($line_input =~ /^>cluster_(\S+)/){
			$cluster = $1;
			$bool=0;
		}		
		if ($bool==1){
			push @{$coreGenes{$cluster}}, $line_input;
		}
		if ($bool == 0){
			$bool = 1;
		}
	}
#	print Dumper(%coreGenes);
	return \%coreGenes;
}

#For the most part, I lifted this code from Justin's clusterGenes.pl file.
#There are some differences though, fyi
sub load_fasta{
	my $fasta_file = shift;
	my %identity;
	my %strains;

	open FASTA, "<", $fasta_file  or  die "ERROR: could not open fasta file $fasta_file!\n";
	BEGIN_REC: while ( my $line_input = <FASTA> ){
		
		# look for the beginning of a record
		if ( $line_input =~ /^\s*\d+>>>/ ){
			
			$line_input =~ s/ctg\d+_/_/i;
			# we found a record.  get strain and locus.
			my ($gene1_name) = $line_input =~ />>>([A-Za-z0-9]+_\d+)/;
			my ($strain, $locus, $length) = $line_input =~ />>>([A-Za-z0-9]+)_(\d+).*\s+-\s+(\d+) nt/; 
			die "Problem parsing the strain, locus, or length from the beginning record on line $.\n" unless ($strain && $locus && $length);

			#added this to account for times we don't want to only examine a subset.
			if (!(exists $strain_hash{$strain})){
				print STDERR "strain $strain doesn't exist in our list, next record\n";
				next BEGIN_REC;
			}
			#Added by Josh, create the list of strains
			if (!exists $strains{$strain}){$strains{$strain}=1}
			
			# search for the start of the matches
			while ( $line_input = <FASTA> ){
				last if ($line_input =~ /^The best scores are:/);
			}

			# gather matches until a blank line is encountered
			MATCH: while ( $line_input = <FASTA> ){
				last if ( $line_input =~ /^>>>/ || $line_input =~ /^\n/ );
				next if ($line_input=~/^\+-/);
				chomp $line_input;

				$line_input =~ s/ctg\d+_/_/i;
				my ($header,$data) = split( /\t/, $line_input );
				die "Something wrong with the pattern match, trying to split on a tab on line $. :\n$line_input\n" unless ($header && $data);
				my ($gene2_name) = $header =~ /^([A-Za-z0-9]+)/;
				#print "The gene match name it is printing is $gene2_name\n";
				die "Something wrong with the pattern match, trying to get beginning gene name on line $. :\n$line_input\n" unless ($gene2_name);
				my ($match_strain) = $header =~ /^([A-Za-z0-9]+)/i;
				#now, trim off the contig ends, if there
				$match_strain =~ s/ctg\d+//i;
				die "Something wrong with the pattern match, trying to get strain and locus on line $. :\n$line_input\n" unless ($match_strain);

				#Again, added this to take care of pesky fasta alignment files, when we only want to look at a subset.
				#print $line_input."\n";
				#print $match_strain."\n";
				if (!(exists $strain_hash{$match_strain})){
					#print STDERR "strain $strain doesn't exist in our list, next record\n";
					next MATCH;
				}
				
				my @match_data = split( /\s+/, $data );
				my $percent_identity = $match_data[0];
				my $alignment_length = $match_data[3];

				#Apparently, you can get two alignments for the same sequence, both the 
				#forward and also the reverse.  Here I don't care to much, just take the largest
				#value of the two for the percent identity.  I do however, write the sequence IDs
				#to a log file. Also, for the fasta files I appear to be finding ortho/paralogs.  If they 
				#exist in the same contig, then the highest scoring one is taken.  If in different
				#contigs they will be averaged in with the rest.
				
				if (!exists $identity{$match_strain}->{$gene1_name}){
					$identity{$match_strain}->{$gene1_name} = $percent_identity;
				}
				else{
					if ($identity{$match_strain}->{$gene1_name}<$percent_identity){
						print W "Identity already exists between $gene1_name and $match_strain, which is $identity{$match_strain}->{$gene1_name}.  New identity is $percent_identity\n";
						$identity{$match_strain}->{$gene1_name} = $percent_identity;
					}
				}
			}
			# we're at the end of the record, start looking for the next record	
		}
	}
	close FASTA;
	my @strains;
	foreach my $s(keys %strains){push @strains, $s;} 
	return \%identity, \@strains;
}

sub build_matrix{
	my $cgRef 		= shift;
	my $idenRef 	= shift;
	my $strainRef 	= shift;
	my %matrix;
	my %geneCounter;		#keep track of the number of core genes seen, use as denominator to normalize
	my $badCounter = 0;		#keep track fo the number of "bad" i.e. genes with no %identity
	my $totCounter = 0;
	
	#Initialize the scoring matrix
	foreach my $a (@{$strainRef}){
		foreach my $b (@{$strainRef}){
			$geneCounter{$a}{$b} = 0;
			$matrix{$a}{$b}=0;
		}
	}
	open IDEN, '>', 'identities' or die "Couldn't open identities: $!\n";
	open BAD, '>', 'badMatches' or die "Couldn't open badMatches: $!\n";
	open IDEN2, '>', 'identities_geneNames' or die "Couldn't open the identities_geneNames";
	print IDEN "strain1,strain2,identity\n";
	foreach my $k (keys %{$cgRef}){
		
		foreach my $j (@{${$cgRef}{$k}}){
			my ($strain1, $loc1) = $j =~ /([a-z,A-Z,0-9]+)_(\d+)/;
				
				foreach my $strain2(@{$strainRef}){
				$totCounter++;
				if($strain1 ne $strain2){
					
					if (exists ${$idenRef}{$strain2}{$j}){
						$geneCounter{$strain1}{$strain2} += 1;	
						$geneCounter{$strain2}{$strain1} += 1;	
						$matrix{$strain1}{$strain2} += ${$idenRef}{$strain2}{$j};
						$matrix{$strain2}{$strain1} += ${$idenRef}{$strain2}{$j};
						print IDEN "$strain2,$strain1,${$idenRef}{$strain2}{$j}\n";
						print IDEN2 "$strain2,$j,${$idenRef}{$strain2}{$j}\n";
					}
					else{
						$badCounter++;
						print BAD "Bad match, can't find $j and $strain2 %identity\n";
					}
				}
			}
		}
	}
	
	print "There were $badCounter missing %identities out of a total of $totCounter\n";

	#normalize matrix
	foreach my $a (@{$strainRef}){
		foreach my $b (@{$strainRef}){
			if ($a ne $b){
				#print "$a and $b are the strains\n";
				$matrix{$a}{$b} = $matrix{$a}{$b}/$geneCounter{$a}{$b};
				$stat->add_data($matrix{$a}{$b});
				#Add current edge to the graph (includes vertices)
				$graph->add_weighted_edge($a, $b,(1-$matrix{$a}{$b}));
				if ((!exists $nearestNeighbor{$a}) || ($nearestNeighbor{$a}[1]>(1-$matrix{$a}{$b}))){
					$nearestNeighbor{$a}[0]= $b;
					$nearestNeighbor{$a}[1]= (1-$matrix{$a}{$b});
			}
			}
		}
	}
	
	#Then calculate the mean and std error
	my $mean = (1-$stat->mean());
	print $stat->count()."\n";
	my $se = $stat->standard_deviation/sqrt($stat->count());
	$validNeighbor = $mean-$se;
	printf "Mean is %7.4f Standard error is %7.4f\n", $mean, $se;
	printf "Valid neighbors are separted by distances less than %7.4f.\n\n",$validNeighbor;
	
#	print Dumper(%matrix);
	return \%matrix;
}

sub print_matrix{
	my $matrix 	= shift;
	my $strains = shift;
	
	open LOGFILE, '>', 'myANGlog' or die "Couldn't open the logfile: $!\n";
	
	#print top lables for matrix.
	foreach my $strain1 (@{$strainRef}){
#		my ($short) = $strain1 =~ /(\S\S\S$)/;
#		print "\t ".$short;
		print "\t ".$strain1;
	}

	print "\n";
	
	#Print out the distance matrix  via the strainlist
	foreach my $strain1 (@{$strainRef}){
		
		#print out the shortened name of the current strain
#		my ($short) = $strain1 =~ /(\S\S\S$)/;
		#print out normal name.
		print $strain1;
		print LOGFILE $strain1;
		
		foreach my $strain2 (@{$strainRef}){
			if ($strain1 ne $strain2){
				print "\t".sprintf("%7.5f",(1-${$matrix}{$strain1}{$strain2}));
				print LOGFILE "\t".sprintf("%7.5f",(1-${$matrix}{$strain1}{$strain2}));
			}
			else{
				print "\t0.00000";
				print LOGFILE "\t0.00000";
			}
		}
		print "\n";
		print LOGFILE "\n";
	}
	print LOGFILE "There are ".($#equiv_sets+1+(keys %updatedSingletsHash))." Groups.\n";
	
	my $groupID = 1;
	foreach my $group (@equiv_sets){
		my $newline = 1;
		print LOGFILE "\nGroup $groupID consists of:\n";
		foreach my $strain (@$group){
			print LOGFILE "\t $strain";
			if (($newline % 5)==0){
				print LOGFILE "\n";
			}
			$newline++;
		}
		print LOGFILE "\n";
		foreach my $nn (@$group){
			printf LOGFILE "\nNearest neighbor of $nn is $nearestNeighbor{$nn}[0] and their distance is %7.5f", $nearestNeighbor{$nn}[1];
		}
		print LOGFILE "\n\n";
		$groupID++;
	}
	
	#Handle Singlets
	foreach my $g1 (keys %updatedSingletsHash){
		print LOGFILE "Group $groupID consists of:\n";
		print LOGFILE "\t$g1\n";
		print LOGFILE "$g1 has no valid neighbors and is thus a singlet\n";
		$groupID++;
	}
	
}

sub print_jgraph{
	
###########################
#All my graphing implementation
#Create a minimum spanning tree from the current graph
#initialize a graphviz object to print out image.
###########################
	my $matrix = shift;
	my $mst = $graph->MST_Kruskal();
	open GRAPH_JAVA, '>graph_for_java.csv' or die "Couldn't open the csv file for the java graph implementation";
	foreach my $edge ($mst->edges){
		my $lab;
		$lab=sprintf("%7.7f", (1-${$matrix}{$edge->[0]}{$edge->[1]}));
		print GRAPH_JAVA "$edge->[0],$edge->[1],$lab\n";
	}
	close GRAPH_JAVA;
}


sub build_groups{
	
#############################################
#Function to build the groups.
#############################################	
	my $strainRef 	= shift;
	my $isGroup = 1;#boolean to keep track of existing groups
	my $groupID = 1;#hash key for groupHash
	foreach my $g1 (@{$strainRef}){
		foreach my $g2 (@{$strainRef}){
			if (($g1 ne $g2) && &isNNeighbor($g1,$g2)){
				#These are neighbors, keep track to later print out singlets
				$singletsHash{$g1} =
				$singletsHash{$g2} = 1;
				#Initialize the groupHash with the first neighbor found
				if (!%groupHash){
					$groupHash{$groupID}{$g1}=1;
				}
				#Since these are neighbors we want to either:
				foreach my $g3 (keys %groupHash){
					#Search and add to a currently existing group (set both hash values to one, removes one logic step)
					if ((exists $groupHash{$g3}{$g1})||(exists $groupHash{$g3}{$g2})){
						$groupHash{$g3}{$g1} =
						$groupHash{$g3}{$g2} = 1;
						$isGroup = 1;#A group already existed, no need to make a new one;
						last;#leave the loop, no need to continue iterating
					}
					$isGroup = 0;#otherwise, it will reach this point and a new group will need to be formed.
				}
				#Create a new group, if one wasn't found
				if (!$isGroup){
					$groupID++;
					$groupHash{$groupID}{$g1} =
					$groupHash{$groupID}{$g2} = 1;
				}
			}
		}
	}
	#Handle the case of overlapping groups
	#Changed data to arrays, since equivalence_sets only works with those
	my @sets = map { [ keys %{$groupHash{$_}} ] } keys %groupHash;
	@equiv_sets = equivalence_sets(\@sets);
	
	#update the singlest hash
	foreach my $g1 (@{$strainRef}){
		if (!exists $singletsHash{$g1}){
			$updatedSingletsHash{$g1}=1;
		}
	}
}

sub isNNeighbor{
	my ($n1, $n2) = @_;
	if(($nearestNeighbor{$n1}[0] eq $n2)&&($nearestNeighbor{$n1}[1]<$validNeighbor)){
		return 1;
	}
	else {
		return 0;
	}
}

