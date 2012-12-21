#!/usr/bin/perl
use warnings;
use strict;
use Statistics::Descriptive;
use Data::Dumper;
use Algorithm::EquivalenceSets;
use Graph::Undirected;
use Graph::Easy;
use Graph::Convert;
use GraphViz;
use Carp;
#use GD::Graph;


#####################################################
#Author: Josh Earl
#Last Date Modified: May 29 2009

# OVERVIEW:
# =========
#This is my implementation of Dr. Barry Hall's Neighbor grouping algorithm. I've added a few outputs
#like graph_for_java which can be used with my program GraphMST2.java (creates a minimum spanning tree)

#Usage
# =========
#The inputs to this program is a file called "MyClusterReport" (no extension).  This is a
#file that is the output of formatCR.pl. formatCR.pl takes as input report_output and clusterGenes
#As long as the files are in the same directory, nothing else is needed to successfully run the scripts
#example:  perl JNG.pl
#####################################################





open CR, 'MyClusterReport' or die "MyClusterReport did not successfully open.";

our %hashOfGenomesByCluster;
our %hashOfBinaryScores;	#Keeps track of whether a cluster is represented in a strain(cluster=>strain=>1 or 0)
our %hashOfGenomes;			#List of Genomes... I like to use hashes even for this for the "exists" function
our %hashOfGenes;			#List of the genes.
our %hashOfFinalScores;		#uses binary scores, keeps track of total # of clusters a strain has
our %distanceBetweenGenomes;#straight up distance between the strains(strain1=>strain2=>.35)
our %groupHash;				#hash that holds all the groups (except singlets)
our %updatedGroupHash;		#Hash that holds the merged groups.
our %singletsHash;			#probably more accurate to call this NOT singlets, since its a list of all the strains in groups
our %updatedSingletsHash;	#Store # of actual singlets
our %nearestNeighbor;		#Hash that holds nearest neighbor, and the score
my $numGenes = 0;			#number of clusters in the supragrenome
my $thisDate = localtime();


my $graph = Graph::Undirected->new;  #Create a graph to hold the matrix

#Start timing the program
my $ElapsedTime = 0;
my $EndTime;
my $BeginTime = time;

#Use the input MyClusterReport file to initialize the hash of the genomes,
#hash of the genomes by cluster, count the number of genes observed and
#place them in the hash of genes and build the hash of binary scores
#which we later use as a kind of boolean test to calculate the "real"
#scores of distances between genomes.
my $currCluster;
while (<CR>){
	chomp;
	if ($_ =~ /(>cluster_)([A-Za-z0-9]+)(_[0-9]+)/){
		$currCluster = $2.$3;
		if (!exists $hashOfGenomes{$2}){
			$hashOfGenomes{$2}=0;
		}
		$numGenes++;
	}
	elsif($_ =~ /([A-Za-z0-9]+)(_[0-9]+)/){
		if (defined $currCluster){
			$hashOfGenomesByCluster{$currCluster}{$1}=1;
		}
		else{
			die "The program is trying to use an undefined cluster at line $. which is: \n$_\nin file $ARGV";
		}
		$hashOfBinaryScores{$currCluster}{$1} = 1;
		if (!exists $hashOfGenes{$_}){
			$hashOfGenes{$_}=1;
		}
		$hashOfGenomes{$1}++;
	}
	elsif($_ eq ""){
		print "There is an empty string in your data at line $. in file $ARGV Continuing normally."
	}
	else{
		die "The regular expressions found a line that was not formatted as expected at line number $. which is:\n$_\nin file $ARGV";
	}
}

print "\nNumber of distributed genes is $numGenes.\n";
print "\nNumber of Genomes is ".(keys %hashOfGenomes)."\n";

#Look at every cluster.  If in that cluster both genome1 and genome2 exist, either create
#a new hash entry for genome1 linked to genome2 with score 1, or increase the score by 1
#of an existing count. Do the same if they both do not exist.
foreach my $cluster (keys (%hashOfGenomesByCluster)){
	foreach my $genome1 (keys (%hashOfGenomes)){
		foreach my $genome2 (keys (%hashOfGenomes)){
			if ((exists $hashOfBinaryScores{$cluster}{$genome1})&&(exists $hashOfBinaryScores{$cluster}{$genome2})){
				if (exists $hashOfFinalScores{$genome1}{$genome2}){
					$hashOfFinalScores{$genome1}{$genome2} += 1;
				}
				else {
					$hashOfFinalScores{$genome1}{$genome2} = 1;
				}
			}
			if ((!exists $hashOfBinaryScores{$cluster}{$genome1})&&(!exists $hashOfBinaryScores{$cluster}{$genome2})){
				if (exists $hashOfFinalScores{$genome1}{$genome2}){
					$hashOfFinalScores{$genome1}{$genome2} += 1;
				}
				else {
					$hashOfFinalScores{$genome1}{$genome2} = 1;
				}
			}
			else{
				if (!exists $hashOfFinalScores{$genome1}{$genome2}){
					$hashOfFinalScores{$genome1}{$genome2} = 0;

				}
			}
		}
	}
}

my %nodes;
my $i = 1;
foreach my $g (keys %hashOfGenomes){
	print $i."\t".$g."\n";
	$nodes{$g}=$i;
	$i++;
}
print "\n";
#print "The total final score of CGSSp9vBS293 and CGSSpBS397 is $hashOfFinalScores{CGSSpBS397}{CGSSp9vBS293}.\n";

my $stat = Statistics::Descriptive::Full->new();

foreach my $g1 (keys (%hashOfGenomes)){
	foreach my $g2 (keys (%hashOfGenomes)){
		#As long as the genomes are not the same (both in the score, or the name) then continue
		#The reason we remove the case of the scores being the same is it would give an error: unclear why at this juncture
		if (!($g1 eq $g2)  && (($hashOfFinalScores{$g1}{$g2}/(keys %hashOfGenomesByCluster))!=1)){
			$distanceBetweenGenomes{$g1}{$g2} = (1-($hashOfFinalScores{$g1}{$g2}/(keys %hashOfGenomesByCluster)))
				or die "There was a problem calculating the score distance between genomes $g1 $g2, error message is: \n$!
				\nThe denominator in this equation is ".(keys %hashOfGenomesByCluster)
				." which is defined if next value is 1: ".(defined (keys %hashOfGenomesByCluster))
				."\nThe numerator is $hashOfFinalScores{$g1}{$g2} which is defined if 1: "
				.defined($hashOfFinalScores{$g1}{$g2})."\n";
			$stat->add_data($distanceBetweenGenomes{$g1}{$g2});
			#Add current edge to the graph (includes vertices)
			$graph->add_weighted_edge($g1, $g2, $distanceBetweenGenomes{$g1}{$g2});
			#I'm changing things up here with a hash of arrays. index 0 holds the name of the
			#Nearest Neighbor, and index 1 holds the distance between them
			#basically, if there is no current nearest neighbor, add this one.  And if the current distance is
			#less then the one previously stored, change the values to the new nearest neighbor.
			if ((!exists $nearestNeighbor{$g1}) || ($nearestNeighbor{$g1}[1]>$distanceBetweenGenomes{$g1}{$g2})){
				$nearestNeighbor{$g1}[0]= $g2;
				$nearestNeighbor{$g1}[1]= $distanceBetweenGenomes{$g1}{$g2};
			}
		}
		else{
			$hashOfFinalScores{$g1}{$g2} = 0;
			$distanceBetweenGenomes{$g1}{$g2} = 0;
		}
	}
}
###########################
#All my graphing implementation
#Create a minimum spanning tree from the current graph
#initialize a graphviz object to print out image.
###########################
my $mst = $graph->MST_Kruskal();
my $graphVizMst = GraphViz->new(
    node 		=>	{
	shape	=> 	'box',
	color	=>	'red'
    },
    layout 		=> 'neato',
    directed	=> 0,
    width 		=> 20,
    height 		=> 15,
    pagewidth 	=> 22,
    pageheight 	=> 17,
    overlap		=>'scale'
    );
my $lab;
my $lab2;
#Open csv for the java graph program, convert $mst to the graphviz object
open GRAPH_JAVA, '>graph_for_java.csv' or die "Couldn't open the csv file for the java graph implementation";
foreach my $edge ($mst->edges){
	$lab=sprintf("%7.5f",$distanceBetweenGenomes{$edge->[0]}{$edge->[1]});
#	$graphVizMst->add_edge($edge->[0], $edge->[1], 'len'=>(10+(log($distanceBetweenGenomes{$edge->[0]}{$edge->[1]}))), 'label'=>$lab);
	$graphVizMst->add_edge($edge->[0], $edge->[1], 'len'=>$distanceBetweenGenomes{$edge->[0]}{$edge->[1]}, 'label'=>$lab);
	print GRAPH_JAVA "$edge->[0],$edge->[1],$lab\n";
	$lab2=sprintf("%7.2f",$distanceBetweenGenomes{$edge->[0]}{$edge->[1]})*100;
	print $nodes{$edge->[0]}."\t".$nodes{$edge->[1]}."\t".$lab2."\n";
}
close GRAPH_JAVA;


#foreach my $edge($graph->edges){
#	$lab=sprintf("%7.5f",$distanceBetweenGenomes{$edge->[0]}{$edge->[1]});
#	$lab *= 100000;
#	print $nodes{$edge->[0]}."\t".$nodes{$edge->[1]}."\t".$lab."\n";
#}

#print out graphviz mst image
open GRAPH, '>graph.ps' or die "Couldn't open the graph file";
$graphVizMst->as_ps(\*GRAPH);
close GRAPH;

#Then calculate the mean and std error
my $mean = $stat->mean();
my $se = $stat->standard_deviation/sqrt($stat->count());
my $validNeighbor = $mean-$se;

#############################################
#Function to build the groups.
#############################################

my $isGroup = 1;#boolean to keep track of existing groups
my $groupID = 1;#hash key for groupHash
foreach my $g1 (keys %hashOfGenomes){
	foreach my $g2 (keys %hashOfGenomes){
		if (($g1 ne $g2) && &isNNeighbor($g1,$g2)){
			#These are neighbors, keep track to later print out singlets
			$singletsHash{$g1} =
			$singletsHash{$g2} = 1;
			#Initialize the groupHash with the first neighbor found
			if (!defined %groupHash){
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
my @equiv_sets = equivalence_sets(\@sets);

#update the singlest hash
foreach my $g1 (keys %hashOfGenomes){
	if (!exists $singletsHash{$g1}){
		$updatedSingletsHash{$g1}=1;
	}
}

##########################################################
#Print all this stuff out
##########################################################

open LOGFILE, '>myNNlogfile' or die "Was unable to open the Logfile";
print LOGFILE "NG log file $thisDate.\n\n";
print LOGFILE "This analysis is based upon $numGenes distributed genes in the supragenome as estimaed by Justin Hogg.\n\n";
print LOGFILE "\n\n\n";
print LOGFILE "The values in the matrix below are distances, i.e. the fraction of distributed genes that are present.\n";
print LOGFILE "in one genome but absent in the other.";
#Print the distance Matrix:
foreach my $d1 (keys (%distanceBetweenGenomes)){
	print LOGFILE "\n".$d1;
	foreach my $d2 (keys %{$distanceBetweenGenomes{$d1}}){
		printf LOGFILE "\t%5.4f",$distanceBetweenGenomes{$d1}{$d2}."  ";
	}
}
printf LOGFILE "\n\nThe mean distance between genomes is %7.5f +/- %7.5f.\n",$mean,$se;
printf LOGFILE "Valid neighbors are separted by distances less than %7.5f.\n\n",$validNeighbor;
print LOGFILE "There are ".($#equiv_sets+1+(keys %updatedSingletsHash))." Groups.\n";

$groupID = 1;
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


print "Genome	#of distributed genes\n";
foreach my $k (keys %hashOfGenomes){
	print "$k	$hashOfGenomes{$k}\n";
}

#close out the timer
$EndTime = time;
$ElapsedTime = $EndTime -$BeginTime;
print LOGFILE "\nThis run required $ElapsedTime seconds.\n";



print "Run completed.\n";

############################################
#Subroutines
############################################

sub isNNeighbor{
	my ($n1, $n2) = @_;
	if(($nearestNeighbor{$n1}[0] eq $n2)&&($nearestNeighbor{$n1}[1]<$validNeighbor)){
		return 1;
	}
	else {
		return 0;
	}
}
