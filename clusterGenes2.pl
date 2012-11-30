#!/usr/bin/perl
#
# USAGE:
# clusterGenes.perl all_fasta_sequences tfasty_output fasta_output > clusters
#
# AUTHOR: justin.s.hogg
# CREATED: 25 april 2006
# UPDATED: 1 june 2006
#
# set tab-stops to 3 spaces!

use strict;
use warnings;
use SequenceList;
use Sequence;



# This program has three inputs:
# (1)  all_fasta_sequences:  a multi-fasta file containing the nucleotide sequences of ORFS from all strains
# (2)  tfasty_output:  the output of a 'tfasty3x' all-against-all 6-frame translation homology search, format -m 9.
# (3)  fasta_output:  the output of a 'fasta3x' all ORFS versus all genomic sequence homology search, format -m 9.

# OVERVIEW:
# =========
# The script generates gene clusters based on pair-wise matches found by ssearch
# and tfasty.  The clustering algorithm is single-linkage with a minimum identity
# parameter and a minimum (assymetric) match length parameter.  Multi-sequence
# alignments are generated for each cluster by the external program 'poa'.  Please
# note that the scoring matrix for poa must be present in the run directory.  Clusters
# whose alignments are shorter than a set length are eliminated. Next, an attempt
# is made to split bad clusters by searching for a partition point that nicely
# divides the aligned sequences into two groups.  Sequences which 'scaffold' two disjoint
# groups are listed as adjuncts to the cluster.  Each cluster is assigned a name
# corresponding to a member gene according to the priority scheme set in the paramters.
# Next, any fasta matches between a cluster and genomic sequence are recorded.  Finally
# the output is written to a file.  Output includes:  the genes in the cluster and their
# sequences, the alignments of the genes, a list of adjuncts and their sequences, and a
# list of strains which have genomic homology to the cluster.




###/============\###
###| PARAMETERS |###
###\============/###

# clustering paramters for single-linkage
# june 15 settings for Hflu redo
our $MIN_AMINO_IDENTITY       = 0.70;     # minimum percent identity of amino acids over the matching region
our $MIN_NUCLEO_IDENTITY      = 0.70;     # minimum percent identity of nucleic acids over the matching region
our $MIN_ALIGNMENT_RATIO      = 0.70;     # minimum passing value of the fraction:  lengthOfMatchRegion / lengthOfSequence
our $MIN_CLUSTER_LENGTH       = 120;      # minimum length of cluster alignment, otherwise cluster is filtered
#$ENV{'PATH'} = '/bin:/usr/bin:/usr/local/bin:/root/poa_pkg/poa_pkg';


open FILE, "strain_list" or die "Could not open filename: $!\n";
my @strain_ar=<FILE>;
# cluster naming scheme
our $STRAIN_PRIORITY_1 = $strain_ar[0];
our $STRAIN_PRIORITY_2 = $strain_ar[1];
our $STRAIN_PRIORITY_3 = $strain_ar[2];
our $STRAIN_PRIORITY_4 = $strain_ar[3];
close FILE;


#Added by Josh, the functionality to read in a strain list from a file (woot)
open (STRAIN, "strain_list");
our @STRAIN_LIST ="";

while (<STRAIN>)
{
	chomp $_;
	push(@STRAIN_LIST, $_);
}
close STRAIN;

# no longer needed to list the strain names here:

# alignment options
our $GENERATE_ALIGNMENTS = 0;                
our $ALIGNMENT_INFILE  = '/var/www/html/supragenome/outputs/sequences.fa';
our $ALIGNMENT_OUTFILE = '/var/www/html/supragenome/outputs/sequences.aln';
our $ALIGNMENT_BINARY  = '/usr/local/poaV2/poa';
#our $ALIGNMENT_OPTIONS = "blosum80.mat -tolower -read_fasta $ALIGNMENT_INFILE -pir $ALIGNMENT_OUTFILE";
our $ALIGNMENT_OPTIONS = "/var/www/cgi-bin/blosum80_poa.mat -tolower -read_fasta $ALIGNMENT_INFILE -pir $ALIGNMENT_OUTFILE";

# cluster splitting options
# june15 settings for Hflu refo
our $MIN_GENES_FOR_SPLITTING = 7;
our $MAX_SCAFFOLDING_FOR_SPLITTING = 3;
our $MIN_CLUSTER_SIZE_FOR_SPLITTING = 3;
our $PARTITION_ENDS = 60;
our $MAX_PARTITION_OVERLAP = 40;
our $GAP_CHARACTER = '\.';





###/==============\###
###| MAIN PROGRAM |###
###\==============/###

my $step = 0;

# get command line parameters
my $sequence_file = shift @ARGV;
my $tfasty_file = shift @ARGV;
my $fasta_file = shift @ARGV;

# load sequences
print STDERR "${step}: loading sequences..\n";  $step++;
our $ALL_SEQUENCES = SequenceList->loadFasta( $sequence_file );

#My ignore list
my %ignore;
open IG, "ignore" or die "Can't open the ignore list; $!\n";
while(<IG>){
    chomp;
    $ignore{$_}=1;
}

# initialize @$GENE_PAIR_MATCHES
our $GENE_PAIR_MATCHES = [];
our $GENE_TO_GENOME_MATCHES = [];

# load tfasty matches
print STDERR "${step}: loading tfasty matches..\n";  $step++;
load_tfasty( $tfasty_file, $GENE_PAIR_MATCHES );

# load fasta matches
print STDERR "${step}: loading fasta matches..\n";  $step++;
load_fasta( $fasta_file, $GENE_TO_GENOME_MATCHES );

# initialize clustering data structures
our @CLUSTERS = ();
our %CLUSTER_INDEX = ();

# initialize cluster index
my $seq_ref;
foreach $seq_ref ( $ALL_SEQUENCES->getAllRefs )
{  $CLUSTER_INDEX{$seq_ref} = undef;  }




# first handle @$GENE_PAIR_MATCHES
# read matches from @matches and perform pairwise clustering 
print STDERR "${step}: performing single-linkage clustering..\n";  $step++;
foreach my $gene_pair ( @$GENE_PAIR_MATCHES )
{  
	# get arguments and parse match data
	my ($gene1_match, $gene2_match ) = @$gene_pair;

	#We don't want stuff from the ignore list ~josh
	if (exists $ignore{$gene1_match}||exists $ignore{$gene2_match}){next}

	# get references to sequences
	my ($gene1_ref, $gene2_ref) = $ALL_SEQUENCES->getRefs( $gene1_match, $gene2_match );

    # check if seq1 and seq2 are already part of a cluster(s).
	my $cluster1_ref = $CLUSTER_INDEX{$gene1_ref};
	my $cluster2_ref = $CLUSTER_INDEX{$gene2_ref};

	# decide how to assign sequences to bigClusters
	if ( defined $cluster1_ref  &&  defined $cluster2_ref )
	{
		if ( $cluster1_ref eq $cluster2_ref  )
		{
			# sequences assigned to same cluster
			# nothing to do
		}
		else
		{
			# merge clusters
			merge_clusters( $cluster1_ref, $cluster2_ref );
		}
	}
	elsif ( defined $cluster1_ref  )
	{
		# gene1 assigned to a cluster, but not gene2
		# add gene2 to $cluster1_ref
		&add_gene_to_cluster( $cluster1_ref, $gene2_ref );
	}
	elsif ( defined $cluster2_ref )
	{
		# gene2 assigned to a cluster, but not gene1
		# add gene1 to $cluster2_ref
		&add_gene_to_cluster( $cluster2_ref, $gene1_ref );
	}
	else 
	{
		# neither sequence is in a cluster
		my $cluster_ref = &new_cluster();
		&add_gene_to_cluster( $cluster_ref, $gene1_ref, $gene2_ref );
	}
}

# we're done with @GENE_PAIR_MATCHES
@$GENE_PAIR_MATCHES = ();
undef $GENE_PAIR_MATCHES;




# create singlet clusters for each sequence with no matches
print STDERR "${step}: creating singlet clusters..\n";  $step++;
foreach my $gene_ref ( $ALL_SEQUENCES->getAllRefs )
{
	unless ( defined $CLUSTER_INDEX{$gene_ref} )
	{
		my $cluster_ref = &new_cluster();
		&add_gene_to_cluster( $cluster_ref, $gene_ref );
	}
}




# generate names for clusters
print STDERR "${step}: generating cluster names..\n";  $step++;
foreach my $cluster_ref ( @CLUSTERS )
{
	next unless ( %$cluster_ref );
	$cluster_ref->{name} = &name_cluster( $cluster_ref );
}




# generate alignments, filter short clusters, split bad clusters
if ( $GENERATE_ALIGNMENTS )
{
	print STDERR "${step}: generating POA multi-alignments..\n";  $step++;
	foreach my $cluster_ref ( @CLUSTERS )
	{
		next unless ( %$cluster_ref );
		$cluster_ref->{alignments} = &generate_alignment( $cluster_ref );
	}


	# filter short clusters
	print STDERR "${step}: filtering short clusters..\n";  $step++;
	foreach my $cluster_ref ( @CLUSTERS )
	{
		next unless ( %$cluster_ref );
		&filter_short_clusters( $cluster_ref );
	}


	# split bad clusters
	print STDERR "${step}: splitting bad clusters..\n";  $step++;
	foreach my $cluster_ref ( @CLUSTERS )
	{
		next unless ( %$cluster_ref );
		&split_cluster( $cluster_ref );
	}
}




# add strains to existing clusters using @$GENE_TO_GENOME_MATCHES
print STDERR "${step}: loading strain matches..\n";  $step++;
foreach my $pair ( @$GENE_TO_GENOME_MATCHES )
{  
	# get arguments and parse match data
	my ($gene_match, $genome_match ) = @$pair;

	# get references to sequences
	my ($gene_ref) = $ALL_SEQUENCES->getRefs( $gene_match );

   # get this gene's cluster
	my $cluster_ref = $CLUSTER_INDEX{$gene_ref};
	next unless (defined $cluster_ref);

	if ( $cluster_ref eq 'adjunct' )
	{
		# do nothing.  this should be improved in future.
	}
	else
	{
		# add strain to the cluster
		&add_strain_to_cluster( $cluster_ref, $genome_match );
	}
}

# done with @$GENE_TO_GENOME_MATCHES, clear up memory
@$GENE_TO_GENOME_MATCHES = ();
undef $GENE_TO_GENOME_MATCHES;




# write cluster output
print STDERR "${step}: writing cluster output..\n";  $step++;
foreach my $cluster_ref ( @CLUSTERS )
{
	next unless ( %$cluster_ref );
	&output_cluster( $cluster_ref );
}




exit 0;



 
###/=============\###
###| SUBROUTINES |###
###\=============/###



sub load_tfasty
{
	my $tfasty_file = shift;
	my $gene_pair_matches = shift;

	open TFASTY, "<", $tfasty_file  or  die "ERROR: could not open tfasty file $tfasty_file!\n";
	while ( my $line_input = <TFASTY> )
	{
		# look for the beginning of a record
		if ( $line_input =~ /^\s*\d+>>>/ )
		{
			# we found a record.  get strain and locus.
		
			my ($gene1_name) = $line_input =~ />>>([A-Za-z0-9]+_\d+)/;
			my ($strain, $locus, $length) = $line_input =~ />>>([A-Za-z0-9]+)_(\d+)\s+-\s+(\d+) aa/; 
			
			# search for the start of the matches
			while ( $line_input = <TFASTY> )
			{
				last if ($line_input =~ /^The best scores are:/);
			}

			# gather matches until a blank line is encountered
			while ( $line_input = <TFASTY> ) 
			{
				last if ( $line_input =~ /^\n/ );
				chomp $line_input;

				my ($header,$data) = split( /\t/, $line_input );
				my ($gene2_name) = $header =~ /^([A-Za-z0-9]+_\d+)/;
				my ($match_strain, $match_locus) = $header =~ /^([A-Za-z0-9]+)_(\d+)/;

				my @match_data = split( /\s+/, $data );

				my $percent_identity = $match_data[0];
				my $alignment_length = $match_data[3];

				if ( $percent_identity >= $MIN_NUCLEO_IDENTITY
						and  ($alignment_length/$length) >= $MIN_ALIGNMENT_RATIO  )
				{
					# save the match
					push @$gene_pair_matches, [$gene1_name, $gene2_name];
				}
			}

			# we're at the end of the record, start looking for the next record	
		}
	}
	close TFASTY;
	return;
}


######
######
######


sub load_fasta
{
	my $fasta_file = shift;
	my $gene_to_genome_matches = shift;

	open FASTA, "<", $fasta_file  or  die "ERROR: could not open fasta file $fasta_file!\n";
	while ( my $line_input = <FASTA> )
	{
		# look for the beginning of a record
		if ( $line_input =~ /^\s*\d+>>>/ )
		{
			# we found a record.  get strain and locus.

		    my ($gene_name) = $line_input =~ />>>([A-Za-z0-9]+_\d+)/;
			my ($strain, $locus, $length) = $line_input =~ />>>([A-Za-z0-9]+)_(\d+)\s+-\s+(\d+) nt/; 
#			print "$length\n";

			# search for the start of the matches
			while ( $line_input = <FASTA> )
			{
				last if ($line_input =~ /^The best scores are:/);
			}

			# gather matches until a blank line is encountered
			while ( $line_input = <FASTA> ) 
			{
				last if ( $line_input =~ /^>>><<</ );
				chomp $line_input;
				if ($line_input=~/^\+-/){next}
				my ($header,$data) = split( /\t/, $line_input );				
				#print "$header\n$data\n";
				my ($match_strain) = $header =~ /^([A-Za-z0-9]+)/;

				my @match_data = split( /\s+/, $data );

				my $percent_identity = $match_data[0];
				my $alignment_length = $match_data[3];
				#print $percent_identity."\n";
				#print $alignment_length."\n";
				if ( $percent_identity >= $MIN_AMINO_IDENTITY
						and  ($alignment_length/$length) >= $MIN_ALIGNMENT_RATIO  )
				{
					# save the match
					push @$gene_to_genome_matches, [$gene_name, $match_strain];
				}
			}
			# we're at the end of the record, start looking for the next record	
		}
	}
	close FASTA;
	return;
}


######
######
######


sub new_cluster
{
	# create new cluster structure
	my $cluster_ref
		= {
			'name'       => undef,
			'genes'      => SequenceList->new(),
			'adjuncts'   => SequenceList->new(),
			'alignments' => undef,
			'strains'    => {},
		  };

	foreach my $strain ( @STRAIN_LIST )
	{
		$cluster_ref->{strains}->{$strain} = 0;
	}

	push @CLUSTERS, $cluster_ref;
	return $cluster_ref;
}


######
######
######


sub add_gene_to_cluster
{
	my $cluster_ref = shift;
	foreach my $gene_ref ( @_ )
	{
		# add gene
		$cluster_ref->{genes}->addSequence( $gene_ref );
	
		# add strain
		my ($strain) = $gene_ref->name =~ /^([A-Za-z0-9]+)_/;
		$cluster_ref->{strains}->{$strain} = 1;

		# update CLUSTER_INDEX
		$CLUSTER_INDEX{$gene_ref} = $cluster_ref;
	}
	return;
}


######
######
######


sub add_adjunct_to_cluster
{
	my $cluster_ref = shift;
	foreach my $gene_ref ( @_ )
	{
		# add gene
		$cluster_ref->{adjuncts}->addSequence( $gene_ref );
	
		# add strain
		my ($strain) = $gene_ref->name =~ /^([A-Za-z0-9]+)_/;
		$cluster_ref->{strains}->{$strain} = 1;

		# update CLUSTER_INDEX
		$CLUSTER_INDEX{$gene_ref} = 'adjunct';
	}
	return;
}



######
######
######


sub add_strain_to_cluster
{
	my $cluster_ref = shift;
	foreach my $strain ( @_ )
	{
		# add strain
		$cluster_ref->{strains}->{$strain} = 1;
	}
	return;
}


######
######
######


sub merge_clusters
{
	my $cluster1_ref = shift;
	my $cluster2_ref = shift;

	# update CLUSTER_INDEX and move genes from cluster2 to cluster1
	foreach my $gene_ref ( $cluster2_ref->{genes}->getAllRefs() )
	{
		$CLUSTER_INDEX{$gene_ref} = $cluster1_ref;
		$cluster1_ref->{genes}->addSequence( $gene_ref );
	}

	# add $cluster2 strains to $cluster1
	foreach my $strain ( keys %{$cluster2_ref->{strains}} )
	{
		if ( $cluster2_ref->{strains}->{$strain} == 1 )
		{ $cluster1_ref->{strains}->{$strain} = 1; }
	}

	# clear $cluster2
	%$cluster2_ref = ();

	return;	
}


######
######
######


sub min
{
	return undef unless (@_);

	my $min = shift;
	foreach my $value (@_) 
	{
		$min = $value if ( $value < $min );
	}
	return $min;
}

sub max
{
	return undef unless (@_);

	my $max = shift;
	foreach my $value (@_) 
	{
		$max = $value if ( $value > $max );
	}
	return $max;
}


######
######
######


sub name_cluster
{
	my $cluster_ref = shift;

	my $cluster_name = '';
	my $score = 0;

	foreach my $gene_ref ( $cluster_ref->{genes}->getAllRefs )
	{
		my $gene_name = $gene_ref->name;
		my ($strain) = $gene_name =~ /^([A-Za-z0-9]+)_/;

		if ( $strain =~ /$STRAIN_PRIORITY_1/i )
		{ 
			if ( $score < 5 ) 
				{ $score=5; $cluster_name = $gene_name; }
			elsif ($score == 5 && $gene_name lt $cluster_name)
				{ $cluster_name = $gene_name; }
		}
		elsif ( $strain =~ /$STRAIN_PRIORITY_2/i )
		{
			if ( $score < 4 ) 
				{ $score=4; $cluster_name = $gene_name; }
			elsif ($score == 4 && $gene_name lt $cluster_name)
				{ $cluster_name = $gene_name; }
		}
		elsif ( $strain =~ /$STRAIN_PRIORITY_3/i )
		{
			if ( $score < 3 ) 
				{ $score=3; $cluster_name = $gene_name; }
			elsif ($score == 3 && $gene_name lt $cluster_name)
				{ $cluster_name = $gene_name; }
		}
		elsif ( $strain =~ /$STRAIN_PRIORITY_4/i )
		{
			if ( $score < 2 ) 
				{ $score=2; $cluster_name = $gene_name; }
			elsif ($score == 2 && $gene_name lt $cluster_name)
				{ $cluster_name = $gene_name; }
		}
		else
		{
			if ( $score < 1 ) 
				{ $score=1; $cluster_name = $gene_name; }
			elsif ($score == 1 && $gene_name lt $cluster_name)
				{ $cluster_name = $gene_name; }
		}
	} 
	return $cluster_name;
}


######
######
######


sub filter_short_clusters
{
	my $cluster_ref = shift;
	return(undef) unless ( defined $cluster_ref->{alignments} );

	my ($align_obj) = $cluster_ref->{alignments}->getAllRefs;
	return(undef) unless ( defined $align_obj );

	if ( $align_obj->len < $MIN_CLUSTER_LENGTH )
	{
		foreach my $gene_ref ( $cluster_ref->{genes}->getAllRefs )
		{
			$CLUSTER_INDEX{$gene_ref} = undef;
		}

		# delete cluster
		%$cluster_ref = ();
	}

	## delete sequences too? ##

	return;
}


######
######
######


sub split_cluster
# split cluster using the partition method
{
	my $cluster_ref = shift;
	return unless ( scalar( $cluster_ref->{genes}->getNames ) >= $MIN_GENES_FOR_SPLITTING );
	return unless ( defined $cluster_ref->{alignments} );

	# make interval list
	my @intervals = ();
	my $align_length;
	foreach my $align_obj ( $cluster_ref->{alignments}->getAllRefs )
	{
		my $alignment = $align_obj->sequence;
		my ($seq_obj) = $ALL_SEQUENCES->getRefs( $align_obj->name );
		
		unless ( defined $align_length )
		{  $align_length = $align_obj->len;  }

		$alignment =~ m/[A-Za-z]/g;
		my $start_pos = pos($alignment) - 1;

		my ($stop_string) = $alignment =~ /(${GAP_CHARACTER}*)$/;
		my $stop_pos = $align_length - 1 - length($stop_string);

		push @intervals, { 'seq_obj' => $seq_obj, 'start' => $start_pos, 'stop' => $stop_pos };
	}

	# iterate through bases and look for a good partition, if any.
	SPLITTING_LOOP: for ( my $base_index = $PARTITION_ENDS - 1;
									$base_index < $align_length - $PARTITION_ENDS;
									$base_index += $MAX_PARTITION_OVERLAP )
	{
		my @left_sequences = ();
		my @right_sequences = ();
		my @split_sequences = ();

		# print STDERR "\nSPLITTING LOOP: $base_index | ";

		# partition sequences
		foreach my $seq ( @intervals )
		{
			if    ( $seq->{start} >= $base_index - $MAX_PARTITION_OVERLAP )
			{
			#	print STDERR "right: $seq->{seq_obj} | ";
				push @right_sequences, $seq->{seq_obj};
			}
			elsif ( $seq->{stop}  <= $base_index + $MAX_PARTITION_OVERLAP )
			{
			#	print STDERR "left: $seq->{seq_obj} | ";
				push @left_sequences,  $seq->{seq_obj};
			}
			else
			{
			#	print STDERR "split: $seq->{seq_obj}\n";
				push @split_sequences, $seq->{seq_obj};
			}
		}

		# evaluate partition
		if (	@split_sequences <= $MAX_SCAFFOLDING_FOR_SPLITTING
					and
				@left_sequences >= $MIN_CLUSTER_SIZE_FOR_SPLITTING
					and
				@right_sequences >= $MIN_CLUSTER_SIZE_FOR_SPLITTING )
 		{
			# this is a good place to split cluster
			print STDERR "## split cluster: $cluster_ref->{name}\n";

			# make new left cluster
			my $cluster_left = &new_cluster();
			foreach my $seq_obj ( @left_sequences )
			{
				&add_gene_to_cluster( $cluster_left, $seq_obj );
			}
			$cluster_left->{name} = &name_cluster( $cluster_left );
			$cluster_left->{alignments} = &generate_alignment( $cluster_left );
			print STDERR "   left cluster: $cluster_left->{name}\n";

			# make new right cluster
			my $cluster_right = &new_cluster();
			foreach my $seq_obj ( @right_sequences )
			{
				&add_gene_to_cluster( $cluster_right, $seq_obj );
			}
			$cluster_right->{name} = &name_cluster( $cluster_right );
			$cluster_right->{alignments} = &generate_alignment( $cluster_right );
			print STDERR "   right cluster: $cluster_right->{name}\n\n";

			# add adjuncts to left and right clusters
			foreach my $seq_obj ( @split_sequences )
			{
				&add_adjunct_to_cluster( $cluster_left, $seq_obj );
				&add_adjunct_to_cluster( $cluster_right, $seq_obj );
			}

			# clear old cluster
			%$cluster_ref = ();

			last SPLITTING_LOOP;
		}
	}
	return;
}


######
######
######


sub output_cluster
{
	my $cluster_ref = shift;

	print STDOUT "<cluster start> $cluster_ref->{name}\n";
	print STDOUT "  <genes start>\n";
	foreach my $gene_ref ( sort $cluster_ref->{genes}->getAllRefs )
	{
		$gene_ref->writeFasta( *STDOUT{IO} );
	}
	print STDOUT "  <genes stop>\n";

	if ( defined $cluster_ref->{alignments} )
	{
		print STDOUT "  <alignments start>\n";
		foreach my $gene_ref ( sort $cluster_ref->{alignments}->getAllRefs )
		{
			$gene_ref->writeFasta( *STDOUT{IO} );
		}
		print STDOUT "  <alignments stop>\n";
	}

	print STDOUT "  <adjuncts start>\n";
	foreach my $gene_ref ( sort $cluster_ref->{adjuncts}->getAllRefs )
	{
		$gene_ref->writeFasta( *STDOUT{IO} );
	}
	print STDOUT "  <adjuncts stop>\n";

	print STDOUT "  <strains start>\n";
	foreach my $strain ( sort keys %{$cluster_ref->{strains}} )
	{
		print STDOUT "    $strain\t$cluster_ref->{strains}->{$strain}\n";
	}
	print STDOUT "  <strains stop>\n";
	print STDOUT "<cluster stop>\n";

	return;
}


######
######
######


sub generate_alignment
{
	my $cluster_ref = shift;			
	my $command;
	my $alignment_list;

	# write cluster to multi-fasta format file
	$cluster_ref->{genes}->writeFasta( $ALIGNMENT_INFILE );

	# call multi-sequence alignment software
	$command = "$ALIGNMENT_BINARY $ALIGNMENT_OPTIONS";
	system( "$command" );

	if ( -e $ALIGNMENT_OUTFILE )
	{
		# load a list of alignment edges
		$alignment_list = SequenceList->loadFasta( $ALIGNMENT_OUTFILE );
		return $alignment_list;
	}

	return undef;
}


######
######
######
