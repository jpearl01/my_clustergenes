package pcr;
use cgs;

######
######
######

# version 04april2006

######
######
######



sub makeEndHash
# given a sequenceList object and a length parameter (bp), creates
# a new sequenceList of sequences generated from the ends of the original sequences.  
# returns 0 and the new sequenceList object
#
# notes: if a contig is short, the entire contig is added to the endHash
#        with the suffix "_Nend"
#
{
   # get pointer to contig array and length parameter
   my $contigs = shift @_;
	my $contigEndLength = shift @_;


	# create contigEnd SequenceList
   my $contigEnds = SequenceList->new;   

	# initialize other local variables
	my $sequence;
	my $contigEndName;
	my $endSequence;
	my $newEnd;

	# process each sequence in %{fastaHash}
   foreach $sequence ( $contigs->getAllRefs )
	{
		# determine length of contig end
      if ( $sequence->len >= 2 * $contigEndLength )
		{
	      # make 5' end of contig
			$newEnd = "";
	      $contigEndName = $sequence->name . "_5end";
	      ( $endSequence = $sequence->getSubsequence(0, $contigEndLength) ) =~ s/X/N/ig;
			$newEnd = Sequence->new( "name" => $contigEndName, "sequence" => $endSequence );
			$contigEnds->addSequence( $newEnd );

	      # make 3' end of contig
			$newEnd = "";
	      $contigEndName = $sequence->name . "_3end";
	      ( $endSequence = $sequence->getSubsequence($sequence->len - $contigEndLength) ) =~ s/X/N/ig;
			$newEnd = Sequence->new( "name" => $contigEndName, "sequence" => $endSequence );
			$contigEnds->addSequence( $newEnd );
      }
		elsif ( $sequence->len >= 3/2 * $contigEndLength )
		{
	      # make 5' end of contig
			$newEnd = "";
	      $contigEndName = $sequence->name . "_5end";
	      ( $endSequence = $sequence->getSubsequence(0, cgs::ceiling($sequence->len/2)) )=~ s/X/N/ig;
			$newEnd = Sequence->new( "name" => $contigEndName, "sequence" => $endSequence );
			$contigEnds->addSequence( $newEnd );

	      # make 3' end of contig
			$newEnd = "";
	      $contigEndName = $sequence->name . "_3end";
	      ( $endSequence = $sequence->getSubsequence($sequence->len - cgs::ceiling($sequence->len/2)) ) =~ s/X/N/ig;
			$newEnd = Sequence->new( "name" => $contigEndName, "sequence" => $endSequence );
			$contigEnds->addSequence( $newEnd );
		}
		else
		{
         # make generic "end" from entire contig
			$newEnd = "";
	      $contigEndName = $sequence->name . "_Nend";

	      ( $endSequence = $sequence->sequence ) =~ s/X/N/ig;
			$newEnd = Sequence->new( "name" => $contigEndName, "sequence" => $endSequence );
			$contigEnds->addSequence( $newEnd );
      }
   }

   return( 0, $contigEnds );

}



######
######
######



sub checkRead
# given a sequence object and a minConsecutiveNonMaskedBases parameter,
# returns 1 if the sequence passes the consecutiveNonMaskedBases test
# and returns 0 otherwise.
#
{
	# get arguments
	my $read = shift;
	my $minConsecutiveNonMaskedBases = shift;

	my $baseIndex;
	my $goodBaseCount = 0;

	for ( $baseIndex = 0; $baseIndex < $read->len; $baseIndex++ )
	{
		if ( $read->getSubsequence($baseIndex, 1) =~ /^[agct]$/i )
			{ $goodBaseCount++; }
		else
			{ $goodBaseCount=0; }
	
		return 1 if ( $goodBaseCount >= $minConsecutiveNonMaskedBases );
	}

	return 0;
}



######
######
######



sub maskLowQuality
# given two sequenceLists with a 1-1 correspondence of names and lengths,
# uses the first list as a guide in transforming uppercase bases to lowercase bases.
# note:  this kludge exists because crossmatch returns all caps.
# 
{
	my $contigEnds = shift;
	my $screenedContigEnds = shift;

	my $endName;
	my $endRef;
	my $screenedEndRef;
	my $newEnd;
	my $baseIndex;


	foreach $endName ( $contigEnds->getNames )
	{
		($endRef) = $contigEnds->getRefs( $endName );
		($screenedEndRef) = $screenedContigEnds->getRefs( $endName );

		$newEnd = "";
		
		for ( $baseIndex = 0; $baseIndex < $endRef->len; $baseIndex++ )
		{
			if ( $endRef->getSubsequence($baseIndex, 1) =~ /[agct]/ )
			{
				# make this base lowercase
				$newEnd .= "N";
			}
			else
			{
				$newEnd .= $screenedEndRef->getSubsequence($baseIndex, 1);
			}
		}
		
		$screenedEndRef->sequence( $newEnd );
	}

	return 0;
}




######
######
######



sub maskBeforeRepeats
# given a sequenceList object and a minMaskDistanceToPrimer parameter (#bp),
# masks any bases in a sequence up to "minMaskDistanceToPrimer" bases before (or after) the repeat
# masked region (after if 5end, before if 3end).
# returns 0.
{
	my $screenedContigEnds = shift;
	my $minMaskDistanceToPrimer = shift;

	my $endName;
	my $endRef;
	my $sequence;

	foreach $endName ( $screenedContigEnds->getNames )
	{
		# get reference to sequence object
		($endRef) = $screenedContigEnds->getRefs( $endName );

		if    ( $endName =~ /_5end$/ )
		{
			# get sequence
			$sequence = $endRef->sequence;

			# mask any bases before repeat-masked bases
			$sequence =~ s/(N+[AGCT]{1,$minMaskDistanceToPrimer})/cgs::tandem("N",length($1))/ige;

			# save changes
			$endRef->sequence( $sequence );
		}
		elsif ( $endName =~ /_3end$/ )
		{
			# get sequence
			$sequence = $endRef->sequence;

			# mask any bases before repeat-masked bases
			$sequence =~ s/([AGCT]{1,$minMaskDistanceToPrimer}N+)/cgs::tandem("N",length($1))/ige;

			# save changes
			$endRef->sequence( $sequence );
		}
	}

	return 0;
}



######
######
######



sub loadShowTilingData
# Given a reference to the moleculeHash, the name of a reference file and
# the filename of "show-tiling" output corresponding to the refFile:
# Reads and parses the show-tiling output and loads the data into
# the %moleculeHash data structure.
# returns 0.


# please note: when a tile is found with a "_Nend" suffix, a _3end tile
#              and a _5end tile are created from this data.  Both tiles
#              are identical aside from the name


# Here's an outline of a typical data structure:
#
# %moleculeHash = "molecule1" => reference to molecule1 dataHash
#                 "molecule2" => reference to molecule2 datahash
#
# %molecule1Hash = "file"	 => name of reference file
#						 "length" => length of molecule in bp
#						 "tiles"  => list of references to tileHashes
#						 "links"  => list of references to linkHashes (not created by this sub)
#
# %tileDataHash = "startRef"			=> start position of tile w.r.t. reference molecule
# 						"endRef"				=> stop position              ``
# 						"nextGapLength"	=> gap between this file and the next (to the right)
#						"contigLength"		=> length of tile
#                 "alignmentCover"  => percentLength of the tile matching the reference
#						"percentIdentity"	=> percentIdentity of tile match to reference
#						"orientation"		=> orientation of the match ("+"=uncomp, "-"=complement)
#						"contigID"			=> name of tile
					
{
	# get arguments
   my $moleculeHash = shift;
	my $refFile = shift;
	my $inputFile = shift;
	my $screenedEnds = shift;
	my $minConsecutiveNonMaskedBases = shift;

   # open fasta file
   open TILING,"<$inputFile" || die "ERROR: could not open $inputFile!";

   # initialize tracking index
   my $working = 0;
   
	# declare other local variables
	my $lineInput;
	my $molecule;
	my $length;
	my @tileData;
	my $tileDataHash;
	my $endRef;


   # process lines in fasta file one at a time
   while ( $lineInput = <TILING> )
	{
      # start next molecule if line begins with ">"        
      if ($lineInput =~ /^>/ )
		{
			# set working index
			$working = 1;
			
			# get molecule name and length
			chomp $lineInput;
			($molecule) = $lineInput =~ /^>(\S+)/;
			($length)   = $lineInput =~ /(\d+)\s+bases/;

         #initialize new molecule data structure
         $moleculeHash->{$molecule} =
			{
				'file'        => $refFile,
				'length'      => $length,
				'tiles'	     => [],
				'screenTiles' => [],
				'links'       => [],
				'type'        => 'reference'
			}; 
		}
      # otherwise continue processing current molecules
      else
		{
          #skip if group record not yet found
          next unless ($working);

          chomp $lineInput;
          @tileData = split /\t/, $lineInput;

			$tileDataHash =  {
									 "startRef"        => $tileData[0],
  	                         "endRef"          => $tileData[1],
  	                         "nextGapLength"   => $tileData[2],
  	                         "contigLength"    => $tileData[3],
  	                         "alignmentCover"  => $tileData[4],
  	                         "percentIdentity" => $tileData[5],
  	                         "orientation"     => $tileData[6],
  	                         "contigID"        => $tileData[7]
								  };

			($endRef) = $screenedEnds->getRefs($tileData[7]);

			if ( &checkRead( $endRef, $minConsecutiveNonMaskedBases ) )
				{ push( @{ $moleculeHash->{$molecule}->{tiles} }, $tileDataHash ); }
			else
				{ push( @{ $moleculeHash->{$molecule}->{screenTiles} }, $tileDataHash ); }
      }
   }

   close TILING;
   return 0;
}




######
######
######



sub findLinks
# given a reference to a list of tiles (all on the same molecule) and
# a reference to the "links" array for the respective molecule:  pairs 
# nieghboring (w.r.t. reference) tiles and passes the pairs to &checkLink 
# (which evaluates the pair and returns true if the pair forms a link).
# If the tiles form a valid link-pair, a linkHash is created and added 
# to the "links" array.
# returns 0.  
#
# note: see &checkLink for definition of a valid link-pair.
{
	my $tileList = shift;
	my $linkList = shift;
 	my $maxLinkGap = shift;
   my $maxLinkOverlap = shift;
	my $circularReference = shift;
	my $masterLinkHash = shift;
	my $contigRegex = shift;
	my $outputPrefix = shift;

   # search through tile pairs for valid links 
   # 
   # NOTE:  tiles are ordered by position on reference sequence,
   #        hence it's only necessary to check adjacent tiles for valid links.
   #
   #  5'=====================================3' refSeq
   #      | | | | |               | | | | |
   #     ---------->  ..link?..  ---------->
   #     leftTile                rightTile 
   #  
   # leftTile must be 3' end in sense orientation or 5' end in antisense orientation.
   # rightTile must be 5' end in sense orientation or 3' end in antisense orientation.
   # gap must be less than $maxLinkGap, and greater than -$maxLinkOverlap

	# declare local variables
   my @groupLinkArray = ();
	my $tile;
	my $leftTile; 	my $rightTile;
	my $linkHash;
	my $maxGap; 	my $maxOverlap;


   for ( $tile = 0; $tile < $#{$tileList}; $tile++ )
	{	
		$leftTile = $tileList->[$tile];
		$rightTile = $tileList->[$tile + 1];

		# add valid links to groupLinkArray
		if ( &checkLink( $leftTile, $rightTile, $maxLinkGap, $maxLinkOverlap) )
		{
         # define template name	
       	$leftTile->{contigID} =~ /$contigRegex/;
			$template  = "${outputPrefix}_$1-$2";
			$rightTile->{contigID} =~ /$contigRegex/;
			$template .= "_$1-$2";

			$linkHash = {
								"template"		=> $template,
								"leftTile"		=> $leftTile,
								"rightTile" 	=> $rightTile,
								"refSeq"			=> "",						
								"gapSequence"	=> SequenceList->new
							};

			push @{ $linkList }, $linkHash; 		
			$masterLinkHash->{$linkHash} = $linkHash;
		}

   }    

   if ( $circularReference )
	{ 
		# test last contig end against first contig end

		$leftTile = $tileList->[$#{$tileList}];
		$rightTile = $tileList->[0];

		# add valid links to groupLinkArray
		if ( &checkLink($leftTile, $rightTile, $maxLinkGap, $maxLinkOverlap) )
		{
         # define template name	
       	$leftTile->{contigID} =~ /$contigRegex/;
			$template  = "${outputPrefix}_$1-$2";
			$rightTile->{contigID} =~ /$contigRegex/;
			$template .= "_$1-$2";

			$linkHash = {
								'template'		=> $template,
								'leftTile'		=> $leftTile,
								'rightTile' 	=> $rightTile,
								'refSeq'			=> '',						
								'gapSequence'	=> SequenceList->new
							};

			push @{ $linkList }, $linkHash;
			$masterLinkHash->{$linkHash} = $linkHash; 
		}
	}

	return 0;
}



######
######
######



sub checkLink
# given two tiles (and two parameters), determines if the tiles
# are properly oriented neighbors which can be "contigged" by
# a simple PCR experiment.  Such tile pairs are hence referred
# to as "links"
# returns 1 if the pair is a valid link, otherwise 0.
#
# three criteria are used to evaluate the link:
#  1. orientation of tiles
#  2. gap length between tiles
#  3. overlap of tiles
#
{     
	# get arguments
   my $leftTile 			= shift;		# pointer to leftTile hash
   my $rightTile 			= shift;		# pointer to rightTile hash
   my $maxLinkGap 		= shift;		# maximum gap between contigs for a valid link
   my $maxLinkOverlap 	= shift;		# maximum overlap between contigs for a valid link

	# declare local variables
	my $leftTileEnd;  my $leftTileOrient;
	my $rightTileEnd; my $rightTileOrient;

   # check leftTile           
	unless ( $leftTile->{contigID} =~ /_Nend$/i )
	{
	   $leftTileEnd =  1 if ( $leftTile->{contigID} =~ /_3end$/i );
	   $leftTileEnd = -1 if ( $leftTile->{contigID} =~ /_5end$/i );

	   $leftTileOrient =  1 if ( $leftTile->{orientation} =~ /\+/ );
  		$leftTileOrient = -1 if ( $leftTile->{orientation} =~ /-/  );

	   return 0 unless ( $leftTileEnd * $leftTileOrient * 1 == 1 );
	}


   # check rightTile           
	unless ( $rightTile->{contigID} =~ /_Nend$/i )
	{
	   $rightTileEnd =  1 if ( $rightTile->{contigID} =~ /_3end$/i );
	   $rightTileEnd = -1 if ( $rightTile->{contigID} =~ /_5end$/i );

	   $rightTileOrient =  1 if ( $rightTile->{orientation} =~ /\+/ );
	   $rightTileOrient = -1 if ( $rightTile->{orientation} =~ /-/ );

	   return 0 unless ( $rightTileEnd * $rightTileOrient * -1 == 1 );
	}      

   # check gap length
   return 0 unless ( $leftTile->{nextGapLength} <= $maxLinkGap );
   return 0 unless ( $leftTile->{nextGapLength} >= -$maxLinkOverlap );
  
   return 1;	# link is valid
}



######
######
######



sub findGenericLinks
# loads generic links from $linkFile
# format of generic links:  contig1_end  contig2_end  gap_length
{
   my $moleculeHash = shift;
 	my $linkFile = shift;
	my $screenedEnds = shift;
	my $maxLinkGap = shift;
   my $maxLinkOverlap = shift;
	my $masterLinkHash = shift;
	my $contigRegex = shift;
	my $outputPrefix = shift;

	# declare local variables
	my $lineInput;
	my @linkData;
	my $left_orient;
	my $right_orient;
	my $endRef;
	my $leftTile;
	my $rightTile;
	my $linkHash;
	my $template;


	#initialize generic molecule data structure
	$moleculeHash->{generic} =
	{
		'file'        => undef,
		'length'      => undef,
		'tiles'	     => [],
		'screenTiles' => [],
		'links'       => [],
		'type'        => 'generic'
	}; 

	open( LINKFILE, "<", $linkFile ) or die "ERROR: could not open link file $linkFile!\n";
	while ($lineInput = <LINKFILE>)
	{
		# skip comment lines
		next if ($lineInput =~ /^#/);

		# parse link data
		chomp $lineInput;
		@linkData = split( /\s+/, $lineInput );

		# check the gap length criterion
		next unless ( $linkData[2] < $maxLinkGap  and $linkData[2] > $maxLinkOverlap );


		# get orientation for left tile
		if ( $linkData[0] =~ /_3end$/ )
			{ $left_orient = '+'; }
		else
			{ $left_orient = '-'; }
	
		# get orientation for right tile
		if ( $linkData[1] =~ /_5end$/ )
			{ $right_orient = '+'; }
		else
			{ $right_orient = '-'; }


		# see if the linked contigs are present
		($endRef) = $screenedEnds->getRefs($linkData[0]);
		unless ( defined $endRef )
		{
			$linkData[0] =~ s/_[35]end$/_Nend/;
			($endRef) = $screenedEnds->getRefs($linkData[0]);
			next unless (defined $contigRef);
		}

		$leftTile = {
							"startRef"        => undef,
							"endRef"          => undef,
                  	"nextGapLength"   => $linkData[2],
                  	"contigLength"    => undef,
                  	"alignmentCover"  => undef,
                  	"percentIdentity" => undef,
                  	"orientation"     => $left_orient,
                     "contigID"        => $linkData[0],
						};

		if ( &checkRead( $endRef, $minConsecutiveNonMaskedBases ) )
			{ push( @{ $moleculeHash->{generic}->{tiles} }, $leftTile ); }
		else
			{ push( @{ $moleculeHash->{generic}->{screenTiles} }, $leftTile ); }


		($endRef) = $screenedEnds->getRefs($linkData[1]);
		unless ( defined $endRef )
		{
			$linkData[1] =~ s/_[35]end$/_Nend/;
			($endRef) = $screenedEnds->getRefs($linkData[1]);
			next unless (defined $endRef);
		}

		$rightTile = {
							"startRef"        => undef,
							"endRef"          => undef,
                  	"nextGapLength"   => undef,
                  	"contigLength"    => undef,
                  	"alignmentCover"  => undef,
                  	"percentIdentity" => undef,
                  	"orientation"     => $right_orient,
                     "contigID"        => $linkData[1],
						};

		if ( &checkRead( $endRef, $minConsecutiveNonMaskedBases ) )
			{ push( @{ $moleculeHash->{generic}->{tiles} }, $rightTile ); }
		else
			{ push( @{ $moleculeHash->{generic}->{screenTiles} }, $rightTile ); }


		
		# add link to groupLinkArray
      # define template name	
     	$linkData[0] =~ /$contigRegex/;
		$template = "${outputPrefix}_$1-$2";
		$linkData[1] =~ /$contigRegex/;
		$template .= "_$1-$2";

		$linkHash = {
							'template'		=> $template,
							'leftTile'		=> $leftTile,
							'rightTile' 	=> $rightTile,
							'refSeq'			=> undef,						
							'gapSequence'	=> SequenceList->new,
						};

		push @{ $moleculeHash->{generic}->{links} }, $linkHash; 		
		$masterLinkHash->{$linkHash} = $linkHash;
		  
	}
	close LINKFILE;

	return 0;
}



######
######
######



sub makeCombinedLinks
# given a reference to %referenceLinkHash, clusters all equivalent links
# and returns a reference to the non-redundant %combinedLinks.
#
# notes:
# 1. %moleculeHash consists of keys which are molecule names, 
#    and values which are references to the associated data hash.
# 2. equivalent links are valid-links formed by the same pair of tiles.
#
{
   my $moleculeHash = shift;
	my $outputPrefix = shift;
   my %combinedLinks = (); 

	my $newLink;
	my $link;
	my $template;		
	my $gap;
	my $item;

   foreach $molecule ( keys %{$moleculeHash} )
	{
		# get reference to molecule data hash
		$moleculeRef = $moleculeHash->{$molecule};

      LINKLOOP: foreach $link ( @{$moleculeRef->{links}} )
		{
         # get link info
         $template = $link->{template};
			$gap = $link->{leftTile}->{nextGapLength};

         #check if link is already in @combined Links
         SEARCHLOOP: foreach $item ( values %combinedLinks )
			{
            # check if contig IDs match
            if    ( $item->{contig1}->{contigID} eq $link->{leftTile}->{contigID} 
							and $item->{contig2}->{contigID} eq $link->{rightTile}->{contigID} )
				{
					next SEARCHLOOP unless ( $link->{leftTile}->{orientation} eq $item->{contig1}->{orientation} );
					next SEARCHLOOP unless ( $link->{rightTile}->{orientation} eq $item->{contig2}->{orientation} );

					# at this point we've found the link is already in @combinedLinks
	  	         # update evidence and move on
	  	         $item->{evidence}++;
	  	         $item->{gap}->{$molecule} = $gap;            
					next LINKLOOP;
				}
				elsif ( $item->{contig1}->{contigID} eq $link->{rightTile}->{contigID}
							and $item->{contig2}->{contigID} eq $link->{leftTile}->{contigID} )
				{
					next SEARCHLOOP unless ( $link->{leftTile}->{orientation} ne $item->{contig2}->{orientation} );
					next SEARCHLOOP unless ( $link->{rightTile}->{orientation} ne $item->{contig1}->{orientation} );

					# at this point we've found the link is already in @combinedLinks
	  	         # update evidence and move on
					push @{$item->{links}}, $link;
	  	         $item->{evidence}++;
	  	         $item->{gap}->{$molecule} = $gap;            
					next LINKLOOP;
				}								 			
         }
 
         # if we've made it this far, the link is not in @combinedLinks
         # add $link to @combinedLinks

			$newLink =
			{ 
				'links'			=> [ $link ],							
				'contig1'		=> $link->{leftTile},
				'contig2'		=> $link->{rightTile},
				'evidence' 		=> 1,
				'gap' 			=> { $molecule => $gap },
				'maxGap'			=> undef,
				'minGap'			=> undef,
				'defined'		=> undef,
				'productSize' 	=> {},
				'maxProduct'	=> undef,
				'minProduct'   => undef,
				'primerLname' 	=> undef,
				'primerRname' 	=> undef,
				'primerLseq' 	=> undef,
				'primerRseq' 	=> undef,
				'primerLtm' 	=>	undef,
				'primerRtm' 	=> undef,
				'primerLgc' 	=> undef,
				'primerRgc' 	=>	undef,
				'primerLpos'	=> undef,
				'primerRpos'	=> undef,
				'plateNumber'	=> undef,
				'pcrPlateName' => undef,
				'plateLname'	=> undef,
				'plateRname'	=> undef,
				'well'			=> undef,
			};

         $combinedLinks{$template} = $newLink;
      } # end LINKLOOP
   } # end foreach loop

	# calculate gap length statistics
	foreach $link ( values %combinedLinks )
	{
		$link->{maxGap} = cgs::max( values %{$link->{gap}} );
		$link->{minGap} = cgs::min( values %{$link->{gap}} );	
	}

   return( \%combinedLinks );
}



######
######
######



sub makePrimerTemplates
# given: reference to %combinedLinksHash
#			reference to $screenedContigEnds (sequenceList object)
#			output filename
#			a filler string (typically a series of N's)
#			a primer3 configuration string
#			an integer corresponding to the length of an excluded primer design region
#				(distance between end of contig and start of valid primer design sequence)
# produces a pcr Primer design file suitable for input to primer3.
# returns 0.  Dies if file cannot be opened for write.						
#
{
	# get arguments
	my $combinedLinksHash 	= shift;		# reference to combinedLinksHash
	my $screenedContigEnds	= shift;		# reference to screenedContigEnds object
	my $primerInputFile		= shift;		# name of primer input file to create
	my $nFillString			= shift;		# string used to separate sequences in design template
	my $primer3configure		= shift;		# configuration data for primer3
	my $excludeRegion			= shift;		# length of region to exclude at end of contigs

	my $template;				# key name for contig link
	my $link;					# holds reference to contig links
	my $sequence1;
	my $sequence2;
	my $endRef1;				# reference to left end sequence
	my $endRef2;				# reference to right end sequence
	my $length1;				# length of sequence1
	my $length2;				# length of sequence2
	my $excludeStart;			# start of primer excluded region
	my $excludeStop;			# length of primer excluded region
	my $targetStart;			# start of pcr target region
	my $targetStop;			# length of pcr target region
	my $templateName;			# name of pcr template
	my $templateSequence;	# sequence

	# initialize counting index
	my $index = 0;

	# open primer design file for output
	open PRIMERIN, ">$primerInputFile" || die "ERROR: could not open primer design file $primerInputFile for output!\n";

	# make design template for each link in the array
	foreach $template ( keys %$combinedLinksHash )
	{
		$link = $combinedLinksHash->{$template};

		# get sequences and lengths
		($endRef1) = $screenedContigEnds->getRefs($link->{contig1}->{contigID});
		($endRef2) = $screenedContigEnds->getRefs($link->{contig2}->{contigID});
	
		$length1 	= $endRef1->len;
		$length2 	= $endRef2->len;	

		# decide whether to output uncomplement or complement
		if ( $link->{contig1}->{orientation} eq '-' )
			{ $sequence1 = $endRef1->complement; }
		else
			{ $sequence1 = $endRef1->sequence; }

		if ( $link->{contig2}->{orientation} eq '-' )
			{ $sequence2 = $endRef2->complement; }
		else
			{ $sequence2 = $endRef2->sequence; }

		# make design template sequence
		$templateSequence = $sequence1 . $nFillString . $sequence2;

		# define target region
		$targetStart 	= $length1 + 1;
		$targetStop	 	= length($nFillString);

		# define excluded region
		$excludeStart = $length1 + 1 - $excludeRegion;
		$excludeStop 	= length($nFillString) + 2 * $excludeRegion;

		# format primer3 entry
	   my $primer3output = "PRIMER_SEQUENCE_ID=$template\n"
   							. "SEQUENCE=$templateSequence\n"
   							. "TARGET=$targetStart,$targetStop\n"
   							. "EXCLUDED_REGION=$excludeStart,$excludeStop\n";

		if ($index)
		{
			# output entry without configuration data
			$index++;
			print PRIMERIN $primer3output;
			print PRIMERIN "=\n";
		}
		else
		{
			# output entry with configuration data
			$index++;
			print PRIMERIN $primer3output;
			print PRIMERIN $primer3configure;
			print PRIMERIN "=\n";
		}

	} # end while

	return 0;
}



######
######
######



sub loadPrimers
# given: primer3 output filename
#			filename prefix string
#        reference to %combinedLinksHash
#        filler string
# retrieves primer design info from a primer3 output file and stores
# the resulting information in the appropriate link in %combinedLinksHash.
# returns 0.  Dies if input file cannot be opened.  No format error checking!
#
{
	# get arguments
	my $primerOutputFile = shift;
	my $outputPrefix = shift;
	my $combinedLinksHash = shift;
	my $nFillString = shift;

	# declare local variables
	my $lineInput; 	my $template;
	my $productSize; 	

	my $targetStart;	my $targetLength;
	my $primerLpos;	my $primerRpos;

	# open primer design file
	open PRIMER, "<$primerOutputFile";

   # process lines in fasta file one at a time
   while ( $lineInput = <PRIMER> )
	{
      # start next record if line begins with "="        
      if ($lineInput =~ /^=/ )
		{
         # save last primer design set if a template was found
			next unless ( defined $template );

			# set defined flag if both primers have been found
			$combinedLinksHash->{$template}->{defined} = 1 
				if ( defined $combinedLinksHash->{$template}->{primerLseq} 
					&& defined $combinedLinksHash->{$template}->{primerRseq} );

			# calculate product size
			if (defined $productSize && defined $nFillString)
			{
				foreach $key ( keys %{$combinedLinksHash->{$template}->{gap}} )
				{
					$combinedLinksHash->{$template}->{productSize}->{$key}
						= $productSize + cgs::max(0,$combinedLinksHash->{$template}->{gap}->{$key}) - length($nFillString);
				}

				$combinedLinksHash->{$template}->{maxProduct} = cgs::max( values %{$combinedLinksHash->{$template}->{productSize}} );
				$combinedLinksHash->{$template}->{minProduct} = cgs::min( values %{$combinedLinksHash->{$template}->{productSize}} );
			}
 
			# calculate primer location w.r.t end of contig
			if ( defined $targetStart && defined $primerLpos )
			{
				$combinedLinksHash->{$template}->{primerLpos} = $targetStart - $primerLpos;

				if ( defined $targetLength && defined $primerRpos )
				{
					$combinedLinksHash->{$template}->{primerRpos} 
						= $primerRpos - ($targetStart+$targetLength);
				}
			}
               
         # reset variables to ward off data mixups;
			$template = undef;		$productSize = undef;      
			$targetStart = undef;	$targetLength = undef;
			$primerLpos  = undef;	$primerRpos = undef;
      } 
		else
		{
         if ($lineInput =~ /PRIMER_SEQUENCE_ID=/)
			{
         	( $template ) = $lineInput =~ /=(\S+)/;
				$combinedLinksHash->{$template}->{primerLname} = "${template}_1F";
				$combinedLinksHash->{$template}->{primerRname} = "${template}_1R";
			}
			elsif ($lineInput =~ /^PRIMER_LEFT_SEQUENCE=/ )
				{ ( $combinedLinksHash->{$template}->{primerLseq} ) = $lineInput =~ /=([agct]+)/i; }
			elsif ($lineInput =~ /^PRIMER_RIGHT_SEQUENCE=/ )
				{ ( $combinedLinksHash->{$template}->{primerRseq} ) = $lineInput =~ /=([agct]+)/i; }
			elsif ($lineInput =~ /^PRIMER_LEFT_TM=/ )
			{
         	( $combinedLinksHash->{$template}->{primerLtm} ) = $lineInput =~ /=([\d\.]+)/;
				$combinedLinksHash->{$template}->{primerLtm} 
					= sprintf("%.1f", $combinedLinksHash->{$template}->{primerLtm});
         }
			elsif ($lineInput =~ /^PRIMER_RIGHT_TM=/ )
			{
         	( $combinedLinksHash->{$template}->{primerRtm} ) = $lineInput =~ /=([\d\.]+)/;
				$combinedLinksHash->{$template}->{primerRtm} 
					= sprintf("%.1f", $combinedLinksHash->{$template}->{primerRtm});
         }
			elsif ($lineInput =~ /^PRIMER_LEFT_GC_PERCENT=/ )
			{
            ( $combinedLinksHash->{$template}->{primerLgc} ) = $lineInput =~ /=([\d\.]+)/;
				$combinedLinksHash->{$template}->{primerLgc}
					= sprintf("%.0f", $combinedLinksHash->{$template}->{primerLgc});
         }
			elsif ($lineInput =~ /^PRIMER_RIGHT_GC_PERCENT=/ )
			{          
         	( $combinedLinksHash->{$template}->{primerRgc} ) = $lineInput =~ /=([\d]+)/;
				$combinedLinksHash->{$template}->{primerRgc} 
					= sprintf("%.0f", $combinedLinksHash->{$template}->{primerRgc});
         }
			elsif ($lineInput =~ /^PRIMER_PRODUCT_SIZE=/ )
				{ ( $productSize ) = $lineInput =~ /=([\d]+)/; }
			elsif ($lineInput =~ /^TARGET=/ )
				{ ( $targetStart, $targetLength ) = $lineInput =~ /=(\d+),(\d+)/; }
			elsif ($lineInput =~ /^PRIMER_LEFT=/ )
				{ ( $primerLpos ) = $lineInput =~ /=(\d+)/; }
			elsif ($lineInput =~ /^PRIMER_RIGHT=/ )
				{ ( $primerRpos ) = $lineInput =~ /=(\d+)/;	}
      }
   }

   return 0;
}



######
######
######



sub getPCRassignment
# Given a pcr reaction number and a total number of reactions,
# determines the plate and well position corresponding the the reaction.
# This data is used for primer ordering and pcr setup.
# Returns plate number and well string.
# 
# notes:
# 1. reactions are evenly divided between the minimum number of 96-well plates.
# 2. plates are filled in across letter columns, then down numbered rows.
#
{
   # get arguments
   my $reaction = shift(@_);
   my $totalReactions = shift(@_);
 
   # declare local variables
	my $plates;
   my $reactionsPerPlate;  
   my $pcrPlate;
	my $pcrWellNumber;
	my $pcrColNumber;
	my $pcrRowNumber;
	my $index;
	my $pcrRowLetter;
	my $pcrWell;

   # define parameters
   my $rows = 12;		# numbered rows
   my $columns = 8;	# lettered columns

   # calculate number of plates
   $plates = cgs::ceiling( $totalReactions / ( $rows * $columns) );  

   # calculate number of reactions per plate
	$reactionsPerPlate = cgs::ceiling( $totalReactions / $plates );
   

   # calculate primer plate and well
   $pcrPlate = cgs::ceiling( $reaction / $reactionsPerPlate );
   $pcrWellNumber = $reaction - ($pcrPlate-1) * $reactionsPerPlate;
   $pcrColNumber = cgs::ceiling( $pcrWellNumber / $columns );
   $pcrRowNumber = $pcrWellNumber - ($pcrColNumber-1)*$columns;

   #calculate column letter
   $pcrRowLetter = "A";
   for ( $index=1; $index < $pcrRowNumber; $index++ )
   {
      $pcrRowLetter++;
   }

   #create pcrWell
   $pcrWell = $pcrRowLetter;
   $pcrWell .= "0" if ( $pcrColNumber < 10 );
   $pcrWell .= $pcrColNumber;

   return( $pcrPlate, $pcrWell );
}



######
######
######



1;
