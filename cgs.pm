package cgs;

sub makeFastaHash
# load a multiFasta file indicated by the first argument into a hash
# ( key value are fasta header names, values hold the associated sequence )
# returns the number of sequences loaded and a reference to the hash
{
	# get arguments
	my $inputFile = shift(@_);
   
   # initialize %fastaHash
   my %fastaHash = ();

   # initialize tracking index
   my $workIndex = 0;
   
	# initialize other local variables
	my $lineInput;
	my $contigName;
	my $sequence;

	# open fasta file
   open(FASTA,"<",$inputFile) || die "ERROR: could not open $inputFile!";

   # process lines in fasta file one at a time
   while ( $lineInput = <FASTA> )
	{
      # start next record if line begins with ">"        
      if ($lineInput =~ /^>/ )
		{
          #save old record first
         if ($workIndex)
			{
            $fastaHash{$contigName} = $sequence;           
         }                

         # set workIndex to true
         $workIndex = 1;

         # set contig length to 0 and sequence to null string
         $sequence = "";

         # get contig name
         ( $contigName ) = $lineInput =~ /^>(\S+)/;
		}
		# otherwise continue process current record
		else
		{
          #skip if record not yet found
          next unless ($workIndex);

          #append line to sequence
          chomp($lineInput);
          $sequence .= $lineInput;
      }

   }

   # save last record
   if ($workIndex)
	{
		$fastaHash{$contigName} = $sequence;          
   } 

   close(FASTA);
   return( $workIndex, \%fastaHash );
}



sub outputFasta
# given a reference to a sequenceHash and a filename,
# writes a multiFasta file to filename containing the
# sequences in sequenceHash.  Returns 0.
{
   # get pointer to contig array
   my $fastaHash = shift(@_);

   # get output file
   my $outputFile = shift(@_);

   # set line length parameter
   my $lineLength = 60;

	# declare other local variables
	my $sequenceName;
	my $headerLine;
	my $offset;
	my $sequence;

   # open output file
   open FASTAOUT, ">$outputFile" || die "ERROR: could not open $outputFile!";

   foreach $sequenceName ( keys %{$fastaHash} )
	{
		# make headerLine
      $headerLine = ">$sequenceName\n";
      print FASTAOUT $headerLine;
      
      $offset = 0;
      while ( $offset < length($fastaHash->{$sequenceName}) )
		{
         $sequence = substr( $fastaHash->{$sequenceName}, $offset, $lineLength ) . "\n";
         print( FASTAOUT $sequence);

         $offset += $lineLength;
      }    
   }

   close(FASTAOUT);
	return 0;
}


sub floor 
# returns the largest integer <= number
{
	my $number = undef;   
	
	if (@_)
	{
		$number = shift(@_);

	   return( $number ) if ( $number == int( $number ) );
	   return( int( $number ) );
	}
	else
	{
		return undef;
	}
}   


sub ceiling 
# returns the smallest integer >= number
{
	my $number = undef;

	if (@_)	
	{ 
		my $number = shift(@_);

	   return( $number ) if ( $number == int( $number ) );
	   return( int($number) + 1 );
	}
	else
	{	
		return undef;
	}
}


sub min
# given a list of numbers, returns the smallest
{
	my $min = undef;

	if (@_)	
	{
		$min = shift @_;
		foreach $value (@_) 
		{
			$min = $value if ( $value < $min );
		}
	}

	return $min;
}


sub max
# given a list of numbers, returns the largest
{
	my $max = undef;

	if (@_)
	{
		$max = shift @_;
		foreach $value (@_) 
		{
			$max = $value if ( $value > $max );
		}
	}

	return $max;
}


sub mean
# given a list of numbers, returns the arithmetic mean
{
	my $mean = 0;
	my $n = 0;
	my $value;

	foreach $value ( @_ )
	{
		$mean += $value;
		$n++;
	}

	if ($n)
		{ $mean = $mean/$n; }
	else
		{ $mean = undef; }

	return $mean;
}




sub variance
# given a list of numbers, returns the variance
{
	# get expectation value
	my $expect = cgs::mean( @_ );
	return undef unless (defined $expect);

	# declare and initialize local variables
	my $value;

	my $n = 0;
	my $variance = 0;

	# sum over squares of (value-expect)
	foreach $value ( @_ )
	{
		$variance += ($value - $expect)*($value - $expect);
		$n++;
	}

	# divide by number of samples
	if ($n)
		{ $variance = $variance/$n; }
	else
		{ $variance = undef; }

	return $variance;
}




sub covariance
# given references to two lists of paired values, returns the covariance
{
	my $X = shift @_;
	my $Y = shift @_;

	return undef unless ( @$X == @$Y );
	return undef unless ( @$X > 0 );

	my $covariance = 0;
	my $n = 0;

	my $meanX = cgs::mean( @$X );
	my $meanY = cgs::mean( @$Y );

	for ($index = 0; $index < @$X; $index++ )
	{
		$covariance += ( $X->[$index] - $meanX ) * ( $Y->[$index] - $meanY );
		$n++;
	}

	if ($n)
		{ $covariance = $covariance/$n; }
	else
		{ $covariance = undef; }

	return $covariance;
}



sub linearRegression
# given references to two lists of paired values, returns the linear regression constants
{
	my $X = shift @_;
	my $Y = shift @_;

	return undef unless ( @$X == @$Y );
	return undef unless ( @$X > 1 );

	my $expectX = cgs::mean( @$X );
	my $expectY = cgs::mean( @$Y );
	
	return undef unless ( $expectX );

	my $varianceX = cgs::variance( @$X );
	my $covarianceXY = cgs::covariance( $X, $Y );

	my $a = $covarianceXY / $varianceX;
	my $b = $expectY - ( $covarianceXY / $varianceX ) * $expectX;

	return( $a, $b );
}



sub median
# given a list of values, returns the median
{
	return undef unless (@_);
	my $values = shift @_;
	return undef unless ( @$values );

	my @sortedValues = ( sort { $a <=> $b} @$values );
	my $median;
	
	if ( @sortedValues % 2 )
	{
		$median = $sortedValues[ &cgs::floor(@sortedValues/2) ];
	}
	else
	{
		$median = ( $sortedValues[@sortedValues/2-1] + $sortedValues[@sortedValues/2] ) / 2;
	} 

	return $median;
}



sub tandem
{
	my $letter = shift;
	my $repeats = shift;

	my $tandem = "";
	for (1..$repeats)
	{
		$tandem .= $letter;
	}

	return $tandem;
}

1;
