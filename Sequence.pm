package Sequence;

# changes: 12apr06 - switched to string based sequence,quality storage

sub new
# create new sequence object
# and optionally define attributes as key-value pairs
{
	my $proto = shift;
	my $class = ref($proto) || $proto;
	my $self = {};

	bless( $self, $class );

	if ( ref($proto) eq 'Sequence' )
	{
		# copy object values from $proto
		$self->{name} = $proto->{name};
		$self->{comment} = $proto->{comment};
		@{$self->{sequence}} = $proto->{sequence};
		@{$self->{quality}} = $proto->{quality};
	}
	else
	{	
		$self->{name} = '';
		$self->{comment} = '';
		$self->{sequence} = '';
		$self->{quality} = '';
	}

	if (@_)
	{
		$self->setAttributes( @_ );
	}
	return $self;
}


sub name
# get or set name attribute
{
	my $self = shift;

	if (@_)
		{ $self->{name} = shift; }

	return $self->{name};
}


sub comment
# get or set comment attribute
{
	my $self = shift;

	if (@_)
		{ $self->{comment} = shift; }

	return $self->{comment};
}


sub sequence
# get or set sequence
{
	my $self = shift;

	# set sequence if another argument is defined
	if ( defined $_[0] )
	{ 
		$self->{sequence} = shift;
	
		# filter any white space
		$self->{sequence} =~ s/\s//g;
		
		#erase quality values to avoid corrupt data
		$self->{quality} = '';
	}
	return $self->{sequence};
}


sub complement
# get complement of sequence
# (this doesn't complement a quality sequence)
{
	my $self = shift;
	
	my $index;
	my $complement = '';
	my $base;

	for ( $index = length($self->{sequence}) - 1; $index >= 0; $index-- )
	{
		$base = substr( $self->{sequence}, $index, 1 );
		if    ( $base eq 'A' )
			{ $base = 'T'; }
		elsif ( $base eq 'C' )
			{ $base = 'G'; }
		elsif ( $base eq 'G' )
			{ $base = 'C'; }
		elsif ( $base eq 'T' )
			{ $base = 'A'; }
		elsif ( $base eq 'a' )
			{ $base = 't'; }
		elsif ( $base eq 'c' )
			{ $base = 'g'; }
		elsif ( $base eq 'g' )
			{ $base = 'c'; }
		elsif ( $base eq 't' )
			{ $base = 'a'; }

		$complement .= $base; 
	}

	return $complement;
}


sub quality
# get or set quality values
# assumes input string is delimited with white space
{
	my $self = shift;
	my $quality;
	my $base;

	# set quality values if another parameter is defined
	if ( defined $_[0] )
	{
		@$quality = ( split /\s+/, shift );
		$self->{$quality} ='';

		# make sure the quality scores correspond to the sequence
		return $self->{quality} unless ( @$quality == $self->len );

		foreach $base (@$quality)
		{
			# pad single digit quality values
			$base = sprintf( "%3s", $base);
			$self->{quality} .= $base;
		}
	}

	return $self->{quality};
}


sub len
# get the sequence length
{
	my $self = shift;
	return length($self->{sequence});
}


sub setAttributes
# set attributes as key-value pairs
{
	my $self = shift;

	my $item;
	my $undef = undef;
	my @output = ();

	while ( $item = shift )
	{
		if    ( $item eq 'name' )
		{ 
			push @output, $self->name(shift);
		}
		elsif ( $item eq 'comment' )
		{ 
			push @output, $self->comment(shift);
		}
		elsif ( $item eq 'sequence'  )
		{
			push @output, $self->sequence(shift);
		}
		elsif ( $item eq 'quality' )
		{ 
			push @output, $self->quality(shift);
		}
		else
		{
			push @output, $undef;
		}
	}

	return @output;
}


sub getAttributes
# get attributes
{
	my $self = shift;

	my $item;
	my $undef = undef;
	my @output = ();

	while ($item = shift)
	{
		if    ( $item eq "name" )
		{ 
			push @output, $self->name;
		}
		elsif ( $item eq "comment" )
		{ 
			push @output, $self->comment;
		}
		elsif ( $item eq "sequence"  )
		{
			push @output, $self->sequence;
		}
		elsif ( $item eq "quality" )
		{ 
			push @output, $self->quality;
		}
		else
		{
			push @output, $undef;
		}
	}

	return @output;
}


sub getSubsequence
# get a subsequence of the sequence
{
	my $self = shift;

	my $start;
	my $length;

	if (@_)
	{
		$start = shift;

		# make sure start is within bounds
		return undef if ( $start >= $self->len );
		$start = 0 if ( $start < 0 );

		if (@_)
		{
			$length = shift;
			# make sure length is within bounds
			if ( $length+$start > $self->len )
				{ $length = $self->len - $start; }
		}
		else
			{ $length = $self->len - $start; }
	}
	else
	{
		$start = 0;
		$length = $self->len;
	}

	return substr( $self->{sequence}, $start, $length );
}


sub getSubquality
# get a subsquence of the quality values
{
	my $self = shift;

	my $start;
	my $length;

	if (@_)
	{
		$start = shift;

		# make sure start is within bounds
		return undef if ( $start >= $self->len );
		$start = 0 if ( $start < 0 );

		if (@_)
		{
			$length = shift;
			# make sure length is within bounds
			if ( $length+$start > $self->len )
				{ $length = $self->len - $start; }
		}
		else
			{ $length = $self->len - $start; }
	}
	else
	{
		$start = 0;
		$length = $self->len;
	}

	# multiple by 3 to transform to quality coordinates
	$start *= 3;
	$length *= 3;

	return substr( $self->{quality}, $start, $length);

}


sub writeFasta
# format sequence as a fasta record and write to a file
{
	my $self = shift;

	my $lucy = 0;

	if ( $_[0] eq 'lucy' )
	{
		shift;
		$lucy = 1;
	}

	my $fastaFile = shift;
	my $qualityFile;
	my $fasta;
	my $quality;

	if ( $lucy )
	{
		($fasta, $quality) = $self->prepareFasta("line" => 60, "lucy");
	}
	else
	{
		($fasta, $quality) = $self->prepareFasta("line" => 60);
	}

	if ( ref($fastaFile) eq 'SCALAR' )
	{
		open FILE, ">$fastaFile";
		print FILE $fasta;
		close FILE;
	}
	elsif ( ref($fastaFile) eq 'IO::Handle' )
	{
		print $fastaFile $fasta;
	}

	if (@_ && $qualityFile)
	{
		$qualityFile = shift;
	
		if ( ref($qualityFile) eq 'SCALAR' )
		{
			open FILE, ">$qualityFile";
			print FILE $quality;
			close FILE;
		}
		elsif ( ref($qualityFile) eq 'IO::Handle' )
		{
			print $qualityFile $quality;
		}	
	}

	return 1;
}
	
sub prepareFasta
# generate sequence and quality fasta format entries
{
	# get reference to sequence object
	my $self = shift;

	# turn off lucy clear range extraction 
	my $lucy = 0;

	# set default lineWidth (entire sequence on one line)
	my $lineWidth = $self->len;

	# set default lucy clear range values
	my $clearStart = 0;
	my $clearStop = $self->len - 1;
	my $clearLength = $self->len;

	# declare local holding variables
	my $parameter;
	my $index;
	my $seqFasta;
	my $qualFasta;
	my $name;
	my $comment;

	# check for additional parameters
	while (@_)
	{
		$parameter = shift;

		if    ($parameter eq 'line')
		# set line width parameter
		{
			$lineWidth = shift if (@_);
			$lineWidth = $self->len unless ($lineWidth > 0);
		}
		elsif ($parameter eq 'lucy')
		# set lucy mode true
		{
			$lucy = 1;
		}
	}


	# setup header
	$name = $self->name;
	$comment = $self->comment;
	if ($lucy)
	{
		$seqFasta  = ">$name trimmed\n";
		$qualFasta = ">$name trimmed\n";
	}
	else
	{
		# setup header
		$seqFasta  = ">$name $comment\n";
		$qualFasta = ">$name $comment\n";
	}


	if ($lucy)
	# try to interpret lucy clear range
	{
		# extract clear range from comment
		$lucy = $self->comment =~ /0 (\d+) (\d+)$/;

		if ($lucy)
		{
			$clearStart = $1;
			$clearStop  = $2;
		
			# decrement to convert to firstIndex=0 coordinates
			$clearStart--;
			$clearStop--;
		
			# make sure clearStart and clearStop are in valid ranges
			$clearStart = 0 if ($clearStart < 0 or $clearStart >= $self->len );
			$clearStop = $self->len - 1 if ($clearStop < $clearStart or $clearStop >= $self->len);
		
			# calculate clearLength
			$clearLength = $clearStop - $clearStart + 1;
		}
	}

	for ( $index = $clearStart; $index <= $clearStop; $index += $lineWidth )
	{
		if ( $index + $lineWidth - 1 > $clearStop)
		{
			$lineWidth = $clearStop - $index + 1;
		}

		$seqFasta .= $self->getSubsequence($index, $lineWidth) . "\n";

		if ( $self->{quality} )
		# return a quality fasta record
			{ $qualFasta .= $self->getSubquality($index, $lineWidth ) . "\n"; }
		else
		# return an empty string
			{ $qualFasta = ''; }
	}

	return $seqFasta, $qualFasta;	
}



sub newFasta
# create a new sequence object using input from fasta format
{
	my $proto = shift;
	my $self;
	if ( ref($proto) )
	{
		$self = $proto->new;
	}
	else
	{
		$self = Sequence->new;
	}

	if (@_)
	{
		$self->loadFastaSequence( shift );

		if (@_) 
		{
			$self->loadFastaQuality( shift );
		}
	}

	return $self;
}


sub loadFastaSequence
# set the sequence attribute using a fast format sequence
# (if name and comment attribs are not set.. these are set from the fasta header)
{
	my $self = shift;

	my $fasta = shift;
	chomp $fasta;
	my @data = split( /\n/, $fasta );
	
	my $name;
	my $comment;

	my $header = shift @data;
	( $name, $comment ) = $header =~ /^>(\S+) ?(.*)/;

	$sequence = join('', @data);

	unless ( $self->name	)
		{ $self->name( $name ); }
	unless ( $self->comment )
		{ $self->comment( $comment ) if ($comment); }
	$self->sequence( $sequence );

	return 1;
}


sub loadFastaQuality
# set the quality attribute from a fasta format quality values
# (if name and comment attribs are not defined.. these are set from the fasta header.)
{
	my $self = shift;

	my $fasta = shift;
	chomp $fasta;
	my	@data = split( /\n/, $fasta );

	my $name;
	my $comment;

	my	$header = shift @data;
		( $name, $comment ) = $header =~ /^>(\S+) (.+)/;

	my	$quality = join(' ', @data);

	return 0 unless ( $self->name	eq $name );
	unless ( $self->comment )
		{ $self->name( $comment );	}
	return 0 unless ( $quality eq $self->quality( $quality ) );
	return 1;
}

1;
