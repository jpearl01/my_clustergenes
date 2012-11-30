package SequenceList;
use Sequence;

# version 04april2006

sub new
# creates a new list hash.
# optionally loads sequences if the method is called with a
# list of references to sequence objects.
{
	my $proto = shift;
	my $class = ref($proto) || $proto;
	my $self = {};
	bless( $self, $class );

	# copy constructor
	if ( ref($proto) eq 'SequenceList' )
	{
		%$self = %$proto;
	}

	my $name;
	my $ref;

	# add sequences to the list
	while ($ref = shift)
	{
		next unless (ref($ref) eq "Sequence");
		$name = $ref->name;
		$self->{name} = $name;
	}

	return $self;
}


sub renameSequence
{
	my $self = shift;

	my $oldName;
	my $newName;
	my $seqRef;

	if (@_)
	{
		$oldName = shift;
		if (@_)
		{
			$newName = shift;
		}
		else
		{
			print STDERR 'SequenceList::renameSequence("newName","oldName") - ERROR: invalid usage!\n';
			return 1;
		}
	}
	else
	{
		print STDERR 'SequenceList::renameSequence("newName","oldName") - ERROR: invalid usage!\n';
		return 1;
	}

	# update Sequence name
	( $seqRef ) = $self->getRefs( $oldName );

	$seqRef->name( $newName );
	
	$self->deleteSequence( $oldName );
	$self->addSequence( $seqRef );

	return $seqRef->name;
}


sub getNames
{
	my $self = shift;
	my @names = ( keys %$self );
	return ( keys %$self );
}


sub getRefs
{
	my $self = shift;

	my @list;
	while ( $name = shift )
	{
		if ( exists $self->{$name}	)
		{
			push @list, $self->{$name};
		}
		else
		{
			push @list, undef;
		}
	}	

	return @list;
}


sub getAttributes
{
	my $self = shift;

	my $attribute;
	my @list;


	(my $ref) = $self->getRefs(shift);

	if ( defined $ref )
	{
		while ( $attribute = shift )
		{
			push @list, $ref->$attribute;
		}	
	}
	else
	{
		while ( $attrbute = shift )
		{
			push @list, undef;
		}
	}

	return @list;
}



sub getAllRefs
{
	my $self = shift;
	return ( values %{$self} );
}


sub addSequence
# adds a sequence to the list
# over-writes old sequence if the name is already in the list
{
	my $self = shift;

	my $name;
	my $seqObject;

	while ($seqObject = shift)
	{
		next unless ( ref($seqObject) eq 'Sequence' );
		$name = $seqObject->name;
		$self->{$name} = $seqObject;
	}
}


sub deleteSequence
# deletes a sequence to the list
{
	my $self = shift;
	my $name;

	while ($seq = shift)
	{
		delete $self->{$seq} if ( exists $self->{$seq} );
	}
}


sub loadFasta
# loads a multiFasta file and creates a new sequenceList
{
	my $proto = shift;
	my $class = ref($proto) || $proto;
	my $self;

	# create new SequenceList unless $proto is a SequenceList
	if ( ref($proto) eq 'SequenceList' )
	{
		$self = $proto;
	}
	else
	{
		$self = {};
		bless( $self, $class );
	}

	my $file;
	my $lineInput;
	my $fasta;
	my $newSeq;

	# load fasta file
	if (@_)
	{ 
		$file = shift;
		open FILE, "<$file";

		$fasta = "";

		while( $lineInput = <FILE> )
		{
			if ( $lineInput =~ /^>/ )
			{
				if ($fasta)
				{				
					$newSeq = Sequence->newFasta($fasta);
					$self->addSequence($newSeq);
				}
		
				$fasta = $lineInput;
			}
			else
			{
				$fasta .= $lineInput;
			}
		}		

		# get that last entry
		if ( $fasta	)
		{
			$newSeq = Sequence->newFasta($fasta);
			$self->addSequence($newSeq);
		}
	}

	# load quality file
	if (@_)
	{ 
		$file = shift;
		open FILE, "<$file";
		
		$fasta = "";

		while( $lineInput = <FILE> )
		{
			if ( $lineInput =~ /^>/ )
			{
				if ($fasta)
				{
					($name) = $fasta =~ /^>(\S+)/;
					$newSeq = $self->{$name};
					$newSeq->loadFastaQuality($fasta);
				}
		
				$fasta = $lineInput;
			}
			else
			{
				$fasta .= $lineInput;
			}
		}		

		# get that last entry
		if ( $fasta	)
		{
			($name) = $fasta =~ /^>(\S+)/;
			$newSeq = $self->{$name};
			$newSeq->loadFastaQuality($fasta);
		}
	}

	return $self;
}


sub writeFasta
# write a sequence list in multi-fasta format
{
	my $self = shift;

	my $fastaFile;
	my $qualityFile;
	my $fasta;
	my $quality;

	my $sort;
	my $seq;
	my $seqRef;

	my $lucy = 0;

	if ( ref($_[0]) eq 'CODE'  )
	{
		$sort = shift;
	}

	if ( $_[0] eq 'lucy' )
	{
		# set lucy parameter
		shift;
		$lucy = 1;
	}

	$fastaFile = shift;
	open FASTAOUT, ">$fastaFile"  or  die "ERROR: could not open $fastaFile for output!\n";

	if (@_)
	{
		$qualityFile = shift;
		open QUALOUT, ">$qualityFile";
	}

	


	if (defined $sort)
	{
		foreach $seq ( sort { &$sort($a,$b) } keys %{$self} )
		{
			print STDERR "$seq\n";
			$seqRef = $self->{$seq};

			if ($lucy)
			{
				($fasta,$quality) = $seqRef->prepareFasta('line' => 60, 'lucy');
			}
			else
			{
				($fasta,$quality) = $seqRef->prepareFasta('line' => 60);
			}

			print FASTAOUT $fasta;

			if ($qualityFile && $quality)
				{ print QUALOUT $quality; }
		}
	}
	else
	{
		foreach $seq ( keys %{$self} )
		{
			$seqRef = $self->{$seq};

			if ($lucy)
			{
				($fasta,$quality) = $seqRef->prepareFasta('line' => 60, 'lucy');
			}
			else
			{
				($fasta,$quality) = $seqRef->prepareFasta('line' => 60);
			}

			print FASTAOUT $fasta;

			if ($qualityFile && $quality)
				{ print QUALOUT $quality; }
		}
	}
		
	close FASTAOUT;
	close QUALOUT if ($qualityFile);	

	return 1;
}

1;
