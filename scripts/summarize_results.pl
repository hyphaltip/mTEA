#! /usr/bin/perl
# Thu 03 Nov 2011 04:55:25 PM PDT summarize the results of a run
use strict;
use Getopt::Long;
require File::Temp;
use File::Temp ();
use File::Temp qw/ :seekable /;

my $RESULTSDIRECTORY; #directory with output data
#my $SIZE_SIMILARITY = 0.9; # threshold below which the reference peptide and discovered peptide are too different in size to be the same
#my $MATCH_THRESHOLD = 0.9; # how many matches the two reference and discovered peptide must be to considered the same
my $PEPTIDEOUT = "pep.fas"; # name of the file with peptide sequences
my $NUCOUT = "nuc.fas"; #name of the file with nucleotide sequences

GetOptions(
	'i:s'     => \$RESULTSDIRECTORY,
#	's:s'	  => \$SIZE_SIMILARITY,
#	'm:s'     => \$MATCH_THRESHOLD,
	'p:s'	  => \$PEPTIDEOUT,
	'n:s'	  => \$NUCOUT
 );

if ($RESULTSDIRECTORY eq "") {
	die "usage: perl summarize.pl -i <directory with data ouput> -r <OPTIONAL: fasta file with known TE peptides and standard names>\nsee script for more options";
}

#### part1 read the result files #########
my %hits; #holds the location in the genome as key and [0] sequence, [1] peptide, [2] sequence name;
my %pepof; #holds the location of all the elements as key and an array with all the peptide names as value
my %allpepseq; #holds the name of the peptide sequence as key and the sequence as value

opendir(DIR, $RESULTSDIRECTORY);
my @filenames = readdir(DIR);
foreach my $filename (@filenames) {
	if ($filename =~ /(\S+)-domainBLAST$/) {
		my $resultfilename = "$RESULTSDIRECTORY/" . $filename . "/" . $1 . ".RESTIRIR";
		open (INPUT, $resultfilename) or die "Cannot open: $resultfilename";
		while (my $line = <INPUT>) {
			if ($line =~ />(\S+)-(\S+?:\d+\.\.\d+)_/) {
				my $peptidename = $1;
				my $pepfile = "$RESULTSDIRECTORY/" . $filename . "/" . $peptidename . ".pep"; #name of assciated peptide file				
				my $location = $2;

				push @{ $pepof{$location} }, $pepfile; #record all the peptide file names at this location

				$line = <INPUT>;
				chomp $line;
				$hits{$location}[0] = $line;
				open (INPUT2, $pepfile) or die "cannot open peptide file $pepfile\n";
				<INPUT2>;
				my $pepseq; #sequence of the peptide
				while (my $line2 = <INPUT2>) {
					chomp $line2;
					$pepseq .= $line2;
				}
				$allpepseq{$pepfile} = $pepseq; 
			}
			else {
				die ("error reading file $resultfilename at line\n$line");
			}
		}
		close INPUT;
	}
}


#### part2 print the results ########
open (OUTPUT1, ">$RESULTSDIRECTORY/$PEPTIDEOUT") or die;
open (OUTPUT2, ">$RESULTSDIRECTORY/$NUCOUT") or die;
foreach my $key (keys %hits) {	
	print OUTPUT2 ">$key\n", "$hits{$key}[0]\n";
}

for my $loc ( keys %pepof ) { 
		print OUTPUT1 "$loc: \n";
		my %uniquepepseq; #holds the sequence of the peptide as key and file location as value, used to remove duplicates
		for my $i ( 0 .. $#{ $pepof{$loc} } ) {
			$uniquepepseq{$allpepseq{$pepof{$loc}[$i]}} = "file $pepof{$loc}[$i]";
		}
		foreach my $key (keys %uniquepepseq) {
			print OUTPUT1 "$uniquepepseq{$key}\n";
			print OUTPUT1 "$key\n";
		}
		print OUTPUT1 "\n";
	}
close OUTPUT1;
close OUTPUT2;
print STDERR "Results written to files $PEPTIDEOUT and $NUCOUT\n";
