#!/usr/bin/perl -w
# Author: Peter Arensburger 
# Summary: Takes a fasta file and optional orf bounds as input, returns possible TIR and TSD locations as well as an
# associated score that indicates how "good" the fit is.  Proceedes in three parts: 1) identify all potential TSD/TIR
# pairs for each sequence provided by the FASTA file, 2) find regions of similarity to the left of the ORF (or middle)
# all the FASTA sequences, calculate how close potential TSD/TIR regions are to the edge of similarity, 3) takes all 
# information for each TSD/TIR location and calculate a "score" indicating how likely it is to be real, output the results. 
#
# This script assumes that the following programs are installed and available by within the local directory:
# BLAST+ (command `blastn ...`)
# BLAT (command `blat ... `)
# CLUSTALW2 (command `clustalw2 ...`)
# TODO: test for the presence of these programs when first running the script

# NOTE: using "left" and "right" for sequences on the left or right of the ORF, also using "1" for left or "2" for right
# Eventually should probably standardize all that

use strict;
use Bio::SeqIO;
require File::Temp;
use File::Temp ();
use File::Temp qw/ :seekable /;
use Getopt::Long;
use List::Util qw(min max);
#use Time::Seconds;

# setup and test constants # 
# TODO: need to test that fasta input names are unique # 
my $TSD_CONSENSUS = "NNNNNNNN"; #"N" for any base, case sensitive except for N
# TODO: test that the TIR_CONSENSUS works in real data
my $TIR_CONSENSUS = "NNNNNNNNNN";
my $TIR_DIFFERENCES = 2; #number of differences allowed between the TIR
my $TSD_DIFFERENCES = 2; #number of differences allowed between the TSD
my $TSD_SEARCH_LENGTH = "low"; #if set to low then searches for TSDs right around the rev-compl. areas

# TODO: test the command line arguments #
GetOptions(
    'c|tsdconsensus:s'      => \$TSD_CONSENSUS,
    't|tirconsensus:s'   => \$TIR_CONSENSUS,
    'd|tirdifferences:i'      => \$TIR_DIFFERENCES,
    's|tsddifferences:i'       => \$TSD_DIFFERENCES,
    'l|tsdsearchlength:s'     => \$TSD_SEARCH_LENGTH,
    );

# Global variables setup

my %hypo; # hypo = 'hypothesis', this holds information for all identified possible transposon ends, holds an
	   # index value as key and as value an array with:
	   # [0] fasta title
	   # [1] position first base of TIR on left side of sequence
	   # [2] position first base of TIR on right side of sequence
	   # [3] percent identity between left and right TSD
	   # [4] percent identity between left and right TIRs
	   # [5] on left side, the closest conserved edge between fasta sequences
	   # [6] on right side, the closest conserved edge between fasta sequences
	   # [7] left TIR sequence
	   # [8] right TIR sequence
	   # [9] number of different TSDs adjacent to left TIR
	   # [10] number of different TSDs adjacent to right TIR
my $hindex = 0; # index for %hypo
my $tsd_regex = cons2regx($TSD_CONSENSUS); # regular expression pattern of the TSD_CONSENSUS;
my $tir_regex = cons2regx($TIR_CONSENSUS); # regular expression pattern of the TIR_CONSENSUS;

my $leftseqs_filename = File::Temp->new( UNLINK => 1, SUFFIX => '.fas' ); # used in PART2 for the BLAT run
my $rightseqs_filename = File::Temp->new( UNLINK => 1, SUFFIX => '.fas' ); # used in PART2 for the BLAT run
my $blatoutput_filename = File::Temp->new( UNLINK => 1, SUFFIX => '.blat' ); # used in PART2 for the BLAT run
my %startright; # holds the name of the fasta file as key and the start position of the right side as value,
	          # used in PART2 to get the position on the right side corrected to account for ORF sequences
open (BLAT1, ">$leftseqs_filename") or die; # used in PART2 for the BLAT run
open (BLAT2, ">$rightseqs_filename") or die; # used in PART2 for the BLAT run


#### PART1. Identify in each fasta file the location of all TSD+TIR combinations that satisfy the criteria specified 
#### by the constants

# read input fasta files, if ORF boundaries are provided in the fasta title read them, otherwise split into to equal
# halves
my %orfbounds; # fasta title as key and array of orf boundaries as value
my %left_tirfreq; # left tir sequences as key and array with all associated TSDs as value
my %right_tirfreq; # right tir sequences as key and array with all associated TSDs as value

my($fastafilename) = @ARGV;
my $infasta  = Bio::SeqIO->new(-file => $fastafilename ,
				  -format => 'fasta');
while (my $seq = $infasta->next_seq) {
	my $fastatitle = $seq->display_id; #used for warnings
	my $b1;  # left bound of ORF
	my $b2;  # right bound of ORF

	# Test fasta input file #
	my $seqlength = $seq->length;
	# check input sequences for other characters than the standard uppercase A, C, G, T, N .  
	# can't we just uppercase the sequence then?
	my $seqstr = $seq->seq;
	my $count_uc_standard_char = ($seqstr =~ tr/ACGTN//); #number of uppercase standard characters
	my $count_lc_standard_char = ($seqstr =~ tr/acgtn//); #number of lowercase standard characters
	unless (($count_uc_standard_char == $seq->length) || ($count_lc_standard_char == $seq->length) ) {
	    warn("WARNING: fasta sequence $fastatitle appears to contain a mixture of upper and lower case letters and/or characters other than A,C,G,T,N.  These will all be treated as separate characters\n");
	}

	# check for ORF bounds, if none split the sequence down the middle
	if($seq->display_id =~ /_(\d+)_(\d+)$/) {
		$b1 = $1;
		$b2 = $2;
	} else {
#	    warn("WARNING: no ORF bounds provided for input fasta sequence $fastatitle, looking for TIRs on either side of the middle of the sequence\n");
	    $b1 = int(($seq->length)/2);
	    $b2 = int(($seq->length)/2) + 1;

	}

	$startright{$seq->display_id} = $b2; # holds the position of the first base of right side used to keep locations straight

	# Identify potential TIR sites by blasting one side to the reverse complent of the other 
	#setup temporary files for running blast
	my $temp_seq1_filename = File::Temp->new( UNLINK => 1, SUFFIX => '.fas' );
	my $temp_seq2_filename = File::Temp->new( UNLINK => 1, SUFFIX => '.fas' );
	my $temp_blast_filename = File::Temp->new( UNLINK => 1, SUFFIX => '.blast' );
	
	open (SEQ1,">$temp_seq1_filename") or die;
	open (SEQ2,">$temp_seq2_filename") or die;

	#get sequences and write them to file
	my $seq1 = $seq->subseq(1,$b1);
	my $seq2 = $seq->subseq($b2,$seq->length);
	print SEQ1 ">seq1\n$seq1\n";
	print SEQ2 ">seq2\n$seq2\n";

	print BLAT1 ">$fastatitle\n$seq1\n"; #update BLAT input file used for PART2
	print BLAT2 ">$fastatitle\n$seq2\n"; #update BLAT input file used for PART2

	# at the moment I don't know a good way to capture errors from blastn, 
	# seems like redirecting the output to 2>&1 is not working
	# This is BLAST+
	`blastn -query $temp_seq1_filename -subject $temp_seq2_filename -word_size 4 -out $temp_blast_filename -gapopen 0 -gapextend 4 -reward 2 -penalty -3 -outfmt 6 -dust no`;

	# go through the blast table results
	open (my $fh => $temp_blast_filename) or die $!;
	while (<$fh>) {

		my($q, $s, $p_identities, $identities, $unk, $gaps, 
		       $q_start, $q_end, $s_start, $s_end, $e, $score) = split;

		next if ($s_start < $s_end); # only look at blast plus/minus hits for tirs
		next if ((abs($s_end - $s_start) < ((length $TIR_CONSENSUS)+1)) ||
			     (abs($q_end - $q_start) < ((length $TIR_CONSENSUS)+1))); #skip if TIR is too short

		my $search_error = $TSD_DIFFERENCES;	# how much to search around the end to account for possible
							# substitutions and insertions

		#determine the bounds of the TSD searches and push the values into @leftp_tsds and @rightp_tsds
		my @leftp_tsds;    # holds the start site of left TSDs
	    	my @rightp_tsds;   # holds the start site of the right TSDs
	   	my $ltsdstart = ($q_start - $search_error); #left tsd start position
		if ($ltsdstart < 0) {
			$ltsdstart = 0;
		}

		my $rtsdend = $s_start + $search_error + $b2;  #right tsd end position
		if ($rtsdend > $seq->length) {
			$rtsdend = $seq->length;
		}

		my $ltsdend; #left tsd end position
		if ($TSD_SEARCH_LENGTH eq "low") {
			$ltsdend = $q_start + (length $TSD_CONSENSUS) + $search_error;
		}
		elsif ($TSD_SEARCH_LENGTH eq "very low") {
			$ltsdend = $q_start + $search_error;
		}
		else {
			$ltsdend = $q_end - (length $TSD_CONSENSUS) - (length $TSD_CONSENSUS) - $search_error; 
		}									
		if ($ltsdend < 0) {
			$ltsdend = 0;
		}

		my $rtsdstart; #right tsd start position
		if ($TSD_SEARCH_LENGTH eq "low") {
			$rtsdstart = $s_start - (length $TSD_CONSENSUS) - $search_error + $b2;
		}
		elsif ($TSD_SEARCH_LENGTH eq "very low") {
			$rtsdstart = $s_start - $search_error + $b2;
		}
		else {
			$rtsdstart = $s_end  + (length $TSD_CONSENSUS) + (length $TSD_CONSENSUS) + $search_error + $b2; 
		}
										
		if ($rtsdstart < $ltsdend) {
			$rtsdstart = $ltsdend;
		}

		#update all the positions that will need to be evaluated
		for (my $i = $ltsdstart; $i <= $ltsdend; $i++) {
			push (@leftp_tsds, $i);
		}
		for (my $i = $rtsdstart; $i <= $rtsdend; $i++) {
			push (@rightp_tsds, $i);
		}
		
	    	#evaluate all TSD pairs
	    	foreach my $ltsd (@leftp_tsds) {
		    foreach my $rtsd (@rightp_tsds) {
			my($test_tsd_outcome, $tsd_identity, $tir_identity, $left_tir_seq, $right_tir_seq)
			    = testsds($ltsd, $rtsd, $seq->seq());
			# test for
			# possible tsd+tir
			# identity between the 
			# left and right ends
			if($test_tsd_outcome) {
			    $hypo{$hindex} = [($fastatitle,$ltsd,$rtsd,$tsd_identity,$tir_identity,undef,undef, #5,6 are blank?
					       $left_tir_seq,$right_tir_seq)];
			    $hindex++;
			}

		    }		
		}		
	}	
	#cleanup 
	close $fh;
	close SEQ1;
	close SEQ2;
}

#record tir frequencies
foreach my $index (keys %hypo) {
    push @{$hypo{$index}}, scalar(@{$left_tirfreq{$hypo{$index}->[7]}});
    push @{$hypo{$index}}, scalar(@{$right_tirfreq{$hypo{$index}->[8]}});	
}

#### PART 2. Identify conserved sequences between the input files and record if possible TSD+TIR locations are near the 
#### edge, update values in %hypo

# TODO: the code below is redundant (running BLAT on left then right) could be tidied up

#run blat and interpret blat output for each line on left side
open(my $fh => "blat $leftseqs_filename $leftseqs_filename stdout |") || die "cannot run BLAT: $!"; #run blat
# skip first 5 lines
for (0..4) {
    <$fh>;
}

while (<$fh>) {
    my ($match, $miss, $rep, $N, $qc, $qb, $tc, $tb, $strand, $q_name, 
	$q_size, $q_start, $q_end, $t_name, $t_size, $t_start, $t_end) 
	= split;
    unless ($q_name eq $t_name) { #avoid looking at self matches
	foreach my $index (keys %hypo) { #run through all possible locations
	    if (($hypo{$index}[0]) eq $q_name) {
		my $distance = min( abs($q_start - $hypo{$index}[1]), abs($t_start - $hypo{$index}[1]));
		# can replace 
		$hypo{$index}->[5] = ! exists $hypo{$index}->[5] ? $distance : min($hypo{$index}->[5],$distance);
	    } 
	}
    }
}
close $fh;

#run blat and interpret blat output for each line on right side
open($fh => "blat $rightseqs_filename $rightseqs_filename stdout |") || die "Cannot run BLAT: $!"; #run blat
#print "$blatoutput_filename\n";
for (0..4) { # skip first 5 lines which are header
    <$fh>;
}
while (<$fh>) {
    my ($match, $miss, $rep, $N, $qc, $qb, $tc, $tb, $strand, $q_name, 
	$q_size, $q_start, $q_end, $t_name, $t_size, $t_start, $t_end) 
	= split;
    unless ($q_name eq $t_name) { #avoid looking at self matches
	foreach my $index (keys %hypo) { #run through all possible locations
	    my $adj_q_end = $q_end + $startright{$hypo{$index}[0]}; #adjusted because right side
	    my $adj_t_end = $t_end + $startright{$hypo{$index}[0]};
	    if (($hypo{$index}[0]) eq $q_name) {
		my $distance = min( abs($adj_q_end - $hypo{$index}[2]), abs($adj_t_end - $hypo{$index}[2]));
		if (exists $hypo{$index}[6]) {
		    if ($distance < $hypo{$index}[6]) {
			$hypo{$index}[6] = $distance;
		    }
		}
		else {
		    $hypo{$index}[6] = $distance;
		}
	    } 
	}
    }
}
close $fh;
close BLAT1;
close BLAT2;

#### STEP 3 Take all the information available for each TSD/TIR pair, calculate likeliness measure, and output the results
# factors entering calculations: 1) tsd identity, 2) tir identity, 3) distance from conserved edge on left, 4) distance from 
# conserved edge on right, 5) number of different TSDs on left, 6) number of different TSDs on right 

my %final_score; # holds the potential element index as key and the final score as value

# how important each fator is, weights should add up to 100
my %factor_weight = ("tsd_id"         => 16.7, 
		     "tir_id"         => 16.7, 
		     "distance_left"  => 16.7, 
		     "distance_right" => 16.7, 
		     "num_tsd_left"   => 16.7, 
		     "num_tsd_right"  => 16.7);

for my $index (keys %hypo) { 
    # calculate each factor out of 100
    my $tsdid = $hypo{$index}[3] * 100;	# tsd similarity
    my $tirid = $hypo{$index}[4] * 100;	# tir similarity
    my $dist_left; # left distance to conserved edge determined in STEP2
    if ($hypo{$index}[5] <= 2) {
	$dist_left = 100;
    } elsif ($hypo{$index}[4] <= 10) {
	$dist_left = 50;
    } else {
	$dist_left = 0;
    }

    my $dist_right; # right distance to conserved edge determined in STEP2
    if ($hypo{$index}[6] <= 2) { 
	$dist_right = 100;
    }
    elsif ($hypo{$index}[5] <= 10) {
	$dist_right = 50;
    }
    else {
	$dist_right = 0;
    }

    my $num_tsd_left;	  # number of different TSDs on the left left 
    if ($hypo{$index}[9] > 1) { 
	$num_tsd_left = 100;
    }
    else {
	$num_tsd_left = 0;
    }

    my $num_tsd_right;		#number of different TSDs on the right
    if ($hypo{$index}[10] > 1) { 
	$num_tsd_right = 100;
    }
    else {
	$num_tsd_right = 0;
    }

    #calculate final score
    $final_score{$index} = (($factor_weight{"tsd_id"}/100) * $tsdid) +
	(($factor_weight{"tir_id"}/100) * $tirid) +
	(($factor_weight{"distance_left"}/100) * $dist_left) +
	(($factor_weight{"distance_right"}/100) * $dist_right) +
	(($factor_weight{"num_tsd_left"}/100) * $num_tsd_left) +
	(($factor_weight{"num_tsd_right"}/100) * $num_tsd_left);
}

#output the results
# sort by score, greatest to least
foreach my $index (sort { $final_score {$b} <=> $final_score {$a}} keys %final_score ) {
    print STDOUT "$final_score{$index}\t$hypo{$index}[0]\t$hypo{$index}[1]\t$hypo{$index}[2]\n";
}



################### SUBROUTINES ###################################
# compares two sequences that each include a combinaton of a TSD and TIR, returns 0 if they don't match or various paramters if they do
sub testsds {
    my($tsd1p, $tsd2p, $seq) = @_; #position of left and right tsds as well as sequence and variation parameters
    my $tsd1_seq = substr($seq,($tsd1p - (length $TSD_CONSENSUS) - 1),length $TSD_CONSENSUS); #sequence of the TSDs
    my $tsd2_seq = substr($seq,$tsd2p,length $TSD_CONSENSUS);
    my $tir1_seq = substr($seq, ($tsd1p - 1), length $TIR_CONSENSUS);
    my $tir2_seq = rc(substr($seq, $tsd2p - (length $TIR_CONSENSUS), length $TIR_CONSENSUS));

    # test that the TSDs and TIR match the required consensus
    unless ($tsd1_seq =~ /$tsd_regex/) {return (0)};
    unless ($tsd2_seq =~ /$tsd_regex/) {return (0)};
    unless ($tir1_seq =~ /$tir_regex/) {return (0)};
    unless ($tir2_seq =~ /$tir_regex/) {return (0)};

    # generate the TIR+TSD combinations to test
    # test using clustalw
    my $tsd_identities = (align2seq($tsd1_seq, $tsd2_seq) =~ tr/*//);
    my $tir_identities = (align2seq($tir1_seq, $tir2_seq) =~ tr/*//);

    if ( ($tsd_identities >= ((length $TSD_CONSENSUS) - $TSD_DIFFERENCES)) && ($tir_identities >= ((length $TIR_CONSENSUS) - $TIR_DIFFERENCES)) ) {

	# update the tir frequency hashes
	unless (grep {$_ eq $tsd1_seq} @{$left_tirfreq{$tir1_seq}}) { # prevents adding the same TSD twice
	    push @{$left_tirfreq{$tir1_seq}}, $tsd1_seq;
	}
	unless (grep {$_ eq $tsd2_seq} @{$right_tirfreq{$tir2_seq}}) { # prevents adding the same TSD twice
	    push @{$right_tirfreq{$tir2_seq}}, $tsd2_seq;
	}

	my $tsd_proportion = $tsd_identities/(length $TSD_CONSENSUS); # measures of identity
	my $tir_proportion = $tir_identities/(length $TIR_CONSENSUS);

	return (1, $tsd_proportion, $tir_proportion, $tir1_seq, $tir2_seq);
    }
    else {
	return (0);
    }
}

# converts a consensus pattern to a regular expression, could be improved by allowing IUPAC uncertainty codes
# there is a fairly easy way to do that by grabbing Bio::Tools::IUPAC and going through the blocks for each uncertainty
sub cons2regx {
    my($consensus) = @_;
    my $pattern;
    unless ($consensus) { warn "WARNING: Consensus string provided is empty"};
    for (my $i=0; $i<length $consensus; $i++) {
	if (substr($consensus, $i, 1) =~ /N/i) {
	    $pattern .= ".";
	}
	else {
	    $pattern .= substr($consensus, $i, 1);
	}
    }
    return $pattern;
}

#aligns two input sequences and the consensus line of clustalw
sub align2seq {
    # why are we using CLUSTALW to find a GLOBAL alignment? 
    # exonerate or others would be cleaner
    my($seq1, $seq2) = @_;

    # run clustalw2
    my $temp_seq3_filename = File::Temp->new( UNLINK => 1, SUFFIX => '.fas' ); # clustalw2 input file
    my $temp_seq4_filename = File::Temp->new( UNLINK => 1, SUFFIX => '.aln' ); # clustalw2 ouptut file
    open (SEQ3,">$temp_seq3_filename") or die;
    print SEQ3 ">s1\n$seq1\n>s2\n$seq2\n";
    `clustalw2 -INFILE=$temp_seq3_filename -outfile=$temp_seq4_filename`;

    # parse the output	
    open (INPUT, $temp_seq4_filename) or die;
    for ( 0..3 ) {
	<INPUT>;
    }

    # determine at what position in the file the alignment starts
    my $startpos;
    if (<INPUT> =~ /s2(\s+)\S+/) {
	$startpos = 2 + length($1);
    }
    else {
	die "error reading clustalw ouptut at line\n";
    }
    return(substr(<INPUT>, $startpos, length $seq1));
}

# reverse complement
# also just grab from Bioperl Bio::PrimarySeqI
sub rc {
    my ($sequence) = @_;
    $sequence = reverse $sequence;
    $sequence =~ tr/ACGTRYMKSWacgtrymksw/TGCAYRKMWStgcayrkmws/;
    return ($sequence);
}
