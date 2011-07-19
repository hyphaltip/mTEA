#!/usr/bin/perl -w
# Author: Peter Arensburger 
# takes a fasta file and orf bounds as input, returns possible TIR positions with TSDs within the bounds specified by the constants

use strict;
use Bio::SeqIO;

require File::Temp;
use File::Temp ();
use File::Temp qw/ :seekable /;

#### Constants, adjust for each run #############
my $TSD_CONSENSUS = "NNNNNNNN"; #"N" for any base, case sensitive except for N
my $TSD_SUBSTITIONS = 2; #maximum number of bp differences allowed between TSDs, not counting insertions symbol
my $TSD_INSERTIONS = 1; #maximum number of insertions allowed in TSDs, at the moment only 0 or 1 are allowed
my $TSD_DELETIONS = 1; #maximum number of deletions allowed in TSDs, at the moment only 0 or 1 are allowed
my $TSD_MAX_INDELS = 1; #maximum number of allowed indels 
my $MIN_TIR_LENGTH = 10; #shortest allowed TIR, the blast parameters influence the TIR length

my $SEARCH_PAST_TIR_END = 2; #number of bp to search for TSDs around the TIR end, changing this value from 2 dramatically quickly increase the amount of computation

GetOptions(
    'c|consensus:s'      => \$TSD_CONSENSUS,
    's|substitution:i'   => \$TSD_SUBSTITUTIONS,
    'i|insertion:i'      => \$TSD_INSERTIONS,
    'd|deletion:s'       => \$TSD_DELETIONS,
    'max|maxindel:i'     => \$TSD_MAX_INDELS,
    'l|tirlength:i'      => \$MIN_TIR_LENGTH,
    );

#tesing constant inputs
if ($TSD_INSERTIONS > 1) { 
	warn "Only 0 or 1 TSD insertion are allowed, defaulting to 1";
	$TSD_INSERTIONS = 1;
}
if ($TSD_DELETIONS > 1) { 
	warn "Only 0 or 1 TSD deletions are allowed, defaulting to 1";
	$TSD_DELETIONS = 1;
}


#Misc variables used throught the program
my $tsd_regex = cons2regx($TSD_CONSENSUS); #regular expression pattern of the TSD_CONSENSUS;
my $tsd_length = length $TSD_CONSENSUS; #not all instances this are replaced in the script so far

#Print output file header
print join("\t",qw(fasta_title TIR1_start TIR1_end TIR2_start TIR2_end TIR1_sequence TIR2_sequence
                   a_possible_TSD1 a_possible_TSD2)),"\n";

# read input fasta files, assuming ORF boundaries are provided in the
# fasta title, if not sequence is split into to equal halves
my %orfbounds; #fasta title as key and array of orf boundaries as value
my($fastafilename) = @ARGV;
my $infasta  = Bio::SeqIO->new(-file => $fastafilename ,
				  -format => 'fasta');
while (my $seq = $infasta->next_seq) {
	my $fastatitle = $seq->display_id; #used for warnings
	my $b1;  # left bound of ORF
	my $b2;  # right bound of ORF
	my %tsd; # holds an index number as key, as value has an array
		 # with the following informating for each possible
		 # tsd pair: [0] sequence on the left side of the
		 # sequence, [1] sequence on the right, [2] start
		 # position on left, [3] start postion on right, [4]
		 # acceptable pair or not, boolean, 0 for not
		 # acceptable

	##### Test fasta input file ##############
	my $seqlength = $seq->length;
	#check the input to make sure a sequence of some length is present
	if ($seqlength < 2) {
		warn("WARNING: sequence length of fasta sequence $fastatitle is too short, ignoring\n");
		next;
	}
	# check input sequences for other characters than the standard uppercase A, C, G, T, N .  
	# can't we just uppercase the sequence then?
	my $seqstr = $seq->seq;
	my $count_uc_standard_char = ($seqstr =~ tr/ACGTN//); #number of uppercase standard characters
	my $count_lc_standard_char = ($seqstr =~ tr/acgtn//); #number of lowercase standard characters
	unless (($count_uc_standard_char == $seq->length) || ($count_lc_standard_char == $seq->length) ) {
	    warn("WARNING: fasta sequence $fastatitle appears to contain a mixture of upper and lower case letters and/or characters other than A,C,G,T,N.  These will all be treated as separate characters\n");
	}

	#check for ORF bounds, if none split the sequence down the middle
	if($seq->display_id =~ /_(\d+)_(\d+)$/) {
		$b1 = $1;
		$b2 = $2;
	} else {
	    warn("WARNING: no ORF bounds provided for input fasta sequence $fastatitle, looking for TIRs on either side of the middle of the sequence\n");
	    $b1 = int($seq->length)/2;
	    $b2 = (int($seq->length)/2) + 1;
	}

	#####Identify potential TIR sites by blasting one side to the reverse complent of the other####
	
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

	# at the moment I don't know a good way to capture errors from blastn, 
	# seems like redirecting the output to 2>&1 is not working

	# This is BLAST+ I believe
	`blastn -query $temp_seq1_filename -subject $temp_seq2_filename -word_size 4 -out $temp_blast_filename -gapopen 0 -gapextend 4 -reward 2 -penalty -5 -outfmt 6 -dust no`;


	# go through the blast table results
	open (my $fh => $temp_blast_filename) or die $!;
	while (<$fh>) {
	    my($q, $s, $p_identities, $identities, $unk, $gaps, 
	       $q_start, $q_end, $s_start, $s_end, $e, $score) = split;
	    next if ($s_start < $s_end); #only look at blast plus/minus hits for tirs
	    next if ((abs($s_end - $s_start) < $MIN_TIR_LENGTH) ||
		     (abs($q_end - $q_start) < $MIN_TIR_LENGTH)); #skip if TIR is too short

	    #get the tsd sequences to examine
	    my @leftsds;   #holds the sequences of potential left TSDs
	    my @rightsds;  #hold the sequences of potential right TSDs
	    my @leftp_tsds;	 #holds the start site of left TSDs
	    my @rightp_tsds;   #holds the start site of the right TSDs

	    if ( $SEARCH_PAST_TIR_END > 0) {
		for (my $i=(- $SEARCH_PAST_TIR_END); $i <  $SEARCH_PAST_TIR_END; $i++) {
		    push @leftp_tsds, ($q_start - $i - $tsd_length - 1);
		    push @rightp_tsds, ($s_start + (length $seq1) + $i);
		}
	    } else {
		push @leftp_tsds, ($q_start - $tsd_length - 1);
		push @rightp_tsds, ($s_start + (length $seq1));
	    }

	    #evaluate every possible TSD pair
	    foreach my $ltsd (@leftp_tsds) {
		foreach my $rtsd (@rightp_tsds) {
		    my($test_tsd_outcome,$ltsd_used, 
		       $rtsd_used) = testsds($ltsd, $rtsd, $seq->seq()); #test for possible identidy between tsds
		    if($test_tsd_outcome) { 
			my $tir1_pos_start = $ltsd+(length $TSD_CONSENSUS)+2;
			my $tir1_pos_end = $q_end;
			my $tir2_pos_start = $s_end + (length $seq1) + 1;
			my $tir2_pos_end = $rtsd;		
			my $tir1_seq = substr($seq->seq(),$ltsd+(length $TSD_CONSENSUS)+1, 
					      $q_end - $ltsd - (length $TSD_CONSENSUS) - 1);
			my $tir2_seq = substr($seq->seq(),($s_end + (length $seq1)), 
					      $rtsd - ($s_end + (length $seq1)));

			print join("\t",$fastatitle,$tir1_pos_start,$tir1_pos_end,$tir2_pos_start,
				   $tir2_pos_end,$tir1_seq,$tir2_seq,$ltsd_used,$rtsd_used),"\n";
		    }
		}		
	    }		
	}
	
	#cleanup 
	close $fh;
# unnecessary if you use the "UNLINK=>1" above
#	`rm $temp_seq1_filename`;
#	`rm $temp_seq2_filename`;
#	`rm $temp_blast_filename`;
}

#compare the tsds given the similarity paramters provided as constants
sub testsds {
	my($tsd1p, $tsd2p, $seq) = @_; #position of left and right tsds as well as sequence and variation parameters

	#test that both tsds match the consensus
	unless (substr($seq,$tsd1p,length $TSD_CONSENSUS) =~ /$tsd_regex/) {return (0)};
	unless (substr($seq,$tsd2p,length $TSD_CONSENSUS) =~ /$tsd_regex/) {return (0)};

	#setup hashes that will hold original tsd sequences and all permutations with indels as key, value is the number of changes (i.e. indels made)
	my %tsd1; #tsds on the left side of the ORF
	my %tsd2; #right side

	#generate TSDs with inserts
	unless ($TSD_INSERTIONS < 1) {
		for (my $i=0; $i<length $TSD_CONSENSUS; $i++) {
		    $tsd1{substr($tsd1_seq, 1, $i) . "-" . substr($tsd1_seq, $i + 1, 
								  (length $TSD_CONSENSUS) - $i - 1)} = 1;
		    $tsd2{substr($tsd2_seq, 1, $i) . "-" . substr($tsd2_seq, $i + 1, 
								  (length $TSD_CONSENSUS) - $i - 1)} = 1;
		}
	}

	#generate TSDs with deletions
	unless (($TSD_DELETIONS < 1) || ($tsd1p < $TSD_DELETIONS) || 
		($tsd2p > (length ($seq) - (length $TSD_CONSENSUS) - $TSD_DELETIONS)))  { 
# ignore if no deletions were requested or if the TSD is too close to the ends of the sequence to make deletion

	    my $ext_tsd1_seq = substr($seq, $tsd1p - $TSD_DELETIONS, (length $TSD_CONSENSUS) + 1); #extended sequence of $tsd1_seq to account for the bps that will be deleted
	    my $ext_tsd2_seq = substr($seq, $tsd2p - $TSD_DELETIONS, (length $TSD_CONSENSUS) + 1); #extended sequence of $tsd2_seq to account for the bps that will be deleted;
	    for (my $i=1; $i<=length $TSD_CONSENSUS; $i++) {
		$tsd1{substr($ext_tsd1_seq, 0, $i) . substr($ext_tsd1_seq,$i+1,(length $TSD_CONSENSUS) - $i)} = 1;
		$tsd2{substr($ext_tsd2_seq, 0, $i) . substr($ext_tsd2_seq,$i+1,(length $TSD_CONSENSUS) - $i)} = 1;
	    }
	}

	#add the unmodified tsd
	$tsd1{$tsd1_seq} = 0;
	$tsd2{$tsd2_seq} = 0;

	#compare all the possible left tsds to the right tsds
	foreach my $seq1 (keys %tsd1) {
	    foreach my $seq2 (keys %tsd2) {
		next if ($tsd1{$seq1} + $tsd2{$seq2}) > $TSD_MAX_INDELS; # don't continue if the two TSDs have more indels combined than allowed
		my $bpchanges = 0; #counter of the number differences between sequences
		for (my $i=0; $i<length $TSD_CONSENSUS; $i++) {
		    my $bp1 = substr($seq1, $i, 1);
		    my $bp2 = substr($seq2, $i, 1);
		    unless (($bp1 eq "-") || ($bp2 eq "-")) { #ignore this position if it's an insert, that has already been accounted for earlier
			unless ($bp1 eq $bp2) {
			    $bpchanges++;
			}
		    }
		    last if ($bpchanges > $TSD_SUBSTITIONS); #leave the loop early if past number of allowed changes
		}
		if ($bpchanges <= $TSD_SUBSTITIONS) { #found a match, woot!
		    return (1, $seq1, $seq2);
		} 
	    }
	}
	return(0); #if got to this point then no matches were found
}

#converts a consensus pattern to a regular expression, could be improved by allow IUPAC uncertainty codes
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
