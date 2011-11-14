#! /usr/bin/perl
# Wed 02 Nov 2011 09:48:36 AM PDT Takes a list of know transposase peptide sequences and directs the other scripts to find 
# potentially recently active TEs (Transposable Elements)

use strict;
use Getopt::Long;
use Bio::SeqIO;

###### PARAMETERS ############
# location of programs and files
my $AB_BLAST_DIR = "/home/peter/Desktop/ab-blast"; #directory of AB-BLAST (because multiple blasts are on this system)
my $GENOME_FILE = "/home3/Genomes/aedes_aegypti_46_1a/aedes_aegypti_46_1a.fas";
my $DOMAIN_BLAST_SCRIPT = "./scripts/domainBLAST_to_TE_loci.pl";
my $ID_TIR_SCRIPT = "./scripts/id_TIR_in_FASTA.pl";
my $SUMMARIZE_SCRIPT = "./scripts/summarize_results.pl";

# general parameters from user
my $INPUT_FASTA; #name of the fasta file with TE peptides
my $OUTPUT_DIR; #output directory
GetOptions(
    'i|inputfile:s'     => \$INPUT_FASTA,
    'o|outdir:s'     => \$OUTPUT_DIR,
 );

# abBLAST parameters
my $ab_E = 1e-8;
my $ab_hspsepsmax = 100;
my $ab_hspmax = 5000;

# domainBLAST parameters
my $domB_blast = "WU";
my $domB_flank = 5000;

# fitering peptides paramters (i.e. how likely is the output of domainBLAST to be a real TE transposase)
my $pep_minlen = 400; #length of peptide
my $pep_maxstop = 2; #maximum number of stop codons in the peptide sequence

# id_TIR_in_FASTA paramters
my $TSD_CONSENSUS = "NNNNNNNN";
my $TIR_CONSENSUS = "NNNNNNNNNNN";
my $TIR_DIFFERENCES = 2;
my $TSD_DIFFERENCES = 2;

#test the inputs
if (($INPUT_FASTA eq "") or ($OUTPUT_DIR eq "") ) { 
	die ("usage: perl run_te_discovery.pl -i <input fasta file with TE peptides> -o <output directory name>\nOther parameters must be set in the run_te_discovery.pl script directly");
}

###### END OF PARAMETERS ############


###### START OF RUN ###### 
mkdir($OUTPUT_DIR) unless -d $OUTPUT_DIR; #make output directory

# run ab-blast on the input file
print STDERR "running ab-blast on file $INPUT_FASTA\n";
`$AB_BLAST_DIR/tblastn -d $GENOME_FILE -i $INPUT_FASTA -links -wordmask seg+xnu -E $ab_E -hspsepsmax $ab_hspsepsmax -hspmax $ab_hspmax -mformat 3 -o $OUTPUT_DIR/$INPUT_FASTA.abTBLASTN`;

# for each input peptide sequence: 
# 1) run the domainBLAST script on it
# 2) select from the output of domainBLAST those peptides that have a good chance of being 'real' (e.g. long, no too many stop codons)
my $infasta  = Bio::SeqIO->new(-file => $INPUT_FASTA ,
				  -format => 'fasta');
while (my $seq = $infasta->next_seq) {

	# write the peptide sequence to a fasta file
	my $fastatitle = $seq->display_id;
	my $sequence = $seq->seq();
	open (OUTPUT, ">$OUTPUT_DIR/$fastatitle.fa") or die;
	print OUTPUT ">$fastatitle\n$sequence\n";
	close OUTPUT;

	# write the results of the ab-blast file
	open (INPUT, "./$OUTPUT_DIR/$INPUT_FASTA.abTBLASTN") or die ("cannot open ./$OUTPUT_DIR/$INPUT_FASTA.abTBLASTN\n");
	open (OUTPUT, ">$OUTPUT_DIR/$fastatitle.abTBLASTN") or die;
	while (my $line = <INPUT>) {
		if ($line =~ /^$fastatitle\s/) {
			print OUTPUT "$line";
		}
	}
	close INPUT;
	close OUTPUT;

	# run the domainBLAST script
	print STDERR "running domainBLAST script on sequence $fastatitle\n";
	`perl $DOMAIN_BLAST_SCRIPT -b $domB_blast -f $domB_flank -o $OUTPUT_DIR/$fastatitle-domainBLAST $GENOME_FILE $OUTPUT_DIR/$fastatitle.fa $OUTPUT_DIR/$fastatitle.abTBLASTN`;

	# open each of .pep files from domainBLAST and select those with acceptable peptides
	# format the selected files for use by the id_TIR_in_FASTA.pl script
	opendir(DIR, "$OUTPUT_DIR/$fastatitle-domainBLAST");
	my @filenames = readdir(DIR);
	open (OUTPUT, ">./$OUTPUT_DIR/$fastatitle-domainBLAST/$fastatitle.PREPFORTIRID") or die;
	foreach my $filename (@filenames) {
		if ($filename =~ /(.+)\.pep$/) {
			my $basename = $1;
			my $inpep  = Bio::SeqIO->new(-file => "$OUTPUT_DIR/$fastatitle-domainBLAST/$filename" ,
				  		     -format => 'fasta');
			my $sequence = $inpep->next_seq;
			my $pepsequence = $sequence->seq();
			my $pepdesc = $sequence->desc();
			my $numstop = ($pepsequence =~ tr/\*//); # number of stop codons
			my $pepbound1;
			my $pepbound2;
			if (((length $pepsequence) >= $pep_minlen) && ($numstop < $pep_maxstop) ) { # conditions for keeping a peptide
				(my $nucseq, my $loc1, my $loc2, my $loc3) = loadnuc("$OUTPUT_DIR/$fastatitle-domainBLAST/$basename.fa");
				if ($pepdesc =~ /^(\d+)\.\.(\d+)\s/) {
					$pepbound1 = $1;
					$pepbound2 = $2;
				}
				else {
					die ("peptide description is wrong $pepdesc\n");
				}

				if ($pepbound1 > $pepbound2) {
					$nucseq = rc($nucseq);
					print OUTPUT ">$basename", "_", "$loc1", "$loc3", "..", "$loc2", "_", "$pepbound2", "_", "$pepbound1", "\n", "$nucseq\n";
				}
				else {
					print OUTPUT ">$basename", "_", "$loc1", "$loc2", "..", "$loc3", "_", "$pepbound1", "_", "$pepbound2", "\n", "$nucseq\n";
				}
			}
		}		
	}
	close OUTPUT;

	# find TIR for each selected sequence using the id_TIR_in_FASTA.pl script
	`perl $ID_TIR_SCRIPT -i "./$OUTPUT_DIR/$fastatitle-domainBLAST/$fastatitle.PREPFORTIRID" -c $TSD_CONSENSUS -t $TIR_CONSENSUS -d $TIR_DIFFERENCES -s $TSD_DIFFERENCES -o ./$OUTPUT_DIR/$fastatitle-domainBLAST/$fastatitle.RESTIRIR`;
}

###### summarize the results
`perl $SUMMARIZE_SCRIPT -i $OUTPUT_DIR`;


###### END OF RUN ###### 


###### SUBROUTINES ###########

sub loadnuc { # read the nucleotide sequence
	(my $filename) = @_;
	my $seq;
	my $loc1;
	my $loc2;
	my $loc3;
	if (open (INPUT, $filename) ) {
		while (my $line = <INPUT>) {
			if ($line =~ />(.+)\s$/) {
				my $name = $1;
				if ($name =~ /^\S+\s(\S+:)(\d+)\.\.(\d+)\s/) {
					$loc1 = $1;
					$loc2 = $2;
					$loc3 = $3;
				}
				else {
					$loc1 = -1;
					$loc2 = -1;
					$loc3 = -1;
				}
			}
			else {
				$line =~ s/\s//g;
				$seq .= $line;
			}
		} 
		return ($seq, $loc1, $loc2, $loc3);
	}
	else {
		return (-1, -1, -1);
	}
}


sub rc { #reverse complement
    my ($sequence) = @_;
    $sequence = reverse $sequence;
    $sequence =~ tr/ACGTRYMKSWacgtrymksw/TGCAYRKMWStgcayrkmws/;
    return ($sequence);
}
