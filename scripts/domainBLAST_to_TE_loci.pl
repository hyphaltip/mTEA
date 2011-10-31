#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Bio::SeqIO;
use List::Util qw(min max);
use Getopt::Long;
use IO::String;
my $flank = 5000;
my $minmatch = 50; # percent 
my $name = 'TE';
my $dir  = "TE_out";
my $debug = 0;
my $Min_evalue = 1e-5;
# NCBIblast
# -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos"

#TBLASTN
# -links -hspspesmax=300 

my $blast_flavor = 'WUBLAST'; # could be NCBI NCBI+ or WU-BLAST
# assumes 


GetOptions(
	   'b|blast:s' => \$blast_flavor,
	   'f|flank:s'   => \$flank,
	   'n|name:s'    => \$name,
	   'o|d|dir:s'   => \$dir,
	   'v|verbose!'  => \$debug,
	   );

# output dir
mkdir($dir) unless -d $dir;

# order of argument is
# Chromosome FASTA database file (genome)
# Proteins searched with (pepdb)
# Table of BLAST results in -mformat 3 from WU-BLAST
my ($genome,$pepdb,$table) = @ARGV;

my $pepdbh = Bio::DB::Fasta->new($pepdb);
my $dbh = Bio::DB::Fasta->new($genome);

open(my $fh => $table) || die $!;
my (%groups,%segments);
my $ct = 0;
# parse BLAST (WU-BLAST mformat=3 and -links), key part is the -links option
while(<$fh>) {
    next if /^\#/;
    if( $blast_flavor =~ /WU/i ) {
	my ($qid,$sid, $evalue, $N, $Sprime,$S,$alnlen,$nident,$npos,$nmism,
	    $pcident,$pcpos,$qgaps,$qgaplen,$sgaps,$sgaplen,$qframe,$qstart,$qend,
	    $sframe,$sstart,$send, $links) = split;

	# next if $pcpos < $minmatch;
	my $num = 1;
	if ($links =~ /\((\d+)\)/) {
	    ($num) = $1;
	}
	$links =~ s/[\(\)]//g;
	$groups{$sid}->{$qid}->{$links} = $evalue;
	$segments{$sid}->{$qid}->[$num] = [$sstart => $send,$pcpos];    
    } elsif( $blast_flavor =~ /NCBI/i ) {
	my ($qid,$sid, $pcident, $alnlen,$nmism,$gapopnes,$qstart,$qend,
	    $sstart,$send, $evalue,$bitscore,$pcpos) = split;
	#next if $pcpos < $minmatch;	
	push @{$groups{$sid}->{$qid}}, [$sstart,$send,$qstart,$qend,$evalue,$pcpos];
    } else {
	die("unknown BLAST flavor $blast_flavor\n");
    }
    $ct++;
}


$ct = 0;
# walk through all the chromosomes with hits
for my $chrom ( keys %groups ) {

    for my $q ( keys %{$groups{$chrom}} ) { # walk through all the proteins with hits
	for my $grp ( keys %{$groups{$chrom}->{$q}} ) {	# walk through each of these logical groups (this is gotten by the -links option in WU-BLAST)

	    my ($min,$max,$evalue);
	    if( $blast_flavor =~ /WU/i ) {
		$evalue = $groups{$chrom}->{$q}->{$grp};
		next if $groups{$chrom}->{$q}->{$grp} > $Min_evalue;
		my ($low_pcident);
		for my $num ( split(/-/,$grp) ) {
		    my $hsp = $segments{$chrom}->{$q}->[$num];		
		    if( ! defined $hsp->[0] ) {
			warn("no value for $grp:$chrom:$q\n");
		    }
		    $min = ! defined $min ? $hsp->[0] : min($hsp->[0],$min);
		    $max = ! defined $max ? $hsp->[1] : max($hsp->[1],$max);
		    if( $hsp->[2] < $minmatch ) {
			$low_pcident = 1;
#			warn("low percent identity\n"); # do something else with this info
		    }
		}
		if( ! defined $min ) {
		    warn("no value for min for grp:$grp in $chrom-$q\n");
		    next;
		}
	    } else {
		warn("not currently supported\n");
		
	    }
	    my ($left,$right) = ($min - $flank, 
				 $max   + $flank);
	    my $chromlen = $dbh->length($chrom);
	    $left = 1 if $left < 1;
	    $right = $chromlen if $right > $chromlen;
	    my $segment = $dbh->seq($chrom,$left => $right);	    
	    my $oname = "$dir/$q\_$name\_$ct";
	    my $out = Bio::SeqIO->new(-format => 'fasta',
				      -file   => ">$oname.fa");
	    $out->write_seq(Bio::Seq->new(-id => "$name\_$ct",
					  -desc => sprintf("%s:%d..%d target=%d..%d domain=%s",$chrom,$left,$right,$min,$max,$q),
					  -seq => $segment));
	    $out->close();
	    my $qid = $q;
	    $qid =~ s/\/(\d+)\-(\d+)$//;
	    if( ! -f "$dir/$qid" && ! -z "$dir/$qid" ) {
		Bio::SeqIO->new(-format => 'fasta',
				-file => ">$dir/$qid")->write_seq
				    (Bio::Seq->new(-id => $qid,
						   -seq => $pepdbh->seq($q)));
	    }
	    # Do a realignment of the protein to genome to handle this better than the BLAST alignment. This is supposed to basically be replacement of the PHI heuristic with 
	    # a full dynamic programming alignment

	    my $cds = &run_exonerate("$dir/$qid",$oname);

	    if( ! $cds ) {		
		$cds = &run_exonerate("$dir/$qid",$oname,1);
	    }
	    if( ! $cds ) {
		warn "no CDS for $left-$right $evalue\n";	    
	    } else {
		my $cds_seq = Bio::SeqIO->new(-fh => IO::String->new($cds),
					      -format => 'fasta')->next_seq;
		$out = Bio::SeqIO->new(-format =>'fasta',-file => ">$oname.pep");
		$out->write_seq($cds_seq->translate);		
	    }
	    $ct++;
	    last if $debug;
	}
	last if $debug;
    }
    last if $debug;
}

sub run_exonerate {
    my ($query,$outname, $sensitive) = @_;
    my $region = '--refine region';
    if( $sensitive ) {
	$region = '--exhaustive';
    }
    my $exe = "exonerate -m protein2genome -q $query -t $outname.fa $region --ryo \"##startCDS\n>%ti %tab..%tae (%tcl nt)\n%tcs\n##endCDS\n\"";
    
    open(my $ofh => "$exe |") || die $!;
    my $cds;
    my $keep = 0;
    my $cache;
    while(<$ofh> ) {
	$cache .= $_;
	if(/^\#\#startCDS/ ) {
	    $keep = 1;
	} elsif( /^\#\#endCDS/ ) {
	    $keep = 0;
	} elsif( $keep ) {		       
	    $cds .= $_;
	}	
    }
    $cds;
}
