#!/usr/bin/perl -w
use strict;

use Bio::SeqIO;
use Bio::DB::Fasta;

my $blast_flavor; # could be NCBI NCBI+ or WU-BLAST
# assumes 
