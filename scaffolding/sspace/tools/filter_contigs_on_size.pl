#!/usr/bin/perl

#filter_contigs_on_size.pl
#Marten Boetzer
#21-05-2010
#BaseClear

#This script reads a .fasta file containing contigs, and removes all contigs
#having a size lower than a user-specified treshold.
#run as follows: perl data/scripts/tmp/marten/scripts/filter_contigs_on_size.pl <INFILE> <MIN-SIZE>
#Input; 1) .fasta file containing contig sequences
#       2) minimum desired contig size
#Output; A .fasta file will be generated in the run-folder folder and is called <INFILE>_min<MIN-SIZE>

use strict;
#use warnings;
use File::Basename;

my $input_file = $ARGV[0];
my $min_contig_size = $ARGV[1];
my $fileBaseName;
my $dirName;
my $fileExtension;

#check on empty input parameters)
if (!$input_file) {
  print "please provide valid input file, such as \"velvet-contigs.fa\"\n";
}

if (!$min_contig_size) {
  print "please provide valid minimum contig size, such as \"300\"\n";
}

($fileBaseName, $dirName, $fileExtension) = fileparse($input_file, ('\.fa') );
my $outputfile = $fileBaseName.$fileExtension.'_min'.$min_contig_size;

open (INFILE, $input_file) || die "Can't open input file. Please provide a valid input filename with .fa extension.\n";
open (OUTFILE, "> $outputfile") || die "Can't open output file.\n";

my $counter=0;
my @line;
my $seq = "";
my $name = "";
my $prevname = "";

while (<INFILE>) {
  chomp;
  if ($_ =~ /^[>]/) {
    $name = $_;
    $counter++;
    if($seq ne "" && $prevname ne "" && (length($seq) >= $min_contig_size)){
       print OUTFILE $prevname, "\n";
       print OUTFILE $seq, "\n";
    }
    $prevname = $name;

       $name = "";
       $seq = "";
  }
  else {
     $seq .= $_;
  }
}
if($seq ne "" && $prevname ne "" && length($seq) >= $min_contig_size){
       print OUTFILE $prevname, "\n";
       print OUTFILE $seq, "\n";
       $name = "";
       $seq = "";
}

close (INFILE);
close OUTFILE;
