#!/usr/bin/perl

$bowtiePath = "";
$bowtiePath = "/mnt/software/mappers/bowtie";

# preprocess of reads and map them using bowtie

if( @ARGV != 4 || $ARGV[ 0 ] eq "-h" )
{
	print "Usage:\n";
	print "\tperl preprocess.pl  <cf> <rf1> <rf2> <output> \n";
	print "\t\t<cf>: Multi-fasta contig file\n";
	print "\t\t<rf1>: Fasta/Fastq file for first reads\n";
	print "\t\t<rf2>: Fasta/Fastq file for second reads\n";
	print "\t\t<output>: Name for output file\n";
	exit( 0 );
}

( $contigFile, $readFile1, $readFile2, $outputFile ) = @ARGV;

if( $bowtiePath !~ "/\$" )
{
	$bowtiePath .= "/";
}

# create output folder
@dir = split( "/", $outputFile );
$folder = "";
for( $i = 0; $i < @dir - 1; $i++ )
{
	$folder .= "$dir[ $i ]/";
}
$outputFile = $dir[ @dir - 1 ];
system("mkdir -p $folder");

print "PREPROCESS:\n";

open( READ, "$readFile1" ) or die $!;
while( $line = <READ> )
{
	if( $line =~ "^>" )
	{ $fasta = "-f"; print "Fasta format is recognized\n"; last; }
	elsif( $line =~ "^@" )
	{ $fasta = "-q"; print "Fastq format is recognized\n"; last; }
	else
	{ print "reads format is not correct\n"; exit( 0 ); }
}

$type = "fs";
$line = <READ>;
if( $line =~ m/(0|1|2|3)/ )
{
    $type = "cs";
}

close READ;


# combine the reads together and change the name
$time = localtime;
print "[$time]\t";
print "Combining read files for bowtie mapping...\n";
$readFile = "${folder}combined.fasta";
open( F, "$readFile1" ) or die $!;
open( S, "$readFile2" ) or die $!;
open( OUTPUT, ">$readFile" ) or die $!;
$num = 1;
while( $first = <F> )
{
    if( $fasta eq "-f" )
    {
	print OUTPUT ">$num".".1\n";
	$first = <F>;
	print OUTPUT $first;
		
	$second = <S>;
	print OUTPUT ">$num".".2\n";
	$second = <S>;
	print OUTPUT $second;
    }
    else
    { # fastaq
	print OUTPUT "\@$num".".1\n";
	$first = <F>;
	print OUTPUT $first;
	$first = <F>;
	print OUTPUT "+$num".".1\n";
	$first = <F>;
	print OUTPUT $first;

	$second = <S>;
	print OUTPUT "\@$num".".2\n";
	$second = <S>;
	print OUTPUT $second;
	$second = <S>;
	print OUTPUT "+$num".".2\n";
	$second = <S>;
	print OUTPUT $second;
    }
    $num++;
}
close S;
close F;
close OUTPUT;


# build the index 
$time = localtime;
print "[$time]\t";
print "Building bowtie index...\n";
@contigName = split( "/", $contigFile );
if( ! -f "$folder$contigName[ -1 ].1.ebwt" )
{
    if( $type eq "cs" )
    {
	$command = "${bowtiePath}bowtie-build -C $contigFile $folder$contigName[ -1 ]";
	`$command` or die "ERROR! Running bowtie-build error.\n";
    }
    else
    {
	$command = "${bowtiePath}bowtie-build $contigFile $folder$contigName[ -1 ]";
	`$command` or die "ERROR! Running bowtie-build error.\n";
    }
}
else
{ # no need to re-build
    print "\t\t\t\tIndex already exist. Skipping building index...\n";
}

# map with bowtie
$time = localtime;
print "[$time]\t";
print "Mapping reads using bowtie...\n";
if( $type eq "cs" )
{
    $command = "${bowtiePath}bowtie -v 3 -a -m 1 -t -C $fasta -p 5 $folder$contigName[ -1 ] $readFile 2>${folder}bowtie.err | sort -n > $folder$outputFile";
    `$command`;
}
else
{
    $command = "${bowtiePath}bowtie -v 3 -a -m 1 -t $fasta -p 5 $folder$contigName[ -1 ] $readFile 2>${folder}bowtie.err | sort -n > $folder$outputFile";
    `$command`;
}

# finalize
$time = localtime;
print "[$time]\t";
print "Preprocessing done!\n";
