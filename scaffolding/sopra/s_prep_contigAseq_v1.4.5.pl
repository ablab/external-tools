#!/usr/bin/perl

#	AUTHOR
#  Adel Dayarian

#	The code is available freely, under the GNU Public License, at http://www.physics.rutgers.edu/~anirvans/SOPRA/

use strict;
use Storable;
use Getopt::Long;

my @contig_file;
my @mate_names;

my @f3_names;
my @r3_names;

my $address="";

my $nonseq=0;

&GetOptions( "contig=s{,}" => \@contig_file,
             "mate=s{,}",  => \@mate_names,
             "f3=s{,}",  => \@f3_names,
     	     "r3=s{,}",  => \@r3_names,
     	     "nonseq=i",  => \$nonseq,
     	     "a=s", => \$address,
           ) or die "Unknown option\n";
#-------------------------------------------------   
if(  (! @contig_file)  || ( (! @mate_names) && (! @f3_names) && ( $nonseq == 0 ) ) || ( $address eq "" ) ){
   print "Usage: $0\n\n";
   print "-contig  File(s) in fasta format containing the contigs\n\n";

   print "-mate  File(s) in fasta format containing reads for a paired-end library\n\n";

   print "-f3  File in fasta format containing F3 reads for a mate pair library\n";
   print "-r3  File in fasta format containing R3 reads for the corresponding mate pair library\n\n";

   print "-nonseq  Set equal to 1 if you are only inputting marker data without any short read library\n\n";

   print "-a  Name of the output directory\n\n";   
   
   #print "example: $0  -mate mate1.fasta -d L1 -mate mate2.fasta -d L2 -frag fragment1.fasta fragment2.fasta -a mydir_sout\n\n";
   exit;
}
if ( scalar(@f3_names)!=scalar(@r3_names) ){
   die "Number of f3 files is not equal to the number of r3 files\n";
}

if( $address =~ /^(.*)\/$/ ){
     $address=$1; 
}
if ( ($address ne "." ) && ($address ne "" ) ){     
	 mkdir("$address", 0777) || print $!;
}
#mkdir("$address"."/mapping", 0777) || print $!; 

print "\nContig files:\n@contig_file\n\n";
if ( $nonseq != 0 ){
	print "Only non short read data\n";
}else{
	print "Paired-end files:\n@mate_names\n\n";
	print "f3 files:\n@f3_names\n\n";
	print "r3 files:\n@r3_names\n\n";
}
print "Output directory:\n$address\n\n";

#------------------------------------------- contigs
my $tig_info;
open (OUT,">$address/contigs_sopra.fasta"); 

my ($tig,$ori_header,$contig); my $start=0;

my $num_contig_file=scalar(@contig_file)-1;
for (my $i=0; $i <= $num_contig_file ; $i++){
    open (IN,"<$contig_file[$i]") or die "Cannot open $contig_file[$i] $!\n";
	while (my $line= <IN>){
	   chomp($line); 		        
	   if ($line=~ /^\>(.*)/){
			if($start==1){           
				  if ( length($contig) > 0 ) {
					  $tig++; 				
					  print OUT ">contig$tig"."|$ori_header\n";
					  print OUT "$contig\n";				
					  $tig_info->{$tig}{'length'}=length($contig);
					  $tig_info->{$tig}{'cov'} = 0;
				  }	  
				  $contig='';
			}else{
				  $start=1;
			}
	   		$ori_header=$1;
	   }else{
		    $contig.=$line;
	   }	
	}
	$tig++; 				
	print OUT ">contig$tig"."|$ori_header\n";
	print OUT "$contig\n";				
	$tig_info->{$tig}{'length'}=length($contig);
	$contig='';

	close IN;	
}
close OUT;
mkdir("$address"."/div", 0777) || print $!; print "\n";

$tig_info->{1}{'seq_l_u'}=0;
$tig_info->{1}{'seq_l_d'}=0;

store \%{$tig_info}, "$address"."/div/tig_info"; 

if ( $nonseq != 0 ){
	print "Finished -- no mater pair/paired-end library\n"; exit;
}

#------------------------------------------- paired-end
my $num_pair; my $num_t;  
my $num_matefile=scalar(@mate_names)-1;
if ( $num_matefile > -1){	
	for (my $i=0; $i <= $num_matefile ; $i++){
		my $num_seq=0;
		open (MATE_FILE,"<$mate_names[$i]") or die "Cannot open $mate_names[$i] $!\n";  

		if( $mate_names[$i] =~ /^(.*)\/([^\/]*)\.[Ff][Aa][Ss][Tt][Aa]$/ ){ 
			open (OUT,">$address/$2"."_sopra.fasta") or die "Cannot open $address/$2"."_sopra.fasta $!\n";  
		}elsif( $mate_names[$i] =~ /^(.*)\/([^\/]*)\.[Ff][Aa]$/ ){
			open (OUT,">$address/$2"."_sopra.fasta") or die "Cannot open $address/$2"."_sopra.fasta $!\n";  			
		}elsif( $mate_names[$i] =~ /^([^\/]*)\.[Ff][Aa][Ss][Tt][Aa]$/ ){
			open (OUT,">$address/$1"."_sopra.fasta") or die "Cannot open $address/$1"."_sopra.fasta $!\n";   
		}elsif( $mate_names[$i] =~ /^([^\/]*)\.[Ff][Aa]$/ ){
			open (OUT,">$address/$1"."_sopra.fasta") or die "Cannot open $address/$1"."_sopra.fasta $!\n";  
		}else{
			die "Incorrect name format $mate_names[$i]\n";  
		}	
				
		while (my $line=<MATE_FILE>){
			if( $line =~ /^\#/ ){
			}elsif( $line =~ /^\>/ ){            
				$num_pair++; $num_seq++; $num_t++; 
				print OUT ">s$num_t"."a\n";            
				$line=<MATE_FILE>; print OUT "$line";
				$line=<MATE_FILE>; 
				if( $line =~ /^\>/ ){
					$line=<MATE_FILE>; 				
					print OUT ">s$num_t"."b\n";			
					print OUT "$line";
				}else{
					die "$mate_names[$i] is not in the right format\n";
				}
			}else{
				die "$mate_names[$i] is not in the right format:\n$line\n";
			}
		}   
		close MATE_FILE; close OUT;
		print "\n$mate_names[$i]        number of pairs:  $num_seq\n";
	}
	print "\nTotal number of pairs (from paired-end libraries):  $num_pair\n\n";
}
#------------------------------------------- mate pair
$num_pair=0;

$num_matefile=scalar(@f3_names)-1;

if ( $num_matefile > -1 ){	
	for (my $i=0; $i <= $num_matefile ; $i++){
		my $num_seq=0;
		
		open (F3_FILE,"<$f3_names[$i]") or die "Cannot open $f3_names[$i] $!\n";
		open (R3_FILE,"<$r3_names[$i]") or die "Cannot open $r3_names[$i] $!\n";
		
		if( $f3_names[$i] =~ /^(.*)\/([^\/]*)\.[Ff][Aa][Ss][Tt][Aa]$/ ){
			open (OUT,">$address/$2"."_and_r3_sopra.fasta") or die "Cannot open $address/$2"."_and_r3_sopra.fasta $!\n";  
		}elsif( $f3_names[$i] =~ /^(.*)\/([^\/]*)\.[Ff][Aa]$/ ){
			open (OUT,">$address/$2"."_and_r3_sopra.fasta") or die "Cannot open $address/$2"."_and_r3_sopra.fasta $!\n";  
		}elsif( $f3_names[$i] =~ /^(.*)\/([^\/]*)\.[Cc][Ss][Ff][Aa][Ss][Tt][Aa]$/ ){ 
			open (OUT,">$address/$2"."_and_r3_sopra.fasta") or die "Cannot open $address/$2"."_and_r3_sopra.fasta $!\n";  
		}elsif( $f3_names[$i] =~ /^([^\/]*)\.[Ff][Aa][Ss][Tt][Aa]$/ ){
			open (OUT,">$address/$1"."_and_r3_sopra.fasta") or die "Cannot open $address/$1"."_and_r3_sopra.fasta $!\n";  
		}elsif( $f3_names[$i] =~ /^([^\/]*)\.[Ff][Aa]$/ ){
			open (OUT,">$address/$1"."_and_r3_sopra.fasta") or die "Cannot open $address/$1"."_and_r3_sopra.fasta $!\n";  
		}elsif( $f3_names[$i] =~ /^([^\/]*)\.[Cc][Ss][Ff][Aa][Ss][Tt][Aa]$/ ){
			open (OUT,">$address/$1"."_and_r3_sopra.fasta") or die "Cannot open $address/$1"."_and_r3_sopra.fasta $!\n";  	
		}else{
			die "Incorrect name format $f3_names[$i]\n";  
		}	
		
		while (my $line_f3=<F3_FILE>){
			if( $line_f3 =~ /^\#/ ){
				my $line_r3=<R3_FILE>;             
				if( $line_r3 =~ /^\#/ ){             
				}else{
					die "$r3_names[$i] is not in the right format\n";
				}                    
			}elsif( $line_f3 =~ /^\>/ ){
				$num_seq++; $num_pair++; $num_t++;
				$line_f3=<F3_FILE>;             
				my $line_r3=<R3_FILE>; 
	
				if( $line_r3 =~ /^\>/ ){
					$line_r3=<R3_FILE>; 	
					print OUT ">s$num_t"."a\n";         
					print OUT "$line_r3";
					print OUT ">s$num_t"."b\n";         
					print OUT "$line_f3";				
				}else{
					die "$r3_names[$i] is not in the right format\n";
				}
			}else{
				die "$f3_names[$i] is not in the right format\n";
			}
		}   
		close F3_FILE; close R3_FILE; close OUT;
		print "\n$f3_names[$i] and $r3_names[$i]        number of pairs: $num_seq\n";
	}
	print "\nTotal number of pairs (from mate pair libraries):  $num_pair\n\n";
}	

#------------------------------------------- 
sub reverseComplement{
   $_ = shift;
   $_ = uc();
   tr/ATGC/TACG/;
   return (reverse());
}
#--------------------------------------------
