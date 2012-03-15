#!/usr/bin/perl

#	AUTHOR
#  Adel Dayarian

#	The code is available freely, under the GNU Public License, at http://www.physics.rutgers.edu/~anirvans/SOPRA/

use strict;
use Storable;
use Getopt::Long;

my @parsed_sam_file;
my @distances;

my @parsed_sam_solid;
my @distances_s;

my $opt_c;
my $opt_e=0;
my $opt_pt;
my $address="";

&GetOptions( "parsed=s{,}" => \@parsed_sam_file,
     	     "d=i{,}",  => \@distances,
             "solid_parsed=s{,}",  => \@parsed_sam_solid,
     	     "solid_d=i{,}",  => \@distances_s,
     	     "c=i", => \$opt_c,
     	     "e=i", => \$opt_e,
     	     "pt=i", => \$opt_pt,
     	     "a=s", => \$address,
           ) or die "Unknown option\n";
           
my ($copy_li,$ign_empi)=(5,0);
#-------------------------------------------------   
if(  ((! @parsed_sam_solid) && (! @parsed_sam_file)) || ( $address eq "" ) ){
   print "Usage: $0\n\n";
   print "-parsed  Parsed SAM/BAM file for a paired-end library, has to be followed by option -d\n";
   print "-d  The insert size for the corresponding paired-end library\n\n"; 

   print "-solid_parsed  Parsed version of SAM file for a mate pair SOLiD library, has to be followed by option -d_solid\n";
   print "-solid_d  The insert size for the corresponding mate pair library\n\n"; 

   print "-c  If the number of times a read and its reverse complement appear in the library is equal to or more than this value, the pairing information from that read will be disregarded  (default -c 5)\n\n"; 
   print "-e  If set equal to 1 (-e 1), the empirical value for the insert size will not be used\n";   
   print "-pt  Total number of parts (in case several libraries are being inputted separately due to lack of memory)\n\n";

   print "-a  The directory that you used for output of prep_contig.pl\n\n";   
   
   #print "example: $0  -mate mate1.fasta -d L1 -mate mate2.fasta -d L2 -frag fragment1.fasta fragment2.fasta -a mydir_sout\n\n";
   exit;
}

if( $address =~ /^(.*)\/$/ ){
     $address=$1; 
}
$copy_li=$opt_c if ($opt_c);
if( $opt_e == 1 ){
     $ign_empi=1;
}

print "\nParsed SAM files:\n@parsed_sam_file\n\n";
print "Inser sizes:\n@distances\n\n";

print "Parsed SAM SOLiD files:\n@parsed_sam_solid\n\n";
print "Insert sizes:\n@distances_s\n\n";

print "Output directory:\n$address\n\n";
print "-c $copy_li\n-e $ign_empi\n\n";

if (  scalar(@parsed_sam_file) != scalar(@distances)  ){
   die "You have not specified valid insert size for each SAM file (use -d option)\n";
}
if (  scalar(@parsed_sam_solid) != scalar(@distances_s)  ){
   die "You have not specified valid insert size for each SAM (from SOLiD) file (use -d_solid option)\n";
}
my $date = `date`; chomp($date); print "Start $date\n"; 
my $tig_info=retrieve("$address"."/div/tig_info") or die "$!\n";


my $tig_info;

if( $opt_pt ){
	$tig_info=retrieve("$address"."/div/tig_info_p1") or die "$!\n";
	for (my $i=2; $i <= $opt_pt ; $i++){
		my $temp=retrieve("$address"."/div/tig_info_p$i") or die "$!\n";
		foreach my $tig (sort { $a <=> $b } keys %{$temp}){  
		     if ( defined $temp->{$tig}{'cov'} ){	
		     	$tig_info->{$tig}{'cov'}+=$temp->{$tig}{'cov'}; 
			}
		} 
	}

	open (CDIS, ">$address"."/div/coverage_distribution.txt");
	foreach my $tig (sort { $a <=> $b } keys %{$tig_info}){  
		 if ( $tig_info->{$tig}{'length'} == 0){
		     print "$tig\n";
		 } 
		 print CDIS "$tig $tig_info->{$tig}{'cov'}\n";
	}
	close CDIS;

	store \%{$tig_info}, "$address"."/div/tig_info_cov"; 
	
}else{
	$tig_info=retrieve("$address"."/div/tig_info") or die "$!\n";
}

#------------------------------------------------- paired-end
my $track; my $min_n=200; my $ORI_DIS;

my $num_parsed_sam_file=scalar(@parsed_sam_file)-1;
for (my $i=0; $i <= $num_parsed_sam_file ; $i++){
    open (IN,"<$parsed_sam_file[$i]") or die "Cannot open $parsed_sam_file[$i] $!\n";
    print "\nProcessing $parsed_sam_file[$i]\n";
    if ( $parsed_sam_file[$i] =~ /^(.*)_parsed$/ ){
	      if( $opt_pt ){
			  die "$parsed_sam_file[$i] format is not consistent witth option -pt $opt_pt\n";	       		  
		  }
    }elsif ( $parsed_sam_file[$i] =~ /^(.*)_parsed_p(\d+)$/ ){
	      if( $opt_pt ){
	      	  if ( $opt_pt == $2 ){	      	  
	      	  }else{
             	  
	      	  }
		  }else{
			  die "specify the option: -pt\n";	       
		  }
    }else{
		die "Invalid format $parsed_sam_file[$i]\n";	       
    }
    my $k; undef $track; my $insert=$distances[$i]; my @sizes; my ($dev_1,$dev_2); my $floor=1;
	  if($insert>5000){
		   $dev_1=1.1; $dev_2=2;
	  }elsif($insert>1000){
		   $dev_1=1.3; $dev_2=2.5;
	  }elsif($insert>500){
		   $dev_1=1.6; $dev_2=3.25;
	  }elsif($insert>250){
		   $dev_1=1.6; $dev_2=3.75;
	  }else{
		   $dev_1=1.6; $dev_2=4;
	  }  
    my $assoc_trim=$parsed_sam_file[$i];        
	if($assoc_trim =~ /(.*)\/([^\/]*)$/){
		 $assoc_trim=$2;          
	}           
    open (DDIS, ">$address"."/div/insertsize_distribution_$assoc_trim"."_$insert".".txt");
    open (RD, ">$address"."/div/read_length_$assoc_trim".".txt");
	while (my $line= <IN>){	
	   chomp($line); $floor++;      
	   if ($line=~ /^[+-]\t(\d+)\t(\d+)\t(\d+)\t(\d+)$/){   
		    my @column_a = split(/\t/,$line);		    
		    $line =<IN>; chomp($line); 		      		    
			if ( $column_a[3] < $copy_li ){			
				my @column_b = split(/\t/,$line);
				if ( $column_b[3] < $copy_li ){					
					print RD "$column_a[4] $column_b[4] ";
					my $tig_a=$column_a[1]; my $tig_b=$column_b[1];			
					my ($A_start,$A_end,$B_start,$B_end);	
					if ( $column_a[0] eq "+" ){
						  $A_start = $column_a[2];							      
						  $A_end=$A_start+$column_a[4];
					}else{
						  $A_end = $column_a[2];							      
						  $A_start=$A_end+$column_a[4];					
					}
					if ( $column_b[0] eq "+" ){
						  $B_start = $column_b[2];							      
						  $B_end=$B_start+$column_b[4];
					}else{
						  $B_end = $column_b[2];							      
						  $B_start=$B_end+$column_b[4];					
					}		
					if ( $tig_a == $tig_b){
						 if( $tig_info->{$tig_a}{'length'} > (2*$insert) ){				  
							  if ( $A_start < $A_end ){
									if ($B_start > $B_end ){#  -> <- 
										  my $imp= $B_start-$A_start;  print DDIS "$imp ";
										  my $d=$B_start-$A_start-$insert;
										  if ( (abs($d)<($dev_2*$insert)) && ($imp>0)){
											   push @sizes, $imp; 
										  }else{
										  }
									}else{     												   
									}
							  }else{
								  if ($B_start < $B_end ){#  ->  <-
										  my $imp= $A_start-$B_start; print DDIS "$imp ";
										  my $d=$A_start-$B_start-$insert;
										  if ( (abs($d)<($dev_2*$insert)) && ($imp>0)){
											   push @sizes, $imp;  
										  }else{														   
										  }
								  }else{     												   
								  }
							  }
						 }
					}else{
						  $k++;		
						  my $read_a=$k."a";
						  my $read_b=$k."b";
						  $track->{$read_a}{'tig'}=$tig_a;
						  $track->{$read_b}{'tig'}=$tig_b;		
						  $track->{$read_a}{'start'}=$A_start;
						  $track->{$read_a}{'end'}=$A_end; 
						  $track->{$read_b}{'start'}=$B_start;
						  $track->{$read_b}{'end'}=$B_end; 
					}
				}
			}
	   }else{
			die "Invalid format $parsed_sam_file[$i] :\n$line\n";	   
	   }
	   if ( $floor % 10000000 == 0 ){
		  print "inputed pairs $floor\n"; 
	   }  	
	}
	close IN; close DDIS; close RD;
	print "Associated library: $parsed_sam_file[$i]\n";
	if ( scalar(@sizes) > 0){
		  @sizes= sort @sizes;   
		  my $L=scalar(@sizes);    
		 my $sum=0;
		 foreach my $t (@sizes){
			   $sum+=$t; 
		  }
		  my $new_i=int ($sum/$L)+1;
		  if ( $L > $min_n ){
				 print "Suggested insert_size: $insert     Emperical insert_size: $new_i   (based on $L pairs)     ";
				 if ($ign_empi==1){
				 }else{
					$insert = $new_i;
				 }					 
				 print "We will use: $insert\n\n";
		  }else{
				 print "Suggested insert_size: $insert      Emperical insert_size: $new_i   (based on $L pairs , minimum required $min_n pairs)\n";
				 print "Not enough data points, keeping the suggested insert_size: $insert\n\n";
		  }		  
	}else{
			print "No data point to make emperical estimation of insert_size (no pair located on the same contig), keeping the suggested insert_size: $insert\n\n";
	}    
	get_oridis($track,$k,$insert);	
}

#------------------------------------------------- mate pair
$num_parsed_sam_file=scalar(@parsed_sam_solid)-1;      

for (my $i=0; $i <= $num_parsed_sam_file ; $i++){
    open (IN,"<$parsed_sam_solid[$i]") or die "Cannot open $parsed_sam_solid[$i] $!\n";
    print "\nProcessing $parsed_sam_solid[$i]\n";
    if ( $parsed_sam_solid[$i] =~ /^(.*)_parsed$/ ){
	      if( $opt_pt ){
			  die "$parsed_sam_file[$i] format is not consistent witth option -pt $opt_pt\n";	       		  
		  }
    }elsif ( $parsed_sam_solid[$i] =~ /^(.*)_parsed_p(\d+)$/ ){
	      if( $opt_pt ){
	      	  if ( $opt_pt == $2 ){	      	  
	      	  }else{
           	  
	      	  }
		  }else{
			  die "specify the option: -pt\n";	       
		  }
    }else{
		die "Invalid format $parsed_sam_solid[$i]\n";	       
    }
    my $k; undef $track; my $insert=$distances_s[$i]; my @sizes; my ($dev_1,$dev_2); my $floor=1;
	  if($insert>5000){
		   $dev_1=1.1; $dev_2=2;
	  }elsif($insert>1000){
		   $dev_1=1.3; $dev_2=2.5;
	  }elsif($insert>500){
		   $dev_1=1.6; $dev_2=3.25;
	  }elsif($insert>250){
		   $dev_1=1.6; $dev_2=3.75;
	  }else{
		   $dev_1=1.6; $dev_2=4;
	  }  
    my $assoc_trim=$parsed_sam_solid[$i];        
	if($assoc_trim =~ /(.*)\/([^\/]*)$/){
		 $assoc_trim=$2;          
	}           
    open (DDIS, ">$address"."/div/insertsize_distribution_$assoc_trim"."_$insert".".txt");
    open (RD, ">$address"."/div/read_length_$assoc_trim".".txt");
	while (my $line= <IN>){	
	   chomp($line); $floor++;      
	   if ($line=~ /^[+-]\t(\d+)\t(\d+)\t(\d+)\t(\d+)$/){   
		    my @column_a = split(/\t/,$line);		    
		    $line =<IN>; chomp($line); 		      		    
			if ( $column_a[3] < $copy_li ){			
				my @column_b = split(/\t/,$line);
				if ( $column_b[3] < $copy_li ){						
					print RD "$column_a[4] $column_b[4] ";
					my $tig_a=$column_a[1]; my $tig_b=$column_b[1];		
					my ($A_start,$A_end,$B_start,$B_end);	
					if ( $column_a[0] eq "+" ){
						  $A_start = $column_a[2];							      
						  $A_end=$A_start+$column_a[4];
					}else{
						  $A_end = $column_a[2];							      
						  $A_start=$A_end+$column_a[4];					
					}
					if ( $column_b[0] eq "+" ){
						  $B_start = $column_b[2];							      
						  $B_end=$B_start+$column_b[4];
					}else{
						  $B_end = $column_b[2];							      
						  $B_start=$B_end+$column_b[4];					
					}		
					if ( $tig_a == $tig_b){
						 if( $tig_info->{$tig_a}{'length'} > (2*$insert) ){				  
							  if ( $A_start < $A_end ){
									if ( $B_start < $B_end ){#  -> -> 	      
										  my $imp= $B_start-$A_start; print DDIS "$imp ";
										  my $d=$B_start-$A_start-$insert;	
										  if ( (abs($d)<($dev_2*$insert)) && ($imp>0)){
											   push @sizes, $imp; 
										  }else{
										  }
									}else{     
									}
							  }else{
								  if ($B_start > $B_end ){#  <-  <-
										  my $imp= $A_start-$B_start; print DDIS "$imp ";
										  my $d=$A_start-$B_start-$insert;									  
										  if ( (abs($d)<($dev_2*$insert)) && ($imp>0)){
											   push @sizes, $imp;  
										  }else{
										  }
								  }else{     
								  }
							  }	  
						 }	
					}else{
						  $k++;		
						  my $read_a=$k."a";
						  my $read_b=$k."b";
						  $track->{$read_a}{'tig'}=$tig_a;
						  $track->{$read_b}{'tig'}=$tig_b;		
						  $track->{$read_a}{'start'}=$A_start;
						  $track->{$read_a}{'end'}=$A_end; 
						  $track->{$read_b}{'start'}=$B_start;
						  $track->{$read_b}{'end'}=$B_end; 
					}
				}
			}
	   }else{
			die "Invalid format $parsed_sam_solid[$i] :\n$line\n";	   
	   }
	   if ( $floor % 10000000 == 0 ){
		  print "inputed pairs $floor\n"; 
	   }  	   
	}
	close IN; close DDIS; close RD;
	print "Associated library: $parsed_sam_solid[$i]\n";
	if ( scalar(@sizes) > 0){
		  @sizes= sort @sizes;   
		  my $L=scalar(@sizes);    
		 my $sum=0;
		 foreach my $t (@sizes){
			   $sum+=$t; 
		  }
		  my $new_i=int ($sum/$L)+1;
		  if ( $L > $min_n ){
				 print "Suggested insert_size: $insert      Emperical insert_size: $new_i   (based on $L pairs)      ";
				 if ($ign_empi==1){
				 }else{
					$insert = $new_i;
				 }					 
				 print "We will use: $insert\n\n";
		  }else{
				 print "Suggested insert_size: $insert      Emperical insert_size: $new_i   (based on $L pairs , minimum required $min_n pairs)\n";
				 print "Not enough data points, keeping the suggested insert_size: $insert\n\n";
		  }		  
	}else{
			print "No data point to make emperical estimation of insert_size (no pair located on the same contig), keeping the suggested insert_size: $insert\n\n";
	}    
	get_oridis_solid($track,$k,$insert);	
}
#------------------------------------------- 

store \%{$ORI_DIS}, "$address"."/orientdistinfo_c$copy_li"; 

#------------------------------------------- 
sub get_oridis{
    my ($track,$k,$insert_size) = @_;  
	print "Number of links connecting contigs $k\n\n";  my $kharab1; my $kharab2; my $kharab3; my $kharab4;
	for(my $j=1; $j<=$k; $j++){
		my $read_a=$j."a"; my $read_b=$j."b";	
		if((defined $track->{$read_a}) && (defined $track->{$read_b})){### both pairs assembled
				 if((!defined $track->{$read_a}{'tig'}) || (!defined $track->{$read_b}{'tig'})){
					print "unexpected problem\n"; 
					print "$read_a $track->{$read_a}{'tig'}\n";
					print "$read_b $track->{$read_b}{'tig'}\n";
					exit;
				 }
				 my $dev_1; my $dev_2;
				 if($insert_size>5000){
					   $dev_1=1.1; $dev_2=.2;
				 }elsif($insert_size>1000){
					   $dev_1=1.3; $dev_2=.5;
				 }else{
					   $dev_1=1.6; $dev_2=1;
				 }          
				 my $tig_a = $track->{$read_a}{'tig'};
				 my $tig_b = $track->{$read_b}{'tig'};    
				 my $A_length = $tig_info->{$tig_a}{'length'};
				 my $A_start = $track->{$read_a}{'start'};
				 my $A_end = $track->{$read_a}{'end'};
				 my $B_length = $tig_info->{$tig_b}{'length'};
				 my $B_start = $track->{$read_b}{'start'} ;
				 my $B_end = $track->{$read_b}{'end'};	  
				 
				 if($tig_a != $tig_b){####paired reads located on <> contigs
						if ($track->{$read_a}{'start'} < $track->{$read_a}{'end'}){
							  if ($track->{$read_b}{'start'} > $track->{$read_b}{'end'} ){#  -> <- 
								 my $d=$insert_size+$A_start-$B_start;
				  
								 my $gap=$d-$tig_info->{$tig_a}{'length'};
								 if(abs($gap)<($dev_1*$insert_size)){
										 if ($tig_a<$tig_b){
											$ORI_DIS->{$tig_a}{$tig_b}{'P'}++;
											$ORI_DIS->{$tig_a}{$tig_b}{'PD'}+=$d;											   
										 }else{    
											$ORI_DIS->{$tig_b}{$tig_a}{'P'}++;
											$ORI_DIS->{$tig_b}{$tig_a}{'PD'}-=$d;				                               
										 }							 
								  }else{
									 $kharab1++;
								  }
							  }else{# -> ->
								 my $d=$insert_size+$A_start+$B_start;
														 
								 my $gap=$d-$tig_info->{$tig_a}{'length'}-$tig_info->{$tig_b}{'length'};
								 if(abs($gap)<($dev_1*$insert_size)){
										 if ($tig_a<$tig_b){
											$ORI_DIS->{$tig_a}{$tig_b}{'N'}++;
											$ORI_DIS->{$tig_a}{$tig_b}{'ND'}+=$d;
										 }else{
											$ORI_DIS->{$tig_b}{$tig_a}{'N'}++;
											$ORI_DIS->{$tig_b}{$tig_a}{'ND'}+=$d;												
										 }
								 }else{                              
									 $kharab2++;
								 }
							  }
						}else{
							  if ($track->{$read_b}{'start'} < $track->{$read_b}{'end'}){#  <-  ->
								 my $d=-($insert_size-$A_start+$B_start);	
				  
								 my $gap=-$d-$tig_info->{$tig_b}{'length'};
								 if(abs($gap)<($dev_1*$insert_size)){
										 if ($tig_a<$tig_b){
											$ORI_DIS->{$tig_a}{$tig_b}{'P'}++;
											$ORI_DIS->{$tig_a}{$tig_b}{'PD'}+=$d;
										 }else{
											$ORI_DIS->{$tig_b}{$tig_a}{'P'}++;
											$ORI_DIS->{$tig_b}{$tig_a}{'PD'}-=$d;												
										 }
								 }else{							 
									 $kharab3++;
								 }
							  }else{    #   <- <-  
								 my $d=-($insert_size-$A_start-$B_start);
								 
								 my $gap=-$d;
								 if(abs($gap)<($dev_1*$insert_size)){
										 if ($tig_a<$tig_b){
											$ORI_DIS->{$tig_a}{$tig_b}{'N'}++;
											$ORI_DIS->{$tig_a}{$tig_b}{'ND'}+=$d;					
										 }else{
											$ORI_DIS->{$tig_b}{$tig_a}{'N'}++;
											$ORI_DIS->{$tig_b}{$tig_a}{'ND'}+=$d;											
										 }
								 }else{							 
									 $kharab4++
								 }									 
							  }
						}
				 }else{ ###Clone, paired reads located on the same contig -- could be used to investigate misassemblies
				 
				 } 
		}else{
			print "o 1\n"; exit;
		}
	} 
	#print "kharab: $kharab1 , $kharab2 , $kharab3 , $kharab4\n";
}	
#------------------------------------------- 
sub get_oridis_solid{
    my ($track,$k,$insert_size) = @_;  
	print "number of links connecting contigs $k\n\n";  my $kharab1; my $kharab2; my $kharab3; my $kharab4;
	for(my $j=1; $j<=$k; $j++){
		my $read_a=$j."a"; my $read_b=$j."b";	
		if((defined $track->{$read_a}) && (defined $track->{$read_b})){### both pairs assembled
				 if((!defined $track->{$read_a}{'tig'}) || (!defined $track->{$read_b}{'tig'})){
					print "unexpected problem\n"; 
					print "$read_a $track->{$read_a}{'tig'}\n";
					print "$read_b $track->{$read_b}{'tig'}\n";
					exit;
				 }
				 my $dev_1; my $dev_2;
				 if($insert_size>5000){
					   $dev_1=1.1; $dev_2=.2;
				 }elsif($insert_size>1000){
					   $dev_1=1.3; $dev_2=.5;
				 }else{
					   $dev_1=1.6; $dev_2=1;
				 }          
				 my $tig_a = $track->{$read_a}{'tig'};
				 my $tig_b = $track->{$read_b}{'tig'};    
				 my $A_length = $tig_info->{$tig_a}{'length'};
				 my $A_start = $track->{$read_a}{'start'};
				 my $A_end = $track->{$read_a}{'end'};
				 my $B_length = $tig_info->{$tig_b}{'length'};
				 my $B_start = $track->{$read_b}{'start'} ;
				 my $B_end = $track->{$read_b}{'end'};	  
				 
				 if($tig_a != $tig_b){####paired reads located on <> contigs
						if ($track->{$read_a}{'start'} < $track->{$read_a}{'end'}){
							  if ($track->{$read_b}{'start'} < $track->{$read_b}{'end'} ){#  -> -> 
								 my $d=$insert_size+$A_start-$B_start;
				  
								 my $gap=$d-$tig_info->{$tig_a}{'length'};
								 if(abs($gap)<($dev_1*$insert_size)){
										 if ($tig_a<$tig_b){
											$ORI_DIS->{$tig_a}{$tig_b}{'P'}++;
											$ORI_DIS->{$tig_a}{$tig_b}{'PD'}+=$d;
										 }else{    
											$ORI_DIS->{$tig_b}{$tig_a}{'P'}++;
											$ORI_DIS->{$tig_b}{$tig_a}{'PD'}-=$d;
										 }									 
								  }else{
									 $kharab1++;
								  }
							  }else{# -> <-
								 my $d=$insert_size+$A_start+$B_start;
														 
								 my $gap=$d-$tig_info->{$tig_a}{'length'}-$tig_info->{$tig_b}{'length'};
								 if(abs($gap)<($dev_1*$insert_size)){
										 if ($tig_a<$tig_b){
											$ORI_DIS->{$tig_a}{$tig_b}{'N'}++;
											$ORI_DIS->{$tig_a}{$tig_b}{'ND'}+=$d;
										 }else{
											$ORI_DIS->{$tig_b}{$tig_a}{'N'}++;
											$ORI_DIS->{$tig_b}{$tig_a}{'ND'}+=$d;
										 }
								 }else{                             
									 $kharab2++;
								 }
							  }
						}else{
							  if ($track->{$read_b}{'start'} > $track->{$read_b}{'end'}){#  <-  <-
								 my $d=-($insert_size-$A_start+$B_start);	
				  
								 my $gap=-$d-$tig_info->{$tig_b}{'length'};
								 if(abs($gap)<($dev_1*$insert_size)){
										 if ($tig_a<$tig_b){
											$ORI_DIS->{$tig_a}{$tig_b}{'P'}++;
											$ORI_DIS->{$tig_a}{$tig_b}{'PD'}+=$d;
										 }else{
											$ORI_DIS->{$tig_b}{$tig_a}{'P'}++;
											$ORI_DIS->{$tig_b}{$tig_a}{'PD'}-=$d;
										 }
								 }else{							 
									 $kharab3++;
								 }
							  }else{    #   <- ->  
								 my $d=-($insert_size-$A_start-$B_start);
								 
								 my $gap=-$d;
								 if(abs($gap)<($dev_1*$insert_size)){
										 if ($tig_a<$tig_b){
											$ORI_DIS->{$tig_a}{$tig_b}{'N'}++;
											$ORI_DIS->{$tig_a}{$tig_b}{'ND'}+=$d;
										 }else{
											$ORI_DIS->{$tig_b}{$tig_a}{'N'}++;
											$ORI_DIS->{$tig_b}{$tig_a}{'ND'}+=$d;
										 }
								 }else{						 
									 $kharab4++
								 }									 
							  }
						}
				 }else{ ###Clone, paired reads located on the same contig -- could be used to investigate misassemblies
				 
				 } 
		}else{
			print "o 1\n"; exit;
		}
	} 
	#print "kharab: $kharab1 , $kharab2 , $kharab3 , $kharab4\n";
}	
#------------------------------------------- 
sub reverseComplement{
   $_ = shift;
   $_ = uc();
   tr/ATGC/TACG/;
   return (reverse());
}
#--------------------------------------------