#!/usr/bin/perl

#	AUTHOR
#  Adel Dayarian

#	The code is available freely, under the GNU Public License, at http://www.physics.rutgers.edu/~anirvans/SOPRA/

use strict;
use Storable;
use Getopt::Long;

my @sam_file;
my @bam_file;
my $sam_tools_address="";
my $opt_p;
my $address="";

&GetOptions( "sam=s{,}" => \@sam_file,
			 "bam=s{,}" => \@bam_file,            
     	     "st_p=s", => \$sam_tools_address,
     	     "p=i", => \$opt_p,
     	     "a=s", => \$address,
           ) or die "Unknown option\n";
           
my ($copy_li,$ign_empi,$s_l_u,$s_l_d)=(5,0,0,0);

#-------------------------------------------------   
if(  (! @sam_file)  || ( $address eq "" ) ){
   print "Usage: $0\n\n";
   print "-sam  File(s) in SAM format containing the result of alignment against contigs\n\n";

   print "-bam  File(s) in BAM format containing the result of alignment against contigs, make sure Samtools is installed and the BAM file is indexed. You also need to specify option -st_p\n";
   print "-st_p  Path to the directory where Samtools is installed (to find out, in the terminal type: which samtools)\n\n";
   
   print "-p  Part number (in case several libraries are being inputted separately due to lack of memory)\n\n";
   print "-a  The directory that you used for output of prep_contig.pl\n\n";   
   
   #print "example: $0  -mate mate1.fasta -d L1 -mate mate2.fasta -d L2 -frag fragment1.fasta fragment2.fasta -a mydir_sout\n\n";
   exit;
}
my $date = `date`; chomp($date); print "\nStart $date\n\n"; 

if( $address =~ /^(.*)\/$/ ){
     $address=$1; 
}

if( $sam_tools_address =~ /^(.*)\/$/ ){
     $sam_tools_address=$1; 
}

if( $opt_p ){
	print "Running: $0\n-sam @sam_file\n-bam @bam_file\n-st_p $sam_tools_address\n-a $address\n-p $opt_p\n\n";
}else{
    print "Running: $0\n-sam @sam_file\n-bam @bam_file\n-st_p $sam_tools_address\n-a $address\n\n";
}

my $tig_info=retrieve("$address"."/div/tig_info") or die "$!\n"; my $id;
#------------------------------------------------- calculate copy number 
my $bin; my $num_re=0; 
my $num_sam_file=scalar(@sam_file)-1;
for (my $i=0; $i <= $num_sam_file ; $i++){
    open (IN,"<$sam_file[$i]") or die "Cannot open $sam_file[$i] $!\n";
	my $date = `date`; chomp($date); print "Processing $sam_file[$i]   $date\n";
    my $current; my $whole; my $floor=1;
	my $line= <IN>;	
	while ( $line=~ /^[\@HD,\@SQ,\@PG]/ ){
	   		$line= <IN>;	   
    }
	chomp($line); 	
	my @column = split(/\t/,$line);
	$current = $column[0] ; #print "first id $current\n";
	$whole=$line;

	while (my $line= <IN>){	   
	   chomp($line); $floor++;
       my @column = split(/\t/,$line);
	   my $next = $column[0];	  
	   my $tedad=1;
	   while ( $current eq $next ){	   
	    		$line= <IN>; chomp($line); $floor++;
	    		my @column = split(/\t/,$line);
	   			$next = $column[0] ;	 
	   			$tedad++;
				if ( $current eq $next ){	    
					if ( $column[9]=~ /^([ATCGN]*)$/){
						  my $read=$column[9];
						  my $tig;
						  if ( $column[2]=~ /^contig(\d+)|/){
								  $tig = $1;
						  }else{	
								die "$sam_file[$i] is not in the right format: $column[2]\n";	   	     
						  }
						  $tig_info->{$tig}{'cov'}+=length($read);
					}else{	
						die "$sam_file[$i] is not in the right format: $column[9]\n";	   	     
					}
				}	
	   }
	   if ( $tedad == 1 ){
	   	     my @column = split(/\t/,$whole);	   	     
			 
			if ( ($column[1] & (1<<2)) !=0 ){
				$column[1]=4;
			}else{
				if ( ( $column[1] & (1<<4)) !=0 ){
					$column[1]=16;
				}else{
					$column[1]=0;
				}	
			}
	   	     if ( ($column[1] == 0 ) || ( $column[1] == 16 ) ){
					if( $column[0] =~ /^s(\d*)a$/i){  						
						@{$id->{$1}{'a'}}=("+" );   		
					}elsif( $column[0] =~ /^s(\d*)b$/i){  	 
						 @{$id->{$1}{'b'}}=("+");   						
					}else{	
		                die "$sam_file[$i] is not in the right format: $column[0]\n";	   	     
					}
					if ( $column[9]=~ /^([ATCGN]*)$/){
					      my $read=$column[9];
						  if( defined $bin->{$read} ){
							   $bin->{$read}++;
						  }else{
							   $bin->{reverseComplement($read)}++;	
						  }	
						   if( $s_l_d == 0 ){
							   $s_l_d=length($read);
						   }elsif( length($read) < $s_l_d ){
							   $s_l_d=length($read);
						   }
						   if( length($read) > $s_l_u ){
							   $s_l_u=length($read);
						   }							   
						   my $tig;
						   if ( $column[2]=~ /^contig(\d+)|/){
								  $tig = $1;
						   }else{	
								die "$sam_file[$i] is not in the right format: $column[2]\n";	   	     
						   }
						   $tig_info->{$tig}{'cov'}+=length($read);
					}else{	
		                die "$sam_file[$i] is not in the right format: $column[9]\n";	   	     
					}	      
			 }
	   }else{ 
	   	     $num_re++;
	   }	
	   $current= $next;
	   $whole=$line;
	   if ( $floor % 5000000 == 0 ){
		  print "line $floor\n"; 
	   }  
	}
	close IN; 
}
#------------------------------------------------- more  copy number 
my $num_bam_file=scalar(@bam_file)-1;
for (my $i=0; $i <= $num_bam_file ; $i++){
	open (IN,"$sam_tools_address view $bam_file[$i] |") or die "Cannot open bam_file[$i] $!\n";
	my $date = `date`; chomp($date); print "Processing $bam_file[$i]   $date\n";
    my $current; my $whole; my $floor=1;
	my $line= <IN>;	
	while ( $line=~ /^[\@HD,\@SQ,\@PG]/ ){
	   		$line= <IN>;	   
    }
	chomp($line); 	
	my @column = split(/\t/,$line);
	$current = $column[0] ; #print "first id $current\n";
	$whole=$line;

	while (my $line= <IN>){	   
	   chomp($line); $floor++;
       my @column = split(/\t/,$line);
	   my $next = $column[0] ;	  
	   my $tedad=1;
	   while ( $current eq $next ){	   
	    		$line= <IN>; chomp($line); $floor++;
	    		my @column = split(/\t/,$line);
	   			$next = $column[0] ;	 
	   			$tedad++;
				if ( $current eq $next ){	    
					if ( $column[9]=~ /^([ATCGN]*)$/){
						  my $read=$column[9];
						  my $tig;
						  if ( $column[2]=~ /^contig(\d+)|/){
								  $tig = $1;
						  }else{	
								die "$bam_file[$i] is not in the right format: $column[2]\n";	   	     
						  }
						  $tig_info->{$tig}{'cov'}+=length($read);
					}else{	
						die "$bam_file[$i] is not in the right format: $column[9]\n";	   	     
					}
				}	
	   }
	   if ( $tedad == 1 ){
	   	     my @column = split(/\t/,$whole);	   	     
			 
			if ( ($column[1] & (1<<2)) !=0 ){
				$column[1]=4;
			}else{
				if ( ( $column[1] & (1<<4)) !=0 ){
					$column[1]=16;
				}else{
					$column[1]=0;
				}	
			}
	   	     if ( ($column[1] == 0 ) || ( $column[1] == 16 ) ){
					if( $column[0] =~ /^s(\d*)a$/i){  						
						@{$id->{$1}{'a'}}=("+" );   		
					}elsif( $column[0] =~ /^s(\d*)b$/i){  	 
						@{$id->{$1}{'b'}}=("+");   						
					}else{	
		                die "$bam_file[$i] is not in the right format: $column[0]\n";	   	     
					}
					if ( $column[9]=~ /^([ATCGN]*)$/){
					      my $read=$column[9];
						  if( defined $bin->{$read} ){
							   $bin->{$read}++;
						  }else{
							   $bin->{reverseComplement($read)}++;	
						  }	
						   if( $s_l_d == 0 ){
							   $s_l_d=length($read);
						   }elsif( length($read) < $s_l_d ){
							   $s_l_d=length($read);
						   }
						   if( length($read) > $s_l_u ){
							   $s_l_u=length($read);
						   }	
						   my $tig;
						   if ( $column[2]=~ /^contig(\d+)|/){
								  $tig = $1;
						   }else{	
								die "$bam_file[$i] is not in the right format: $column[2]\n";	   	     
						   }
						   $tig_info->{$tig}{'cov'}+=length($read);
					}else{	
		                die "$bam_file[$i] is not in the right format: $column[9]\n";	   	     
					}	      
			 }
	   }else{ 
	   	     $num_re++;
	   }	
	   $current= $next;
	   $whole=$line;
	   if ( $floor % 5000000 == 0 ){
		  print "line $floor\n"; 
	   }  
	}
	close IN; 
}
print "\nmax read length $s_l_u bp      min read length $s_l_d bp\n\n";
$tig_info->{1}{'seq_l_u'}=$s_l_u;
$tig_info->{1}{'seq_l_d'}=$s_l_d;

#-------------------------------------------------  
if( $opt_p ){
	open (CN, ">$address"."/div/copynumber_p$opt_p".".txt" ) or die "$!";
}else{
	open (CN, ">$address"."/div/copynumber.txt" ) or die "$!";
}
foreach my $r ( sort { $bin->{$b} <=> $bin->{$a} }  keys %{$bin} ){
     print CN "$bin->{$r} ";
}
close CN; 

if( $opt_p ){
	open (CDIS, ">$address"."/div/coverage_distribution_p$opt_p".".txt");
	foreach my $tig (sort { $a <=> $b } keys %{$tig_info}){  
		 $tig_info->{$tig}{'cov'}= .001 * int(1000 * $tig_info->{$tig}{'cov'}/$tig_info->{$tig}{'length'});
		 print CDIS "$tig $tig_info->{$tig}{'cov'}\n";
	}
	close CDIS;
	store \%{$tig_info}, "$address"."/div/tig_info_p$opt_p"; 
	
}else{
	open (CDIS, ">$address"."/div/coverage_distribution.txt");
	foreach my $tig (sort { $a <=> $b } keys %{$tig_info}){  
		 if ( $tig_info->{$tig}{'length'} == 0){
		     print "$tig\n";
		 } 
		 $tig_info->{$tig}{'cov'}= .01 * int(100 * $tig_info->{$tig}{'cov'}/$tig_info->{$tig}{'length'});
		 print CDIS "$tig $tig_info->{$tig}{'cov'}\n";
	}
	close CDIS;
	
	store \%{$tig_info}, "$address"."/div/tig_info_cov"; 
}
#-------------------------------------------------  
my $tot;
$num_sam_file=scalar(@sam_file)-1;
for (my $i=0; $i <= $num_sam_file ; $i++){
	if ( ( $sam_file[$i]=~ /.zip$/) || ( $sam_file[$i]=~ /.ZIP$/) ){ 
    	open (IN,"unzip -p $sam_file[$i] |") or die "Cannot open unzip -p sam_file[$i] $!\n";	
    }elsif (  $sam_file[$i]=~ /.gz$/ ){ 
    	open (IN,"gunzip -c $sam_file[$i] |") or die "Cannot open gunzip -c sam_file[$i] $!\n";		
	}else{
    	open (IN,"<$sam_file[$i]") or die "Cannot open sam_file[$i] $!\n";
    }
	my $date = `date`; chomp($date); print "Sorting pairs, $sam_file[$i]   $date\n";
    my $current; my $whole; my $floor=1; my $local_id; my $one_leg=0;
	my $line= <IN>;	
	while ( $line=~ /^[\@HD,\@SQ,\@PG]/ ){
	   		$line= <IN>;	   
    }
	chomp($line); 	
	my @column = split(/\t/,$line);
	$current = $column[0] ; #print "first id $current\n";
	$whole=$line;
	while (my $line= <IN>){	   
	   chomp($line); $floor++;
       my @column = split(/\t/,$line);
	   my $next = $column[0] ;	  
	   my $tedad=1;
	   while ( $current eq $next ){	   
	    		$line= <IN>; chomp($line); 	$floor++;
	    		my @column = split(/\t/,$line);
	   			$next = $column[0] ;	 
	   			$tedad++;
	   }	   
	   @column = split(/\t/,$whole); my $tag;	   
	   if( $column[0] =~ /^s(\d*)a$/i){  						
			$tag=$1;
	   }elsif( $column[0] =~ /^s(\d*)b$/i){  	
			$tag=$1;					
	   }else{
	   		die "$sam_file[$i] is not in the right format: $column[0]\n";	   	   
	   }
	   if ( (defined $id->{$tag}{'a'}) && (defined $id->{$tag}{'b'}) ) { 
			if ( ($column[1] & (1<<2)) !=0 ){
				$column[1]=4;
			}else{
				if ( ( $column[1] & (1<<4)) !=0 ){
					$column[1]=16;
				}else{
					$column[1]=0;
				}	
			}	   	     
	   	    if ( ($column[1] == 0 ) || ( $column[1] == 16 ) ){
					my $tig; my $start; my $read; my $tcount;
					if ( $column[2]=~ /^contig(\d+)|/){
					      $tig = $1;
					}else{	
		                die "$sam_file[$i] is not in the right format: $column[2]\n";	   	     
					}
					if ( $column[3]=~ /^(\d+)$/){
						  $start = $1;		
					}else{	
		                die "$sam_file[$i] is not in the right format: $column[3]\n";	   	     
					}	
					if ( $column[9]=~ /^([ATCGN]*)$/){
					      $read=$column[9];
					}else{	
		                die "$sam_file[$i] is not in the right format: $column[9]\n";	   	     
					}	      
					if( defined $bin->{$read} ){
						$tcount=$bin->{$read};
					}else{
						$tcount=$bin->{reverseComplement($read)};
					}				
					my $L=length($read); 
					if( $column[0] =~ /^s(\d*)a$/i){  						
						  if ( $column[1] == 0 ) {
						  	  @{$local_id->{$1}{'a'}}=("+" , "$tig" , "$start" , "$tcount", "$L" );   		
						  }elsif ( $column[1] == 16 ){ 
						  	  @{$local_id->{$1}{'a'}}=("-" , "$tig" , "$start" , "$tcount", "$L" );   								  
						  }
					}elsif( $column[0] =~ /^s(\d*)b$/i){  	 
						  if ( $column[1] == 0 ) {
						  	  @{$local_id->{$1}{'b'}}=("+" , "$tig" , "$start" , "$tcount", "$L" );   		
						  }elsif ( $column[1] == 16 ) {
						  	  @{$local_id->{$1}{'b'}}=("-" , "$tig" , "$start" , "$tcount", "$L" );   								  
						  }					
					}else{	
		                die "$sam_file[$i] is not in the right format: $column[0]\n";	   	     
					}					
			}else{
				 print "Unexpected error\n"; exit;
			}
	   }else{
			if ( (defined $id->{$tag}{'a'}) || (defined $id->{$tag}{'b'}) ) { 
		  		$one_leg++;
			}
	   }	
	   $current= $next ;
	   $whole=$line;
	   if ( $floor % 5000000 == 0 ){
		  print "line $floor\n"; 
	   }  
	}
	close IN; 
	if( $opt_p ){
	    open (OUT,">$sam_file[$i]"."_parsed_p$opt_p") or die "$!\n"; 
    }else{
	    open (OUT,">$sam_file[$i]"."_parsed") or die "$!\n"; 
    }
    my $loc_tot; 
	foreach my $h (keys %{$local_id}){
		  if ( (defined $local_id->{$h}{'a'}) && (defined $local_id->{$h}{'b'}) ) { 
		        $loc_tot++; $tot++;
		        print OUT "${$local_id->{$h}{'a'}}[0]\t${$local_id->{$h}{'a'}}[1]\t${$local_id->{$h}{'a'}}[2]\t${$local_id->{$h}{'a'}}[3]\t${$local_id->{$h}{'a'}}[4]\n";
		        print OUT "${$local_id->{$h}{'b'}}[0]\t${$local_id->{$h}{'b'}}[1]\t${$local_id->{$h}{'b'}}[2]\t${$local_id->{$h}{'b'}}[3]\t${$local_id->{$h}{'a'}}[4]\n";
		  }
	}
	close OUT; $one_leg=int(.5*$one_leg);	
	#print "Number of pairs for which only one leg mapped $one_leg\n";	
	print "Number of useful pairs $loc_tot\n\n";
}
#------------------------------------------- 
$num_bam_file=scalar(@bam_file)-1;
for (my $i=0; $i <= $num_bam_file ; $i++){	
	open (IN,"$sam_tools_address view $bam_file[$i] |") or die "Cannot open bam_file[$i] $!\n";
	my $date = `date`; chomp($date); print "Sorting pairs, $bam_file[$i]   $date\n";
    my $current; my $whole; my $floor=1; my $local_id; my $one_leg=0;
	my $line= <IN>;	
	while ( $line=~ /^[\@HD,\@SQ,\@PG]/ ){
	   		$line= <IN>;	   
    }
	chomp($line); 	
	my @column = split(/\t/,$line);
	$current = $column[0] ; #print "first id $current\n";
	$whole=$line;
	while (my $line= <IN>){	   
	   chomp($line); $floor++;
       my @column = split(/\t/,$line);
	   my $next = $column[0] ;	  
	   my $tedad=1;
	   while ( $current eq $next ){	   
	    		$line= <IN>; chomp($line); 	$floor++;
	    		my @column = split(/\t/,$line);
	   			$next = $column[0] ;	 
	   			$tedad++;
	   }	   
	   @column = split(/\t/,$whole); my $tag;	   
	   if( $column[0] =~ /^s(\d*)a$/i){  						
			$tag=$1;
	   }elsif( $column[0] =~ /^s(\d*)b$/i){  	
			$tag=$1;					
	   }else{
	   		die "$sam_file[$i] is not in the right format: $column[0]\n";	   	   
	   }	   
	   if ( (defined $id->{$tag}{'a'}) && (defined $id->{$tag}{'b'}) ) { 	   
	   	     my @column = split(/\t/,$whole);
			if ( ($column[1] & (1<<2)) !=0 ){
				$column[1]=4;
			}else{
				if ( ( $column[1] & (1<<4)) !=0 ){
					$column[1]=16;
				}else{
					$column[1]=0;
				}	
			}	   	     
	   	     if ( ($column[1] == 0 ) || ( $column[1] == 16 ) ){
					my $tig; my $start; my $read; my $tcount;
					if ( $column[2]=~ /^contig(\d+)|/){
					      $tig = $1;
					}else{	
		                die "$bam_file[$i] is not in the right format: $column[2]\n";	   	     
					}
					if ( $column[3]=~ /^(\d+)$/){
						  $start = $1;		
					}else{	
		                die "$bam_file[$i] is not in the right format: $column[3]\n";	   	     
					}	
					if ( $column[9]=~ /^([ATCGN]*)$/){
					      $read=$column[9];
					}else{	
		                die "$bam_file[$i] is not in the right format: $column[9]\n";	   	     
					}	      
					if( defined $bin->{$read} ){
						$tcount=$bin->{$read};
					}else{
						$tcount=$bin->{reverseComplement($read)};
					}	
					my $L=length($read);
					if( $column[0] =~ /^s(\d*)a$/i){  						
						  if ( $column[1] == 0 ) {
						  	  @{$local_id->{$1}{'a'}}=("+" , "$tig" , "$start" , "$tcount", "$L" );   		
						  }elsif ( $column[1] == 16 ){ 
						  	  @{$local_id->{$1}{'a'}}=("-" , "$tig" , "$start" , "$tcount", "$L" );   								  
						  }
					}elsif( $column[0] =~ /^s(\d*)b$/i){  	 
						  if ( $column[1] == 0 ) {
						  	  @{$local_id->{$1}{'b'}}=("+" , "$tig" , "$start" , "$tcount", "$L" );   		
						  }elsif ( $column[1] == 16 ) {
						  	  @{$local_id->{$1}{'b'}}=("-" , "$tig" , "$start" , "$tcount", "$L" );   								  
						  }								
					}else{	
		                die "$bam_file[$i] is not in the right format: $column[0]\n";	   	     
					}					
			 }else{
				 print "Unexpected error\n"; exit;
			}
	   }else{
			if ( (defined $id->{$tag}{'a'}) || (defined $id->{$tag}{'b'}) ) { 
		  		$one_leg++;
			}
	   }	
	   $current= $next ;
	   $whole=$line;
	   if ( $floor % 5000000 == 0 ){
		  print "line $floor\n"; 
	   }  
	}
	close IN; 
	if( $opt_p ){
	    open (OUT,">$bam_file[$i]"."_parsed_p$opt_p") or die "$!\n"; 
    }else{
	    open (OUT,">$bam_file[$i]"."_parsed") or die "$!\n"; 
    }
    my $loc_tot;    
	foreach my $h (keys %{$local_id}){
		  if ( (defined $local_id->{$h}{'a'}) && (defined $local_id->{$h}{'b'}) ) { 
		        $loc_tot++; $tot++;
		        print OUT "${$local_id->{$h}{'a'}}[0]\t${$local_id->{$h}{'a'}}[1]\t${$local_id->{$h}{'a'}}[2]\t${$local_id->{$h}{'a'}}[3]\t${$local_id->{$h}{'a'}}[4]\n";
		        print OUT "${$local_id->{$h}{'b'}}[0]\t${$local_id->{$h}{'b'}}[1]\t${$local_id->{$h}{'b'}}[2]\t${$local_id->{$h}{'b'}}[3]\t${$local_id->{$h}{'a'}}[4]\n";
		  }
	}
	close OUT;	
	#print "Number of pairs for which only one leg mapped $one_leg\n";	
	print "Number of useful pairs $loc_tot\n\n";
}

print "Total number of useful pairs $tot\n\n";
print "Removed $num_re repetitive reads\n";

$date = `date`; chomp($date); print "\nFinished $date\n\n"; 

#------------------------------------------- 
sub reverseComplement{
   $_ = shift;
   $_ = uc();
   tr/ATGC/TACG/;
   return (reverse());
}
#--------------------------------------------