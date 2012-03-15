#!/usr/bin/perl

#	AUTHOR
#  Adel Dayarian

#	The code is available freely, under the GNU Public License, at http://www.physics.rutgers.edu/~anirvans/SOPRA/

use Storable;
use strict; 
use Data::Dumper;
use Getopt::Long;

my $opt_w;
my $opt_L;
my $opt_u;
my $opt_s;
my $opt_d;
my $opt_k;
my $opt_o="";
my $opt_v;
my $opt_h;
my $opt_a;
my $non_seq;
my $opt_t;

my $address="";

&GetOptions( "weight=i",  => \$opt_w,
     	     "length=i",  => \$opt_L,
     	     "t=i",  => \$opt_t,
     	     "u=i",  => \$opt_u,
     	     "s=i",  => \$opt_s,
     	     "d=i",  => \$opt_d,
     	     "k=i",  => \$opt_k,
     	     "o=s",  => \$opt_o,
     	     "v=i",  => \$opt_v,
     	     "h=f",  => \$opt_h,
     	     "n=s", => \$non_seq,
     	     "a=s", => \$opt_a,
           ) or die "Unknown option\n";
#-------------------------------------------------   
my ($minlink,$minlength,$high_cov)=(4,150,2.2); my $stretch_thr=1500;            
my ($tig_info, $Ori_Dis, $verbs, @lengths, $tah); my $strength=1;

$minlink = $opt_w if ($opt_w);  $minlength = $opt_L if ($opt_L); $verbs = 1 if ($opt_v); 
$high_cov = $opt_h if ($opt_h); $strength = $opt_t if ($opt_t);

if( ($opt_o=~/(.*)orientdistinfo_c(\d*)$/) && ( $opt_a) ){
       $tah=$2;
       $address= $opt_a;
       if( $address =~ /^(.*)\/$/ ){
           $address=$1; 
       }

       $tig_info=retrieve("$address"."/div/tig_info_cov") or die;  
       $Ori_Dis=retrieve("$opt_o") or die;
       if( $non_seq ){
       		$Ori_Dis =input_non_seq($non_seq,$Ori_Dis);   
       }
}elsif( ( $non_seq ) && ( $opt_a) ){       
	   print "No mater pair/paired-end library\n\n";	
       $tah=-3.1;       
       
       $address= $opt_a;
       if( $address =~ /^(.*)\/$/ ){
           $address=$1; 
       }
       $tig_info=retrieve("$address"."/div/tig_info") or die;  
       $Ori_Dis =input_non_seq($non_seq,$Ori_Dis);          
}else{
   print "Usage: $0\n";
   print "-o  The file named orientdistinfo_c(some number), outputted by read_sam.pl\n";
   print "-w  Minimum number of mate pair/paired-end links between two contigs (default -w $minlink".")\n";
   print "-L  Minimum length of contigs to be used in scaffold assembly (default -L $minlength".")\n";
   print "-h  High coverage contigs (above mean_coverage + h * std_coverage) are not considered in the scaffold assembly mainly to exclude reads from repetitive regions (default -h $high_cov".")\n";
   print "-n  Tab delimited file containing marker information (from other sources such as shared synteny, physical map, etc.)\n";
   print "-t  Trust in marker information relative to short read data (default -t $strength".")\n";

   print "-a  The directory that you used for output of the formating script (format_base.pl or format_col.pl)\n\n";   
   
   print "example: $0 -o orientdistinfo_c(some number) -a mydir_sout\n\n";
   
   exit;
} 

#--------------------------------    
my $date = `date`; chomp($date); print "start $date\n";	

my $paakon;

if ( $tah != -3.1 ){ 
	if ( $non_seq ){ 
		print "\nRunning: $0\n-h $high_cov\n-L $minlength\n-w $minlink\n-n $non_seq\n-t $strength\n-a $address\n\n";	
	}else{
		print "\nRunning: $0\n-h $high_cov\n-L $minlength\n-w $minlink\n-a $address\n\n";		
	}
	my $rem_cov; my $cut1=0; my $cut2=0;
	for (my $t=0; $t<1; $t++){
		my $mean=0; my $dus=0; my $N=0;
		foreach my $tig (sort { $tig_info->{$b}{'cov'} <=> $tig_info->{$a}{'cov'} } keys %{$tig_info}){  
			if($tig_info->{$tig}{'length'}>150){
				$mean+=$tig_info->{$tig}{'cov'}; $dus+=$tig_info->{$tig}{'cov'}**2; $N++;
			}    
		}
		$mean=$mean/$N; $dus=$dus/$N; my $std=($dus-$mean**2)**.5; 
		my $cut=$mean+$high_cov*$std; $cut1+=(1-$t)*$cut; $cut2+=$t*$cut;  
		#print "\nmean $mean std $std cut1 $cut1 cut2 $cut2 n $N\n";
		$mean=.01* int(100*$mean); $std=.01* int(100*$std); $cut1=.01* int(100*$cut1);
		print "\naverage coverage: $mean       standard deviation: $std       coverage cutoff: $cut1\n";
		print "coverage cutoff $cut1\n" if ($verbs);
		foreach my $tig (sort { $tig_info->{$b}{'cov'} <=> $tig_info->{$a}{'cov'} } keys %{$tig_info}){  
			  if($tig_info->{$tig}{'cov'}>$cut){
					$paakon->{$tig}{'highcov'}=$tig_info->{$tig}{'cov'};  
					if($tig_info->{$tig}{'length'}>$minlength){
						 $rem_cov++; print "removed coverage  $paakon->{$tig}{'highcov'}   -   length  $tig_info->{$tig}{'length'} , contig number $tig\n";
					}
					delete $tig_info->{$tig};  
			  }
		}
		foreach my $tig (keys %$tig_info){
			if (!defined $tig_info->{$tig}{'cov'}){    
			   delete $tig_info->{$tig}; print "eee\n";
			}
		}
	}
	print "\nRemoved $rem_cov contigs from scaffold assembly because of their high coverage\n"; 
}else{
	print "\nRunning: $0\n-L $minlength\n-w $minlink\n-n $non_seq\n-t $strength\n-a $address\n\n";	
}

if ( $tah != -3.1 ){
	$tig_info=retrieve("$address"."/div/tig_info_cov") or die;   
}else{
	$tig_info=retrieve("$address"."/div/tig_info") or die;   
}

my $seq_l=$tig_info->{1}{'seq_l_d'};  

if ( $seq_l == 0 ){
	 $seq_l = 50;
}
#--------------------------------  
my $problem_str=1; my $cyc=0; my $tnum_stretch=0;  my $R_fanar; my $R_garden; my $R_tree; my $ign_mem; 
#--------------------------------     
my ($du_avej,$du_nume,$du_max)=get_ave($Ori_Dis,$paakon,$ign_mem);   
if ( $tah != -3.1 ){
	if ( $non_seq ){ 
		print "Average number of links between two contigs using minlength $minlength, minlink 2 and marker trust $strength is $du_avej , max number of edges: $du_max , number of edges: $du_nume\n";
	}else{
		print "Average number of links between two contigs using minlength $minlength and minlink 2 is $du_avej , max number of edges: $du_max , number of edges: $du_nume\n";	
	}
}else{
	print "Average number of links between two contigs using minlength $minlength, minlink 2 and marker trust $strength is $du_avej , max number of edges: $du_max , number of edges: $du_nume\n";
}
#--------------------------------  
my $Z_combthre=6;	
while($problem_str==1){
    $problem_str=0; my $start_ori;
    $cyc++; my $date = `date`; chomp($date); 
    print "\nStarting cycle $cyc of orientation assignment  $date\n"; 
    my ($j,$avej,$nume)=gen_JforIsing($Ori_Dis,$paakon,$ign_mem);   
	if ( $tah != -3.1 ){
	    if ( $non_seq ){ 
		    print "Average number of links between two contigs using minlength $minlength, minlink $minlink and marker trust $strength is $avej , number of edges: $nume\n";
		}else{
	    	print "Average number of links between two contigs using minlength $minlength and minlink $minlink is $avej , number of edges: $nume\n";		
		}
	}else{
	    print "Average number of links between two contigs using minlength $minlength, minlink $minlink and marker trust $strength is $avej , number of edges: $nume\n";
	}
    my $garden=find_treeforopt($j,$cyc); 
	
	my $neighbor; 
	
	foreach my $m (keys %$j){  
		 $neighbor->{$m}= scalar(keys %{$j->{$m}});
	}	
	
    foreach my $tree (sort { $a <=> $b } keys %{$garden}){  
          my $braking=1; my $floor=0; my $articulate; my $not_articulate;  my $Seq;  
          #print "tree $tree\n"; 
          while($braking==1) {
				  if ( $floor % 11 == 10 ){
					  print "f $floor "; 
				  }          
                $braking=0;    my $n_sub=0; 
                foreach my $sub (sort keys %{$garden->{$tree}{'cycle'}{$floor}} ){ 
                     $garden->{$tree}{'cycle'}{$floor}{$sub}{'final_block'}=1;
					 
					 my $gg=gen_Jforblockfrustration($Ori_Dis,$paakon,$ign_mem,$garden,$tree,$floor,$sub,$articulate,$cyc);    
					 my $Barray=assigntemporaryorient($gg,$cyc,$tree,$floor,$sub);    
					 my ($frust,$e_f)=test_frustration($Barray,$gg);   
					 if ($frust == 0) {                
					        my $du=scalar(@{$garden->{$tree}{'cycle'}{$floor}{$sub}{'member'}});
                            #print "\ntree $tree floor $floor sub $sub no search for artic since frust is zero size $du\n";                     
                     }else{                                         
							 if ( (scalar(@{$garden->{$tree}{'cycle'}{$floor}{$sub}{'member'}}) > (1.8*$Z_combthre)) && (scalar(@{$garden->{$tree}{'cycle'}{$floor}{$sub}{'member'}}) < 50000 ) ){
									 FORTS: foreach my $m (@{$garden->{$tree}{'cycle'}{$floor}{$sub}{'member'}}){ 
											if( $neighbor->{$m}>2 ){
												   if( (! defined $articulate->{$m}) && (!defined $not_articulate->{$m}) ){
														#print "\ntree $tree floor $floor sub $sub consider point $m\n";
														my @r_ne;
														my ($g,$aa,$du)=gen_JforIsing($Ori_Dis,$paakon,$ign_mem);    
														delete $g->{$m};
														foreach my $n (keys %{$g}){    
															 if (!defined $garden->{$tree}{'cycle'}{$floor}{$sub}{'o'}{$n}){
																  delete $g->{$n}; 
															 }
														}   
														foreach my $n (keys %{$g}){    
															 if (defined $g->{$n}{$m} ){
																 delete $g->{$n}{$m};
																 push @r_ne , $n;
															 }
															 foreach my $k (keys %{$g->{$n}}){    
																	if (! defined $garden->{$tree}{'cycle'}{$floor}{$sub}{'o'}{$k}){
																		 if ( defined $articulate->{$n} ){
																			 delete $g->{$n}{$k};  #print "\n$n is articulate and connected to $k which is not here\n\n";
																		 }else{
																			  print "$n $k tree $tree floor $floor $sub\n"; exit;
																		 }
																	}
															 }
														}  	
														my $approved=1; my ($basket,$Barray);
														foreach my $k (keys %$g){ 
															   if (defined $g->{$k}){
																  $basket++; push @{$Barray->{$basket}{'member'}}, $k;   
																  $Barray->{$basket}{'o'}{$k}=1; 
																  my $Extension=1;
																  while ($Extension==1){
																	 $Extension=0; 
																	 foreach my $n (@{$Barray->{$basket}{'member'}}){ 
																		  foreach my $s (keys %{$g->{$n}}){  
																			  if (!defined $Barray->{$basket}{'o'}{$s}){
																				   $Extension=1;
																				   push @{$Barray->{$basket}{'member'}}, $s;  
																				   $Barray->{$basket}{'o'}{$s}=1; 
																				   delete $g->{$n}{$s}; delete $g->{$s}{$n}
																			  }else{
																			  }
																		  }
																	 }
																  } 
																  my $du=scalar(@{$Barray->{$basket}{'member'}}); 
																  my $Lm=scalar(@{$garden->{$tree}{'cycle'}{$floor}{$sub}{'member'}});
																  if($du>($Lm-2)){
																	  $approved=0;
																  }
																  foreach my $n (@{$Barray->{$basket}{'member'}}){ 
																	  delete $g->{$n};
																  }
																  @{$Barray->{$basket}{'member'}}=sort { $a <=> $b } @{$Barray->{$basket}{'member'}};
															   }
														}
														if ( $basket > 1)  {				                              
															  #print "$m is articulate, basket $basket ,  approved $approved\n";                          
															  $articulate->{$m}=1;
															  $garden->{$tree}{'cycle'}{$floor}{$sub}{'final_block'}=0;
															  $garden->{$tree}{'cycle_sharedpoint'}{$floor+1}{$m}{'dad_sub'}=$sub;
															  
															  foreach my $r (@r_ne){
																   $neighbor->{$r}--;  
															  }
															  
															  foreach my $bas ( keys %{$Barray} ){
																   $n_sub++;
																   push @{$garden->{$tree}{'cycle'}{$floor+1}{$n_sub}{'member'}}, $m;
																   $garden->{$tree}{'cycle'}{$floor+1}{$n_sub}{'o'}{$m}=1;
																   foreach my $s (  @{$Barray->{$bas}{'member'}} ){ 												   
																		 push @{$garden->{$tree}{'cycle'}{$floor+1}{$n_sub}{'member'}}, $s;
																		 $garden->{$tree}{'cycle'}{$floor+1}{$n_sub}{'o'}{$s}=1;
																   }				                   
																   push @{$garden->{$tree}{'cycle_sharedpoint'}{$floor+1}{$m}{'new_subs'}} , $n_sub;
																   #print "$floor+1 nsub $n_sub : @{$garden->{$tree}{'cycle'}{$floor+1}{$n_sub}{'member'}}\n";
															  }
															  
															 $braking=1;
															 last FORTS;
														}else{
															  $not_articulate->{$m}=1;
														}
												   }   
											}  # if nei > 2     
									 }   # for this sub
							 }else{
								   #my $du=scalar(@{$garden->{$tree}{'cycle'}{$floor}{$sub}{'member'}}); print "size du $du\n";
								   #print "tree $tree floor $floor sub $sub is small : @{$garden->{$tree}{'cycle'}{$floor}{$sub}{'member'}}\n";
							 }
							 #print "tree $tree floor $floor sub $sub final? $garden->{$tree}{'cycle'}{$floor}{$sub}{'final_block'}\n\n";
							 
					 }	 
					 
                }
                $floor++;
          } # while
          
          my $tot_floor=$floor-1;   undef  $floor;  
		  my $du=scalar(@{$garden->{$tree}{'cycle'}{0}{1}{'member'}});
          print "tree $tree tot_floor $tot_floor , $du members\n"; 
          for (my $floor=$tot_floor ; $floor >= 0 ; $floor--){    
                 foreach my $sub (keys %{$garden->{$tree}{'cycle'}{$floor}} ){ 
                       if ( $garden->{$tree}{'cycle'}{$floor}{$sub}{'final_block'}==1 ){    
                                my $g=gen_Jforblockfrustration($Ori_Dis,$paakon,$ign_mem,$garden,$tree,$floor,$sub,$articulate,$cyc);    
								my $Barray=assigntemporaryorient($g,$cyc,$tree,$floor,$sub);    
								my ($frust,$e_f)=test_frustration($Barray,$g);   
								#print "tree $tree floor $floor sub $sub frust $frust e_f $e_f\n"; 
                                if ($frust == 0) {
                                      foreach my $s (keys %{$Barray->{1}{'o'}}){  
                                             $garden->{$tree}{'cycle'}{$floor}{$sub}{'orient'}{$s}=$Barray->{1}{'o'}{$s};
                                      }                               
                                }else{
                                       #print "\nfind zstructure for tree $tree floor $floor sub $sub frust $frust\n\n"; 
                                       
                                      my ($Z,$big)=find_Zstructure($g); 
                                      print "floor $floor max Z size $big frust $frust\n";
                                      if ($big==0){
                                                my ($mat_old,$mat_cur,$Score,$L_old,$L_cur); 
                                             	for (my $current_lay=1; $current_lay < 2; $current_lay++){ 
                                             	      $L_cur=scalar(keys %{$Z->{$current_lay}{'o'}} );
													  my $dum=2**($L_cur-1);    
													  for (my $j=1; $j <= $L_cur; $j++){  
														   my $num_block=2**($j-1);	   
														   for (my $i=1; $i <= $num_block; $i+=2){  
																for (my $z=(1+($i-1)*(2**($L_cur-$j))); $z <= ($i*(2**($L_cur-$j))); $z++){  
																	  $mat_cur->{$z}{$j}=1;
																	  if ( $j>1 ){
																		  $mat_cur->{$z+$dum}{$j}=1;
																	  }    
																}                       
														   }
													  }
													  $mat_cur->{2**$L_cur}{1}=0;
													  for (my $c=1; $c <= 2**$L_cur; $c++){  
													         my $spin_cur;
													         for (my $m=1; $m <= $L_cur; $m++){
													               if ($mat_cur->{$c}{$m}==1){
													                     $spin_cur->{${$Z->{$current_lay}{'mem'}}[$m-1]}=1;
													               }else{
													                     $spin_cur->{${$Z->{$current_lay}{'mem'}}[$m-1]}=-1;
													               }  
													         }
													         my $inside_c; 
															 foreach my $mem_c (keys %{$Z->{$current_lay}{'o'}} ){    
																	foreach my $mem_cc (keys %{$Z->{$current_lay}{'o'}} ){    
																		   if (defined $g->{$mem_c}{$mem_cc}){
																				 $inside_c-=$g->{$mem_c}{$mem_cc}*$spin_cur->{$mem_c}*$spin_cur->{$mem_cc};
																		   }
																		   
																	}
															 }													              
															 $inside_c*=.5;  
															 my $transfer_o_c;    
															 my $mem_o=${$Z->{0}{'mem'}}[0];
															 foreach my $mem_c (keys %{$Z->{$current_lay}{'o'}} ){    
																 if (defined $g->{$mem_o}{$mem_c}){
																	  $transfer_o_c-=$g->{$mem_o}{$mem_c}*$spin_cur->{$mem_c};
																 }
															 }							
													          $Seq->{$current_lay}{$c}=1;
													          $Score->{$c}=$transfer_o_c+$inside_c;
													  }										  
                                             	}
                                                
                                             	for (my $current_lay=2; $current_lay < scalar(keys %{$Z}); $current_lay++){ 
                                             	      undef $mat_old;
                                             	      foreach my $s ( keys %{$mat_cur} ){  
                                             	            foreach my $k ( keys %{$mat_cur->{$s}} ){  
                                             	                  $mat_old->{$s}{$k}=$mat_cur->{$s}{$k};
                                             	            }
                                             	      }
                                             	      $L_old=$L_cur;
                                             	      $L_cur=scalar(keys %{$Z->{$current_lay}{'o'}} );
                                             	      undef $mat_cur;
													  my $dum=2**($L_cur-1);    
													  for (my $j=1; $j <= $L_cur; $j++){  
														   my $num_block=2**($j-1);	   
														   for (my $i=1; $i <= $num_block; $i+=2){  
																for (my $z=(1+($i-1)*(2**($L_cur-$j))); $z <= ($i*(2**($L_cur-$j))); $z++){  
																	  $mat_cur->{$z}{$j}=1;
																	  if ( $j>1 ){
																		  $mat_cur->{$z+$dum}{$j}=1;
																	  }    
																}                       
														   }
													  }
													  $mat_cur->{2**$L_cur}{1}=0; 
													  my $TScore; 
													  for (my $c=1; $c <= 2**$L_cur; $c++){  
													         my $spin_cur;
													         for (my $m=1; $m <= $L_cur; $m++){
													               if ($mat_cur->{$c}{$m}==1){
													                     $spin_cur->{${$Z->{$current_lay}{'mem'}}[$m-1]}=1;
													               }else{
													                     $spin_cur->{${$Z->{$current_lay}{'mem'}}[$m-1]}=-1;
													               }  
													         }
													         my $inside_c; 
															 foreach my $mem_c (keys %{$Z->{$current_lay}{'o'}} ){    
																	foreach my $mem_cc (keys %{$Z->{$current_lay}{'o'}} ){    
																		   if (defined $g->{$mem_c}{$mem_cc}){
																				 $inside_c-=$g->{$mem_c}{$mem_cc}*$spin_cur->{$mem_c}*$spin_cur->{$mem_cc};
																		   }
																		   
																	}
															 }													              
															 $inside_c*=.5;  	
															 my $TempScore; 
													          for (my $o=1; $o <= 2**$L_old; $o++){
																	 my $spin_old;
																	 for (my $m=1; $m <= $L_old; $m++){
																		   if ($mat_old->{$o}{$m}==1){
																				 $spin_old->{${$Z->{$current_lay-1}{'mem'}}[$m-1]}=1;
																		   }else{
																				 $spin_old->{${$Z->{$current_lay-1}{'mem'}}[$m-1]}=-1;
																		   }  
																	 }	
																	 my $transfer_o_c;
													                 foreach my $mem_o (keys %{$Z->{$current_lay-1}{'o'}} ){    
													                        foreach my $mem_c (keys %{$Z->{$current_lay}{'o'}} ){    
                     															   if (defined $g->{$mem_o}{$mem_c}){
                     															         $transfer_o_c-=$g->{$mem_o}{$mem_c}*$spin_old->{$mem_o}*$spin_cur->{$mem_c};
                     															   }
													                               
													                        }
													                 }
													                 $TempScore->{$o}=$Score->{$o} + $transfer_o_c;
													          }
													          foreach my $z (sort { $TempScore->{$a} <=> $TempScore->{$b} } keys %{$TempScore}){  
													                 $TScore->{$c}=$TempScore->{$z}+$inside_c;  
													                 $Seq->{$current_lay}{$c}=$z;
													                 last;
													          }													  
													  }
													  undef $Score;
													  foreach my $c (keys %{$TScore}){
														  $Score->{$c}=$TScore->{$c};
													  }   					
                                             	}
                                       
		                                        my $c;
											    foreach my $z (sort { $Score->{$a} <=> $Score->{$b} } keys %{$Score}){  
													  $c=$z;													  
													  my $L_cur=scalar(keys %{$Z->{scalar(keys %{$Z})-1}{'o'}} );
													  my $spin_cur;
 													  for (my $m=1; $m <= $L_cur; $m++){
															  my $block_size=2**($L_cur-$m);
															  my $which_block=$c/$block_size;
															  if ( $which_block==int($which_block) ){
															  
															  }else{
																   $which_block=int($which_block)+1;
															  }
															   if ( ($which_block % 2)==1 ){
																	$garden->{$tree}{'cycle'}{$floor}{$sub}{'orient'}{${$Z->{scalar(keys %{$Z})-1}{'mem'}}[$m-1]}=1;
															  }else{
																	$garden->{$tree}{'cycle'}{$floor}{$sub}{'orient'}{${$Z->{scalar(keys %{$Z})-1}{'mem'}}[$m-1]}=-1;
															  }
													  }
													  last;
											    }                                       
                                                
                                                for (my $current_lay=scalar(keys %{$Z})-2;  0 < $current_lay ; $current_lay--){ 
													  $c=$Seq->{$current_lay+1}{$c};  												  
													  my $L_cur=scalar(keys %{$Z->{$current_lay}{'o'}});
													  my $spin_cur;
													  for (my $m=1; $m <= $L_cur; $m++){
															  my $block_size=2**($L_cur-$m);
															  my $which_block=$c/$block_size;
															  if ( $which_block==int($which_block) ){
															  
															  }else{
																   $which_block=int($which_block)+1;
															  }
															   if ( ($which_block % 2)==1 ){
																	$garden->{$tree}{'cycle'}{$floor}{$sub}{'orient'}{${$Z->{$current_lay}{'mem'}}[$m-1]}=1;
															  }else{
																	$garden->{$tree}{'cycle'}{$floor}{$sub}{'orient'}{${$Z->{$current_lay}{'mem'}}[$m-1]}=-1;
															  }
													  }    
                                                }	    
                                                $garden->{$tree}{'cycle'}{$floor}{$sub}{'orient'}{${$Z->{0}{'mem'}}[0]}=1;
                                                 #draw_graph_1($g,$garden,$cyc,$tree,$floor,$sub);
												  my @B; my $minE; my $E;
												  foreach my $m (keys %$g){  
														$B[$m]=0;
														foreach my $n (keys %{$g->{$m}}){ 
															$B[$m]+= $g->{$m}{$n}*$garden->{$tree}{'cycle'}{$floor}{$sub}{'orient'}{$n};  
															$minE-=abs($g->{$m}{$n});
														}
														$E+=- $B[$m]*$garden->{$tree}{'cycle'}{$floor}{$sub}{'orient'}{$m};  
												  }
												  $minE*=.5; $E*=.5; 
												  #print "\nZstructure ::: tree $tree floor $floor sub $sub frust $frust transferM minE $minE E $E\n\n" ;
                                      }else{                 
                                                $garden=simulated_annealing($g,$garden,$tree,$floor,$sub);
 
                                                    #draw_graph_1($g,$garden,$cyc,$tree,$floor,$sub);
                                                	  my @B; my $minE; my $E;
													  foreach my $m (keys %$g){  
															$B[$m]=0;
															foreach my $n (keys %{$g->{$m}}){ 
																$B[$m]+= $g->{$m}{$n}*$garden->{$tree}{'cycle'}{$floor}{$sub}{'orient'}{$n};  
																$minE-=abs($g->{$m}{$n});
															}
															$E+=- $B[$m]*$garden->{$tree}{'cycle'}{$floor}{$sub}{'orient'}{$m};  
													  }
													  $minE*=.5; $E*=.5; 
													  #print "tree $tree floor $floor sub $sub frust $frust minE $minE E $E\n";    
													  #$minE*=(20/$avej); $E*=(20/$avej); $minE=int($minE); $E=int($E); print "minE $minE final E $E\n"; 
                                      }      
                                }
                       }
                 }  
          }
            #print "putting things back together\n";
            if ( $tot_floor > 0 ) {
				for (my $floor=$tot_floor-1 ; $floor >= 0 ; $floor--){            # putting things back together
						  foreach my $art (keys %{$garden->{$tree}{'cycle_sharedpoint'}{$floor+1}} ) {
							   my $dad_sub=$garden->{$tree}{'cycle_sharedpoint'}{$floor+1}{$art}{'dad_sub'};
							   #print "floor $floor art $art dad_sub $dad_sub\n";
							   foreach my $n_sub (  @{$garden->{$tree}{'cycle_sharedpoint'}{$floor+1}{$art}{'new_subs'}} ){ 	
									  if ( $garden->{$tree}{'cycle'}{$floor+1}{$n_sub}{'orient'}{$art}==-1 ) {
											foreach my $s (keys %{$garden->{$tree}{'cycle'}{$floor+1}{$n_sub}{'orient'}} ){
												  $garden->{$tree}{'cycle'}{$floor+1}{$n_sub}{'orient'}{$s}*=-1;
											}
									  }elsif ( $garden->{$tree}{'cycle'}{$floor+1}{$n_sub}{'orient'}{$art}==1 ) {
									  
									  }else{
										  print "error $garden->{$tree}{'cycle'}{$floor+1}{$n_sub}{'orient'}{$art}\n"; exit;
									  }
									  foreach my $s (keys %{$garden->{$tree}{'cycle'}{$floor+1}{$n_sub}{'orient'}} ){
											$garden->{$tree}{'cycle'}{$floor}{$dad_sub}{'orient'}{$s}=$garden->{$tree}{'cycle'}{$floor+1}{$n_sub}{'orient'}{$s};
									  }
							   }
						  }
				 }
		     }else{
		     
		     }
			 foreach my $s (keys %{$garden->{$tree}{'cycle'}{0}{1}{'o'}} ){
			      if ( ! defined $garden->{$tree}{'cycle'}{0}{1}{'orient'}{$s}){
			          print "s not defined $s\n";  exit;
			      }
			 }
			 foreach my $s (keys %{$garden->{$tree}{'cycle'}{0}{1}{'orient'}} ){
						 $start_ori->{$s}=$garden->{$tree}{'cycle'}{0}{1}{'orient'}{$s};
   			 }          
			 #my $g=gen_Jforblockfrustration($Ori_Dis,$paakon,$ign_mem,$garden,$tree,0,1,$articulate,$cyc);    
			 #draw_graph_1($g,$garden,$cyc,$tree,0,1);
    }

	  my @B; my $minE; my $E;
	  foreach my $m (keys %$j){  
			$B[$m]=0;
			foreach my $n (keys %{$j->{$m}}){ 
				$B[$m]+= $j->{$m}{$n}*$start_ori->{$n};  
				$minE-=abs($j->{$m}{$n});
			}
			$E+=- $B[$m]*$start_ori->{$m};  
	  }
	  $minE*=.5; $E*=.5; 
	  print "\nminE $minE final E $E\n";  # $minE*=(20/$avej); $E*=(20/$avej); $minE=int($minE); $E=int($E); print "minE $minE final E $E\n"; 
    
    #-----------------------------------------   #after start_ori
	my ($J,$rededge)=gen_Jforspring($Ori_Dis,$start_ori);	
	
	$garden=find_treeforspring($J,$cyc,$rededge,$start_ori);      	
	
	my $dfanar;  ($J,$dfanar)=gen_fanar($Ori_Dis,$start_ori);
	#-----------------------------------------
	foreach my $tree (sort { $a <=> $b } keys %{$garden}){ 
       $garden=linear_eq($garden,$tree,$dfanar,$J);
        my $nnum_stretch=0;  my $kesh=0; my $candide_kesh;
		foreach my $m (keys %$dfanar){  
			if(defined ($garden->{$tree}{'cycle'}{0}{1}{'position'}{$m})){
				 foreach my $n (keys %{$dfanar->{$m}}){     
					if(defined ($garden->{$tree}{'cycle'}{0}{1}{'position'}{$n})){		
						 if($m<$n){             
							   my $S=$garden->{$tree}{'cycle'}{0}{1}{'position'}{$n}-$garden->{$tree}{'cycle'}{0}{1}{'position'}{$m}-$dfanar->{$m}{$n};  
							   if(abs($S)>$stretch_thr){
									$problem_str=1; 
									$kesh++;
									if(abs($S)>(1.5*$stretch_thr) ){
										if ($tig_info->{$m}{'length'} < $tig_info->{$n}{'length'}){
											$candide_kesh->{$m}=$tig_info->{$m}{'length'};
										}else{
											$candide_kesh->{$n}=$tig_info->{$n}{'length'};
										} 
									}	
							   }else{
 							   }
						  }
					}	  
				 }
			} 
		}
		my $cand_size=4;
		if( $kesh > 0 ){
			if( scalar(keys %{$candide_kesh}) > $cand_size ){	
				my @cand=sort { $candide_kesh->{$a} <=> $candide_kesh->{$b} } keys %{$candide_kesh};
				for (my $i=0; $i < $cand_size; $i++){  
					$paakon->{$cand[$i]}{'stretch'}++;
					$nnum_stretch+=1; 		
					#print "removed L $candide_kesh->{$cand[$i]}\n";
				}
			}elsif( scalar(keys %{$candide_kesh}) > 0 ){	
				foreach my $tig (keys %{$candide_kesh}){  
					$paakon->{$tig}{'stretch'}++;
					$nnum_stretch+=1; 
				}			
			}
		}

        if( ($kesh > 0) && ($nnum_stretch==0) ){
            print "kesh $kesh nnum_stretch $nnum_stretch\n";
			foreach my $m (keys %$dfanar){  
				if(defined ($garden->{$tree}{'cycle'}{0}{1}{'position'}{$m})){
					 foreach my $n (keys %{$dfanar->{$m}}){     
						if(defined ($garden->{$tree}{'cycle'}{0}{1}{'position'}{$n})){		
							 if($m<$n){             
								   my $S=$garden->{$tree}{'cycle'}{0}{1}{'position'}{$n}-$garden->{$tree}{'cycle'}{0}{1}{'position'}{$m}-$dfanar->{$m}{$n};  
								   if(abs($S)>$stretch_thr){
										#print "tree $tree  stretch $m $n\n";
										if ($tig_info->{$m}{'length'} < $tig_info->{$n}{'length'}){
											$candide_kesh->{$m}=$tig_info->{$m}{'length'};
										}else{
											$candide_kesh->{$n}=$tig_info->{$n}{'length'};
										} 
										$problem_str=1;  
								   }else{	
								   }
							  }
						}	  
					 }
				} 
			}
			
			if( scalar(keys %{$candide_kesh}) > $cand_size ){	
				my @cand=sort { $candide_kesh->{$a} <=> $candide_kesh->{$b} } keys %{$candide_kesh};
				for (my $i=0; $i < $cand_size; $i++){  
					$paakon->{$cand[$i]}{'stretch'}++;
					$nnum_stretch+=1; 		
					#print "-> removed L $candide_kesh->{$cand[$i]}\n";
				}
			}elsif( scalar(keys %{$candide_kesh}) > 0 ){	
				foreach my $tig (keys %{$candide_kesh}){  
					$paakon->{$tig}{'stretch'}++;
					$nnum_stretch+=1; 
				}			
			}else{
				print "problem cutting strings\n"; exit;
			}
			
            print "-> kesh $kesh nnum_stretch $nnum_stretch\n";
        }elsif( ($kesh > 0) && ($nnum_stretch>0)){
            print "kesh $kesh nnum_stretch $nnum_stretch\n";
        }
		$tnum_stretch+=$nnum_stretch;   

		#print "Number of stretched springs $tnum_stretch\n";	
		if($nnum_stretch==0){
			 $R_tree++;	     #print " R_tree $R_tree\n";	
			 foreach my $m (keys %{$garden->{$tree}{'cycle'}{0}{1}{'position'}}){
				 $ign_mem->{$m}=1;
				 push @{$R_garden->{$R_tree}{'cycle'}{0}{1}{'member'}}, $m;  
				 $R_garden->{$R_tree}{'cycle'}{0}{1}{'position'}{$m}=$garden->{$tree}{'cycle'}{0}{1}{'position'}{$m};
				 $R_garden->{$R_tree}{'cycle'}{0}{1}{'orientation'}{$m}=$start_ori->{$m};
				 foreach my $n (keys %{$dfanar->{$m}}){
					 if(defined $garden->{$tree}{'cycle'}{0}{1}{'position'}{$n}){ 
						  $R_fanar->{$m}{$n}=$dfanar->{$m}{$n}; $R_fanar->{$n}{$m}=$dfanar->{$m}{$n};
					 }else{
	
					 }
				 }
				 
			 }
		}else{
		    my $du=scalar(keys %{$garden->{$tree}{'cycle'}{0}{1}{'position'}});
		    print "tree $tree members $du number of stretched springs $nnum_stretch    tot number so far $tnum_stretch\n\n";	
		}	
	}		
}

 
store \%{$R_fanar}, "$address"."/div/fanar"."_h$high_cov"."_L$minlength"."_w$minlink";   
store \%{$R_garden}, "$address"."/div/garden"."_h$high_cov"."_L$minlength"."_w$minlink";  
store \%{$paakon}, "$address"."/div/paakon"."_h$high_cov"."_L$minlength"."_w$minlink";     


$date = `date`; chomp($date); print "Finished orientation assignment\n"; 

undef $paakon; undef $R_fanar; undef $R_garden; undef $R_tree; undef $ign_mem;     
#--------------------------------------------------------

my $garden; 

if ( $tah != -3.1 ){
	$tig_info=retrieve("$address"."/div/tig_info_cov") or die "$!";  
	$Ori_Dis=retrieve("$address"."/orientdistinfo_c$tah") or die "$!"; 
}else{
	$tig_info=retrieve("$address"."/div/tig_info") or die "$!";  
	undef $Ori_Dis;
}
if( $non_seq ){
	$Ori_Dis =input_non_seq($non_seq,$Ori_Dis);   
}
$garden=retrieve("$address"."/div/garden"."_h$high_cov"."_L$minlength"."_w$minlink" ) or die "$!"; 
$paakon=retrieve("$address"."/div/paakon"."_h$high_cov"."_L$minlength"."_w$minlink" ) or die "$!";

my ($max_walk,$min_walk,$do_walk)=($seq_l,14,1);

my ($combin_thr,$thr_bad_den)=(8,30);

my ($talkw,$talkrc,$talksc)=(0,0,0);


$max_walk = $opt_u if ($opt_u);
$min_walk = $opt_s if ($opt_s);
$do_walk = 0 if ($opt_d);
$combin_thr = $opt_k if ($opt_k);

print "\nminlength $minlength stretch_thr $stretch_thr  min_walk $min_walk thr_bad_den $thr_bad_den minlink $minlink \n\n";   

open (SCFD, ">$address"."/scaffolds"."_h$high_cov"."_L$minlength"."_w$minlink" .".fasta") or die "$!";

my ($walk,$Pbaagh); my $daarbast=0; 
foreach my $tree (sort { $a <=> $b } keys %{$garden}){  
     my $sub=1; my $cycle=0; $|=1; ###clear buffer

     foreach my $m (@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'member'}}){  
         $garden->{$tree}{'cycle'}{$cycle}{$sub}{'position'}{$m}=int($garden->{$tree}{'cycle'}{$cycle}{$sub}{'position'}{$m});
     }
     my $ontop;
	 ($garden,$ontop)=find_den($garden,$tree,$cycle,$sub);	  
 
	 if(($garden->{$tree}{'cycle'}{$cycle}{$sub}{'Bad_den'}<=$thr_bad_den)&&($do_walk==1)){  
			 #print "good tree $tree find walk\n";
			 foreach my $m (@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'member'}}){  
				foreach my $n (@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'member'}}){  
				   if($m<$n){
						 if((($garden->{$tree}{'cycle'}{$cycle}{$sub}{'Rstart'}{$n}-$garden->{$tree}{'cycle'}{$cycle}{$sub}{'Rend'}{$m})<15)&&(($garden->{$tree}{'cycle'}{$cycle}{$sub}{'Rstart'}{$n}-$garden->{$tree}{'cycle'}{$cycle}{$sub}{'Rend'}{$m})>-150)){   
							      ($walk,$garden)=findwalk($walk,$m,$n,$garden,$tree,$cycle,$sub);
						 }elsif((($garden->{$tree}{'cycle'}{$cycle}{$sub}{'Rstart'}{$m}-$garden->{$tree}{'cycle'}{$cycle}{$sub}{'Rend'}{$n})<15)&&(($garden->{$tree}{'cycle'}{$cycle}{$sub}{'Rstart'}{$m}-$garden->{$tree}{'cycle'}{$cycle}{$sub}{'Rend'}{$n})>-150)){   
							      ($walk,$garden)=findwalk($walk,$n,$m,$garden,$tree,$cycle,$sub);					      
						 }  
				   } 
				}  
			 }     
			 my $problem_walk; my @walk_slist;
			 ($walk,$garden,$problem_walk,@walk_slist)=processwalk($walk,$garden,$tree,$cycle,$sub);
			 #print "end find walk problem_walk $problem_walk\n";  
			 if($problem_walk==1){
                   push @{$Pbaagh->{$cycle}{$tree}{1}},@walk_slist;
                   $garden->{$tree}{'cycle'}{$cycle}{$sub}{'status'}='problem_walk';
                   delete $walk->{$tree}{'cycle'}{$cycle}{$sub}; 
                   delete $garden->{$tree}{'cycle'}{$cycle}{$sub}{'village'};
			 }else{
			       buildScaffold($garden,$tree,$cycle,$sub,$walk);  
			       $garden->{$tree}{'cycle'}{$cycle}{$sub}{'status'}='good';
			 }
     }
     if($garden->{$tree}{'cycle'}{$cycle}{$sub}{'Bad_den'}>$thr_bad_den){   
 			$garden->{$tree}{'cycle'}{$cycle}{$sub}{'status'}='problem_density';
			$garden=findhump($garden,$tree,$cycle,$sub);
			
			my ($cut_edge, $still_over)=color_nodes($garden,$tree,$cycle,$sub,$ontop);
			
# 			my ($cut_edge1, $still_over1)=color_nodes($garden,$tree,$cycle,$sub,$ontop);
# 			my $cut_edge; my @n1=keys %{$cut_edge1};
# 			if(scalar(@n1)<4){
# 				foreach my $s1 (keys %{$cut_edge1}){
# 					   $cut_edge->{$s1}=1;		   
# 				}
# 			}else{
# 			    my ($cut_edge2, $still_over2)=color_nodes($garden,$tree,$cycle,$sub,$ontop); my @n2=keys %{$cut_edge2};  print "cut1: @n1  cut2 @n2\n";
# 				foreach my $s1 (keys %{$cut_edge1}){
# 				   if(defined $cut_edge2->{$s1}){
# 					   $cut_edge->{$s1}=1;
# 				   }else{			   
# 				   }
# 				}						
# 				if(scalar(keys %{$cut_edge})==0){
# 					if((scalar(@n1)>0)  && (scalar(@n1)<=scalar(@n2))){
# 						foreach my $s1 (keys %{$cut_edge1}){
# 							   $cut_edge->{$s1}=1;
# 						}				    
# 					}elsif((scalar(@n2)>0)  && (scalar(@n2)<=scalar(@n1))){
# 						foreach my $s1 (keys %{$cut_edge2}){
# 							   $cut_edge->{$s1}=1;
# 						}				  
# 					}else{
# 					
# 					}
# 				}			
# 			}

	        if(scalar(keys %{$cut_edge})==0){	
				 foreach my $m (keys %{$ontop}){    
					 foreach my $n (keys %{$ontop->{$m}}){
						  if($ontop->{$m}{$n}>300){ 
                                 $cut_edge->{$m}=1; $cut_edge->{$n}=1; 
						  }elsif($ontop->{$m}{$n}>100){ 
							   if(($tig_info->{$m}{'length'}<300) || ($tig_info->{$n}{'length'}<300)){
		       			           $cut_edge->{$m}=1; $cut_edge->{$n}=1;
							   }
						  }
					 }   
				  }
	        }
			
			$garden=makeSlist($cut_edge,$garden,$tree,$cycle,$sub); 			
	# 		foreach my $m (keys %{$still_over}){
	# 		     foreach my $n (keys %{$still_over->{$m}}){  
	# 		           if ($m<$n){
	# 		                print "still overlap $m and $n. $ontop->{$m}{$n}\n";   
	# 		           }
	# 		     }
	# 		}
			my $hal_s=0;	
			KOHAN: foreach my $camel (keys %{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}}){
				if (defined @{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Slist'}}){
					 my $du=scalar(@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Slist'}});
					 if($du<$combin_thr){
						  push @{$Pbaagh->{$cycle}{$tree}{1}},@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Slist'}};	
						  #print "for slist chosed camel $camel num_snodes $du rn1 @{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Slist'}}\n";
						  
						  $garden->{$tree}{'cycle'}{$cycle}{$sub}{'colors'}=$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'summit'};
						  
						  $hal_s=1;
						  last KOHAN;
					 }
				}else{
					#print "\nwhy not defined tree $tree hump $camel\n";
				}
			}
			if($hal_s==0){
				 GHOLE:foreach my $camel (keys %{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}}){
					if (defined @{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Sumit_Slist'}}){
						 my $du=scalar(@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Sumit_Slist'}});
						 if($du<$combin_thr){
							  push @{$Pbaagh->{$cycle}{$tree}{1}},@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Sumit_Slist'}};	 
							  #print "for slist chosed summit camel $camel num_snodes $du rn2 @{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Sumit_Slist'}}\n";
							  
							  $garden->{$tree}{'cycle'}{$cycle}{$sub}{'colors'}=$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'summit'};
							  
							  $hal_s=1;
							  last GHOLE;
						 }
					} 
				 }	
			}
			if($hal_s==0){
				 AVALGHOL:foreach my $camel (keys %{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}}){
					if (defined @{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Sumit_Slist'}}){
						  push @{$Pbaagh->{$cycle}{$tree}{1}},@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Sumit_Slist'}};
						  my $du=scalar(@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Sumit_Slist'}}); #print "for slist chosed summit camel $camel num_snodes $du rn3 @{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Sumit_Slist'}}\n";
						  
						  $garden->{$tree}{'cycle'}{$cycle}{$sub}{'colors'}=$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'summit'};
						  
						  $hal_s=1;
						  last AVALGHOL;
					}     
				 }	
			}
			if($hal_s==0){
				 AVALKOH:foreach my $camel (keys %{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}}){
					if (defined @{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Slist'}}){
						  push @{$Pbaagh->{$cycle}{$tree}{1}},@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Slist'}};
						  my $du=scalar(@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Slist'}}); #print "for slist chosed camel $camel num_snodes $du rn4 @{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Slist'}}\n";
						  
						  $garden->{$tree}{'cycle'}{$cycle}{$sub}{'colors'}=$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'summit'};
						  
						  $hal_s=1;
						  last AVALKOH;
					}
				 }	
			}		
			if($hal_s==0){
				#print "\n\n$tree che khabareee??\n\n";
				$garden->{$tree}{'cycle'}{$cycle}{$sub}{'colors'}=1.5;
				foreach my $s1 (keys %{$cut_edge}){
					  push @{$Pbaagh->{$cycle}{$tree}{$sub}}, $s1;	
				}	
				#print "for slist chosed rn @{$Pbaagh->{$cycle}{$tree}{$sub}}\n";					
			}			 
     }      
}

my $cycle=0;  
while(defined $Pbaagh->{$cycle}){
    my $Old_cycle=$cycle; $cycle++;  print "\n\nstarting cycle $cycle of density resolving\n\n"; $|=1; ###clear buffer
	foreach my $tree (sort { $a <=> $b } keys %{$Pbaagh->{$Old_cycle}}){  
	        my $nsub=0;  
	        foreach my $Osub (keys %{$Pbaagh->{$Old_cycle}{$tree}}){
	             my $Barray; my $Solved_den=0; my $Lm=scalar(@{$garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'member'}});
	             my @all_Snodes=@{$Pbaagh->{$Old_cycle}{$tree}{$Osub}}; my $Ln=scalar(@all_Snodes); my $orL=$Ln; my $rcounter=0; #print "\nstarting cycle $cycle tree $tree osub $Osub num_remo $Ln : @{$Pbaagh->{$Old_cycle}{$tree}{$Osub}}\n"; 
	             my $wemade=0;
	             if($Ln>=$combin_thr){
	                 @all_Snodes=@all_Snodes[0..($combin_thr-2)]; $Ln=scalar(@all_Snodes); $wemade=1;
	                 #print "\nwe made new rn num_remo $Ln : @all_Snodes\n"; 
	             }
	             print "t $tree osub $Osub n_label $garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'colors'} ori_n_rem $orL n_rem $Ln  "; 	             
				 my $mat_b;  my $only_more_than4s=0;
				 if((scalar(@all_Snodes)>2) && ($garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'colors'}>2)){
				     # $only_more_than4s=1;  #print "\n\nonly more than 2 is on:  colors: $garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'colors'}\n\n";
				 }
                 if( (($garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'colors'}<5) || ($orL<40)) && ($garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'status'} ne 'problem_walk') ){ 
					  for (my $j=1; $j <= $Ln; $j++){  
						   my $num_block=2**($j-1);	   
						   for (my $i=1; $i <= $num_block; $i+=2){  
								for (my $z=(1+($i-1)*(2**($Ln-$j))); $z <= ($i*(2**($Ln-$j))); $z++){  
									  $mat_b->{$z}{$j}=1;
								}                       
						   }
					  }
					  for (my $i=1; $i <= 2**($Ln-1); $i++){  
							my @S0; my @S1;
							for (my $j=1; $j <= $Ln; $j++){ 
								   if ( defined $mat_b->{$i}{$j} ){
									   push @S1,$all_Snodes[$j-1];               
								   }else{
									   push @S0,$all_Snodes[$j-1];
								   } 
							}
							if(($only_more_than4s==0) || (scalar(@S0)>2) ){
								my $g=generate_g($garden,$tree,$Old_cycle,$Osub,@S0); 
								foreach my $m (keys %{$g}){    
									 if (!defined $garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'position'}{$m}){
										  delete $g->{$m};  #print "dddd\n"; exit;
									 }
								}  
								foreach my $m (keys %{$g}){    
									 foreach my $n (keys %{$g->{$m}}){    
											if (!defined $garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'position'}{$n}){
												#print "$m $n $tree $Old_cycle $Osub\n"; exit;
											}
									 }
								}  
								my $basket=0; $rcounter++; my $approved=1; @{$Barray->{$rcounter}{'removednodes'}}=@S0; $Barray->{$rcounter}{'num_removed'}=scalar(@S0); $Barray->{$rcounter}{'score'}=0; my $duh;
								foreach my $m (keys %$g){ 
									   if (defined $g->{$m}){
										  $basket++; push @{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}, $m;   $Barray->{$rcounter}{'basket'}{$basket}{'orientation'}{$m}=$garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'orientation'}{$m}; $Barray->{$rcounter}{'basket'}{$basket}{'position'}{$m}=$garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'position'}{$m}; 
										  my $Extension=1;
										  while ($Extension==1){
											 $Extension=0; 
											 foreach my $n (@{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}){ 
												  foreach my $s (keys %{$g->{$n}}){  
													  if (!defined $Barray->{$rcounter}{'basket'}{$basket}{'position'}{$s}){
														   $Extension=1;
														   push @{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}, $s;  
														   $Barray->{$rcounter}{'basket'}{$basket}{'orientation'}{$s}=$garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'orientation'}{$s}; $Barray->{$rcounter}{'basket'}{$basket}{'position'}{$s}=$garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'position'}{$s}; 	
														   delete $g->{$n}{$s}; delete $g->{$s}{$n}; $duh->{$s}=1; $duh->{$n}=1;
													  }else{
														 #print "$n $s vay vayyy 1\n";
													  }
												  }
											 }
										  } 
										  my $du=scalar(@{$Barray->{$rcounter}{'basket'}{$basket}{'member'}});  
										  if($du>($Lm-scalar(@S0)-2)){
											  $approved=0;
											  #print "rcounter $rcounter not approved $Lm - @S0 - $du \n" if($talkrc==1);
										  }
										  foreach my $n (@{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}){ 
											  delete $g->{$n};
										  }
										  @{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}=sort { $a <=> $b } @{$Barray->{$rcounter}{'basket'}{$basket}{'member'}};
									   }else{
										  if(!defined $duh->{$m}){
											 #print "na baba\n";
										  }  
									   }
								}	
								if($approved==0){
									delete $Barray->{$rcounter}; $rcounter--; 
								}
							}
							if(($only_more_than4s==0) || (scalar(@S1)>2) ){						
								my $g=generate_g($garden,$tree,$Old_cycle,$Osub,@S1); 
								foreach my $m (keys %{$g}){    
									 if (!defined $garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'position'}{$m}){
										  delete $g->{$m}; #print "dddd33\n"; exit;
									 }
								}  
								foreach my $m (keys %{$g}){    
									 foreach my $n (keys %{$g->{$m}}){    
											if (!defined $garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'position'}{$n}){
												#print "$m $n $tree $Old_cycle $Osub\n"; exit;
											}
									 }
								}  
								my $basket=0; $rcounter++; my $approved=1; @{$Barray->{$rcounter}{'removednodes'}}=@S1; $Barray->{$rcounter}{'num_removed'}=scalar(@S1); $Barray->{$rcounter}{'score'}=0;
								foreach my $m (keys %$g){ 
									   if (defined $g->{$m}){
										  $basket++; push @{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}, $m;   $Barray->{$rcounter}{'basket'}{$basket}{'orientation'}{$m}=$garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'orientation'}{$m}; $Barray->{$rcounter}{'basket'}{$basket}{'position'}{$m}=$garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'position'}{$m}; 
										  my $Extension=1;
										  while ($Extension==1){
											 $Extension=0; 
											 foreach my $n (@{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}){ 
												  foreach my $s (keys %{$g->{$n}}){  
													  if (!defined $Barray->{$rcounter}{'basket'}{$basket}{'position'}{$s}){
														   $Extension=1;
														   push @{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}, $s;  
														   $Barray->{$rcounter}{'basket'}{$basket}{'orientation'}{$s}=$garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'orientation'}{$s}; $Barray->{$rcounter}{'basket'}{$basket}{'position'}{$s}=$garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'position'}{$s}; 	
														   delete $g->{$n}{$s}; delete $g->{$s}{$n}
													  }else{
													   # print "$n $s vay vayyy\n";
													  }
												  }
											 }
										  } 
										  my $du=scalar(@{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}); 
										  if($du>($Lm-scalar(@S1)-2)){
											  $approved=0;
											  print "rcounter $rcounter not approved $Lm - @S1 - $du \n" if($talkrc==1);
										  }
										  foreach my $n (@{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}){ 
											  delete $g->{$n};
										  }
										  @{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}=sort { $a <=> $b } @{$Barray->{$rcounter}{'basket'}{$basket}{'member'}};
									   }
								}	
								if($approved==0){
									delete $Barray->{$rcounter}; $rcounter--; 
								}	
							}	
					  }
				 }else{
				     print "include all nodes ";
				     if ( $garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'status'} eq 'problem_walk' ){ 
				        		print "(cause walk) ";
				     }
				 }
				 print ", "; 
				 undef $mat_b;
	             if($rcounter==1){ 
	                 if(scalar(keys %{$Barray->{$rcounter}{'basket'}})==1){
	                      if(scalar(@{$Barray->{$rcounter}{'basket'}{1}{'member'}})>100){
	                            if(scalar(@{$Barray->{$rcounter}{'removednodes'}})<scalar(@all_Snodes)){
									if($wemade==1){
									    #print "\n\nonly rcounter $rcounter only one basket\n\n";
										delete $Barray->{$rcounter};
										$rcounter=0;
									}	
	                            }
	                      }
	                 }
	             }
				 print "rcounter $rcounter ";
	             if($rcounter==0){
	                    $rcounter++; @{$Barray->{$rcounter}{'removednodes'}}=@all_Snodes; $Barray->{$rcounter}{'num_removed'}=$Ln; $Barray->{$rcounter}{'score'}=0; my $basket=0;    
	                    my $g=generate_g($garden,$tree,$Old_cycle,$Osub,@all_Snodes); 
						foreach my $m (keys %{$g}){    
							 if (!defined $garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'position'}{$m}){
								  delete $g->{$m}; #print "dddd33\n"; exit;
							 }
						}  
						foreach my $m (keys %{$g}){    
							 foreach my $n (keys %{$g->{$m}}){    
									if (!defined $garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'position'}{$n}){
										#print "$m $n $tree $Old_cycle $Osub\n"; exit;
									}
							 }
						}  
						foreach my $m (keys %$g){ 
							   if (defined $g->{$m}){
								  $basket++; push @{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}, $m;   $Barray->{$rcounter}{'basket'}{$basket}{'orientation'}{$m}=$garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'orientation'}{$m}; $Barray->{$rcounter}{'basket'}{$basket}{'position'}{$m}=$garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'position'}{$m}; 
								  my $Extension=1;
								  while ($Extension==1){
									 $Extension=0; 
									 foreach my $n (@{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}){ 
										  foreach my $s (keys %{$g->{$n}}){  
											  if (!defined $Barray->{$rcounter}{'basket'}{$basket}{'position'}{$s}){
												   $Extension=1;
												   push @{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}, $s;  
												   $Barray->{$rcounter}{'basket'}{$basket}{'orientation'}{$s}=$garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'orientation'}{$s}; $Barray->{$rcounter}{'basket'}{$basket}{'position'}{$s}=$garden->{$tree}{'cycle'}{$Old_cycle}{$Osub}{'position'}{$s}; 	
												   delete $g->{$n}{$s}; delete $g->{$s}{$n}
											  }else{
											   # print "$n $s vay vayyy\n";
											  }
										  }
									 }
								  } 
								  my $du=scalar(@{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}); 
								  foreach my $n (@{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}){ 
									  delete $g->{$n};
								  }
								  @{$Barray->{$rcounter}{'basket'}{$basket}{'member'}}=sort { $a <=> $b } @{$Barray->{$rcounter}{'basket'}{$basket}{'member'}};
							   }
						}				
	             }
	             if ( $rcounter != 1 ){	 
	             	 undef $rcounter;
					 S_LIST:foreach my $rcounter (sort {$Barray->{$a}{'num_removed'} <=> $Barray->{$b}{'num_removed'}} keys %{$Barray}){
							 print "\ntree $tree Osub $Osub rcounter $rcounter removed_n @{$Barray->{$rcounter}{'removednodes'}}\n" if($talkrc==1); 
							 $Barray->{$rcounter}{'status'}='good';
							 foreach my $basket (sort keys %{$Barray->{$rcounter}{'basket'}}){
								   my $Lb=scalar(@{$Barray->{$rcounter}{'basket'}{$basket}{'member'}});  
								   print "basket $basket has $Lb members\n" if($talkrc==1);
								   my $min; my $max; 
								   ($Barray)=B_find_den($Barray,$rcounter,$basket); 
								   if($Barray->{$rcounter}{'basket'}{$basket}{'Bad_den'}>$thr_bad_den){  
										  $Barray->{$rcounter}{'status'}='problem_density';
										  print "rcounter $rcounter basket $basket high density $Barray->{$rcounter}{'basket'}{$basket}{'Bad_den'}\n" if($talkrc==1);
			#	                          my $x=$L-($Barray->{$rcounter}{'basket'}{$basket}{'Bad_den'}/$thr_bad_den)**2;
			#						      if ($x>0){
			#								   $Barray->{$rcounter}{'score'}+=$x**2;
			#							  }else{
			#								   $Barray->{$rcounter}{'score'}-=$x**2;  
			#							  }
			#							   print "x $x\n";
								   }else{
										   my $x=0; 
										   if ($Lb>1){
											   #$x=$L-(4*$Barray->{$rcounter}{'basket'}{$basket}{'Bad_den'}/$thr_bad_den)**4;
											   #$x=$Lb-($Barray->{$rcounter}{'basket'}{$basket}{'Bad_den'});
											   $x=$Lb;
											   if ($x>0){
												   $Barray->{$rcounter}{'score'}+=$x**1.5;
											   }else{
												   $Barray->{$rcounter}{'score'}-=(-$x)**1.5; 
												   print "\nakhe chara $Lb - $Barray->{$rcounter}{'basket'}{$basket}{'Bad_den'} \n" if($talkrc==1);  
											   }
										   }	 
										   print "rcounter $rcounter basket $basket low density $Barray->{$rcounter}{'basket'}{$basket}{'Bad_den'} , $x total score $Barray->{$rcounter}{'score'}\n\n" if($talkrc==1);  
								   }  	                           
							 }
							 print "rcounter $rcounter score $Barray->{$rcounter}{'score'}\n" if($talkrc==1);
							 $Barray->{$rcounter}{'score'}-=50*scalar(@{$Barray->{$rcounter}{'removednodes'}});
							 print "rcounter $rcounter score $Barray->{$rcounter}{'score'}\n" if($talkrc==1);
							 if($Barray->{$rcounter}{'status'} eq 'good'){
									print "Osub $Osub - rcounter $rcounter worked\n" if($talkrc==1);
									$Solved_den=1;
							   }else{
									print "Osub $Osub -  bad rcounter $rcounter\n" if($talkrc==1);
							   } 				   
					 } # foreach my $rcounter
					 if($Solved_den==1){
						 #print "Osub $Osub solved, choose the best rcounter\n";
						 my $max_score; my $best_counter;
						 FORMAX: foreach my $rcounter (keys %{$Barray}){ 
							   if($Barray->{$rcounter}{'status'} eq 'good'){
									  $max_score=$Barray->{$rcounter}{'score'};
									  $best_counter=$rcounter;
									  last FORMAX;
							   }		  
						  } 
						 foreach my $rcounter (keys %{$Barray}){ 
							   if($Barray->{$rcounter}{'status'} eq 'yes'){
									  if($Barray->{$rcounter}{'score'} > $max_score){
										  $max_score=$Barray->{$rcounter}{'score'};
										  $best_counter=$rcounter;
									  }
							   }		  
						  } 
						  #print "best_counter $best_counter tree $tree\n";
						  #print "removed nodes @{$Barray->{$best_counter}{'removednodes'}}\n";
						  foreach my $r (@{$Barray->{$best_counter}{'removednodes'}}){ 
										$paakon->{$r}{'tree'}=$tree;
										$paakon->{$r}{'cycle'}=$Old_cycle;
										$paakon->{$r}{'Osub'}=$Osub;
						  }
						  my $du=scalar(@{$Barray->{$best_counter}{'removednodes'}}); #
						  print " , solved, n_rem $du  new subs: "; #
						  foreach my $basket (keys %{$Barray->{$best_counter}{'basket'}}){    
								  $nsub++;
								  @{$garden->{$tree}{'cycle'}{$cycle}{$nsub}{'member'}}=@{$Barray->{$best_counter}{'basket'}{$basket}{'member'}};
								  #@{$array->{$tree}{'cycle'}{$cycle}{$sub}{'density_win'}}=@{$Barray->{$best_counter}{'basket'}{$basket}{'density_win'}};
								  foreach my $m (@{$garden->{$tree}{'cycle'}{$cycle}{$nsub}{'member'}}){    
									  $garden->{$tree}{'cycle'}{$cycle}{$nsub}{'orientation'}{$m}=$Barray->{$best_counter}{'basket'}{$basket}{'orientation'}{$m};
									  $garden->{$tree}{'cycle'}{$cycle}{$nsub}{'position'}{$m}=$Barray->{$best_counter}{'basket'}{$basket}{'position'}{$m};
									#  $array->{$tree}{'cycle'}{$cycle}{$sub}{'Bad_den'}=$Barray->{$best_counter}{'basket'}{$basket}{'Bad_den'}; 						
								  }
								  my $du=scalar(@{$garden->{$tree}{'cycle'}{$cycle}{$nsub}{'member'}});
								  print "$nsub m $du - "; 
								  
						  } 
						  print "\n"; 
					 }else{
						  #print "Osub $Osub failed, therefore tree $tree is still problematic\n";
						  my $max_score=0; my $best_counter=0;
						  foreach my $rcounter (keys %{$Barray}){   
									  $max_score=$Barray->{$rcounter}{'score'};
									  $best_counter=$rcounter;
									  last;
						  }  
						  foreach my $rcounter (keys %{$Barray}){   
								  if($Barray->{$rcounter}{'score'} > $max_score){
									  $max_score=$Barray->{$rcounter}{'score'};
									  $best_counter=$rcounter;
								  }
						  }  
						  #print "best_counter $best_counter tree $tree removed nodes @{$Barray->{$best_counter}{'removednodes'}}\n";
						  foreach my $r (@{$Barray->{$best_counter}{'removednodes'}}){ 
										$paakon->{$r}{'tree'}=$tree;
										$paakon->{$r}{'cycle'}=$cycle-1;
										$paakon->{$r}{'Osub'}=$Osub;
						  }
						  my $du=scalar(@{$Barray->{$best_counter}{'removednodes'}}); #
						  print " , not solved, n_rem $du  new subs: "; 
						  foreach my $basket (keys %{$Barray->{$best_counter}{'basket'}}){    
								  $nsub++;
								  @{$garden->{$tree}{'cycle'}{$cycle}{$nsub}{'member'}}=@{$Barray->{$best_counter}{'basket'}{$basket}{'member'}};
								  #@{$array->{$tree}{'cycle'}{$cycle}{$sub}{'density_win'}}=@{$Barray->{$best_counter}{'basket'}{$basket}{'density_win'}};
								  foreach my $m (@{$garden->{$tree}{'cycle'}{$cycle}{$nsub}{'member'}}){    
									  $garden->{$tree}{'cycle'}{$cycle}{$nsub}{'orientation'}{$m}=$Barray->{$best_counter}{'basket'}{$basket}{'orientation'}{$m};
									  $garden->{$tree}{'cycle'}{$cycle}{$nsub}{'position'}{$m}=$Barray->{$best_counter}{'basket'}{$basket}{'position'}{$m};
									#  $array->{$tree}{'cycle'}{$cycle}{$sub}{'Bad_den'}=$Barray->{$best_counter}{'basket'}{$basket}{'Bad_den'}; 						
								  }
								  my $du=scalar(@{$garden->{$tree}{'cycle'}{$cycle}{$nsub}{'member'}});
								  if ( defined $Barray->{$best_counter}{'basket'}{$basket}{'Bad_den'} ){
									   print "$nsub den $Barray->{$best_counter}{'basket'}{$basket}{'Bad_den'} m $du , "; 								  
								  }else{
									   print "$nsub m $du , "; 								  								  
								  }
						  }     
						  print "\n"; #
					 } 
				 }else{ 
				     my $best_counter=1;
					  foreach my $r (@{$Barray->{$best_counter}{'removednodes'}}){ 
									$paakon->{$r}{'tree'}=$tree;
									$paakon->{$r}{'cycle'}=$cycle-1;
									$paakon->{$r}{'Osub'}=$Osub;
					  }
					  print " , one opt, new sub: "; 
					  foreach my $basket (keys %{$Barray->{$best_counter}{'basket'}}){    
							  $nsub++;  print "$nsub "; #
							  @{$garden->{$tree}{'cycle'}{$cycle}{$nsub}{'member'}}=@{$Barray->{$best_counter}{'basket'}{$basket}{'member'}};
							  #@{$array->{$tree}{'cycle'}{$cycle}{$sub}{'density_win'}}=@{$Barray->{$best_counter}{'basket'}{$basket}{'density_win'}};
							  foreach my $m (@{$garden->{$tree}{'cycle'}{$cycle}{$nsub}{'member'}}){    
								  $garden->{$tree}{'cycle'}{$cycle}{$nsub}{'orientation'}{$m}=$Barray->{$best_counter}{'basket'}{$basket}{'orientation'}{$m};
								  $garden->{$tree}{'cycle'}{$cycle}{$nsub}{'position'}{$m}=$Barray->{$best_counter}{'basket'}{$basket}{'position'}{$m};
								#  $array->{$tree}{'cycle'}{$cycle}{$sub}{'Bad_den'}=$Barray->{$best_counter}{'basket'}{$basket}{'Bad_den'}; 						
							  }
					  }    	 			 	 
					  print "\n"; #
				 }
	        } # foreach $Osub   
	}
	print "\n"; #
	foreach my $tree (sort { $a <=> $b } keys %{$Pbaagh->{$Old_cycle}}){  
	     foreach my $sub (keys %{$garden->{$tree}{'cycle'}{$cycle}}){
				 #my $L=scalar(@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'member'}});    
				 #print "tree $tree cycle $cycle sub $sub has $L members\n"; undef $L;  
				 my $ontop;
				 ($garden,$ontop)=find_den($garden,$tree,$cycle,$sub);	 # ($garden,$tree,$cycle,$sub);  

	 			 if(($garden->{$tree}{'cycle'}{$cycle}{$sub}{'Bad_den'}<=$thr_bad_den)&&($do_walk==1)){  
						 #print "good tree $tree c $cycle sub $sub  find walk\n";
						 foreach my $m (@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'member'}}){  
							foreach my $n (@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'member'}}){  
							   if($m<$n){
									 if((($garden->{$tree}{'cycle'}{$cycle}{$sub}{'Rstart'}{$n}-$garden->{$tree}{'cycle'}{$cycle}{$sub}{'Rend'}{$m})<15)&&(($garden->{$tree}{'cycle'}{$cycle}{$sub}{'Rstart'}{$n}-$garden->{$tree}{'cycle'}{$cycle}{$sub}{'Rend'}{$m})>-200)){   
											  ($walk,$garden)=findwalk($walk,$m,$n,$garden,$tree,$cycle,$sub);
									 }elsif((($garden->{$tree}{'cycle'}{$cycle}{$sub}{'Rstart'}{$m}-$garden->{$tree}{'cycle'}{$cycle}{$sub}{'Rend'}{$n})<15)&&(($garden->{$tree}{'cycle'}{$cycle}{$sub}{'Rstart'}{$m}-$garden->{$tree}{'cycle'}{$cycle}{$sub}{'Rend'}{$n})>-200)){   
											  ($walk,$garden)=findwalk($walk,$n,$m,$garden,$tree,$cycle,$sub);					      
									 }  									  
							   } 
							}  
						 }     
						 my $problem_walk; my @walk_slist;
						 ($walk,$garden,$problem_walk,@walk_slist)=processwalk($walk,$garden,$tree,$cycle,$sub);
						 #print "end find walk problem_walk $problem_walk\n";  
						 if($problem_walk==1){
							   push @{$Pbaagh->{$cycle}{$tree}{$sub}},@walk_slist;
							   $garden->{$tree}{'cycle'}{$cycle}{$sub}{'status'}='problem_walk';
							   delete $walk->{$tree}{'cycle'}{$cycle}{$sub}; 
							   delete $garden->{$tree}{'cycle'}{$cycle}{$sub}{'village'};
						 }else{
							   buildScaffold($garden,$tree,$cycle,$sub,$walk);  
							   $garden->{$tree}{'cycle'}{$cycle}{$sub}{'status'}='good';
						 }
				 }			
				 if($garden->{$tree}{'cycle'}{$cycle}{$sub}{'Bad_den'}>$thr_bad_den){   
 					$garden->{$tree}{'cycle'}{$cycle}{$sub}{'status'}='problem_density';
					$garden=findhump($garden,$tree,$cycle,$sub);       				
					
					my ($cut_edge, $still_over)=color_nodes($garden,$tree,$cycle,$sub,$ontop);
					
					if(scalar(keys %{$cut_edge})==0){		 
						 foreach my $m (keys %{$ontop}){    
							 foreach my $n (keys %{$ontop->{$m}}){
								  if($ontop->{$m}{$n}>300){ 
										 $cut_edge->{$m}=1; $cut_edge->{$n}=1; #print "cause cutedge was empty, we added $m and $n\n";
								  }elsif($ontop->{$m}{$n}>100){ 
									   if(($tig_info->{$m}{'length'}<300) || ($tig_info->{$n}{'length'}<300)){
										   $cut_edge->{$m}=1; $cut_edge->{$n}=1; #print "cause cutedge was empty, we added $m and $n\n";
									   }
								  }
							 }   
						  }
					}					
					$garden=makeSlist($cut_edge,$garden,$tree,$cycle,$sub);  
			# 		foreach my $m (keys %{$still_over}){
			# 		     foreach my $n (keys %{$still_over->{$m}}){  
			# 		           if ($m<$n){
			# 		                print "still overlap $m and $n. $ontop->{$m}{$n}\n";   
			# 		           }
			# 		     }
			# 		}
					my $hal_s=0;	
					KOHAN: foreach my $camel (keys %{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}}){
						if (defined @{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Slist'}}){
							 my $du=scalar(@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Slist'}});
							 if($du<$combin_thr){
								  push @{$Pbaagh->{$cycle}{$tree}{$sub}},@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Slist'}};	
								  #print "for slist chosed camel $camel num_snodes $du rn1 @{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Slist'}}\n";
								  
								  $garden->{$tree}{'cycle'}{$cycle}{$sub}{'colors'}=$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'summit'};
								  
								  $hal_s=1;
								  last KOHAN;
							 }
						}else{
							#print "\nwhy not defined tree $tree hump $camel\n";  
						}
					}
					if($hal_s==0){
						 GHOLE:foreach my $camel (keys %{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}}){
							if (defined @{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Sumit_Slist'}}){
								 my $du=scalar(@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Sumit_Slist'}});
								 if($du<$combin_thr){
									  push @{$Pbaagh->{$cycle}{$tree}{$sub}},@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Sumit_Slist'}};	 
									  #print "for slist chosed summit camel $camel num_snodes $du rn2 @{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Sumit_Slist'}}\n";
									  
									  $garden->{$tree}{'cycle'}{$cycle}{$sub}{'colors'}=$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'summit'};
									  
									  $hal_s=1;
									  last GHOLE;
								 }
							} 
						 }	
					}
					if($hal_s==0){
						 AVALGHOL:foreach my $camel (keys %{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}}){
							if (defined @{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Sumit_Slist'}}){
								  push @{$Pbaagh->{$cycle}{$tree}{$sub}},@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Sumit_Slist'}};
								  my $du=scalar(@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Sumit_Slist'}}); #print "for slist chosed summit camel $camel num_snodes $du rn3 @{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Sumit_Slist'}}\n";
								  
								  $garden->{$tree}{'cycle'}{$cycle}{$sub}{'colors'}=$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'summit'};
								  
								  $hal_s=1;
								  last AVALGHOL;
							} 
						 }	
					}
					if($hal_s==0){
						 AVALKOH:foreach my $camel (keys %{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}}){
							if (defined @{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Slist'}}){
								  push @{$Pbaagh->{$cycle}{$tree}{$sub}},@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Slist'}};
								  my $du=scalar(@{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Slist'}}); #print "for slist chosed camel $camel num_snodes $du rn4 @{$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'Slist'}}\n";
								  
								  $garden->{$tree}{'cycle'}{$cycle}{$sub}{'colors'}=$garden->{$tree}{'cycle'}{$cycle}{$sub}{'hump'}{$camel}{'summit'};
								  
								  $hal_s=1;
								  last AVALKOH;
							}   
						 }	
					}		
					if($hal_s==0){
						#print "\n\n$tree che khabare?\n\n";
						foreach my $s1 (keys %{$cut_edge}){
							  push @{$Pbaagh->{$cycle}{$tree}{$sub}}, $s1;	
						}	
						$garden->{$tree}{'cycle'}{$cycle}{$sub}{'colors'}=1.5;
						#print "for slist chosed rn @{$Pbaagh->{$cycle}{$tree}{$sub}}\n";     					
					}			 
				 }      
	     }
	}
} 
#-------------------------------------------
my $numd=$daarbast;
#print "scaffold $numd\n";

my $sum;
foreach my $a (@lengths){
     $sum+=$a;
}	

@lengths=sort {$b <=> $a} @lengths;  

my $larges=0;my $N50;
SBASSE: for (my $i=0; $i < scalar(@lengths); $i++){        
      $larges+=$lengths[$i];
      if($larges>($sum/2)){
            $N50=$lengths[$i];
            last SBASSE;
      }
}
print "\nWithout single contigs:\n";  
print "\nN50 $N50       sum $sum\n\n\nAdding single contigs:\n";  

open (INC, "<$address"."/contigs_sopra.fasta") or die "$!";

while (<INC>){
	if ($_=~ /^\>/){
		 my @string = split(/\|/,$_);
		 my $m=substr($string[0],7);
		 if(!defined $tig_info->{$m}{'used'}){
			   my $contig=<INC>; chomp($contig); $numd++; my $Ld=length($contig);
			   
			   printf SCFD ">scaf%i|size%i|single_contig\n", ($numd,$Ld);    
			   push @lengths, $Ld; 
			   $sum+=$Ld;
			   print SCFD "$contig\n"; 
		 }		   
	}   
}
close INC;
print "scaffold $numd\n";
close SCFD;

@lengths=sort {$b <=> $a} @lengths;  

$larges=0; $N50=0;
SBASSE: for (my $i=0; $i < scalar(@lengths); $i++){        
      $larges+=$lengths[$i];
      if($larges>($sum/2)){
            $N50=$lengths[$i];
            last SBASSE;
      }
}
print "\nN50 $N50       sum $sum\n\n";
 
$date = `date`; chomp($date); print "\nFinished $date\n\n"; 
#----------------------------------------------------
sub color_nodes{
	  my ($gar,$derakht,$dayere,$shakhe,$OnTop) =  @_; 
	  my $L= scalar(@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}});
	  my $num_color=${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[0];  
	  for (my $i=1; $i < $L; $i++){   
		 my $pelle=${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i]-${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[($i-1)]; 
		 if ($pelle>0){
			 $num_color+=$pelle;
		 }
	  } undef $L;
	  #print "tree $derakht c $dayere sub $shakhe num_color $num_color\n";  
	  if($num_color==0){
	      print "unexpected error\n"; exit;
	  }	  
	  my $Ln=scalar(@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'member'}});   
	  my $fat=0;#2/($Ln**2);
	  my $C=generate_C($gar,$derakht,$dayere,$shakhe);
	  foreach my $m (keys %{$C}){    
		 if (!defined $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'position'}{$m}){
			  delete $C->{$m};
		 }
	  }  
	  foreach my $m (keys %{$C}){    
		 foreach my $n (keys %{$C->{$m}}){    
				if (!defined $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'position'}{$n}){
					print "unexpected error\n"; exit;
				}
		 }
	  }  
	  my $vasl;
	  foreach my $m (keys %{$C}){    
		 foreach my $n (keys %{$C->{$m}}){  
			 $vasl->{$m}{$n}=$C->{$m}{$n};
		 }
		 
	  }   
	  my $minE_edges=0;
	  foreach my $m (keys %{$C}){   
		 foreach my $n (keys %{$C->{$m}}){
			  $minE_edges-=1;
		 }
	  }
	  $minE_edges/=2; 
# 	      my $MfatE=$fat*$L**2; print "MaxfatE $MfatE\n";
	  my $repul;
	  foreach my $m (keys %{$OnTop}){    
		 foreach my $n (keys %{$OnTop->{$m}}){
			  if($OnTop->{$m}{$n}>300){ 
				 my $du=($OnTop->{$m}{$n})/250;
				 if($du>5){
					 $du=5;
				 }
				 $C->{$m}{$n}-=$du;
				 $repul+=$du;
			  }elsif($OnTop->{$m}{$n}>50){ 
			       if(($tig_info->{$m}{'length'}<300) || ($tig_info->{$n}{'length'}<300)){
			             my $bad=0;
			             if (defined $vasl->{$m}){
			                  if (defined $vasl->{$m}{$n}){
			                     $bad=1;
			                  }
			             }
                         if ($bad==0){
							 my $du=1.3; 
							 $C->{$m}{$n}-=$du;
							 $repul+=$du; 
							 if($m<$n){
								 #print "\nfor color repulsion: $m length $tig_info->{$m}{'length'} - $n length $tig_info->{$n}{'length'} - top $OnTop->{$m}{$n}\n\n";
							 }			   
						 }	 
			       }
			  }
		 }   
	  }  
	  $repul/=2; 

	  my $col;  
	  foreach my $m (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'member'}}){ 
		 $col->{'nodes'}{$m}=0;
		 $col->{'num_mem'}{0}++;
	  }
	  $repul=0;
	  foreach my $m (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'member'}}){ 
		  foreach my $n (keys %{$C->{$m}}){    
			  if($C->{$m}{$n}<0){
				  $repul+=-$C->{$m}{$n};
			  }    
		  }
	  } 
	  $repul/=2; 
	  my $E=0;
	  foreach my $m (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'member'}}){ 
		  foreach my $n (keys %{$C->{$m}}){    
			  if($col->{'nodes'}{$m}==$col->{'nodes'}{$n}){
				  $E+=-$C->{$m}{$n};
			  }    
		  }
	  }
	  $E/=2;  
	  my $du=.01* int(100*$E); $minE_edges=.01* int(100*$minE_edges); my $date = `date`; chomp($date); 
	  print "Current label energy $du     minE_edges $minE_edges    number of labels $num_color   tree $derakht   $date\n";
	  #open (ENR, ">coloring/Cenergy_T$derakht"."_c$dayere"."_sub$shakhe".".txt");   
	  #print ENR "$E ";
	  my $dama=.1; my $time=0; my $saturated=0; my @Memo_E;  my $cut_sat=1+$num_color+int(-$minE_edges/175); #print "minE_edges $minE_edges cut_sat $cut_sat\n"; 
	  for (my $i=0; $i<20; $i++){
		 $Memo_E[$i]=$E;
	  }
	  #my $dama_step=10**(-7);
	  
	  $L=scalar(keys %{$col->{'nodes'}});
	  
	  #my $es=1500*(1+.1* int($L/200)); print "1500 es $es L $L\n";
	  
	  my $es=600*(1+.05* int($L/10)); # print "500 es $es L $L\n";
	  
	  while($saturated==0){    
          if (($time % $es)==0){
              $dama*= 1.001;
          }        		  
		  $time++;
		  my $change=0;
		  my $i = int(rand($Ln));
		  my $element = ${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'member'}}[$i];
		  my $old_col=$col->{'nodes'}{$element};
		  my $new_col= ($old_col+ 1 + int(rand($num_color-1))) % $num_color;
		  my $delta_1=0;
		  #my $delta_1=$fat*(($col->{$tree}{'num_mem'}{$old_col}-1)**2+($col->{$tree}{'num_mem'}{$new_col}+1)**2-($col->{$tree}{'num_mem'}{$old_col})**2-($col->{$tree}{'num_mem'}{$new_col})**2);
		  my $delta_2=0;
		  foreach my $n (keys %{$C->{$element}}){    
			  if($col->{'nodes'}{$n}==$old_col){
				  $delta_2+=$C->{$element}{$n};
			  }elsif($col->{'nodes'}{$n}==$new_col){
				  $delta_2+=-$C->{$element}{$n};
			  }
		  }
		  my $delta=$delta_1+$delta_2;
		  if ($delta>0){   
			   if (rand()<exp(-$dama*$delta)){  
				  $change=1;
			   }  
		  }else{       
			   $change=1;
		  }
		  if ($change){
			 $col->{'nodes'}{$element}=$new_col; 
			 $col->{'num_mem'}{$old_col}-=1;
			 $col->{'num_mem'}{$new_col}+=1;
			 $E+=$delta;               
		  }
		  if (($time % 50000)==0){
			 #print ENR "$E ";
			 push @Memo_E, $E;
			 shift(@Memo_E);
			 if (($time % 1000000)==0){
				$saturated=1;
				for (my $i=0; $i < 15; $i++){
					 if( (($Memo_E[$i]-$E)>$cut_sat) || (($Memo_E[$i]-$E)<0) ){ 
						  $saturated=0;
					 }
				}
			 }							    
		  }
		  if (($time % 3000000)==0){
		     my $du=.01* int(100*$E);
			 print  "Current label energy $du\n";
		  }
	  } #print ENR "$E ";
	  foreach my $m (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'member'}}){    # for flip
		  my $delta; my $best_color;
		  for (my $i=0; $i <  $num_color; $i++){   
			  $delta->{$i}=0;
			  foreach my $n (keys %{$C->{$m}}){    
				  if($col->{'nodes'}{$n}==$i){
					  $delta->{$i}+=-$C->{$m}{$n};
				  }    
			  }
		  }		
		  foreach my $rang (sort { $delta->{$a} <=> $delta->{$b} } keys %{$delta}) {
				  $best_color=$rang;
				  last;
		  }
		  if(($col->{'nodes'}{$m}!=$best_color)&&($delta->{$best_color}<$delta->{$col->{'nodes'}{$m}})){
			  #print "\nflip colorrrrrrrrrrrrrrrrrrrrrrrrrr  of $m new E $delta->{$best_color} old $delta->{$col->{'nodes'}{$m}}\n\n";
			  $E+=$delta->{$best_color}-$delta->{$col->{'nodes'}{$m}};
			  $col->{'num_mem'}{$col->{'nodes'}{$m}}-=1;
			  $col->{'num_mem'}{$best_color}+=1;
			  $col->{'nodes'}{$m}=$best_color;
		  }
	  }  #print ENR "$E ";
# 	  for (my $i=0; $i <  $num_color; $i++){  
# 		 open (DEN, ">coloring/color_T$derakht"."c$dayere"."_sub$shakhe"."_col$i".".txt");   
# 		 foreach my $m (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'member'}}){ 
# 			 if($col->{'nodes'}{$m}==$i){
# 				print DEN "$m\n";   
# 			 }
# 		 } 
# 		 close DEN;
# 	  }
	  my $NE=0;
	  foreach my $m (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'member'}}){ 
		  foreach my $n (keys %{$C->{$m}}){    
			  if($col->{'nodes'}{$m}==$col->{'nodes'}{$n}){
				  $NE+=-$C->{$m}{$n};
			  }    
		  }
	  }
	  $NE/=2; 
	  $du=.01* int(100*$E); $date = `date`; chomp($date); 
	  print "Final label energy $du     $date\n\n";
	  if($E!=$NE){
		  #print "NE $NE E $E\n"; 
	  }
# 	      for (my $i=0; $i < $num_color; $i++){
# 	          $E+=$fat*($col->{'num_mem'}{$i})**2; 
# 	      }
# 	      print "tot E $E\n";
	  my $c_edge; $du=0; my $Still_over;
	  foreach my $m (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'member'}}){ 
		  foreach my $n (keys %{$vasl->{$m}}){    
			  if(defined $vasl->{$m}{$n}){ 
				  if($col->{'nodes'}{$m}!=$col->{'nodes'}{$n}){
					  $c_edge->{$m}{$n}=1;
					  $c_edge->{$n}{$m}=1;
				  } 
			  }
		  }
	  }
	  foreach my $m (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'member'}}){ 
		  foreach my $n (keys %{$C->{$m}}){    
			  if($C->{$m}{$n}<0){ 
				  if($col->{'nodes'}{$m}==$col->{'nodes'}{$n}){
					  $Still_over->{$m}{$n}=1;
					  $Still_over->{$n}{$m}=1;
				  } 
			  }
		  }
	  }
	  return $c_edge, $Still_over;
}
#----------------------------------------------------
sub makeSlist{
    my ($c_edge,$gar,$derakht,$dayere,$shakhe) = @_;
    my $used_s;
#     foreach my $S (keys %{$c_edge}){
#         print "s $S $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rstart'}{$S} $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rend'}{$S}\n";
#     } print "\n";
	my (@S0,$neig); my $g=generate_g($gar,$derakht,$dayere,$shakhe,@S0); 
	
	foreach my $m (keys %$g){  
			$neig->{$m}= scalar(keys %{$g->{$m}});
	}	

	foreach my $camel (keys %{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}}){
    #print "camel $camel\n";
		 foreach my $S (keys %{$c_edge}){
			  if ((($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit_start'}<$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rstart'}{$S})&&($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit_end'}>$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rstart'}{$S})) || (($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit_start'}<$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rend'}{$S})&&($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit_end'}>$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rend'}{$S}))){
				   if ( $neig->{$S} > 2 ){   
			    	   push(@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'Sumit_Slist'}},$S); #$used_s->{$S}=1;
			       }
			  }
			  if ((($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'start'}<$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rstart'}{$S})&&($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'end'}>$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rstart'}{$S})) || (($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'start'}<$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rend'}{$S})&&($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'end'}>$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rend'}{$S}))){
				    if ( $neig->{$S} > 2 ){   
				       push(@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'Slist'}},$S); #$used_s->{$S}=1;
			        }
			  }
		 }
	}
	return $gar;	
}		
#----------------------------------------------------		 
sub findwalk{
       my ($W,$a,$b,$gar,$derakht,$dayere,$shakhe) = @_;
       my $read_a; my $read_b;

       if ($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'orientation'}{$a}==-1){
			 $read_a =reverseComplement(findDnascontig($a));
	   }else{
			 $read_a = findDnascontig($a);
       }
       if ($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'orientation'}{$b}==-1){
			 $read_b =reverseComplement(findDnascontig($b));
	   }else{
			 $read_b = findDnascontig($b);
       }  
	   my @read_a_end=split(//,substr($read_a,-$max_walk));
       my @read_b_start=split(//,substr($read_b,0,$max_walk));
       BASEDIGE: for (my $i=$max_walk; $i > $min_walk; $i--){
           my $yes_m=0; my $no_m=0;
           for (my $j=0; $j < $i; $j++){
               if($read_a_end[$j] eq $read_b_start[$j]){
                   $yes_m++;
               }else{
                   $no_m++;
               }
           }
           if(($yes_m/$i)>.95){
               $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$a}{$b}=1;
               push(@{$W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'num_m'}{$a}{$b}}, $i);
               push(@{$W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'score'}{$a}{$b}}, $yes_m);
               last BASEDIGE;
           }
           shift(@read_a_end);
           pop(@read_b_start);
       }

#        if(defined $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$a}){
# 		   if(defined $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$a}{$b}){
# 			 #print "$gar->{$tree}{'orientation'}{$a} $a $garden->{$tree}{'orientation'}{$b} $b walk loc @{$W->{'num_m'}{$a}{$b}} score @{$W->{'score'}{$a}{$b}}\n";
# 			 #my $ad=$gar->{$tree}{'Rstart'}{$b}-$gar->{$tree}{'Rend'}{$a};
# 			 #print "-@{$W->{'num_m'}{$a}{$b}} $ad\n";
# 		   }
# 	   }	   
# 	   if(defined $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$b}){
# 		   if(defined $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$b}{$a}){
# 			 #print "$gar->{$tree}{'orientation'}{$b} $b $gar->{$tree}{'orientation'}{$a} $a walk loc @{$W->{'num_m'}{$b}{$a}} score @{$W->{'score'}{$b}{$a}}\n";
# 			 #my $ad=$gar->{$tree}{'Rstart'}{$a}-$gar->{$tree}{'Rend'}{$b};
# 			 #print "-@{$W->{'num_m'}{$b}{$a}} $ad\n";
# 		   }
# 	   }
	   if((defined $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$a})&&(defined $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$b})){
		   if((defined $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$a}{$b})&&(defined $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$b}{$a})){
			 print "\n$b $a both walk\n" if ($talkw==1);
			 delete $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$a};
			 delete $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$b};
			 $W->{'deleted'}{$a}{$b}=1;
			 $W->{'deleted'}{$a}{$b}=1;
		   }
	   }	   
       return $W, $gar;
} 
#----------------------- 
sub processwalk{     
     my ($W,$gar,$derakht,$dayere,$shakhe) = @_;
	 my $entry; 
	 foreach my $m (keys %{$W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}}){
		   my @exit_node=keys %{$W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$m}};
		   if (scalar(@exit_node)>1){
			#	print "delete double exit for $m -> @exit_node\n";
				foreach my $n (keys %{$W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$m}}){
					  $W->{'deleted'}{$m}{$n}=1;
				}
				delete $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$m};
		   }elsif(scalar(@exit_node)==1){
			   $entry->{$exit_node[0]}++;
		   }else{
			   print "@exit_node m $m\n"; exit;
		   }
	 }
	 foreach my $n (keys %{$entry}){
		 if ($entry->{$n}==0){
			 print "kessh\n"; exit;
		 }elsif($entry->{$n}>1){
			 foreach my $m (keys %{$W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}}){
				 if (defined $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$m}{$n}){
					  delete $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$m};
					  $W->{'deleted'}{$m}{$n}=1; #	  print "delete double entry for $n, deleted $m\n";
				 }
			 }
		 }
	 }
	 my $problem_w=1;
	 while($problem_w==1){
		 $problem_w=0;
		 delete $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'};	
		 my $starter;
		 my $entry;  # second round
		 foreach my $m (keys %{$W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}}){
			   my @exit_node=keys %{$W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$m}};
			   if (scalar(@exit_node)>1){
					print "hala vayy bay\n"; exit;	
			   }elsif(scalar(@exit_node)==1){
				   $entry->{$exit_node[0]}++;
			   }else{
				   print "@exit_node m $m\n"; exit;
			   }
		 }
		 my $num_starter=0;
		 foreach my $m (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'member'}}){  
			 if ( !defined $entry->{$m} ){
				$starter->{$m}=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rstart'}{$m};
				$num_starter++; # print "starter $m\n";
			 }elsif($entry->{$m}==1){
				  
			 }elsif($entry->{$m}>1){
			 
			 }		 
		 }   
		 print "num_starter $num_starter\n" if ($talkw==1);
		 my $roosta=0;
		 foreach my $S (sort { $starter->{$a} <=> $starter->{$b} } keys %{$starter}){  
				$roosta++; 
				print "starter $S Rstart $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rstart'}{$S} village $roosta\n" if ($talkw==1);
				$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$roosta}{'start'}=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rstart'}{$S};
				push (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$roosta}{'members'}},$S);  
		    	if ($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'orientation'}{$S}==-1){
					 $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$roosta}{'Dstring'}=reverseComplement(findDnascontig($S));  
				}else{
					 $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$roosta}{'Dstring'}=findDnascontig($S);  
				}
				my $extension=1;
				my $nok=$S;
				while($extension==1){
					  $extension=0;
					  if(defined $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$nok}){
						  my @exit_node=keys %{$W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$nok}};
						  if (scalar(@exit_node)>1){
								print "aha\n"; exit;
						  }elsif(scalar(@exit_node)==1){
							   $extension=1;
							   my $overlap=${$W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'num_m'}{$nok}{$exit_node[0]}}[0];
							   my $du=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rstart'}{$exit_node[0]}-$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rend'}{$nok};
							   if($du>50){
							       #print "\n\n\n\ngap $du -- $nok to $exit_node[0]\n\n\n";
							   }							   
							   $nok=$exit_node[0]; 
							   push (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$roosta}{'members'}},$nok);  
							   if ($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'orientation'}{$nok}==-1){ 
									$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$roosta}{'Dstring'}.=substr(reverseComplement(findDnascontig($nok)),$overlap);
							   }else{
									$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$roosta}{'Dstring'}.=substr(findDnascontig($nok),$overlap);
							   }  
							   
						  }
					  }		  
				}
				$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$roosta}{'end'}=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$roosta}{'start'}+length($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$roosta}{'Dstring'})-1;
				     
		 }
		 for (my $i=1; $i < $num_starter; $i++){
			  for (my $j=$i+1; $j <= $num_starter; $j++){
					 if(($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$j}{'end'}<$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$i}{'end'})&& ($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$j}{'start'}<($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$i}{'end'}-2000))){
                    		 print "dehat $j is in the middle of $i\n" if ($talkw==1);
							 foreach my $m (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$i}{'members'}}){  
								 if (defined $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$m}){
										foreach my $n (keys %{$W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$m}}){
											  if(($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rstart'}{$m}<=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$j}{'end'}) && ($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rend'}{$n}>=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$j}{'start'})){
												   #print "$m chap $n raast - \n";
												   $W->{'deleted'}{$m}{$n}=1; 
												   delete $W->{$derakht}{'cycle'}{$dayere}{$shakhe}{'jump'}{$m};
												   $problem_w=1;
											  }
										}
								 }		
							 }
					 }
			  }
		 }        
	 }

	 my $num_village=scalar(keys %{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}}); 

	 for (my $i=1; $i < $num_village; $i++){
		  for (my $j=$i+1; $j <= $num_village; $j++){
				 if(($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$j}{'end'}<$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$i}{'end'})&& ($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$j}{'start'}<($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$i}{'end'}-2000))){			   		   
					   $problem_w=1;
				 }
		  }	 
	 }  
	 my @W_slist;
     if ( $problem_w==1 ){
		 for (my $i=1; $i < $num_village; $i++){
			  for (my $j=$i+1; $j <= $num_village; $j++){
					 if(($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$j}{'end'}<$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$i}{'end'})&& ($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$j}{'start'}<($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$i}{'end'}-2000))){
							 foreach my $m (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$i}{'members'}}){  
								   foreach my $n (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$j}{'members'}}){  
											  if(( $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rend'}{$n} < $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rend'}{$m})&& ( $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rstart'}{$n}<($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rend'}{$m}-2000))){			   		   
													push @W_slist , $n;																	
											  }
								   }		
							 }
					 }
			  }
		 } 
	 }	 		 
	 if ( ($problem_w==1) && ( scalar(@W_slist) == 0) ){
	 	 print "cause it is not defined\n\n\n\n\n\n\n\n";	                
		 for (my $i=1; $i < $num_village; $i++){
			  for (my $j=$i+1; $j <= $num_village; $j++){
					 if(($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$j}{'end'}<$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$i}{'end'})&& ($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$j}{'start'}<($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$i}{'end'}-2000))){			   		   
						   push @W_slist , @{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$i}{'members'}};
						   push @W_slist , @{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$j}{'members'}};
						   $problem_w=1;
					 }
			  }	 
		 }                  
	 }
     return $W, $gar, $problem_w, @W_slist;
}       
#----------------------- -AD      
sub buildScaffold{
     my ($gar,$derakht,$dayere,$shakhe,$W) = @_;
	 $daarbast++;
	 my ($num_Dgood, $num_Dbad, $totL_Dgood, $totL_Dbad, $tot_Lgap, $num_goodGap, $num_badGap, $num_verybadGap, $num_2verybadGap, $num_3verybadGap)=(0,0,0,0,0,0,0,0,0,0);
	 my $num_village=scalar(keys %{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}});
	 my $current_vil=1;  		
	 foreach my $m (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{1}{'members'}}){
	      $tig_info->{$m}{'used'}=$daarbast;
	 }
	 print "shoroo darbast $daarbast ba $current_vil numbervil $num_village tree $derakht c $dayere sub $shakhe\n" if ($talksc==1);    
     
	 my $bast_Dstring=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$current_vil}{'Dstring'};       
	 for (my $i=2; $i <= $num_village; $i++){
			  if ($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$current_vil}{'end'}<$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$i}{'start'}){
				   my $du=($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$i}{'start'}-$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$current_vil}{'end'});
				    
				   $du=5*(1+int($du/5)); my $gap='N' x $du;				   
				   $bast_Dstring.=$gap;
				   $bast_Dstring.=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$i}{'Dstring'};     
 
		  		   #print "$i added to $current_vil after $Lgap gap R_gap $R_gap\n" if ($talksc==1);
			  }else{
			  	   my $gap='N' x 5;
				   $bast_Dstring.=$gap;   				   
				   $bast_Dstring.=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$i}{'Dstring'};     

		  		   #print "$i added to $current_vil after X gap R_gap $R_gap\n" if ($talksc==1);
			  }
			  foreach my $m (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'}{$i}{'members'}}){
	               $tig_info->{$m}{'used'}=$daarbast;
	          }
			  $current_vil=$i;		  
	 }
	 delete $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'village'};
	 
	 my $Ld=length($bast_Dstring);  	 
	 printf SCFD ">scaf%i|size%i\n", ($daarbast,$Ld);    
	 print SCFD "$bast_Dstring\n";  
	 push @lengths, $Ld;
	 print "scaffold $daarbast length $Ld\n";           
}       
#----------------------------------------------------
sub findhump{
     my ($gar,$derakht,$dayere,$shakhe) = @_;
	 #print "min $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'min'} max $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'max'}\n";
	 my $camel=0;
	 if (${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[0]>1){
		   $camel++;
		   $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'start'}=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'min'};
		   $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit'}=${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[0];
		   $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit_start'}=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'min'};
	 }
	 my $L= scalar(@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}});
	 for (my $i=1; $i < $L; $i++){   
		 my $pelle=${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i]-${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[($i-1)]; 
		 if (($pelle>0)&&(${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[($i-1)]==1)){
			 $camel++;
			 $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'start'}=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'min'}+($i-1)*1000;  
			 $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit'}=${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i];
			 $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit_start'}=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'min'}+($i-1)*1000;
		 }
		 if (($pelle<0)&&(${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[($i)]==1)){
			 $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'end'}=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'min'}+$i*1000+2000;  
		 }
		 if($camel>0){
			 if(${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i]>$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit'}){
					$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit'}=${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i];
					$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit_start'}=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'min'}+($i-1)*1000;
			 }
		 }		 
		 if (($pelle<0)&&(${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[($i-1)]==$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit'})){
				$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit_end'}=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'min'}+$i*1000+2000;  
		 }
	 }
	 if (${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[($L-1)]>1){
	   $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'end'}=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'max'};  
	   if (${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[($L-1)]==$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit'}){
		   $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit_end'}=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'max'};  
	   }
	}
	#print "tree $derakht c $dayere sub $shakhe num_hump $camel\n";  
	foreach my $camel (keys %{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}}){   
		#print "hump $camel $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'start'} $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'end'}\n";	
		#print "summit $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit'} start $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit_start'} end $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'hump'}{$camel}{'summit_end'}\n\n";
	}       
	return $gar;
}       
#-----------------------      
sub find_den{    
   my ($gar,$derakht,$dayere,$shakhe) =  @_;    
   my $OnTop; 	my $e1; my $e2; $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Bad_den'}=0;
   foreach my $z (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'member'}}){  
	   if ($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'orientation'}{$z}==1){
	       $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rstart'}{$z}=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'position'}{$z}; 
           $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rend'}{$z}=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'position'}{$z}+$tig_info->{$z}{'length'}-1; 
           
           $e1=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'position'}{$z}; $e2=($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'position'}{$z}+$tig_info->{$z}{'length'});
           
           $e1=round_to_10($e1); $e2=round_to_10($e2);		   
           
		   for (my $i=$e1; $i <= $e2; $i+=10){   
			   push (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'whoishere'}{$i}}, $z);
		   }
	   }else{
	       $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rend'}{$z}=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'position'}{$z}; 
           $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Rstart'}{$z}=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'position'}{$z}-$tig_info->{$z}{'length'}+1;  
           
           $e1=$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'position'}{$z}; $e2=($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'position'}{$z}-$tig_info->{$z}{'length'});
           
           $e1=round_to_10($e1); $e2=round_to_10($e2);		  
           
		   for (my $i=$e1; $i >= $e2 ; $i-=10){    
			   push (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'whoishere'}{$i}}, $z);
		   }       
	   } 
   }         
     my $Max; my $Min; 
	 my @Sposition= keys %{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'whoishere'}};
	 @Sposition= sort {$a <=> $b} @Sposition;
	 $Max=pop(@Sposition); $Min=shift(@Sposition); my $interval=$Max-$Min; 

	 $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'max'}=$Max; $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'min'}=$Min; 
	 
	 my $num_window=int($interval/1000);
	 for (my $k=0; $k < ($num_window-1); $k++){  
		   my $local_den=0; 
		   for (my $r=($Min+$k*1000); $r < ($Min+$k*1000+2000); $r+=10){
			  if ( defined $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'whoishere'}{$r} ){
			  	  $local_den+=10* scalar(@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'whoishere'}{$r}});
			  }				  
		   }
		   $local_den/=2000;
		   
		   $local_den=.001* int(1000*$local_den);
		   
		   push (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}},$local_den);
		   if ($local_den>2){
		        $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Bad_den'}+=$thr_bad_den;
		   }elsif($local_den>1.7){ 
		        $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Bad_den'}+=$thr_bad_den;
		   }elsif($local_den>1.4){
		        $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Bad_den'}+=$thr_bad_den; 
		   }elsif($local_den>1.2){
		        $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Bad_den'}+=.5*$thr_bad_den;
# 		   }elsif($local_den>1.05){
# 		        $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Bad_den'}+=1;
		   }
	 }
	 my $local_den=0; 
	 for (my $r=($Min+($num_window-1)*1000); $r < $Max; $r+=10){
			  if ( defined $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'whoishere'}{$r} ){
			  	  $local_den+=10* scalar(@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'whoishere'}{$r}});
			  }		
	 }
	 $local_den/=($interval-($num_window-1)*1000); $local_den=.001* int(1000*$local_den);
	 push (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}},$local_den);
	   if ($local_den>2){
			$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Bad_den'}+=$thr_bad_den;
	   }elsif($local_den>1.7){ 
			$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Bad_den'}+=$thr_bad_den;
	   }elsif($local_den>1.4){
			$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Bad_den'}+=$thr_bad_den; 
	   }elsif($local_den>1.2){
			$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Bad_den'}+=.5*$thr_bad_den;
	# 		   }elsif($local_den>1.05){
	# 		        $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Bad_den'}+=1;
	   }
	 for (my $i=$Min; $i <= $Max; $i+=10){
		  foreach my $m (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'whoishere'}{$i}}){
			   foreach my $n (@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'whoishere'}{$i}}){
				   if ($m!=$n){
						$OnTop->{$m}{$n}+=10;
						$OnTop->{$n}{$m}+=10;
				   }
			   }
		  }
	 } 
	 
	 if($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Bad_den'}<=$thr_bad_den){
	     #print "min $Min max $Max  Tree $derakht low density $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Bad_den'}\n\n";
	 }
	 if($gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Bad_den'}>$thr_bad_den){   
	     print "\nmin $Min max $Max  tree $derakht sub $shakhe cycle $dayere high density $gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'Bad_den'}\n\n";
	     my $L= scalar(@{$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}});
		 for (my $i=0; $i < $L; $i++){   
			 if ( (${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i]-int(${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i]))<.2  ){
				 ${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i]=int(${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i]);
			 }else{
				 ${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i]=int(${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i])+1;
			 }
			 if (${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i]==0){
				${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i]=1;
			 }
		 }
		 for (my $i=1; $i < ($L-1); $i++){   
			 my $Left_D=${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i-1]-${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i];
			 my $Right_D=${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i+1]-${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i];
			 if (($Left_D!=0)||($Right_D!=0)){
				   if (($Left_D==1)&&($Right_D==1)){
					   ${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i]+=1;
				   }
				   if (($Left_D==-1)&&($Right_D==-1)){
					   #${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i]-=1;
				   }
				   if ((($Left_D==1) && ($Right_D==0)) || (($Left_D==0) && ($Right_D==1)) ){
					  if ((($i-1)>0)&&(($i+2)<$L)){   
						   my $LLeft_D=${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i-2]-${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i];
						   my $RRight_D=${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i+2]-${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i];  
						   if (($LLeft_D==1)&&($RRight_D==1)){
							   ${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i]+=1;
							   #print "fill the gap\n";
						   }
						  # if (($LLeft_D==-1)&&($RRight_D==-1)){
						  #	   ${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[$i]-=1;
						  # }
					  }		   
				   }
			 }
		 } 
		 if ( ${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[0] < ${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[1] ){
		       ${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[0]=${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[1];
		 }
		 if ( ${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[($L-1)] < ${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[($L-2)] ){
		       ${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[($L-1)]=${$gar->{$derakht}{'cycle'}{$dayere}{$shakhe}{'density_win'}}[($L-2)];
		 }		 
	 }     
	 return $gar, $OnTop;
}     
#------------------------------------------- 
sub B_find_den{    
   my ($Ba,$rc,$bas) =  @_;  
   	my $e1; my $e2; $Ba->{$rc}{'basket'}{$bas}{'Bad_den'}=0;
   foreach my $z (@{$Ba->{$rc}{'basket'}{$bas}{'member'}}){  
	   if ($Ba->{$rc}{'basket'}{$bas}{'orientation'}{$z}==1){
	       
	       $e1=$Ba->{$rc}{'basket'}{$bas}{'position'}{$z}; $e2=$Ba->{$rc}{'basket'}{$bas}{'position'}{$z}+$tig_info->{$z}{'length'}; 

           $e1=round_to_10($e1); $e2=round_to_10($e2);		   	  

		   for (my $i=$e1; $i <= $e2 ; $i+=10){   
			   $Ba->{$rc}{'basket'}{$bas}{'density'}{$i}+=10;
		   }
	   }else{
	   	   $e1=$Ba->{$rc}{'basket'}{$bas}{'position'}{$z};  $e2=$Ba->{$rc}{'basket'}{$bas}{'position'}{$z}-$tig_info->{$z}{'length'};		

           $e1=round_to_10($e1); $e2=round_to_10($e2);		   

		   for (my $i=$e1; $i >= $e2 ; $i-=10){   
			   $Ba->{$rc}{'basket'}{$bas}{'density'}{$i}+=10;
		   }       
	   } 
   }         
     my $Max; my $Min;
	 my @Sposition= keys %{$Ba->{$rc}{'basket'}{$bas}{'density'}};
	 @Sposition= sort {$a <=> $b} @Sposition;
	 $Max=pop(@Sposition); $Min=shift(@Sposition); my $interval=$Max-$Min; 
	 my $num_window=int($interval/1000); 
	 for (my $k=0; $k < ($num_window-1); $k++){  
		   my $local_den=0; 
		   for (my $r=($Min+$k*1000); $r < ($Min+$k*1000+2000); $r+=10){
		   	  if ( defined $Ba->{$rc}{'basket'}{$bas}{'density'}{$r} ){
			  	   $local_den+=$Ba->{$rc}{'basket'}{$bas}{'density'}{$r};
			  }	
		   }
		   $local_den/=2000; $local_den=.001* int(1000*$local_den);
		   if ($local_den>2){
		        $Ba->{$rc}{'basket'}{$bas}{'Bad_den'}+=$thr_bad_den;
		   }elsif($local_den>1.7){ 
		        $Ba->{$rc}{'basket'}{$bas}{'Bad_den'}+=$thr_bad_den;
		   }elsif($local_den>1.4){
		        $Ba->{$rc}{'basket'}{$bas}{'Bad_den'}+=$thr_bad_den; 
		   }elsif($local_den>1.2){
		        $Ba->{$rc}{'basket'}{$bas}{'Bad_den'}+=.5*$thr_bad_den;
# 		   }elsif($local_den>1.05){
# 		        $Ba->{$rc}{'basket'}{$bas}{'Bad_den'}+=1;
		   }
	 }
		   
	 my $local_den=0; 
	 for (my $r=($Min+($num_window-1)*1000); $r < $Max; $r+=10){
		  if ( defined $Ba->{$rc}{'basket'}{$bas}{'density'}{$r} ){
			   $local_den+=$Ba->{$rc}{'basket'}{$bas}{'density'}{$r};
		  }	
	 }
	 $local_den/=($interval-($num_window-1)*1000); $local_den=.001* int(1000*$local_den);
	 if ($local_den>2){
		$Ba->{$rc}{'basket'}{$bas}{'Bad_den'}+=$thr_bad_den;
	 }elsif($local_den>1.7){
		$Ba->{$rc}{'basket'}{$bas}{'Bad_den'}+=$thr_bad_den;
	 }elsif($local_den>1.4){
		$Ba->{$rc}{'basket'}{$bas}{'Bad_den'}+=$thr_bad_den; 
	 }elsif($local_den>1.2){
		$Ba->{$rc}{'basket'}{$bas}{'Bad_den'}+=.5*$thr_bad_den;
# 	 }elsif($local_den>1.05){
# 		$Ba->{$rc}{'basket'}{$bas}{'Bad_den'}+=1;
	 }
	 	 
	 delete $Ba->{$rc}{'basket'}{$bas}{'density'}; 
	 return $Ba;
}     
#-------------------------------------------       
sub generate_C{
    my ($gaa,$derakh,$dayer,$shakh) =  @_;
    my $G;
	foreach my $m (keys %$Ori_Dis){ 
	   if(defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'position'}{$m}){
				 foreach my $n (keys %{$Ori_Dis->{$m}}){ 
					 if(defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'position'}{$n}){
						 if ((defined $Ori_Dis->{$m}{$n}{'P'}) && (! defined $Ori_Dis->{$m}{$n}{'N'})){
							if ($Ori_Dis->{$m}{$n}{'P'}>=$minlink ){
									if((defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$m})&&(defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$n})){
											 if (($gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$m}*$gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$n})==1){               
												$G->{$m}{$n}=1;
												$G->{$n}{$m}=1;
											 }else{
				
											 }
									}else{
									       print "ajabaa $m n $n"; exit;
												$G->{$m}{$n}=1;
												$G->{$n}{$m}=1;						
									}
							}
						 }elsif ((! defined $Ori_Dis->{$m}{$n}{'P'}) && (defined $Ori_Dis->{$m}{$n}{'N'})){
							if ($Ori_Dis->{$m}{$n}{'N'}>=$minlink ){ 
								if((defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$m})&&(defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$n})){
											 if (($gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$m}*$gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$n})==-1){               
												$G->{$m}{$n}=1;
												$G->{$n}{$m}=1;
											 }else{
				
											 }
								}else{
								     print "ajabaa $m n $n"; exit;
												$G->{$m}{$n}=1;
												$G->{$n}{$m}=1; 										
								}
							}		 
						 }elsif((defined $Ori_Dis->{$m}{$n}{'P'}) && (defined $Ori_Dis->{$m}{$n}{'N'})){
							if ((.4*$Ori_Dis->{$m}{$n}{'P'})>$Ori_Dis->{$m}{$n}{'N'}){
								if ($Ori_Dis->{$m}{$n}{'P'}>=$minlink ){	
									 if((defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$m})&&(defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$n})){
											if (($gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$m}*$gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$n})==1){                      
												$G->{$m}{$n}=1;
												$G->{$n}{$m}=1;
											}
									 }else{
									      print "ajabaa $m n $n"; exit;
												$G->{$m}{$n}=1;
												$G->{$n}{$m}=1;	 								 
									 }
								}		
							}elsif ((.4*$Ori_Dis->{$m}{$n}{'N'})>$Ori_Dis->{$m}{$n}{'P'}){
								if ($Ori_Dis->{$m}{$n}{'N'}>=$minlink ){	
									 if((defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$m})&&(defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$n})){
											if (($gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$m}*$gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$n})==-1){               
												$G->{$m}{$n}=1;
												$G->{$n}{$m}=1;
											}    
									 }else{
									      print "ajabaa $m n $n"; exit;
												$G->{$m}{$n}=1;
												$G->{$n}{$m}=1;								 
									 }
								}	
							}else{
	
							}
							
						 }
					 }		 
				 }
	   }		 
	}
	return $G;
} 
#------------------------------------------- 
sub generate_g{
    my ($gaa,$derakh,$dayer,$shakh,@rn)= @_; 
    my $Lpaak;
    foreach my $m (@rn){
       $Lpaak->{$m}=1;
    }
    my $G;
	foreach my $m (keys %$Ori_Dis){ 
	   if((defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'position'}{$m})&&(!defined $Lpaak->{$m})){
				 foreach my $n (keys %{$Ori_Dis->{$m}}){ 
					 if((defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'position'}{$n})&&(!defined $Lpaak->{$n})){
						 if ((defined $Ori_Dis->{$m}{$n}{'P'}) && (! defined $Ori_Dis->{$m}{$n}{'N'})){
							if ($Ori_Dis->{$m}{$n}{'P'}>=$minlink ){
									if((defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$m})&&(defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$n})){
											 if (($gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$m}*$gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$n})==1){         
												$G->{$m}{$n}=1;
												$G->{$n}{$m}=1;
											 }else{
				
											 }
									}else{
									      print "ajabaa $m n $n"; exit;
												$G->{$m}{$n}=1;
												$G->{$n}{$m}=1;						
									}
							}
						 }elsif ((! defined $Ori_Dis->{$m}{$n}{'P'}) && (defined $Ori_Dis->{$m}{$n}{'N'})){
							if ($Ori_Dis->{$m}{$n}{'N'}>=$minlink ){ 
								if((defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$m})&&(defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$n})){
											 if (($gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$m}*$gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$n})==-1){       
												$G->{$m}{$n}=-1;
												$G->{$n}{$m}=-1;
											 }else{
				
											 }
								}else{  
								        print "ajabaa $m n $n"; exit;
												$G->{$m}{$n}=-1;
												$G->{$n}{$m}=-1; 										
								}
							}		 
						 }elsif((defined $Ori_Dis->{$m}{$n}{'P'}) && (defined $Ori_Dis->{$m}{$n}{'N'})){
							if ((.4*$Ori_Dis->{$m}{$n}{'P'})>$Ori_Dis->{$m}{$n}{'N'}){
								if ($Ori_Dis->{$m}{$n}{'P'}>=$minlink ){	
									 if((defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$m})&&(defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$n})){
											if (($gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$m}*$gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$n})==1){ 
												$G->{$m}{$n}=1;
												$G->{$n}{$m}=1;
											}
									 }else{
									       print "ajabaa $m n $n"; exit;
												$G->{$m}{$n}=1;
												$G->{$n}{$m}=1;	 								 
									 }
								}		
							}elsif ((.4*$Ori_Dis->{$m}{$n}{'N'})>$Ori_Dis->{$m}{$n}{'P'}){
								if ($Ori_Dis->{$m}{$n}{'N'}>=$minlink ){	
									 if((defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$m})&&(defined $gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$n})){
											if (($gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$m}*$gaa->{$derakh}{'cycle'}{$dayer}{$shakh}{'orientation'}{$n})==-1){ 
												$G->{$m}{$n}=-1;
												$G->{$n}{$m}=-1;
											}    
									 }else{
									      print "ajabaa $m n $n"; exit;
												$G->{$m}{$n}=-1;
												$G->{$n}{$m}=-1;								 
									 }
								}	
							}else{
	
							}
							
						 }
					 }		 
				 }
	   }		 
	}
	return $G;
} 
#------------------------------------------- 
sub reverseComplement{
   $_ = shift;
   $_ = uc();
   tr/ATGC/TACG/;
   return (reverse());
}
#--------------------------------------------
sub findDnascontig{      
	  my $num = shift;
	  open (INC, "<$address"."/contigs_sopra.fasta") or die "$!";
	  my $contig;
	  LINE: while (<INC>){
		if ($_=~ /^\>/){
		   my @string = split(/\|/,$_);
		   if (substr($string[0],7)==$num){
			   $contig=<INC>;
			   chomp($contig);
			   last LINE;
		   }
		}
	  }
	  close INC;
	  return $contig;
}
#--------------------------------------------
sub find_treeforopt{
    my ($inj,$cy) =  @_;  
    my $int;
    foreach my $m (keys %$inj){ 
		 foreach my $n (keys %{$inj->{$m}}){  
	         $int->{$m}{$n}=$inj->{$m}{$n};
         }	
	}			  
	my $derakht; my $gar;      	
	foreach my $m (keys %$int){ 
	   if (defined $int->{$m}){
		  $derakht++; my $con; push @{$gar->{$derakht}{'cycle'}{0}{1}{'member'}}, $m;   
		  $gar->{$derakht}{'cycle'}{0}{1}{'o'}{$m}=1;
		  my $Extension=1;
		  while ($Extension==1){
			 $Extension=0; 
			 foreach my $n (@{$gar->{$derakht}{'cycle'}{0}{1}{'member'}}){ 
				  foreach my $s (keys %{$int->{$n}}){  

				      $con->{$n}{$s}=$int->{$n}{$s};   
					  $con->{$s}{$n}=$int->{$n}{$s};   
					  
					  if (!defined $gar->{$derakht}{'cycle'}{0}{1}{'o'}{$s}){
						   $Extension=1;
						   push @{$gar->{$derakht}{'cycle'}{0}{1}{'member'}}, $s;  
						   $gar->{$derakht}{'cycle'}{0}{1}{'o'}{$s}=1;		
						   delete $int->{$n}{$s}; delete $int->{$s}{$n}
					  }else{

					  }
				  }
			 }
		  } 
		   foreach my $n (@{$gar->{$derakht}{'cycle'}{0}{1}{'member'}}){ 
			  delete $int->{$n};
		   }
		   @{$gar->{$derakht}{'cycle'}{0}{1}{'member'}}=sort { $a <=> $b } @{$gar->{$derakht}{'cycle'}{0}{1}{'member'}};

		   #my $du=scalar(@{$gar->{$derakht}{'cycle'}{0}{1}{'member'}});
		   #print "tree $derakht  num memeber $du\n";
		   
# 			  open (DT, ">d_startorcyc$cy"."_t$derakht"."_$du".".txt"); 
# 			  print DT "graph simple_hierarchy {\n"; 
# 			  print DT "   rankdir=LR;\n"; 
# 			  print DT "   fixedsize=true;\n";
# 			  print DT "   overlap=scale;\n";
# 			  print DT "   splines=true;\n"; 
# 			  print DT "   node[shape=box,fontsize=8,height=.05];\n"; 
# 			
# 			  foreach my $m (keys %$con){ 
# 					 my $w;
# 					 if ($tig_info->{$m}{'length'}<3000){
# 						  $w=.03+int($tig_info->{$m}{'length'}/100)*.03
# 					 }else{
# 						  $w=1.3;
# 					 }
# 					 print DT "   $m [width=$w"."]\n";
# 			  }					  
# 			  foreach my $m (keys %$con){        
# 					   foreach my $n (keys %{$con->{$m}}){ 
# 							  if( $m < $n){
# 								 print DT "   $m"."--$n". " [fontsize=8,label=$con->{$m}{$n}" ."]\n";	  
# 							  }
# 					   }
# 			  }		 
# 			  print DT "}";
# 			  close DT;	
	   }
	}
    return $gar;
}	
#--------------------------------------------
sub input_non_seq{
   my ($nfile,$OD) =  @_; 
   my $tot_n; my $temp_od; my $dis_prev; 
   if ( $strength < 1){
   		$strength=1;
   }
   print "Reading $nfile\n";
    open (IN,"<$nfile") or die "Cannot open $nfile $!\n";
	while (my $line= <IN>){
		chomp($line);       	    
		if ($line=~ /^(\d+)\t(\d+)\t([+-])\t(\d+)\t(\d+)\t([+-])\t(\d+)/){   
			my $tig_a=$1; my $A_start;	my $A_end; my $insert_size=$7;

			 my $dev_1; my $dev_2;
			 if($insert_size>5000){
				   $dev_1=1.1; $dev_2=.2;
			 }elsif($insert_size>1000){
				   $dev_1=1.3; $dev_2=.5;
			 }else{
				   $dev_1=1.6; $dev_2=1;
			 } 
			 
			if ( $3 eq "+" ){					
				$A_start=$2; $A_end=$A_start+20;
			}else{
				$A_end=$2; $A_start=$A_end+20;
			}				
			
			my $tig_b=$4; my $B_start;	my $B_end;
			
			if ( $6 eq "+" ){					
				$B_start=$5; $B_end=$B_start+20;
			}else{
				$B_end=$5; $B_start=$B_end+20;
			}				
			
			if( ($tig_info->{$tig_a}{'length'}>$minlength) && ($tig_info->{$tig_b}{'length'}>$minlength) && ($tig_a != $tig_b) ){
			
				 my $incorpor=1; my $m; my $n; my $x=1;
				
				 if( defined $OD->{$tig_a} ){
					 if( defined $OD->{$tig_a}{$tig_b} ){
							$m=$tig_a; $n=$tig_b;
					 }
				 }
				 if( defined $OD->{$tig_b} ){
					 if( defined $OD->{$tig_b}{$tig_a} ){
							$m=$tig_b; $n=$tig_a;
					 }
				 }		 		 
				 if ( ($m > 0) && ( $n > 0 ) ) {
					 if ( ( defined $OD->{$m}{$n}{'P'} ) && (! defined $OD->{$m}{$n}{'N'} ) ){
						if ($OD->{$m}{$n}{'P'}>=$minlink ){
								$x=int($OD->{$m}{$n}{'PD'}/$OD->{$m}{$n}{'P'}); $incorpor=0;							
								if ( $x==0 ){
									 if ($OD->{$m}{$n}{'PD'} > 0) {
										  $x=1;  
									 }else{
										  $x=-1; 
									 }
								}									
						}
					 }elsif ((! defined $OD->{$m}{$n}{'P'}) && (defined $OD->{$m}{$n}{'N'})){
						if ($OD->{$m}{$n}{'N'}>=$minlink ){ 		 
								$x=int($OD->{$m}{$n}{'ND'}/$OD->{$m}{$n}{'N'}); $incorpor=0;	
								if ( $x==0 ){
									 if ($OD->{$m}{$n}{'ND'} > 0) {
										  $x=1; 
									 }else{
										  $x=-1; 
									 }
								}		
						}		 
					 }elsif((defined $OD->{$m}{$n}{'P'}) && (defined $OD->{$m}{$n}{'N'})){
						if ((.4*$OD->{$m}{$n}{'P'})>$OD->{$m}{$n}{'N'}){
							if ($OD->{$m}{$n}{'P'}>=$minlink ){	
									$x=int($OD->{$m}{$n}{'PD'}/$OD->{$m}{$n}{'P'}); $incorpor=0;	
									if ( $x==0 ){
										 if ($OD->{$m}{$n}{'PD'} > 0) {
											  $x=1; 
										 }else{
											  $x=-1;  
										 }
									}													
							}		
						}elsif ((.4*$OD->{$m}{$n}{'N'})>$OD->{$m}{$n}{'P'}){
							if ($OD->{$m}{$n}{'N'}>=$minlink ){	
								$x=int($OD->{$m}{$n}{'ND'}/$OD->{$m}{$n}{'N'}); $incorpor=0;	
								if ( $x==0 ){
									 if ($OD->{$m}{$n}{'ND'} > 0) {
										  $x=1;  
									 }else{
										  $x=-1;  
									 }
								}													
							}	
						}
					 }					 
				 }else{
				 }
				 
				 if($tig_a != $tig_b){ ####paired reads located on <> contigs
						if ( $A_start < $A_end ){
							  if ( $B_start < $B_end ){#  -> -> 
								 my $d=$insert_size+$A_start-$B_start;					  
								 my $gap=$d-$tig_info->{$tig_a}{'length'};
								 if(abs($gap)<($dev_1*$insert_size)){
										 if ($tig_a<$tig_b){
											if ( ( $incorpor==1 ) || (  abs(($x-$d)/$x)< .1  ) ){
												$temp_od->{$tig_a}{$tig_b}{'P'}+=$strength;
												$temp_od->{$tig_a}{$tig_b}{'PD'}+=($strength*$d); 
											}else{
												$dis_prev++;
											}
										 }else{    
											if ( ( $incorpor==1 ) || (  abs(($x+$d)/$x)< .1  ) ){
												$temp_od->{$tig_b}{$tig_a}{'P'}+=$strength;
												$temp_od->{$tig_b}{$tig_a}{'PD'}-=($strength*$d); 
											}else{
												$dis_prev++;
											}
										 }									 
								  }else{
								  }
							  }else{# -> <-
								 my $d=$insert_size+$A_start+$B_start;															 
								 my $gap=$d-$tig_info->{$tig_a}{'length'}-$tig_info->{$tig_b}{'length'};
								 if(abs($gap)<($dev_1*$insert_size)){
										 if ($tig_a<$tig_b){
											if ( ( $incorpor==1 ) || (  abs(($x-$d)/$x)< .1  ) ){		
												$temp_od->{$tig_a}{$tig_b}{'N'}+=$strength;
												$temp_od->{$tig_a}{$tig_b}{'ND'}+=($strength*$d); 
											}else{
												$dis_prev++;
											}
										 }else{
											if ( ( $incorpor==1 ) || (  abs(($x-$d)/$x)< .1  ) ){
												$temp_od->{$tig_b}{$tig_a}{'N'}+=$strength;
												$temp_od->{$tig_b}{$tig_a}{'ND'}+=($strength*$d); 
											}else{
												$dis_prev++;
											}
										 }
								 }else{                             
								 }
							  }
						}else{
							  if ( $B_start > $B_end ){#  <-  <-
								 my $d=-($insert_size-$A_start+$B_start);						  
								 my $gap=-$d-$tig_info->{$tig_b}{'length'};
								 if(abs($gap)<($dev_1*$insert_size)){
										 if ($tig_a<$tig_b){
											if ( ( $incorpor==1 ) || (  abs(($x-$d)/$x)< .1  ) ){	
												$temp_od->{$tig_a}{$tig_b}{'P'}+=$strength;
												$temp_od->{$tig_a}{$tig_b}{'PD'}+=($strength*$d); 
											}else{
												$dis_prev++;
											}	
										 }else{
											if ( ( $incorpor==1 ) || (  abs(($x+$d)/$x)< .1  ) ){
												$temp_od->{$tig_b}{$tig_a}{'P'}+=$strength;
												$temp_od->{$tig_b}{$tig_a}{'PD'}-=($strength*$d); 
											}else{
												$dis_prev++;
											}	
										 }
								 }else{							 
								 }
							  }else{    #   <- ->  
								 my $d=-($insert_size-$A_start-$B_start);									 
								 my $gap=-$d;
								 if(abs($gap)<($dev_1*$insert_size)){
										 if ($tig_a<$tig_b){
											if ( ( $incorpor==1 ) || (  abs(($x-$d)/$x)< .1  ) ){
												$temp_od->{$tig_a}{$tig_b}{'N'}+=$strength;
												$temp_od->{$tig_a}{$tig_b}{'ND'}+=($strength*$d); 
											}else{
												$dis_prev++;
											}
										 }else{
											if ( ( $incorpor==1 ) || (  abs(($x-$d)/$x)< .1  ) ){
												$temp_od->{$tig_b}{$tig_a}{'N'}+=$strength;
												$temp_od->{$tig_b}{$tig_a}{'ND'}+=($strength*$d);  
											}else{
												$dis_prev++;
											}
										 }
								 }else{						 
								 }									 
							  }
						}
				 }else{ ###Clone, paired reads located on the same contig -- could be used to investigate misassemblies
				 } 								 
			}				
		}else{
			die "Invalid format $nfile :\n$line\n";	   
	   } 	
	} 		
	foreach my $m (keys %$temp_od){  
		foreach my $n (keys %{$temp_od->{$m}}){ 
			if (defined $temp_od->{$m}{$n}{'P'} ){
					$OD->{$m}{$n}{'P'}+=$temp_od->{$m}{$n}{'P'};
					$OD->{$m}{$n}{'PD'}+=$temp_od->{$m}{$n}{'PD'}; $tot_n++;
			}
			if (defined $temp_od->{$m}{$n}{'N'} ){
					$OD->{$m}{$n}{'N'}+=$temp_od->{$m}{$n}{'N'};
					$OD->{$m}{$n}{'ND'}+=$temp_od->{$m}{$n}{'ND'}; $tot_n++;
			}	
			
			
		}
	}
	if ( $dis_prev > 0 ){
		print "Number of discarded pairs due to conflict with short read data $dis_prev\n";
	}
	print "Total number of useful pairs $tot_n\n\n";	close IN;
    return $OD;
}	
#-------------------------------------------- did u register for the new class
sub gen_JforIsing{
   my ($OD,$paak,$ign_m) =  @_; 
   my $int; my $maxJ=0; my $numJ=0; my $sumJ=0;   

	foreach my $m (keys %$OD){  
		 if(($tig_info->{$m}{'length'}>$minlength)&&(!defined $paak->{$m})&&(!defined $ign_m->{$m})){
				 foreach my $n (keys %{$OD->{$m}}){ 
						  if(($tig_info->{$n}{'length'}>$minlength)&&(!defined $paak->{$n})&&(!defined $ign_m->{$n})){
								if ((defined $OD->{$m}{$n}{'P'}) && (! defined $OD->{$m}{$n}{'N'})){
									 if ($OD->{$m}{$n}{'P'}>=$minlink){
										 $int->{$m}{$n}=$OD->{$m}{$n}{'P'};
										 $int->{$n}{$m}=$OD->{$m}{$n}{'P'};
										  $numJ++; $sumJ+=$OD->{$m}{$n}{'P'};
										  if($maxJ<$OD->{$m}{$n}{'P'}){
										       $maxJ=$OD->{$m}{$n}{'P'};
										  }
									 }     
								}elsif ((! defined $OD->{$m}{$n}{'P'}) && (defined $OD->{$m}{$n}{'N'})){
									 if ($OD->{$m}{$n}{'N'}>=$minlink){
										 $int->{$m}{$n}=-$OD->{$m}{$n}{'N'};
										 $int->{$n}{$m}=-$OD->{$m}{$n}{'N'};
										 $numJ++; $sumJ+=$OD->{$m}{$n}{'N'};
										 if($maxJ<$OD->{$m}{$n}{'N'}){
										       $maxJ=$OD->{$m}{$n}{'N'};
										  }
									 }     
								}elsif((defined $OD->{$m}{$n}{'P'}) && (defined $OD->{$m}{$n}{'N'})){   
								   if ((.4*$OD->{$m}{$n}{'P'})>$OD->{$m}{$n}{'N'}){
									   if ($OD->{$m}{$n}{'P'}>=$minlink){  
										  $int->{$m}{$n}=$OD->{$m}{$n}{'P'};
										  $int->{$n}{$m}=$OD->{$m}{$n}{'P'};
										  $numJ++; $sumJ+=$OD->{$m}{$n}{'P'};
										  if($maxJ<$OD->{$m}{$n}{'P'}){
										       $maxJ=$OD->{$m}{$n}{'P'};
										  }
									   }  
								   }elsif ((.4*$OD->{$m}{$n}{'N'})>$OD->{$m}{$n}{'P'}){
									   if ($OD->{$m}{$n}{'N'}>=$minlink){   
										  $int->{$m}{$n}=-$OD->{$m}{$n}{'N'};
										  $int->{$n}{$m}=-$OD->{$m}{$n}{'N'};
										  $numJ++; $sumJ+=$OD->{$m}{$n}{'N'};
										  if($maxJ<$OD->{$m}{$n}{'N'}){
										       $maxJ=$OD->{$m}{$n}{'N'};
										  }
									   }   
								   }  
								}
						  }	
				 }
		 }		 
	} 
	my $avej=0;
	if ( $numJ > 0){
		 $avej=$sumJ/$numJ;  $avej=.01* int(100*$avej); 	
	}else{
	   print "\n\n\neee numJ $numJ\n\n";
	}
	#print "Average number of links between two contigs using minlength $minlength and minlink $minlink is $avej\n";    
    return $int,$avej,$numJ;
}	
#--------------------------------------------
sub get_ave{
   my ($OD,$paak,$ign_m) =  @_; 
   my $maxJ=0; my $numJ=0; my $sumJ=0;   
	foreach my $m (keys %$OD){  
		 if(($tig_info->{$m}{'length'}>$minlength)&&(!defined $paak->{$m})&&(!defined $ign_m->{$m})){
				 foreach my $n (keys %{$OD->{$m}}){ 
						  if(($tig_info->{$n}{'length'}>$minlength)&&(!defined $paak->{$n})&&(!defined $ign_m->{$n})){
								if ((defined $OD->{$m}{$n}{'P'}) && (! defined $OD->{$m}{$n}{'N'})){
									 if ($OD->{$m}{$n}{'P'}>= 2){
										  $numJ++; $sumJ+=$OD->{$m}{$n}{'P'};
										  if($maxJ<$OD->{$m}{$n}{'P'}){
										       $maxJ=$OD->{$m}{$n}{'P'};
										  }
									 }     
								}elsif ((! defined $OD->{$m}{$n}{'P'}) && (defined $OD->{$m}{$n}{'N'})){
									 if ($OD->{$m}{$n}{'N'}>= 2 ){
										 $numJ++; $sumJ+=$OD->{$m}{$n}{'N'};
										 if($maxJ<$OD->{$m}{$n}{'N'}){
										       $maxJ=$OD->{$m}{$n}{'N'};
										  }
									 }     
								}elsif((defined $OD->{$m}{$n}{'P'}) && (defined $OD->{$m}{$n}{'N'})){   
								   if ((.4*$OD->{$m}{$n}{'P'})>$OD->{$m}{$n}{'N'}){
									   if ($OD->{$m}{$n}{'P'}>= 2){  
										  $numJ++; $sumJ+=$OD->{$m}{$n}{'P'};
										  if($maxJ<$OD->{$m}{$n}{'P'}){
										       $maxJ=$OD->{$m}{$n}{'P'};
										  }
									   }  
								   }elsif ((.4*$OD->{$m}{$n}{'N'})>$OD->{$m}{$n}{'P'}){
									   if ($OD->{$m}{$n}{'N'}>= 2){   
										  $numJ++; $sumJ+=$OD->{$m}{$n}{'N'};
										  if($maxJ<$OD->{$m}{$n}{'N'}){
										       $maxJ=$OD->{$m}{$n}{'N'};
										  }
									   }   
								   }  
								}
						  }	
				 }
		 }		 
	} 
	my $avej=0;
	if ( $numJ > 0 ){
     	$avej=$sumJ/$numJ; $avej=.01* int(100*$avej); 	
	}
	#print "Average number of links between two contigs using minlength $minlength and minlink $minlink is $avej\n";    
    return $avej,$numJ,$maxJ;
}	
#--------------------------------------------
sub gen_Jforblockfrustration{
   my ($OD,$paak,$ign_m,$gar,$derakht,$tabaghe,$shakhe,$articul,$cy) =  @_; 
	my ($int,$avej,$du)=gen_JforIsing($OD,$paak,$ign_m);
	foreach my $n (keys %{$int}){    
		 if (!defined $gar->{$derakht}{'cycle'}{$tabaghe}{$shakhe}{'o'}{$n}){
			  delete $int->{$n}; 
		 }
	}   
	foreach my $n (keys %{$int}){    
		 foreach my $k (keys %{$int->{$n}}){    
				if (! defined $gar->{$derakht}{'cycle'}{$tabaghe}{$shakhe}{'o'}{$k}){
					 if ( defined $articul->{$n} ){
						 delete $int->{$n}{$k};  
					 }else{
						  print "$n $k $derakht floor $tabaghe $shakhe\n"; exit;
					 }
				}
		 }
	}      
    return $int;
}	
#--------------------------------------------
sub assigntemporaryorient{
    my ($inj,$cy,$derakht,$tabaghe,$shakhe) = @_;
    my $int;
    foreach my $m (keys %$inj){ 
		 foreach my $n (keys %{$inj->{$m}}){  
	         $int->{$m}{$n}=$inj->{$m}{$n};
         }	
	}	
    my ($B,$bas,$con);
	foreach my $k (keys %$int){ 
		   if (defined $int->{$k}){
			  $bas++; push @{$B->{$bas}{'member'}}, $k;   
			  $B->{$bas}{'o'}{$k}=1; 
			  my $Extension=1;
			  while ($Extension==1){
				 $Extension=0; 
				 foreach my $n (@{$B->{$bas}{'member'}}){ 
					  foreach my $s (keys %{$int->{$n}}){  
					      #$con->{$n}{$s}=$int->{$n}{$s};   
					      #$con->{$s}{$n}=$int->{$n}{$s};   
						  if (!defined $B->{$bas}{'o'}{$s}){
							   $Extension=1;
							   push @{$B->{$bas}{'member'}}, $s;  
							   if ( $int->{$n}{$s}>0 ){
									$B->{$bas}{'o'}{$s}= $B->{$bas}{'o'}{$n}; 
							   }else{
									 $B->{$bas}{'o'}{$s}=-$B->{$bas}{'o'}{$n}; 
							   }
							   delete $int->{$n}{$s}; delete $int->{$s}{$n}
						  }else{
						  }
					  }
				 }
			  } 
			  foreach my $n (@{$B->{$bas}{'member'}}){ 
				  delete $int->{$n};
			  }
			  @{$B->{$bas}{'member'}}=sort { $a <=> $b } @{$B->{$bas}{'member'}};
		   }
	}								
	if ( $bas >1 ){
	     print "basket $bas (should be 1 )\n"; exit;
	}
	
#	  my $du=scalar(@{$B->{1}{'member'}});
# 	  open (DT, ">d_blockfrustration$cy"."_t$derakht"."_f$tabaghe"."_s$shakhe"."_L$du".".txt"); 
# 	  print DT "graph simple_hierarchy {\n"; 
# 	  print DT "   rankdir=LR;\n"; 
# 	  print DT "   fixedsize=true;\n";
# 	  print DT "   overlap=scale;\n";
# 	  print DT "   splines=true;\n"; 
# 	  print DT "   node[shape=box,fontsize=8,height=.05];\n"; 
# 	
# 	  foreach my $m (keys %$con){ 
# 			 my $w;
# 			 if ($tig_info->{$m}{'length'}<3000){
# 				  $w=.03+int($tig_info->{$m}{'length'}/100)*.03
# 			 }else{
# 				  $w=1.3;
# 			 }
# 			 if ( $B->{1}{'o'}{$m} == 1){
# 			     print DT "   $m [style=filled, width=$w".",color=deepskyblue1"."]\n";
# 			 }elsif( $B->{1}{'o'}{$m} == -1){
# 			     print DT "   $m [style=filled, width=$w".",color=green"."]\n";
# 			 }else{
# 			     print "ori assign nist\n"; exit;
# 			 }
# 	  }					  
# 	  foreach my $m (keys %$con){        
# 			   foreach my $n (keys %{$con->{$m}}){ 
# 					  if( $m < $n){
# 						  if ( $con->{$n}{$m}>0 ){
# 							   if ( ($B->{1}{'o'}{$m}*$B->{1}{'o'}{$n}) == 1 ){
# 							       print DT "   $m"."--$n". " [fontsize=8,label=$con->{$m}{$n}" ."]\n"; 
# 							   }else{
# 								   print DT "   $m"."--$n". " [fontsize=8, color=red, label=$con->{$m}{$n}" ."]\n";
# 							   }
# 						  }else{
# 							   if ( ($B->{1}{'o'}{$m}*$B->{1}{'o'}{$n}) == -1 ){
# 							       print DT "   $m"."--$n". " [fontsize=8,label=$con->{$m}{$n}" ."]\n";
# 							   }else{
# 								   print DT "   $m"."--$n". " [color=red, fontsize=8, label=$con->{$m}{$n}" ."]\n"; 
# 							   }              
# 						  }	
# 					  }
# 			   }
# 	  }		 
# 	  print DT "}";    
# 	  close DT;	 
    return $B;
}	
#--------------------------------------------
sub test_frustration{
    my ($B,$int) = @_; 
    my $frus=0; my $e_fr=0;
    foreach my $n (keys %$int){ 
         if (!defined $B->{1}{'o'}{$n}){
             print "chera too barray nist $n\n"; exit;
         }
         foreach my $s (keys %{$int->{$n}} ){
			  if (!defined $B->{1}{'o'}{$s}){
				   print "chera too barray nist $s\n"; exit;
			  }
              if ( $int->{$n}{$s}>0 ){
                   if ( ($B->{1}{'o'}{$s}*$B->{1}{'o'}{$n}) == 1 ){
                   
                   }else{
                       $frus++;  $e_fr-= $int->{$n}{$s}*$B->{1}{'o'}{$s}*$B->{1}{'o'}{$n};
                       #print "frust $n $s    $int->{$n}{$s}\n";
                   }
              }else{
                   if ( ($B->{1}{'o'}{$s}*$B->{1}{'o'}{$n}) == -1 ){
                   
                   }else{
                       $frus++;  $e_fr-= $int->{$n}{$s}*$B->{1}{'o'}{$s}*$B->{1}{'o'}{$n};
                       #print "frust $n $s    $int->{$n}{$s}\n"
                   }              
              }
         
         }
   } 
   $frus*=.5; $e_fr*=.5;
   return $frus,$e_fr;
}	
#--------------------------------------------
sub find_Zstructure{
   my $inj =  shift;
   my $int;
    foreach my $m (keys %$inj){ 
		 foreach my $n (keys %{$inj->{$m}}){  
	         $int->{$m}{$n}=$inj->{$m}{$n};
         }	
	}
   #open(MT, ">conect.txt");
   my $layer=0; my $ze;
	my $bas; my $gar;
	foreach my $m (keys %$int){ 
	   if (defined $int->{$m}){ 
		  $ze->{$layer}{'o'}{$m}=1;  push @{$ze->{$layer}{'mem'}}, $m;  
		  $bas++; push @{$gar->{$bas}{'member'}}, $m;  $gar->{$bas}{'o'}{$m}=1;
		  my $Extension=1;
		  while ($Extension==1){
			 $Extension=0;  $layer++;
			 my @aa=@{$gar->{$bas}{'member'}};
			 foreach my $n (@aa){ 
				  foreach my $s (keys %{$int->{$n}}){  
				      #print MT "$n $s 1\n"; print MT "$s $n 1\n";
				      #print "n $n      s $s\n"; 
					  if (!defined $gar->{$bas}{'o'}{$s}){
					       $ze->{$layer}{'o'}{$s} = $int->{$n}{$s};   
					       push @{$ze->{$layer}{'mem'}}, $s;  
						   $Extension=1;
						   push @{$gar->{$bas}{'member'}}, $s;  
						   $gar->{$bas}{'o'}{$s}=1;		
						   delete $int->{$n}{$s}; delete $int->{$s}{$n};
					  }else{
                            delete $int->{$n}{$s}; delete $int->{$s}{$n};
					  }
				  }
				  delete $int->{$n};
			 }
			 if($Extension==1){
				 # print "Zlayer $layer: @{$ze->{$layer}{'mem'}}\n";
			 }	 
		  } 
		   foreach my $n (@{$gar->{$bas}{'member'}}){ 
			  delete $int->{$n};
		   }
	   }
	}	
	#close MT;
	if ( $bas==1){
	
	}else{
	    print "bas? $bas\n"; exit;
	}
	my $bigg=0;  
	foreach my $lay (keys %{$ze} ){     
	    my $du=scalar( keys %{$ze->{$lay}{'o'}} );
	    if ( $du> $Z_combthre ){
	           if ( $bigg < $du ){
	                  $bigg= $du; 
	           }
	           #print "du $du bigg $bigg\n";
	    }
	}
    return $ze,$bigg;
}	
#--------------------------------------------
sub draw_graph_1{
  my ($con,$gar,$cy,$derakht,$tabaghe,$shakhe) =  @_; 
  my $du=scalar(keys %{$con});
  open (DT, ">d_blockfinalori_c$cy"."_t$derakht"."_f$tabaghe"."_s$shakhe"."_L$du".".txt"); 
  print DT "graph simple_hierarchy {\n"; 
  print DT "   rankdir=LR;\n"; 
  print DT "   fixedsize=true;\n";
  print DT "   overlap=scale;\n";
  print DT "   splines=true;\n"; 
  print DT "   node[shape=box,fontsize=8,height=.05];\n"; 

  foreach my $m (keys %$con){ 
		 my $w;
		 if ($tig_info->{$m}{'length'}<3000){
			  $w=.03+int($tig_info->{$m}{'length'}/100)*.03
		 }else{
			  $w=1.3;
		 }
		 if ( $gar->{$derakht}{'cycle'}{$tabaghe}{$shakhe}{'orient'}{$m} == 1){
			 print DT "   $m [style=filled, width=$w".",color=deepskyblue1"."]\n";
		 }elsif( $gar->{$derakht}{'cycle'}{$tabaghe}{$shakhe}{'orient'}{$m}  == -1){
			 print DT "   $m [style=filled, width=$w".",color=green"."]\n";
		 }else{
			 print "ori assign nist t $derakht floor $tabaghe sub $shakhe\n"; exit;
		 }
  }					  
  foreach my $m (keys %$con){        
		   foreach my $n (keys %{$con->{$m}}){ 
				  if( $m < $n){
					  if ( $con->{$n}{$m}>0 ){
						   if ( ($gar->{$derakht}{'cycle'}{$tabaghe}{$shakhe}{'orient'} {$m}*$gar->{$derakht}{'cycle'}{$tabaghe}{$shakhe}{'orient'} {$n}) == 1 ){
							   print DT "   $m"."--$n". " [fontsize=8,label=$con->{$m}{$n}" ."]\n"; 
						   }else{
							   print DT "   $m"."--$n". " [color=red, fontsize=8, label=$con->{$m}{$n}" ."]\n";
						   }
					  }else{
						   if ( ($gar->{$derakht}{'cycle'}{$tabaghe}{$shakhe}{'orient'} {$m}*$gar->{$derakht}{'cycle'}{$tabaghe}{$shakhe}{'orient'} {$n})  == -1 ){
							   print DT "   $m"."--$n". " [fontsize=8, label=$con->{$m}{$n}" ."]\n";
						   }else{
							   print DT "   $m"."--$n". " [fontsize=8,color=red, label=$con->{$m}{$n}" ."]\n"; 
						   }              
					  }	
				  }
		   }
  }		 
  print DT "}";    
  close DT;	
	  
}	
#--------------------------------------------
sub simulated_annealing{
    my ($inj,$gar,$derakht,$tabaghe,$shakhe)= @_; 
    my $j;
    foreach my $m (keys %$inj){ 
		 foreach my $n (keys %{$inj->{$m}}){  
	         $j->{$m}{$n}=$inj->{$m}{$n};
         }	
	}	
    my $J; my $L=0; my $Num; my $minE=0; my @B; my @S; my $E;
    
    #print "sim anneal: tree $derakht floor $tabaghe sub $shakhe\n";    
    
    my $numJ=0; my $sumJ=0;   
    foreach my $m (keys %$j){  
			foreach my $n (keys %{$j->{$m}}){ 
			      $numJ++; $sumJ+=abs($j->{$m}{$n});  
			}
    } 
    my $avj=$sumJ/$numJ; $avj=int($avj); 
    #print "Average links t $derakht f $tabaghe s $shakhe is $avj\n";    
    
    if($avj>20){ 
		foreach my $m (keys %$j){  
			foreach my $n (keys %{$j->{$m}}){ 
				 $j->{$m}{$n}*=(20/$avj);
			}
		}  
	}	
	foreach my $m (keys %$j){  
		 if(! defined $Num->{$m}){
			$Num->{$m}=$L;
			$S[$L]=1;
			$L++;
		 }
		 foreach my $n (keys %{$j->{$m}}){  
			 if(! defined $Num->{$n}){
				$Num->{$n}=$L;
				$S[$L]=1;
				$L++;
			 }
			 $J->{$Num->{$m}}{$Num->{$n}}=$j->{$m}{$n};
		 }  
	}
	undef $j; #print "L $L\n";
	foreach my $m (keys %$J){  
		 $B[$m]=0;
		 foreach my $n (keys %{$J->{$m}}){ 
			$B[$m]+= $J->{$m}{$n}*$S[$n];
			$minE-=abs($J->{$m}{$n});
		 }
		 $E+=- $B[$m]*$S[$m];   
	}
	$E=$E/2;
	$minE=$minE/2; my $cut_sat=40+int(-$minE/4000);  
	$minE=int($minE)-1; print "Min energy $minE\n";  
	#open (EN, ">Energy_c$cyc"."_minE$minE".".txt");

    my $du=int($E); 
	#print EN "$du ";
	my $T=.01; my $time=0; my $saturated=0; my @Memo_E;
	for (my $i=0; $i<20; $i++){
	     $Memo_E[$i]=$E;
	}
	#my $T_step=3*(10**(-10))/(.75+$L/2100);	
	#my $es=2000*(.75+.1 * int($L/300) ); print "2000 es $es L $L\n";
     
     my $es=500* (1 +.01* int($L/5) ); #print "500 es $es L $L\n";
     
	while($saturated==0){   
        if (($time %  $es)==0){
              $T*= 1.0001;
         }
		 $time++;
		 my $change=0;
		 my $i = int(rand($L));
		 my $element = $S[$i];
		 my $delta=$B[$i]*$S[$i];
		 if ($delta>0){   # need to check wheter or not change $S[i]          
			 if (rand()<exp(-$T*$delta)){  #will change $S[i] 
				$change=1;
			 }  
		 }else{  # will change $S[i]       
			 $change=1;
		 }
		 if ($change){
		   $E+=2*$delta;
		   $S[$i]*=-1;
		   foreach my $n (keys %{$J->{$i}}){ # update other B 
			$B[$n]+= 2*$J->{$i}{$n}*$S[$i];
			}      
		 }
		 if (($time % 2000000)==0){
		    my $du=int($E);  
			#print EN "$du ";
			push @Memo_E, $E;
			shift(@Memo_E);
			if ($time > 42000000){
				$saturated=1;
				for (my $i=0; $i < 16; $i++){
					 if( (($Memo_E[$i]-$E)>$cut_sat) || (($Memo_E[$i]-$E)<0) ){ 
						  $saturated=0;
					 }
				}
				if(abs($minE-$E)<.001){
				    $saturated=1;  
				}
		    }		
		 }
		 if (($time % 50000000)==0){
		    my $du=int($E);
			my $date = `date`;	chomp($date);	print "Current energy $du     $date\n";
		 }	
	}
	#----------------------------------
	my $flip; 
	foreach my $i (keys %{$J}){
	   if(($B[$i]*$S[$i])<0){
		   $flip++;
		   $E+=2*$B[$i]*$S[$i];
		   $S[$i]*=-1;
		   foreach my $n (keys %{$J->{$i}}){ # update other B 
			$B[$n]+= 2*$J->{$i}{$n}*$S[$i];
			} 
	   }
	}
	#print EN "$E"; close EN; 
	undef $J; $E=int($E); print "Final energy $E\n";  
    my $s_ori;  
	foreach my $m (keys %$Num){       
	     $gar->{$derakht}{'cycle'}{$tabaghe}{$shakhe}{'orient'}{$m}=$S[$Num->{$m}];
	}  
    return $gar;
}	
#--------------------------------------------
sub gen_Jforspring{
   my ($OD,$s_ori) =  @_; 
   my $int; my $rededge;
	foreach my $m (keys %$OD){ 
	   if(defined $s_ori->{$m}){
				 foreach my $n (keys %{$OD->{$m}}){ 
					 if(defined $s_ori->{$n}){
						 if ((defined $OD->{$m}{$n}{'P'}) && (! defined $OD->{$m}{$n}{'N'})){
							if ($OD->{$m}{$n}{'P'}>=$minlink ){
								 if (($s_ori->{$m}*$s_ori->{$n})==1){        
									$int->{$m}{$n}=$OD->{$m}{$n}{'P'};
									$int->{$n}{$m}=$OD->{$m}{$n}{'P'};
								 }else{
	                                $rededge->{$m}{$n}=$Ori_Dis->{$m}{$n}{'P'}; 
								 }
							}
						 }elsif ((! defined $OD->{$m}{$n}{'P'}) && (defined $OD->{$m}{$n}{'N'})){
							if ($OD->{$m}{$n}{'N'}>=$minlink ){ 		 
								 if (($s_ori->{$m}*$s_ori->{$n})==-1){      
									$int->{$m}{$n}=$OD->{$m}{$n}{'N'};
									$int->{$n}{$m}=$OD->{$m}{$n}{'N'};
								 }else{
	           						$rededge->{$m}{$n}=-$Ori_Dis->{$m}{$n}{'N'};	
								 }
							}		 
						 }elsif((defined $OD->{$m}{$n}{'P'}) && (defined $OD->{$m}{$n}{'N'})){
							if ((.4*$OD->{$m}{$n}{'P'})>$OD->{$m}{$n}{'N'}){
								if ($OD->{$m}{$n}{'P'}>=$minlink ){	
									if (($s_ori->{$m}*$s_ori->{$n})==1){              
										$int->{$m}{$n}=abs($OD->{$m}{$n}{'P'}-$OD->{$m}{$n}{'N'});
										$int->{$n}{$m}=abs($OD->{$m}{$n}{'P'}-$OD->{$m}{$n}{'N'});
									}else{
									    $rededge->{$m}{$n}=abs($Ori_Dis->{$m}{$n}{'P'}-$Ori_Dis->{$m}{$n}{'N'});
									}
								}		
							}elsif ((.4*$OD->{$m}{$n}{'N'})>$OD->{$m}{$n}{'P'}){
								if ($OD->{$m}{$n}{'N'}>=$minlink ){	
									if (($s_ori->{$m}*$s_ori->{$n})==-1){               
										$int->{$m}{$n}=abs($OD->{$m}{$n}{'P'}-$OD->{$m}{$n}{'N'});
										$int->{$n}{$m}=abs($OD->{$m}{$n}{'P'}-$OD->{$m}{$n}{'N'});
									}else{
									    $rededge->{$m}{$n}=-abs($Ori_Dis->{$m}{$n}{'P'}-$Ori_Dis->{$m}{$n}{'N'});
									}     
								}	
							}else{
	
							}
							
						 }
					 }		 
				 }
	   }		 
	} 
    return $int,$rededge;
}	
#--------------------------------------------
sub find_treeforspring{
    my ($inj,$cy,$rededge,$s_ori) =  @_;    
    my $int;
    foreach my $m (keys %$inj){ 
		 foreach my $n (keys %{$inj->{$m}}){  
	         $int->{$m}{$n}=$inj->{$m}{$n};
         }	
	}	
	my $derakht; my $gar;  my $con;    	
	foreach my $m (keys %$int){ 
	   if (defined $int->{$m}){
		  $derakht++; my $con; push @{$gar->{$derakht}{'cycle'}{0}{1}{'member'}}, $m;   
		  $gar->{$derakht}{'cycle'}{0}{1}{'o'}{$m}=1;
		  my $Extension=1;
		  while ($Extension==1){
			 $Extension=0; 
			 foreach my $n (@{$gar->{$derakht}{'cycle'}{0}{1}{'member'}}){ 
				  foreach my $s (keys %{$int->{$n}}){  

				      $con->{$n}{$s}=$int->{$n}{$s};   
					  $con->{$s}{$n}=$int->{$n}{$s};   
					  
					  if (!defined $gar->{$derakht}{'cycle'}{0}{1}{'o'}{$s}){
						   $Extension=1;
						   push @{$gar->{$derakht}{'cycle'}{0}{1}{'member'}}, $s;  
						   $gar->{$derakht}{'cycle'}{0}{1}{'o'}{$s}=1;		
						   delete $int->{$n}{$s}; delete $int->{$s}{$n}
					  }else{

					  }
				  }
			 }
		  } 
		   my $du=scalar(@{$gar->{$derakht}{'cycle'}{0}{1}{'member'}});
		   #print "tree $derakht  num memeber $du\n";
		   foreach my $n (@{$gar->{$derakht}{'cycle'}{0}{1}{'member'}}){ 
			  delete $int->{$n};
		   }
		   
		   @{$gar->{$derakht}{'cycle'}{0}{1}{'member'}}=sort { $a <=> $b } @{$gar->{$derakht}{'cycle'}{0}{1}{'member'}};
# 		
# 			  open (DT, ">d_oricyc$cy"."_t$derakht"."_red_L$du".".txt"); 
# 			  print DT "graph simple_hierarchy {\n"; 
# 			  print DT "   rankdir=LR;\n"; 
# 			  print DT "   fixedsize=true;\n";
# 			  print DT "   overlap=scale;\n";
# 			  print DT "   splines=true;\n"; 
# 			  print DT "   node[shape=box,fontsize=8,height=.05];\n"; 
# 			  
# 			  	  foreach my $m (keys %$con){ 
# 						 my $w;
# 						 if ($tig_info->{$m}{'length'}<3000){
# 							  $w=.03+int($tig_info->{$m}{'length'}/100)*.03
# 						 }else{
# 							  $w=1.3;
# 						 }
# 						 if ( $s_ori->{$m} == 1){
# 							 print DT "   $m [style=filled, width=$w".",color=orange"."]\n";
# 						 }elsif( $s_ori->{$m} == -1){
# 							 print DT "   $m [style=filled, width=$w".",color=green"."]\n";
# 						 }else{
# 							 print "ori assign nist inja t $derakht m $m\n"; exit;
# 						 }
# 				  }		
# 	  
# 			  foreach my $m (keys %$con){                 #black edge between mates
# 					   foreach my $n (keys %{$con->{$m}}){ 
# 							  if( $m < $n){
# 								 print DT "   $m"."--$n". " [fontsize=8,label=$con->{$m}{$n}"."]\n";
# 							  }		   
# 					   }
# 			  }		  
# 			  foreach my $m (keys %$rededge){                 #red edge between mates
# 			      if (defined $con->{$m}){
# 						   foreach my $n (keys %{$rededge->{$m}}){ 
# 						       if (defined $con->{$n}){
# 								  if( $m < $n){
# 									   print DT  "   $m"."--$n". " [fontsize=8,color=red, label=$rededge->{$m}{$n}" ."]\n";	
# 								  }		   
# 								}else{
# 								     print DT  "   $m"."--$n". " [fontsize=8,color=red, label=$rededge->{$m}{$n}" ."]\n";	
# 								     print DT  "   $n [style=filled".",color=deepskyblue1"."]\n";
# 								}
# 						   }
# 				  } 	   
# 			  }				 		  
#     		  print DT "}";
# 	    	  close DT;			  
	   }
	}
    return $gar;
}	
#--------------------------------------------
sub linear_eq{
   my ($gar,$derakht,$fan,$int) =  @_; 
	my $A; my $b; my $C; my $r; my $z; my $d;
	foreach my $m (@{$gar->{$derakht}{'cycle'}{0}{1}{'member'}}){ 
		 foreach my $n (keys %{$int->{$m}}){ 
			  $A->{$m}{$n}=-$int->{$m}{$n};
			  $A->{$m}{$m}+=$int->{$m}{$n};
			  if ($m<$n){
				  $b->{$m}-=$int->{$m}{$n}*$fan->{$m}{$n};
			  }else{
				  $b->{$m}+=$int->{$m}{$n}*$fan->{$m}{$n}; 
			  }
		 }  
	}
	my $m=${$gar->{$derakht}{'cycle'}{0}{1}{'member'}}[0];
	$gar->{$derakht}{'cycle'}{0}{1}{'position'}{$m}=0;
	delete $b->{$m}; delete $A->{$m};
	foreach my $n (keys %{$A}){ 
	   if (defined  $A->{$n}{$m}){
		   delete $A->{$n}{$m};
	   }
	}
	undef $m;
	my $e0=0;
	foreach my $m (keys %{$b}){
		  $e0+=$b->{$m}**2;
	}
	$e0=$e0**.5;
	foreach my $m (keys %{$A}){ 
		$C->{$m}{$m}=1/$A->{$m}{$m};
		$gar->{$derakht}{'cycle'}{0}{1}{'position'}{$m}=0;
		$r->{$m}=$b->{$m};
		$z->{$m}=$C->{$m}{$m}*$r->{$m};
		$d->{$m}=$z->{$m};
	}
	if ( $e0 != 0 ){
	    my $dafe=scalar(@{$gar->{$derakht}{'cycle'}{0}{1}{'member'}});  
	    if( $dafe < 6000){
	         $dafe=6000;
	    }
		ITER: for (my $k=0; $k <= $dafe; $k++){  
			my $alpha=0;
			foreach my $m (keys %{$z}){ 
			   $alpha+=$z->{$m}*$r->{$m};
			}
			my $forbeta=$alpha;
			my $du=0;  #update alpha
			foreach my $m (keys %{$A}){        
				  my $du2=0; 
				  foreach my $n (keys %{$A->{$m}}){ 
					   $du2+=$A->{$m}{$n}*$d->{$n};             
				  } 
				  $du+=$d->{$m}*$du2;
			}  
			if($du==0){
				foreach my $m (keys %{$A}){        
					  my $du2=0; 
					  foreach my $n (keys %{$A->{$m}}){ 
						   $du2+=$A->{$m}{$n}*$d->{$n}; 
					  } 
					  $du+=$d->{$m}*$du2;
				}  
			}
			$alpha/=$du; undef $du;
			foreach my $m (keys %{$gar->{$derakht}{'cycle'}{0}{1}{'position'}}){  
				 if ( defined $d->{$m} ){
				 	$gar->{$derakht}{'cycle'}{0}{1}{'position'}{$m}+= ($alpha*$d->{$m});
				 }
			}
			foreach my $m (keys %{$r}){    #update r
				  my $du=0; 
				  foreach my $n (keys %{$A->{$m}}){ 
					   $du+=$A->{$m}{$n}*$d->{$n};             
				  }   
				  $r->{$m}-=$alpha*$du;
			}
			foreach my $m (keys %{$z}){    #update z
				  $z->{$m}=$C->{$m}{$m}*$r->{$m};
			}
			my $beta=0; 
			foreach my $m (keys %{$z}){     #update beta
			   $beta+=$z->{$m}*$r->{$m}; 
			}
			$beta/=$forbeta;
			foreach my $m (keys %{$d}){     #update d
				 $d->{$m}=$z->{$m}+$beta*$d->{$m};
			}
			my $res; my $e=0;
			foreach my $m (keys %{$b}){     
				  my $du=0; 
				  foreach my $n (keys %{$A->{$m}}){ 
					   $du+=$A->{$m}{$n}*$gar->{$derakht}{'cycle'}{0}{1}{'position'}{$n};             
				  }   
				  $res->{$m}=$b->{$m}-$du;
				  $e+=$res->{$m}**2;
			}
			if (($k % 2000)==0){
				if ( $k > 0) {  
					print "iteration $k   e $e e0 $e0\n"; $|=1; ###clear buffer
				}
			}
			$e=$e**.5;
			if(($e/$e0)<10**(-11)){
				last ITER;   
			}    
		}   
	}else{
	     foreach my $m (@{$gar->{$derakht}{'cycle'}{0}{1}{'member'}}){ 
      		print "$m $gar->{$derakht}{'cycle'}{0}{1}{'position'}{$m}\n";
         }		
	}
    return $gar;
}	
#--------------------------------------------
sub gen_fanar{
   my ($OD,$s_ori) = @_; 
   my ($int,$fan);
	foreach my $m (keys %$OD){ 
	   if(defined $s_ori->{$m}){
				 foreach my $n (keys %{$OD->{$m}}){ 
					 if(defined $s_ori->{$n}){
						 if ((defined $OD->{$m}{$n}{'P'}) && (! defined $OD->{$m}{$n}{'N'})){
							if ($OD->{$m}{$n}{'P'}>=$minlink ){
								 if (($s_ori->{$m}*$s_ori->{$n})==1){        
									$int->{$m}{$n}=$OD->{$m}{$n}{'P'};
									$int->{$n}{$m}=$OD->{$m}{$n}{'P'};
									my $x=int($OD->{$m}{$n}{'PD'}/$OD->{$m}{$n}{'P'}); 
									if ( $x==0 ){
									     if ($OD->{$m}{$n}{'PD'} > 0) {
									          $x=1;  #print "m $m n $n , x $x\n"; 
									     }else{
									          $x=-1;  #print "m $m n $n , x $x\n"; 
									     }
									}
									if ($s_ori->{$m}==1){
										$fan->{$m}{$n}=$x;  
										$fan->{$n}{$m}=$x;  
									}else{
										$fan->{$m}{$n}=-$x;
										$fan->{$n}{$m}=-$x;
									}
								 }else{
	
								 }
							}
						 }elsif ((! defined $OD->{$m}{$n}{'P'}) && (defined $OD->{$m}{$n}{'N'})){
							if ($OD->{$m}{$n}{'N'}>=$minlink ){ 		 
								 if (($s_ori->{$m}*$s_ori->{$n})==-1){      
									$int->{$m}{$n}=$OD->{$m}{$n}{'N'};
									$int->{$n}{$m}=$OD->{$m}{$n}{'N'};
									my $x=int($OD->{$m}{$n}{'ND'}/$OD->{$m}{$n}{'N'});
									if ( $x==0 ){
									     if ($OD->{$m}{$n}{'ND'} > 0) {
									          $x=1;  #print "m $m n $n , x $x\n"; 
									     }else{
									          $x=-1;  #print "m $m n $n , x $x\n"; 
									     }
									}									
									if ($s_ori->{$m}==1){
										$fan->{$m}{$n}=$x; 
										$fan->{$n}{$m}=$x;  
									}else{
										$fan->{$m}{$n}=-$x;
										$fan->{$n}{$m}=-$x;
									}
								 }else{
	
								 }
							}		 
						 }elsif((defined $OD->{$m}{$n}{'P'}) && (defined $OD->{$m}{$n}{'N'})){
							if ((.4*$OD->{$m}{$n}{'P'})>$OD->{$m}{$n}{'N'}){
								if ($OD->{$m}{$n}{'P'}>=$minlink ){	
									if (($s_ori->{$m}*$s_ori->{$n})==1){              
										$int->{$m}{$n}=abs($OD->{$m}{$n}{'P'}-$OD->{$m}{$n}{'N'});
										$int->{$n}{$m}=abs($OD->{$m}{$n}{'P'}-$OD->{$m}{$n}{'N'});
									    my $x=int($OD->{$m}{$n}{'PD'}/$OD->{$m}{$n}{'P'}); 
										if ( $x==0 ){
											 if ($OD->{$m}{$n}{'PD'} > 0) {
												  $x=1; #print "m $m n $n , x $x\n"; 
											 }else{
												  $x=-1; #print "m $m n $n , x $x\n"; 
											 }
										}									    
										if ($s_ori->{$m}==1){
										   $fan->{$m}{$n}=$x;  
										   $fan->{$n}{$m}=$x;  
										}else{
										   $fan->{$m}{$n}=-$x;
										   $fan->{$n}{$m}=-$x;
										} 
									}
								}		
							}elsif ((.4*$OD->{$m}{$n}{'N'})>$OD->{$m}{$n}{'P'}){
								if ($OD->{$m}{$n}{'N'}>=$minlink ){	
									if (($s_ori->{$m}*$s_ori->{$n})==-1){               
										$int->{$m}{$n}=abs($OD->{$m}{$n}{'P'}-$OD->{$m}{$n}{'N'});
										$int->{$n}{$m}=abs($OD->{$m}{$n}{'P'}-$OD->{$m}{$n}{'N'});
										my $x=int($OD->{$m}{$n}{'ND'}/$OD->{$m}{$n}{'N'});
										if ( $x==0 ){
											 if ($OD->{$m}{$n}{'ND'} > 0) {
												  $x=1;  #print "m $m n $n , x $x\n"; 
											 }else{
												  $x=-1;  #print "m $m n $n , x $x\n"; 
											 }
										}									    
										if ($s_ori->{$m}==1){
										   $fan->{$m}{$n}=$x;  
										   $fan->{$n}{$m}=$x;  
										}else{
										   $fan->{$m}{$n}=-$x;
										   $fan->{$n}{$m}=-$x;
										} 
									}     
								}	
							}else{
	
							}
						 }
					 }		 
				 }
	   }		 
	}     
    return $int,$fan;
}	
#------------------------------------------- 
sub round_to_10{    																 
	my $e = shift;
	#print "$e ";
	   if ( $e> 0 ){           		
			if ( ($e - (10 * int($e/10))) >5 ){
				  $e= 10+10 * int($e/10); 	
			}else{
				  $e= 10 * int($e/10); 	
			}
	   }else{
			if ( ( (10 * int($e/10)) - $e ) >5 ){
				  $e= -10+10 * int($e/10); 	
			}else{
				  $e= 10 * int($e/10); 	
			}           
	   }
	#print "round to $e\n";	   
 	return $e;          
}           
#------------------------------------------- 
sub round_to_20{    																 
	my $e = shift;
	#print "$e ";
	   if ( $e> 0 ){           		
			if ( ($e - (20 * int($e/20))) > 10 ){
				  $e= 20+20 * int($e/20); 	
			}else{
				  $e= 20 * int($e/20); 	
			}
	   }else{
			if ( ( (20 * int($e/20)) - $e ) >10 ){
				  $e= -20+20 * int($e/20); 	
			}else{
				  $e= 20 * int($e/20); 	
			}           
	   }
	#print "round to $e\n";	   
 	return $e;          
}           
#------------------------------------------- 