##	Script to apply to trasnfer domain and site protein coordinates to ENSEMBL protein coordinates.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;


my $input=$ARGV[0];
my $input2=$ARGV[1];
my $input3=$ARGV[2];
my $input4=$ARGV[3];
my $input5=$ARGV[4];
my $output1=$ARGV[5];
my $output2=$ARGV[6];

my %hash1=();
my %hash_sites=();

# Here we parse the field containing the coordinates of each functional site referred to its UniProt identifier

my $time='['. timestamp(). ']'."\n";
print "Tiempo de incio de carga_hash1:$time\n";
if(open(INPUT,$input))
{
	# Input=  site_coordinates.txt
	
	#~ P31946  MOD_RES 1       1       N-acetylmethionine; in 14-3-3 protein
	#~ P31946  MOD_RES 106     106     Nitrated tyrosine (By similarity).
	#~ P31946  MOD_RES 117     117     N6-acetyllysine.
	#~ P31946  MOD_RES 186     186     Phosphoserine (By similarity).
	#~ P31946  MOD_RES 2       2       N-acetylthreonine; in 14-3-3 protein
	#~ P31946  MOD_RES 60      60      Phosphoserine (By similarity).
	#~ P31946  MOD_RES 70      70      N6-acetyllysine.
	#~ P31946  MOD_RES 84      84      Nitrated tyrosine (By similarity).
	#~ P62258  MOD_RES 1       1       N-acetylmethionine.
	#~ P62258  MOD_RES 118     118     N6-acetyllysine.
	#~ P62258  MOD_RES 123     123     N6-acetyllysine.
	#~ P62258  MOD_RES 210     210     Phosphoserine.
	#~ P62258  MOD_RES 50      50      N6-acetyllysine.
	#~ P62258  MOD_RES 69      69      N6-acetyllysine.
	
	while(my $line=<INPUT>)
	{
		chomp $line;
		#~ print "$line\n";
		if($line =~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t.+/)
		{
			#~ print "Hello_world:$line\n";
			my $AC=$1;
			my $SITE=$2;
			my $begin=$3;
			my $end=$4;
			
			$hash_sites{$AC}{'SITE'}{$SITE}{$begin}{$end}=1;
			#~ print "$AC\t'SITE'\t$SITE\t$begin\t$end\n";
		}
		
	}	
}else{print "Unable to open $input\n";}

$time='['. timestamp(). ']'."\n";
print "Processing_hash1:$time\n";

# Here we number the sites as there no overlaping problems; they all come from UniProt and show no overlap with sites of the same type

foreach my $AC_tok(sort keys %hash_sites)
{
	foreach my $SITE_tok(sort keys %{$hash_sites{$AC_tok}{'SITE'}})
	{
		my @tmp_begin=sort{$a<=>$b}keys%{$hash_sites{$AC_tok}{'SITE'}{$SITE_tok}};
		for(my $i=0;$i<scalar(@tmp_begin);$i++)
		{
			my $index=$i+1;
			my $string=join("**",$SITE_tok,$index);
			foreach my $end_tok(sort keys %{$hash_sites{$AC_tok}{'SITE'}{$SITE_tok}{$tmp_begin[$i]}})
			{
				$hash1{$AC_tok}{'SITE'}{$string}{$tmp_begin[$i]}{$end_tok}=1;
				#~ print "$AC_tok\t'SITES'\t$string\t$tmp_begin[$i]\t$end_tok\n";
			}
			
		}
	}
}

# Here we open the file with the domain coordinates we have processed before

%hash_sites=();

$time='['. timestamp(). ']'."\n";
print "Tiempo de incio de carga_hash2:$time\n";

if(open(INPUT2,$input2))
{
	# Input=  domain_coordinates.txt
	
	#~ A0A5B9  IPR003597__1    16      105
	#~ A0A5B9  IPR007110__1    8       117
	#~ A0A5B9  IPR013783__1    7       39
	#~ A0A5B9  IPR013783__2    73      128
	#~ A0AUZ9  IPR026180__1    1       987
	#~ A0AUZ9  IPR029332__1    795     915
	#~ A0AV02  IPR004841__1    44      411
	#~ A0AV96  IPR000504__1    71      145
	#~ A0AV96  IPR000504__2    151     233
	
	while(my $line=<INPUT2>)
	{
		chomp $line;
		my @tmp=split(/\t/,$line);
		my $AC=$tmp[0];
		my $IPR=$tmp[1];
		$IPR=~s/\_\_/\*\*/g;
		my $begin=$tmp[2];
		my $end=$tmp[3];
		$hash1{$AC}{'IPR'}{$IPR}{$begin}{$end}=1;
		#~ print "$AC\t'IPR'\t$IPR\t$begin\t$end\n";
	}	
}else{print "Unable to open $input2\n";}

$time='['. timestamp(). ']'."\n";
print "Tiempo de incio de carga_hash3:$time\n";

# Here we open the file to obtain the ENSG-SYMBOL-strand equivalence

my %hash2=();

if(open(INPUT3,$input3))
{
	# Input=  gtf_output_ENSG.txt
	
	#~ ENSG00000186092 OR4F5   +       69091   70008
	#~ ENSG00000237683 AL627309.1      -       134901  139379
	#~ ENSG00000235249 OR4F29  +       367640  368634
	#~ ENSG00000185097 OR4F16  -       621059  622053
	
	while(my $line=<INPUT3>)
	{
		chomp $line;
		#~ print "$line\n";
		my @tmp=split(/\t/,$line);
		#~ print "El array es:@tmp\n";
		my $ENSG=$tmp[0];
		#~ print "************$ENSG\n";
		my $SYMBOL=$tmp[1];
		
		$hash2{$ENSG}{$SYMBOL}=1;
	}	
}else{print "Unable to open $input3\n";}


# Here we open the file for the >18.100 instances where ENSEMBL matches and contains the UniProt displayed isoform


$time='['. timestamp(). ']'."\n";
print "Tiempo de incio de carga_hash4:$time\n";

if(open(INPUT4,$input4) && open(OUTPUT,'>'.$output1))
{
	# Input=  ENSEMBL_includes_UNIPROT.txt
	
	#~ A0A183  ENSG00000235942 ENST00000431011 UPSTREAM=0      DOWNSTREAM=0
	#~ A0AUZ9  ENSG00000144445 ENST00000281772 UPSTREAM=0      DOWNSTREAM=0
	#~ A0AV02  ENSG00000221955 ENST00000393469 UPSTREAM=0      DOWNSTREAM=0
	#~ A0AV96  ENSG00000163694 ENST00000295971 UPSTREAM=0      DOWNSTREAM=0
	#~ A0AVF1  ENSG00000105948 ENST00000464848 UPSTREAM=0      DOWNSTREAM=0
	#~ A0AVI4  ENSG00000168936 ENST00000382936 UPSTREAM=0      DOWNSTREAM=0
	#~ A0AVK6  ENSG00000129173 ENST00000250024 UPSTREAM=0      DOWNSTREAM=0
	#~ A0AVT1  ENSG00000033178 ENST00000322244 UPSTREAM=0      DOWNSTREAM=0
	#~ A0FGR9  ENSG00000158220 ENST00000389567 UPSTREAM=0      DOWNSTREAM=0

	
	while(my $line=<INPUT4>)
	{
		chomp $line;
		#~ print "$line\n";
		if($line =~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\tUPSTREAM=(.+)\tDOWNSTREAM=(.+)/)
		{
			#~ print "Hello_world:$line\n";
			my $AC=$1;
			my $ENSG=$2;
			my $ENST=$3;
			my $upstream=$4;
			my $downstream=$5;
			if(exists($hash2{$ENSG}))
			{
			
				print OUTPUT "$AC\t$ENSG\t$ENST\t";
				
				foreach my $SITE_tok(sort keys %{$hash1{$AC}{'SITE'}})
				{
				foreach my $begin_SITE_tok(sort {$a<=>$b}keys %{$hash1{$AC}{'SITE'}{$SITE_tok}})
				{
				foreach my $end_SITE_tok(sort {$a<=>$b}keys %{$hash1{$AC}{'SITE'}{$SITE_tok}{$begin_SITE_tok}})
				{
					# In the cases where the  N-ter of ENSEMBL is bigger we have to add the offset to every site
					
					my $realocated_begin=$begin_SITE_tok+$upstream;
					my $realocated_end=$end_SITE_tok+$upstream;
					
					my $distance_begin=($realocated_begin-1)*3;
					my $distance_feature=($realocated_end-$realocated_begin +1)*3;
					
					print OUTPUT "FEATURE:$SITE_tok"."__"."$distance_begin"."__"."$distance_feature\t";
				}	
				}	
				}
				foreach my $IPR_tok(sort keys %{$hash1{$AC}{'IPR'}})
				{
				foreach my $begin_IPR_tok(sort {$a<=>$b}keys %{$hash1{$AC}{'IPR'}{$IPR_tok}})
				{
				foreach my $end_IPR_tok(sort {$a<=>$b}keys %{$hash1{$AC}{'IPR'}{$IPR_tok}{$begin_IPR_tok}})
				{
					# In the cases where the  N-ter of ENSEMBL is bigger we have to add the offset to every protein domain
					
					my $realocated_begin=$begin_IPR_tok+$upstream;
					my $realocated_end=$end_IPR_tok+$upstream;
					
					my $distance_begin=($realocated_begin-1)*3;
					my $distance_feature=($realocated_end-$realocated_begin +1)*3;
					
					print OUTPUT "FEATURE:$IPR_tok"."__"."$distance_begin"."__"."$distance_feature\t";
				}	
				}	
				}
				
				print OUTPUT "\n";
			}
		}
		
	}	
}else{print "Unable to open $input4\n";}
$time='['. timestamp(). ']'."\n";
print "Tiempo de incio de carga_hash5:$time\n";

# Here we open the file for the 712 instances where ENSEMBL doesn't match the UniProt displayed isoform for which we have calculated
# the no-display intervals (intervals present in UniProt but absent in ENSEMBL), the partial offsets of intrasequence loops and the differences in
# N-ter and C-ter length condensed in the upstream and downstream offsets.

if(open(INPUT5,$input5) && open(OUTPUT,'>>'.$output1))
{
	# Input=  aligned_no_display_offset.txt
	
	#~ O60763  ENSG00000138768 UPSTREAM=LIMIT:-72      DOWNSTREAM=0    NO_DISPLAY=0    OFFSET=0        
	#~ O60774  ENSG00000117507 UPSTREAM=0      DOWNSTREAM=LIMIT:418    NO_DISPLAY=0    OFFSET=0        
	#~ O75140  ENSG00000100150 UPSTREAM=0      DOWNSTREAM=0    NO_DISPLAY=725__733     OFFSET=734__-9  
	#~ O75949  ENSG00000130054 UPSTREAM=0      DOWNSTREAM=0    NO_DISPLAY=369__369     OFFSET=370__-1  
	#~ O95072  ENSG00000100918 UPSTREAM=0      DOWNSTREAM=0    NO_DISPLAY=236__236     OFFSET=237__-1  

	
	while(my $line=<INPUT5>)
	{
		chomp $line;
		#~ print "$line\n";
		if($line =~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\tUPSTREAM=(.+)\tDOWNSTREAM=(.+)\tNO_DISPLAY=(.+)\tOFFSET=(.+)/)
		{
			#~ print "Hello_world:$line\n";
			my $AC=$1;
			my $ENSG=$2;
			my $ENST=$3;
			my $upstream=$4;
			my $downstream=$5;
			my $no_display=$6;
			my $offset=$7;
			
			# Here we discard UniProt dispalyed isoforms for which there was no alignment (upstream eq XX)
			
			if(exists($hash2{$ENSG}) && $upstream ne 'XX')
			{
				# LIMITS regions in N-ter or C-ter absent in ENSEMBL but present in UniProt
				
				my $upstream_limit=0;
				my $upstream_def=0;
				
				# If there is a N-ter limit the substitute the original value set to 0
				
				if($upstream =~/LIMIT:(.+)/)
				{
					$upstream_limit=$1;
					$upstream_def=$1;
					$upstream_limit =~ s/\-//g;
				}
				else
				{
					$upstream_def=$upstream;
					$upstream_def =~ s/LIMIT\://g;
				}
				
				# If there is a N-ter limit the substitute the original value set to 10.000 (more residues than any known protein)
				
				my $downstream_limit=100000;
				my $downstream_def=0;
				
				if($downstream =~/LIMIT:(.+)/)
				{
					$downstream_limit=$1;
					$downstream_def=$1;
					$downstream_limit =~ s/\-//g;
				}
				else
				{
					$downstream_def=$downstream;
					$downstream_def =~ s/LIMIT\://g;
				}
				
				#~ print "Hello_worldII:$AC\t$ENSG\t$ENST\t$upstream_def\tLIMIT:\t$upstream_limit\t$downstream_def\tLIMIT:\t$downstream_limit\t$no_display\t$offset\n";
				print OUTPUT "$AC\t$ENSG\t$ENST\t";
				
				my %forbidden_positions=();
				
				# Here we create the hash of forbidden positions; for every transcript all the features lying completely in these postions will not be transferred
				# and those lying partially will be adjusted
				
				if($no_display =~ /\_\_/)
				{
					my @no_display_tmp=split(";",$no_display);
					#~ print "El array es:@no_display_tmp\n";
					foreach my $no_display_tmp_tok(@no_display_tmp)
					{
						if($no_display_tmp_tok=~ /([^\_]+)\_\_(.+)/)
						{
							my $begin=$1;
							my $end=$2;
							for(my$i=$begin;$i<=$end;$i++)
							{
								$forbidden_positions{$i}=1;
							}
						}
					}
				}
				
				# Here we parse the begin and extension of the partial intrasequence offsets
				
				my %offset_positions=();
				
				if($offset=~ /\_\_/)
				{
					my @offset_tmp=split(";",$offset);
					#~ print "El array es:@offset_tmp\n";
					foreach my $offset_tmp_tok(@offset_tmp)
					{
						if($offset_tmp_tok=~ /([^\_]+)\_\_(.+)/)
						{
							my $begin=$1;
							my $expand=$2;
							$offset_positions{$begin}{$expand}=1;
						}
					}
				}
				
				
				foreach my $SITE_tok(sort keys %{$hash1{$AC}{'SITE'}})
				{
					my %allowed_positions=();
					my %allowed_offset_positions=();
					
				foreach my $begin_SITE_tok(sort {$a<=>$b}keys %{$hash1{$AC}{'SITE'}{$SITE_tok}})
				{
				foreach my $end_SITE_tok(sort {$a<=>$b}keys %{$hash1{$AC}{'SITE'}{$SITE_tok}{$begin_SITE_tok}})
				{
					if($end_SITE_tok > $upstream_limit && $begin_SITE_tok < $downstream_limit)
					{
						# First we exclude features that are outside limits of the maximal alignment. We redo the boundaries of features that lie in the middle of these limits.
						
						if($begin_SITE_tok > $upstream_limit && $end_SITE_tok < $downstream_limit)
						{
							# Do nothing
						}
						elsif($begin_SITE_tok <= $upstream_limit && $end_SITE_tok < $downstream_limit)
						{
							$begin_SITE_tok=$upstream_limit+1;
						}
						elsif($begin_SITE_tok > $upstream_limit && $end_SITE_tok > $downstream_limit)
						{
							$end_SITE_tok=$downstream_limit-1;
						}
						elsif($begin_SITE_tok <= $upstream_limit && $end_SITE_tok > $downstream_limit)
						{
							$begin_SITE_tok=$upstream_limit+1;
							$end_SITE_tok=$downstream_limit-1;
						}
						
						# Second we exclude feature positions that lie in loops; parts of the UNIPROT that are not present in ENSEMBL
						
						for(my $j=$begin_SITE_tok;$j<=$end_SITE_tok;$j++)
						{
							if(exists($forbidden_positions{$j}))
							{
								# Do nothing
							}
							else
							{
								$allowed_positions{$j}=1;
							}
						}
						
						# We skip all the sites which lie completely in loops
						
						my @threshold=sort{$a<=>$b}keys%allowed_positions;
						
						unless(scalar(@threshold) == 0)
						{
							# Third the offset created by these loops in downstream residues and the offset of reidues present in ENSEMBL but not in UNIPROT (GAP_UNIPROT) is applied to each position.
								
							foreach my $POS_allowed(@threshold)
							{
								#~ print "FEATURE:$SITE_tok\tPOS_allowed=$POS_allowed\n";
								my @offset_apply_tmp=();
								
							foreach my $POS_offset(sort{$a<=>$b}keys%offset_positions)
							{
							foreach my $expand_offset(sort{$a<=>$b}keys%{$offset_positions{$POS_offset}})
							{
								# Here we calculate all the offsets to be applied on a POS_allowed
								
								if($POS_allowed >= $POS_offset)
								{
									push(@offset_apply_tmp,$expand_offset)
								}
							}		
							}
							
							# Now we apply the total offset on POS_allowed
							
								my $total_amount_offset_POS = 0;
								for ( @offset_apply_tmp )
								{
									$total_amount_offset_POS += $_;
								}
							
									my $POS_def=$POS_allowed+$total_amount_offset_POS;
									$allowed_offset_positions{$POS_def}=1;
							}
		
							my @tmp_allowed=sort{$a<=>$b}keys%allowed_offset_positions;
							
							#~ print "FEATURE:$SITE_tok\tEl array ALLOWED_OFFSETED ES: es:@tmp_allowed\n";
							
							my $new_begin=$tmp_allowed[0];
							my $new_end=$tmp_allowed[scalar(@tmp_allowed)-1];
							
							# Fourth and finally we apply initial (upstream_def) offset to each position
								
							my $realocated_begin=$new_begin+$upstream_def;
							my $realocated_end=$new_end+$upstream_def;
							
							# Finally we calculate distance begin and ditance feature in 3-base cordinates (turning residues length to codon length)
							
							my $distance_begin=($realocated_begin-1)*3;
							my $distance_feature=($realocated_end-$realocated_begin +1)*3;
									
							#~ print  "FEATURE:$SITE_tok"."__"."$new_begin"."__"."$new_end\t";
							print  OUTPUT "FEATURE:$SITE_tok"."__"."$distance_begin"."__"."$distance_feature\t";
						} #unless
					}
				}	
				}	
				}
				foreach my $IPR_tok(sort keys %{$hash1{$AC}{'IPR'}})
				{
					my %allowed_positions=();
					my %allowed_offset_positions=();
					
				foreach my $begin_IPR_tok(sort {$a<=>$b}keys %{$hash1{$AC}{'IPR'}{$IPR_tok}})
				{
				foreach my $end_IPR_tok(sort {$a<=>$b}keys %{$hash1{$AC}{'IPR'}{$IPR_tok}{$begin_IPR_tok}})
				{
					if($end_IPR_tok > $upstream_limit && $begin_IPR_tok < $downstream_limit)
					{
						# First we exclude features that are outside limits of the maximal alignment. We redo the boundaries of features that lie in the middle of these limits.
						
						if($begin_IPR_tok > $upstream_limit && $end_IPR_tok < $downstream_limit)
						{
							# Do nothing
						}
						elsif($begin_IPR_tok <= $upstream_limit && $end_IPR_tok < $downstream_limit)
						{
							$begin_IPR_tok=$upstream_limit+1;
						}
						elsif($begin_IPR_tok > $upstream_limit && $end_IPR_tok > $downstream_limit)
						{
							$end_IPR_tok=$downstream_limit-1;
						}
						elsif($begin_IPR_tok <= $upstream_limit && $end_IPR_tok > $downstream_limit)
						{
							$begin_IPR_tok=$upstream_limit+1;
							$end_IPR_tok=$downstream_limit-1;
						}
						
						# Second we exclude feature positions that lie in loops; parts of the UNIPROT that are not present in ENSEMBL
						
						for(my $j=$begin_IPR_tok;$j<=$end_IPR_tok;$j++)
						{
							if(exists($forbidden_positions{$j}))
							{
								# Do nothing
							}
							else
							{
								$allowed_positions{$j}=1;
							}
						}
						
						# We skip all the sites which lie completely in loops
						
						my @threshold=sort{$a<=>$b}keys%allowed_positions;
						
						unless(scalar(@threshold) == 0)
						{
							# Third the offset created by these loops in downstream residues and the offset of reidues present in ENSEMBL but not in UNIPROT (GAP_UNIPROT) is applied to each position.
								
							foreach my $POS_allowed(@threshold)
							{
								#~ print "FEATURE:$IPR_tok\tPOS_allowed=$POS_allowed\n";
								my @offset_apply_tmp=();
								
							foreach my $POS_offset(sort{$a<=>$b}keys%offset_positions)
							{
							foreach my $expand_offset(sort{$a<=>$b}keys%{$offset_positions{$POS_offset}})
							{
								# Here we calculate all the offsets to be applied on a POS_allowed
								
								if($POS_allowed >= $POS_offset)
								{
									push(@offset_apply_tmp,$expand_offset)
								}
							}		
							}
							
							# Now we apply the total offset on POS_allowed
							
								my $total_amount_offset_POS = 0;
								for ( @offset_apply_tmp )
								{
									$total_amount_offset_POS += $_;
								}
							
									my $POS_def=$POS_allowed+$total_amount_offset_POS;
									$allowed_offset_positions{$POS_def}=1;
							}
		
							my @tmp_allowed=sort{$a<=>$b}keys%allowed_offset_positions;
							
							#~ print "FEATURE:$IPR_tok\tEl array ALLOWED_OFFSETED ES: es:@tmp_allowed\n";
							
							my $new_begin=$tmp_allowed[0];
							my $new_end=$tmp_allowed[scalar(@tmp_allowed)-1];
							
							# Fourth and finally we apply initial (upstream_def) offset to each position
								
							my $realocated_begin=$new_begin+$upstream_def;
							my $realocated_end=$new_end+$upstream_def;
							
							# Finally we calculate distance begin and ditance feature in 3-base cordinates (turning residues length to codon length)
							
							my $distance_begin=($realocated_begin-1)*3;
							my $distance_feature=($realocated_end-$realocated_begin +1)*3;
									
							#~ print  "FEATURE:$IPR_tok"."__"."$new_begin"."__"."$new_end\t";
							print  OUTPUT "FEATURE:$IPR_tok"."__"."$distance_begin"."__"."$distance_feature\t";
						} #unless
					}
				}	
				}	
				}			
				print  OUTPUT "\n";
			}
		}
	}	
}else{print "Unable to open $input5\n";}


sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
