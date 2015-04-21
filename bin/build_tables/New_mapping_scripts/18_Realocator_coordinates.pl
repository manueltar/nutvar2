##	Script to calculate the offsets and regions no displayed from the BLAST alignment.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;


my $input=$ARGV[0];
my $output1=$ARGV[1];
my $output2=$ARGV[2];

my $AC="NaN";
my $ENSG="NaN";
my $ENST="NaN";
my %inner_alignment_hash=();
my %hash1=();

# Here we open the file with the overall information of the blast alignment in the last line of each register

my $time='['. timestamp(). ']'."\n";
	print "Tiempo de incio de carga_hash1:$time\n";
if(open(INPUT,$input) && open(OUTPUT,'>'.$output1))
{
	# Input= aligned_features_coordinates.txt
	
			#~ >sp|A0FGR8|ENSG00000117868_ENST00000251527
		#~ Missense:60     R__G
		#~ Missense:62     R__A
		#~ Missense:64     K__R
		#~ Missense:65     T__H
		#~ Missense:66     A__C
		#~ Missense:67     R__G
		#~ Missense:68     G__A
		#~ Missense:69     L__M
		#~ Missense:70     R__S
		#~ Missense:72     H__A
		#~ Missense:74     Q__G
		#~ Missense:75     R__E
		#~ Missense:77     A__P
		#~ Missense:78     G__E
		#~ Missense:81     L__A
		#~ Missense:82     S__G
		#~ Missense:83     R__G
		#~ Missense:84     P__A
		#~ Missense:86     S__G
		#~ Missense:87     A__R
		#~ Missense:88     R__A
		#~ GAP_ENST:91     S
		#~ GAP_ENST:92     P
		#~ GAP_ENST:93     P
		#~ GAP_ENST:94     R
		#~ GAP_ENST:95     P
		#~ GAP_ENST:96     G
		#~ GAP_ENST:97     G
		#~ GAP_ENST:98     P
		#~ **Length:       query:863(893;39|893)   subject:863(921;59|921)**

	
	
	while(my $line=<INPUT>)
	{
		chomp $line;
		
		# In the last line we unload the hashes and reinitialize hashes and variables
		
		if($line =~ /^Length:\tquery:([^\(]+)\(([^;]+);([^\|]+)\|([^\)]+)\)\tsubject:([^\(]+)\(([^;]+);([^\|]+)\|([^\)]+)\)/)
		{
			#~ print ">$AC\t$ENSG\t$ENST\n";
			
			# We parse the fields we need; from the query (ENSEMBL polipetide): the aligned residues; the total residues; the begin and the end of the alignment
			
			my $query_aligned=$1;
			my $query_total=$2;
			my $query_begin=$3;
			my $query_end=$4;
			
			# We parse the fields we need; from the subject (UniProt displayed isoform): the aligned residues; the total residues; the begin and the end of the alignment
			
			my $subject_aligned=$5;
			my $subject_total=$6;
			my $subject_begin=$7;
			my $subject_end=$8;
			
			if($ENST eq 'ENST00000525985')
			{
				print "FIN:$query_aligned\($query_total\;$query_begin\|$query_end\)\t$subject_aligned\($subject_total\;$subject_begin\|$subject_end\)\n";
			}
			
			my $FLAG1=0;
			my $FLAG2=0;
			my $FLAG3=0;
			
			# If the aligned start in residue 1 in both ->condition 1
			
			if($query_begin == $subject_begin && $query_begin ==1)
			{
				$FLAG1=1;
			}
			
			# If the aligned end in the same residue and comprises the whole of both proteins -> condition 2
			
			if($query_end == $subject_end && $query_end == $query_total && $subject_end == $subject_total)
			{
				$FLAG2=1;
			}
			
			# If the amount of aligned sequence is equal between query and subject (remember gaps) and is equal or bigger (gaps) than the
			# total length of the query -> condition 3
			
			if($query_aligned == $subject_aligned && $query_aligned >= $query_total)
			{
				$FLAG3=1;
			}
			
			if($ENST eq 'ENST00000525985')
			{
				print "FLAGS:$FLAG1\t$FLAG2\t$FLAG3\n";
			}
				# This is just to check the missense changes between UniProt and ENSEMBL, I am not using the information.
			
				foreach my $missense_POS(sort {$a<=>$b}keys %{$inner_alignment_hash{'missense'}})
				{
				foreach my $missense_change(sort keys %{$inner_alignment_hash{'missense'}{$missense_POS}})
				{
					#~ print "MISSENSE:$missense_POS\t$missense_change\t";
				}	
				}
				#~ print "\n";
				
				# We unload the GAPS in ENSEMBL. This internal gaps represent regions present in UniProt and absent in ENSEMBL
				# They give rise to offsets that need to be deducted for the protein coordinates of features lying c-terminal of them. 
				# Besides all the features lying totally in these regions won't be transferred to the ENSEMBL coordinates and those laying partially
				# have to be adapted so their boundaries lie within regions present in the ENSEMBL polipetide.
				# Here we collapse this gaps in intervals of protein coordinates
				
				my @GAP_ENST_POS=sort {$a<=>$b}keys %{$inner_alignment_hash{'GAP_ENST'}};
				my $amount_GAP_ENST=scalar(@GAP_ENST_POS);
				my @positions=();
				my $counter=0;
				my $max=scalar(@GAP_ENST_POS);
				
				while(scalar(@GAP_ENST_POS) !=0)
				{
					my $contender=shift(@GAP_ENST_POS);
					push(@positions,$contender);
					$counter++;
					#~ print  "GAP_ENST\tCounter:$counter\tEl array es:@positions\n";
					
					# For the first position we open the first interval
					
					if($counter == 1)
					{
						print OUTPUT  ">$AC\t$ENSG\t$ENST\tGAP_ENST\t\[$positions[$counter-1]\-";
					}
					
					# If the gap extends for a consecutive position do nothing
					
					elsif($counter != 1 && $counter != $max)
					{
						if($positions[$counter-1] == 1 + $positions[$counter-2])
						{
							# Do nothing
						}
					# If the following postion is no consecutive close the previous interval and open the new one
					
						elsif($positions[$counter-1] != 1 + $positions[$counter-2])
						{
							print OUTPUT  "$positions[$counter-2]\]\n";
							print OUTPUT  ">$AC\t$ENSG\t$ENST\tGAP_ENST\t\[$positions[$counter-1]\-";
						}
					}
					
					# If we have reached the last position of all the gap array
					
					elsif($counter == $max)
					{
						# If it is consecutive to the previous one, close the array
						
						if($positions[$counter-1] == 1 + $positions[$counter-2])
						{
							print OUTPUT  "$positions[$counter-1]\]\n";
						}
						
						# If it is not consecutive to the previous one then close the previous array and create a las one with the last position
						# as beginning and end of the interval
						
						elsif($positions[$counter-1] != 1 + $positions[$counter-2])
						{
							print OUTPUT  "$positions[$counter-2]\]\n";
							print OUTPUT  ">$AC\t$ENSG\t$ENST\tGAP_ENST\t\[$positions[$counter-1]\-$positions[$counter-1]\]\n";
						}
					}
				}
				
				# This strategy is new; if there was only one position for ENSEMBL gaps then the position array has only
				# one element and therefore we can close the initial interval using the same element.
				
				if(scalar(@positions) == 1)
				{
					print OUTPUT  "$positions[$counter-1]\]\n";
				}
				
				my @GAP_UNIPROT_POS_tmp=sort {$a<=>$b}keys %{$inner_alignment_hash{'GAP_UNIPROT'}};
				my @amount_GAP_UNIPROT=();
				
				# The gaps in UniProt are regions that are absent in the UniProt polypetide. Therefore from the extension of the gap has
				# to be added to UniProt features lying c-ter of the initial position of the GAP.
				
				foreach my $GAP_UNIPROT_POS_tok(@GAP_UNIPROT_POS_tmp)
				{
					my @tmp_residues=sort keys %{$inner_alignment_hash{'GAP_UNIPROT'}{$GAP_UNIPROT_POS_tok}};
					my $amount_GAP_UNIPROT=scalar(@tmp_residues);
					push(@amount_GAP_UNIPROT,$amount_GAP_UNIPROT);
					print OUTPUT  ">$AC\t$ENSG\t$ENST\tGAP_UNIPROT\t\[$GAP_UNIPROT_POS_tok\]\{$amount_GAP_UNIPROT\}\n";
				}
				
				# Here we measure the total expansion of the gap
				
				my $total_amount_GAP_UNIPROT = 0;
				for ( @amount_GAP_UNIPROT )
				{
					$total_amount_GAP_UNIPROT += $_;
				}
				
				#~ print OUTPUT "FLAGS:\t$FLAG1\t$FLAG2\t$FLAG3\n";
			
				# If conditions 1,2 and 3 apply then the offsets are 0
				
				if($FLAG1 ==1 && $FLAG2==1 && $FLAG3==1)
				{
					print OUTPUT "ENDLINE:$AC\t$ENSG\t$ENST\tUPSTREAM=0"."__"."DOWNSTREAM=0\n";
				}
				
				# If the origin in query and subject is one but the end is different from the total length of both
				
				elsif($FLAG1 ==1 && $FLAG2 !=1)
				{
					
					if($ENST eq 'ENST00000525985')
					{
						print "**********************************$query_end\t$amount_GAP_ENST\t$query_total\t$subject_end\t$total_amount_GAP_UNIPROT\t$subject_total\n";
					}
					
					# If the UniProt GAPS can account for the difference in alignment length for subject and query the offset is 0 (with a UniProt GAP)
					
					if($query_aligned + $amount_GAP_ENST >= $query_total && $subject_aligned + $total_amount_GAP_UNIPROT >= $subject_total)
					{
						print OUTPUT "ENDLINE:$AC\t$ENSG\t$ENST\tUPSTREAM=0"."__"."DOWNSTREAM=0\n";
						if($ENST eq 'ENST00000525985')
						{
							print "ENDLINE:$AC\t$ENSG\t$ENST\tUPSTREAM=0"."__"."DOWNSTREAM=0\n";
						}
					}
					
					# If not
					
					else
					{
						# Then the limit for features to be displayed is the length of the query. Beyond that the positions in UniProt are absent
						# in ENSEMBL
						
						print OUTPUT "ENDLINE:$AC\t$ENSG\t$ENST\tUPSTREAM=0"."__"."DOWNSTREAM=LIMIT:$query_total\n";
						if($ENST eq 'ENST00000525985')
						{
							print "ENDLINE:$AC\t$ENSG\t$ENST\tUPSTREAM=0"."__"."DOWNSTREAM=LIMIT:$query_total\n";
						}
					}
				}
				
				# If alignment does not begin in residue 1 of query and subject
				
				elsif($FLAG1 !=1)
				{
					#~ **Length:       query:863(893;39|893)   subject:863(921;59|921)**
					# We calculate the N-ter offset of the alignment, it can be positive or negative
					
					my $difference_upstream=$query_begin-$subject_begin;
					
					# If the initial offset plus the posible UniProt GAPs plus the end of the query equals the total subject
			
					if($query_end - $difference_upstream + $amount_GAP_ENST == $subject_end)
					{
						print OUTPUT "ENDLINE:$AC\t$ENSG\t$ENST\tUPSTREAM=LIMIT:$difference_upstream"."__"."DOWNSTREAM=0\n";
					}
					else
					{
						# If there  is no BLAST alignment
						
						if($query_aligned == 0)
						{
							print OUTPUT "ENDLINE:$AC\t$ENSG\t$ENST\tUPSTREAM=XX"."__"."DOWNSTREAM=XX\n";
						}
						
						# if there is alignment but C-ter and N-ter of the subjects are bigger (without counting on possible UniProt intra sequence gaps)
						
						else
						{
							print OUTPUT "ENDLINE:$AC\t$ENSG\t$ENST\tUPSTREAM=LIMIT:$difference_upstream"."__"."DOWNSTREAM=LIMIT:$query_total\n";
						}
					}
				}
			 $AC="NaN";
			 $ENSG="NaN";
			 $ENST="NaN";
			 %inner_alignment_hash=();
		}
		
		# We obtain the identifiers AC number, ENSG and ENST for the query
		
		elsif ($line=~/^>/)
		{
			$line=~s/>//ig;
			#~ print OUTPUT "INICIO:$line\n";
			my @tmp=split(/\|/,$line);
			$AC=$tmp[1];
			my $ENSG_ENST=$tmp[2];
			my @ENSG_ENST_tmp=split("_",$ENSG_ENST);
			$ENSG=$ENSG_ENST_tmp[0];
			$ENST=$ENSG_ENST_tmp[1];
		}
		
		# We load the missense hash 
		
		elsif($line !~ /^>/  || $line !~ /^Length:/)
		{
			if($line =~ /^Missense:([^\t]+)\t(.+)/)
			{
				my $POS=$1;
				my $change=$2;
				$inner_alignment_hash{'missense'}{$POS}{$change}=1;
				if($ENST eq 'ENST00000525985')
				{
					print "'missense'\t$POS\t$change\n";
				}
			}
			
			# We load the hash of ENSEMBL gaps
			
			elsif($line =~ /^GAP_ENST:([^\t]+)\t(.+)/)
			{
				my $POS=$1;
				my $change=$2;
				$inner_alignment_hash{'GAP_ENST'}{$POS}{$change}=1;
				if($ENST eq 'ENST00000525985')
				{
					print "'GAP_ENST'\t$POS\t$change\n";
				}
			}
			
			# We load the hash of UniProt gaps
			
			elsif($line =~ /^GAP_UNIPROT:([^\t]+)\t(.+)/)
			{
				my $POS=$1;
				my $change=$2;
				$inner_alignment_hash{'GAP_UNIPROT'}{$POS}{$change}=1;
				if($ENST eq 'ENST00000525985')
				{
					print "'GAP_UNIPROT'\t$POS\t$change\n";
				}
			}
		}
		
	}	
}else{print "Unable to open $input\n";}

# Here we open the file we have just created with condensed intervals for ENSEMBL gaps and Position an extension for UniProt Gaps
# and the outer offfsets due to N and C ter longer in UniProt calculated.

my %offset_hash=();

if(open(OUTPUT, $output1) && open(OUTPUT2, '>'. $output2))
{
	#~ Input=aligned_features_coordinates_midC.txt
	#~ 
	#~ >Q9Y679 ENSG00000115307 ENST00000377526 GAP_ENST        [23-35]
	#~ >Q9Y679 ENSG00000115307 ENST00000377526 GAP_ENST        [53-59]
	#~ >Q9Y679 ENSG00000115307 ENST00000377526 GAP_ENST        [73-79]
	#~ >Q9Y679 ENSG00000115307 ENST00000377526 GAP_ENST        [95-109]
	#~ >Q9Y679 ENSG00000115307 ENST00000377526 GAP_ENST        [145-153]
	#~ >Q9Y679 ENSG00000115307 ENST00000377526 GAP_ENST        [159-177]
	#~ >Q9Y679 ENSG00000115307 ENST00000377526 GAP_UNIPROT     [65]{4}
	#~ ENDLINE:Q9Y679  ENSG00000115307 ENST00000377526 UPSTREAM=0__DOWNSTREAM=0


	#~ A0A183  ENSG00000235942 ENST00000431011 UPSTREAM=0      DOWNSTREAM=0
		#~ A0AUZ9  ENSG00000144445 ENST00000281772 UPSTREAM=0      DOWNSTREAM=0

	while(my $line=<OUTPUT>)
	{
		chomp $line;
		
		# If we are in the last line we unload the hashes and re-initialize the variables
		
		if($line =~ /^ENDLINE:([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\_\_]+)\_\_(.+)/)
		{
			my $AC=$1;
			my $ENSG=$2;
			my $ENST=$3;
			my $UPSTREAM=$4;
			my $DOWNSTREAM=$5;
			
			# Here we undo "no display regions" as begin__end intervals and we add them to the category No Display
			
			print OUTPUT2"$AC\t$ENSG\t$ENST\t$UPSTREAM\t$DOWNSTREAM\t";
			my @no_display=();
			foreach my $begin_tok(sort{$a<=>$b} keys %{$offset_hash{'NO_DISPLAY'}})
			{
			foreach my $end_tok(sort keys %{$offset_hash{'NO_DISPLAY'}{$begin_tok}})
			{
				my $string=join("__",$begin_tok,$end_tok);
				push(@no_display,$string);
			}	
			}
			my $no_display=join(";",@no_display);
			
			if (scalar(@no_display) >0)
			{
				print OUTPUT2"NO_DISPLAY=$no_display\t";
			}
			else
			{
				print OUTPUT2"NO_DISPLAY=0\t";
			}
			
			# For the same GAPs in ENSEMBL polipetide we the partial offsets (negative or ositive) we have to apply for features lying C-ter
			# to the position of the ENSEMBL (negative) or UniProt (positive) GAP.
			
			my @offset=();
			foreach my $POS_tok(sort{$a<=>$b} keys %{$offset_hash{'OFFSET'}})
			{
			foreach my $expand_tok(sort keys %{$offset_hash{'OFFSET'}{$POS_tok}})
			{
				my $string=join("__",$POS_tok,$expand_tok);
				push(@offset,$string);
			}	
			}
			my $offset=join(";",@offset);
			
			if (scalar(@offset) >0)
			{
				print OUTPUT2"OFFSET=$offset\t";
			}
			else
			{
				print OUTPUT2"OFFSET=0\t";
			}
			
			print OUTPUT2"\n";
			%offset_hash=();
		}
		
		# We load the hash for ENSEMBL GAPS
		
		elsif($line =~ /^>([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			my $interval=$5;
			if($interval =~ /\[(.+)\-(.+)\]/)
			{
				# We calculate begin/end intervals of no display
				
				my $begin=$1;
				my $end=$2;
				$offset_hash{'NO_DISPLAY'}{$begin}{$end}=1;
				
				# We calculate the position and extension of the GAP we have to deduct from UniProt features lying c-ter of the GAP.
				# It is a deduction hence the - sign in the distance formula.
				
				my $POS_offset=$end+1;
				my $distance=-1*($end-$begin+1);
				$offset_hash{'OFFSET'}{$POS_offset}{$distance}=1;
			}
			
			# Here we calculate the offset we have to add for UniProt GAPS; ENSEMBL regions that are absent in UniProt.
			# It is a positive amount.
			
			elsif($interval =~ /\[(.+)\]\{(.+)\}/)
			{
				my $begin=$1;
				my $expand=$2;
				my $POS_offset=$begin;
				$offset_hash{'OFFSET'}{$POS_offset}{$expand}=1;
			}
			
		}
	}	
}

$time='['. timestamp(). ']'."\n";
print "Fin del script:$time\n";

sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
