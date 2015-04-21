##	Script to calculate the percentage of the sequence affected by the truncation. 2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne) 


use strict;
use warnings;
use Time::localtime;
use Memory::Usage;
    
my $mu = Memory::Usage->new();

    # Record amount of memory used by current process
$mu->record('starting work');

my %hash1=();
my %hash2=();
my %hash3=();
my %hash_splice=();
my %RESULT_hash=();

my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $output=$ARGV[2];

my $time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

if (open (INPUT1, $input1))
{
	#INPUT1=XXX_out_vep_parsed.vcf
	
	#	13      32906565        .       CA      C       33.73   .       frameshift_variant&feature_truncation;BRCA2;ENST00000380152;protein_coding;;0.5;;;;0.5;;;;NaN;NaN;NaN;NaN;PIRSF_domain:PIRSF002397;;
	#	13      32906565        .       CA      C       33.73   .       frameshift_variant&feature_truncation;BRCA2;ENST00000530893;protein_coding;;0.5;;;;0.5;;;;NaN;NaN;NaN;NaN;NaN;;
	#	13      32906565        .       CA      C       33.73   .       frameshift_variant&feature_truncation;BRCA2;ENST00000544455;protein_coding;;0.5;;;;0.5;;;;NaN;NaN;NaN;NaN;PIRSF_domain:PIRSF002397;;
	#	13      32906576        .       CA      C       36.73   .       frameshift_variant&feature_truncation;BRCA2;ENST00000380152;protein_co
	
	#INPUT1b=XXX-eff_parsed.vcf
	
	#13      32906565        .       CA      C       33.73   .       frameshift_variant;BRCA2;ENST00000380152;protein_coding;T317;10;;CODING;aca/;3418;HIGH;1;0.5;NaN;NaN    LOF_positive_snpEff;BRCA2;ENSG00000139618;6;0.50        
	#13      32906565        .       CA      C       33.73   .       frameshift_variant;BRCA2;ENST00000530893;protein_coding;T194;10;;CODING;aca/;480;HIGH;1;0.5;WARNING_TRANSCRIPT_INCOMPLETE;NaN   LOF_positive_snpEff;BRCA2;ENSG00000139618;6;0.50        
	#13      32906565        .       CA      C       33.73   .       frameshift_variant;BRCA2;ENST00000544455;protein_coding;T317;10;;CODING;aca/;3418;HIGH;1;0.5;NaN;NaN    LOF_positive_snpEff;BRCA2;ENSG00000139618;6;0.50        
	#13      32906576        .       CA      C       36.73   .       frameshift_variant;BRCA2;ENST00000380152;protein_coding;Q321;10;;CODING;caa/;3418;HIGH;1;0.5;NaN;NaN    LOF_positive_snpEff;BRCA2;ENSG00000139618;6;0.50        

while(my $line=<INPUT1>)
	{
		chomp $line;
		#print "$line\n";
		if($line=~/^([^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)\t[^\t]+\t[^\t]+\t(.+)/)
		{
			my $CHROM=$1;
			my $POS=$2;
			my $REF=$3;
			my $ALT=$4;
			my $fields=$5;
			#print "hello_world:$CHROM\t$POS\t$fields\t$REF\t$ALT\n";
			my @fields_tmp=split(";",$fields);
			my $Effect=$fields_tmp[0];
			my $SYMBOL=$fields_tmp[1];
			my $ENST=$fields_tmp[2];
			
			unless($Effect =~/splice_donor_variant/ || $Effect =~/splice_acceptor_variant/)
			{
				$hash1{$CHROM}{$SYMBOL}{$ENST}{$POS}{$REF}{$ALT}{$Effect}=1;
				#print "hello_world:$CHROM\t$POS\t$REF\t$ALT\t$Effect\t$SYMBOL\t$ENST\n";
			}
			elsif($Effect=~/splice_acceptor_variant/ || $Effect=~/splice_donor_variant/)
			{
				unless($Effect=~/nc_transcript_variant/)
				{
					#print"***********************$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t$Effect\n";
					$hash_splice{$CHROM}{$SYMBOL}{$ENST}{$POS}{$REF}{$ALT}{$Effect}=1;
				}
			}
		}
	}
}else{print "Unable to open INPUT1\n";}

    # Record amount in use afterwards
$mu->record();

    # Spit out a report
$mu->dump();


$time='['. timestamp(). ']'."\n";
print "Start charging hash2:$time\n";

if (open(INPUT2, $input2))
{    
	while(my $line=<INPUT2>)
	{
		chomp($line);
		# Input=ENST_table_full_condensed.txt
		
		#~ 1       A3GALT2 -       ENST00000330379 [33777791-33777822]     
		#~ 1       A3GALT2 -       ENST00000330379 ENST00000442999 [33772370-33773054]     [33777653-33777790]     
		#~ 1       A3GALT2 -       ENST00000442999 [33778102-33778191]     [33778408-33778491]     [33786677-33786699]     
		#~ 1       AADACL3 +       ENST00000332530 [12776344-12776347]     
		#~ 1       AADACL3 +       ENST00000332530 ENST00000359318 [12780885-12780948]     [12785189-12785960]     
		#~ 1       AADACL3 +       ENST00000359318 [12779480-12779693]
		
		#~ print "Hello worldI:$line**\n";
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)\t$/)
		{
			my %order_hash=();
			#~ print "Hello worldII:$line**\n";
			my $CHROM=$1;
			my $SYMBOL=$2;
			my $strand=$3;
			my $fields=$4;
			#~ print "$CHROM\t$strand\t$SYMBOL\t**$fields**\n";
			
			my @fields_tmp=split(/\t/,$fields);
			foreach my $fields_tmp_tok(@fields_tmp)
			{
				if($fields_tmp_tok=~/^ENST\d+/)
				{
					#~ print "AA:$CHROM\t$strand\t$SYMBOL\t$fields_tmp_tok:AA\n";
					$order_hash{'ENST'}{$fields_tmp_tok}=1;
				}
				else
				{
					#~ print "CC:$CHROM\t$strand\t$SYMBOL\t$fields_tmp_tok:CC\n";
					$order_hash{'INTERVAL'}{$fields_tmp_tok}=1;
				}
			}
			foreach my $ENST_tok(sort keys%{$order_hash{'ENST'}})
			{
			foreach my $interval_tok(sort keys%{$order_hash{'INTERVAL'}})
			{
				
				$hash2{$CHROM}{$strand}{$SYMBOL}{$ENST_tok}{$interval_tok}=1;
				#~ print "$CHROM\t$strand\t$SYMBOL\t$ENST_tok\t$interval_tok**\n";
			}	
			}
		}
	}
}else {print "impossible to open INPUT2\n";die;}




    # Record amount in use afterwards
$mu->record();

    # Spit out a report
$mu->dump();

$time='['. timestamp(). ']'."\n";
print "Start processing hash2:$time\n";

foreach my $CHROM_tok(sort keys%hash2)
{
foreach my $strand_tok(sort keys %{$hash2{$CHROM_tok}})
{
foreach my $SYMBOL_tok(sort keys %{$hash2{$CHROM_tok}{$strand_tok}})
{
foreach my $ENST_tok(sort keys %{$hash2{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}})
{
foreach my $interval_tok(sort keys %{$hash2{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENST_tok}})
{
	if($interval_tok=~/\[(.+)\-(.+)\]/)
	{
		my $begin=$1;
		my $end=$2;
		my $distance=$end-$begin+1;
		# As there might be repeated distance values total distances should be kept in an array rahter than in a hash
		push(@{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'DISTANCE'}},$distance);
		$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'BEGIN'}{$begin}{$end}=1;
		$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'END'}{$end}{$begin}=1;
	}
}
}	
}	
}
}

%hash2=();

$time='['. timestamp(). ']'."\n";
print "Start PRINTING:$time\n";

if(open(OUTPUT, '>'.$output))
{
	foreach my $CHROM_tok(sort keys%hash1)
	{
	foreach my $SYMBOL_tok(sort keys%{$hash1{$CHROM_tok}})
	{
	foreach my $ENST_tok(sort keys%{$hash1{$CHROM_tok}{$SYMBOL_tok}})
	{
	foreach my $POS_tok(sort keys%{$hash1{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}})
	{
	foreach my $REF_tok(sort keys%{$hash1{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}})
	{
	foreach my $ALT_tok(sort keys%{$hash1{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}{$REF_tok}})
	{
	foreach my $Effect_tok(sort keys%{$hash1{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}{$REF_tok}{$ALT_tok}})
	{
		my $string=join("\t",$CHROM_tok,$SYMBOL_tok,$POS_tok,$REF_tok,$ALT_tok,$Effect_tok,$ENST_tok);
		foreach my $strand_tok(sort keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}})
		{
			#~ print "INICIO:$string\t$strand_tok\n";
			my @distance_tmp=sort{$a<=>$b}@{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'DISTANCE'}};
			#~ print "El array es:$ENST_tok\t@distance_tmp\n";
			my $TOTAL_distance_ENST = 0;
			for ( @distance_tmp )
			{
				$TOTAL_distance_ENST += $_;
			}
			#~ print "La distancia total es:$ENST_tok\t$TOTAL_distance_ENST\n";
				if($strand_tok eq '+')
				{
					my @positions=();
					my @downstream_distances=();
					my @Begin_tmp=sort{$a<=>$b}keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'BEGIN'}};
					my @End_tmp=sort{$a<=>$b}keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'END'}};
					if($POS_tok >= $Begin_tmp[0] && $POS_tok <= $End_tmp[scalar(@End_tmp)-1])
					{
						my $contender=0;
						my $distance_affected_initial=0;
						my $counter=-1;
						if(scalar(@Begin_tmp) == 1)
						{
							#~ print "El array es:$ENST_tok\t@distance_tmp\n";
							#~ print "$string\tEl array begin es:@Begin_tmp\n";
							foreach my $end_tok(sort{$a<=>$b} keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'BEGIN'}{$Begin_tmp[0]}})
							{
								my $distance_affected=$end_tok-$POS_tok+1;
								#~ print "La distancia afectada es:$end_tok-$POS_tok+1:$distance_affected\n";
								my $Percentage_affected=100*($distance_affected/$TOTAL_distance_ENST);
								#~ $RESULT_hash{$string}{$Percentage_affected}=1;
								print OUTPUT "$string\t$Percentage_affected\n";
							}
						}
						elsif(scalar(@Begin_tmp) > 1)
						{
							my @Begin_shifted=@Begin_tmp;
							#~ print "Array begin inicial:@Begin_tmp\n";
							while(scalar(@Begin_shifted) !=0)
							{
								$contender=shift(@Begin_shifted);
								push(@positions,$contender);
								if($contender <= $POS_tok){$counter++;}
								#~ print "El array begin es:@Begin_shifted\t$counter\n";
								#~ print "El array contender es:@positions\t$counter\n";
							}
							foreach my $end_tok_affected(sort keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'BEGIN'}{$Begin_tmp[$counter]}})
							{
								#~ print "***$Begin_tmp[$counter]\t$POS_tok\t$end_tok_affected****\n";
								unless($POS_tok > $end_tok_affected)
								{
									$distance_affected_initial=$end_tok_affected-$POS_tok+1;
								}
								#~ print "distance_affected_initial:$distance_affected_initial\t$end_tok_affected\t$POS_tok\n";
							}
							if($counter < scalar(@Begin_tmp)-1)
							{
								for (my $i=$counter+1;$i < scalar(@Begin_tmp);$i++)
								{
									foreach my $end_downstream_tok(sort keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'BEGIN'}{$Begin_tmp[$i]}})
									{
										my $distance=1+$end_downstream_tok-$Begin_tmp[$i];
										#~ print "$distance\t$end_downstream_tok\t$Begin_tmp[$i]\n";
										push(@downstream_distances,$distance);
									}
								}
								my @sorted=sort@downstream_distances;
								#~ print "Las distancias downstream son:@sorted\n";
								my $Affected_downstream_distance = 0;
								for ( @downstream_distances )
								{
									$Affected_downstream_distance += $_;
								}
								#~ print "El total downstream afectado es:$Affected_downstream_distance\n";
								my $distance_affected=$distance_affected_initial+$Affected_downstream_distance;
								#~ print "La distancia afectada total es: $distance_affected_initial\t$Affected_downstream_distance\t$distance_affected\n";
								my $Percentage_affected=100*$distance_affected/$TOTAL_distance_ENST;
								#~ $RESULT_hash{$string}{$Percentage_affected}=1;
								print OUTPUT "$string\t$Percentage_affected\n";
							}
							elsif($counter == scalar(@Begin_tmp)-1)
							{
								my $distance_affected=$distance_affected_initial;
								my $Percentage_affected=100*($distance_affected_initial/$TOTAL_distance_ENST);
								#~ $RESULT_hash{$string}{$Percentage_affected}=1;
								print OUTPUT "$string\t$Percentage_affected\n";
							}
						}
					}
					elsif($POS_tok < $Begin_tmp[0])
					{
						my $Percentage_affected="NaN";
						print OUTPUT "$string\t$Percentage_affected\n";
						
					}
					elsif($POS_tok > $End_tmp[scalar(@End_tmp)-1])
					{
						my $Percentage_affected="NaN";
						print OUTPUT "$string\t$Percentage_affected\n";
					}
				}
				elsif($strand_tok eq '-')
				{
					my @positions=();
					my @upstream_distances=();
					my @Begin_tmp=sort{$a<=>$b}keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'BEGIN'}};
					my @End_tmp=reverse sort{$a<=>$b}keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'END'}};
					my $contender=0;
					my $distance_affected_initial=0;
					my $counter=-1;
					if($POS_tok >= $Begin_tmp[0] && $POS_tok <= $End_tmp[0])
					{
					
						if(scalar(@End_tmp) == 1)
						{
							#~ print "El array es:$ENST_tok\t@distance_tmp\n";
							#~ print "$string\tEl array end es:@End_tmp\n";
							foreach my $begin_tok(sort{$a<=>$b} keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'END'}{$End_tmp[0]}})
							{
								my $distance_affected=$POS_tok-$begin_tok+1;
								#~ print "La distancia afectada es:$begin_tok-$POS_tok+1:$distance_affected\n";
								my $Percentage_affected=100*($distance_affected/$TOTAL_distance_ENST);
								#~ $RESULT_hash{$string}{$Percentage_affected}=1;
								print OUTPUT "$string\t$Percentage_affected\n";
							}
						}
						elsif(scalar(@End_tmp) > 1)
						{
							my @End_shifted=@End_tmp;
							#~ print "Array end inicial:@End_tmp"."\t".scalar(@End_tmp)."\n";
							while(scalar(@End_shifted) !=0)
							{
								$contender=shift(@End_shifted);
								push(@positions,$contender);
								if($contender >= $POS_tok){$counter++;}
								#~ print "El array end es:@End_shifted\t$counter\n";
								#~ print "El array contender es:@positions\t$counter\n";
							}
							foreach my $begin_tok_affected(sort keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'END'}{$End_tmp[$counter]}})
							{
								#~ print "***$End_tmp[$counter]\t$POS_tok\t$begin_tok_affected****\n";
								unless($POS_tok < $begin_tok_affected)
								{
									$distance_affected_initial=$POS_tok-$begin_tok_affected+1;
								}
							}
							if($counter < scalar(@End_tmp)-1)
							{
								for (my $i=$counter+1;$i < scalar(@End_tmp);$i++)
								{
									foreach my $begin_upstream_tok(sort keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'END'}{$End_tmp[$i]}})
									{
										my $distance=1+$End_tmp[$i]-$begin_upstream_tok;
										#~ print "$distance\t$End_tmp[$i]\t$begin_upstream_tok\n";
										push(@upstream_distances,$distance);
									}
								}
								my @sorted=sort@upstream_distances;
								#~ print "Las distancias upstream son:@sorted\n";
								my $Affected_upstream_distance = 0;
								for ( @upstream_distances )
								{
									$Affected_upstream_distance += $_;
								}
								#~ print "El total upstream afectado es:$Affected_upstream_distance\n";
								my $distance_affected=$distance_affected_initial+$Affected_upstream_distance;
								#~ print "La distancia afectada total es: $distance_affected_initial\t$Affected_upstream_distance\t$distance_affected\n";
								my $Percentage_affected=100*$distance_affected/$TOTAL_distance_ENST;
								#~ $RESULT_hash{$string}{$Percentage_affected}=1;
								print OUTPUT "$string\t$Percentage_affected\n";
							}
							elsif($counter == scalar(@End_tmp)-1)
							{
								my $distance_affected=$distance_affected_initial;
								my $Percentage_affected=100*($distance_affected/$TOTAL_distance_ENST);
								#~ $RESULT_hash{$string}{$Percentage_affected}=1;
								print OUTPUT "$string\t$Percentage_affected\n";
							}
						}
					}
					elsif($POS_tok < $Begin_tmp[0])
					{
						my $Percentage_affected="NaN";
						print OUTPUT "$string\t$Percentage_affected\n";
						
					}
					elsif($POS_tok > $End_tmp[0])
					{
						my $Percentage_affected="NaN";
						print OUTPUT "$string\t$Percentage_affected\n";
					}
				}
		}
	}	
	}	
	}	
	}#	
	}	
	}
	}
	
	foreach my $CHROM_tok(sort keys%hash_splice)
	{
	foreach my $SYMBOL_tok(sort keys%{$hash_splice{$CHROM_tok}})
	{
	foreach my $ENST_tok(sort keys%{$hash_splice{$CHROM_tok}{$SYMBOL_tok}})
	{
	foreach my $POS_tok(sort keys%{$hash_splice{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}})
	{
	foreach my $REF_tok(sort keys%{$hash_splice{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}})
	{
	foreach my $ALT_tok(sort keys%{$hash_splice{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}{$REF_tok}})
	{
	foreach my $Effect_tok(sort keys%{$hash_splice{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}{$REF_tok}{$ALT_tok}})
	{
		my $string=join("\t",$CHROM_tok,$SYMBOL_tok,$POS_tok,$REF_tok,$ALT_tok,$Effect_tok,$ENST_tok);
		foreach my $strand_tok(sort keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}})
		{
			#~ print "INICIO:$string\t$strand_tok\n";
			my @distance_tmp=sort{$a<=>$b}@{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'DISTANCE'}};
			#~ print "El array es:$ENST_tok\t@distance_tmp\n";
			my $TOTAL_distance_ENST = 0;
			for ( @distance_tmp )
			{
				$TOTAL_distance_ENST += $_;
			}
			#~ print "La distancia total es:$ENST_tok\t$TOTAL_distance_ENST\n";
				if($strand_tok eq '+')
				{
					my @positions=();
					my @upstream_distances=();
					my @Begin_tmp=sort{$a<=>$b}keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'BEGIN'}};
					my @End_tmp=sort{$a<=>$b}keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'END'}};
					my $contender=0;
					my $distance_affected_initial=0;
					my $counter=-1;

					# Conditions for splices inside coding region
					
					if($POS_tok > $Begin_tmp[0] && $POS_tok < $Begin_tmp[scalar@Begin_tmp-1])
					{
						if(scalar(@Begin_tmp) == 1)
						{
							print "ERROR!:$string\tsplice_in_monoexonic_gene\n";
						}
						elsif(scalar(@Begin_tmp) > 1)
						{
							my @Begin_shifted=@Begin_tmp;
							#~ print "Array begin inicial:@Begin_tmp\n";
							while(scalar(@Begin_shifted) !=0)
							{
								$contender=shift(@Begin_shifted);
								push(@positions,$contender);
								if($contender <= $POS_tok){$counter++;}
								#~ print "El array begin es:@Begin_shifted\t$counter\n";
								#~ print "El array contender es:@positions\t$counter\n";
							}
							if($counter < scalar(@Begin_tmp)-1)
							{
								for (my $i=$counter+1;$i < scalar(@Begin_tmp);$i++)
								{
									foreach my $end_upstream_tok(sort keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'BEGIN'}{$Begin_tmp[$i]}})
									{
										my $distance=1+$end_upstream_tok-$Begin_tmp[$i];
										#~ print "$distance\t$end_upstream_tok\t$Begin_tmp[$i]\n";
										push(@upstream_distances,$distance);
									}
								}
								my @sorted=sort@upstream_distances;
								#~ print "Las distancias upstream son:@sorted\n";
								my $Affected_upstream_distance = 0;
								for ( @upstream_distances )
								{
									$Affected_upstream_distance += $_;
								}
								#~ print "El total upstream afectado es:$Affected_upstream_distance\n";
								my $distance_affected=$distance_affected_initial+$Affected_upstream_distance;
								#~ print "La distancia afectada total es: $distance_affected_initial\t$Affected_upstream_distance\t$distance_affected\n";
								my $Percentage_affected=100*$distance_affected/$TOTAL_distance_ENST;
								#~ $RESULT_hash{$string}{$Percentage_affected}=1;
								print OUTPUT "$string\t$Percentage_affected\n";
							}
							elsif($counter == scalar(@Begin_tmp)-1)
							{
								print "ERROR!:$string\tsplice_in_last_exon\n";
							}
						}
					}
					
					# Special conditions for splice positions outside coding region
					elsif($POS_tok < $Begin_tmp[0])
					{
						my $Percentage_affected="NaN";
						#~ $RESULT_hash{$string}{$Percentage_affected}=1;
						print OUTPUT "$string\t$Percentage_affected\n";
					}
					elsif($POS_tok > $Begin_tmp[scalar@Begin_tmp-1])
					{
						my $Percentage_affected="NaN";
						#~ $RESULT_hash{$string}{$Percentage_affected}=1;
						print OUTPUT "$string\t$Percentage_affected\n";
					}
				}
				elsif($strand_tok eq '-')
				{
					my @positions=();
					my @upstream_distances=();
					my @End_tmp=reverse sort{$a<=>$b}keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'END'}};
					my $contender=0;
					my $distance_affected_initial=0;
					my $counter=-1;
					
					# Conditions for splices inside coding region
					
					if($POS_tok < $End_tmp[0] && $POS_tok > $End_tmp[scalar@End_tmp-1])
					{
						if(scalar(@End_tmp) == 1)
						{
							print "ERROR!:$string\tsplice_in_monoexonic_gene\n";
						}
						elsif(scalar(@End_tmp) > 1)
						{
							my @End_shifted=@End_tmp;
							#~ print "Array end inicial:@End_tmp\n";
							while(scalar(@End_shifted) !=0)
							{
								$contender=shift(@End_shifted);
								push(@positions,$contender);
								if($contender >= $POS_tok){$counter++;}
								#~ print "El array end es:@End_shifted\t$counter\n";
								#~ print "El array contender es:@positions\t$counter\n";
							}
							if($counter < scalar(@End_tmp)-1)
							{
								for (my $i=$counter+1;$i < scalar(@End_tmp);$i++)
								{
									foreach my $begin_upstream_tok(sort keys%{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'END'}{$End_tmp[$i]}})
									{
										my $distance=1+$End_tmp[$i]-$begin_upstream_tok;
										#~ print "$distance\t$begin_upstream_tok\t$End_tmp[$i]\n";
										push(@upstream_distances,$distance);
									}
								}
								my @sorted=sort@upstream_distances;
								#~ print "Las distancias upstream son:@sorted\n";
								my $Affected_upstream_distance = 0;
								for ( @upstream_distances )
								{
									$Affected_upstream_distance += $_;
								}
								#~ print "El total upstream afectado es:$Affected_upstream_distance\n";
								my $distance_affected=$distance_affected_initial+$Affected_upstream_distance;
								#~ print "La distancia afectada total es: $distance_affected_initial\t$Affected_upstream_distance\t$distance_affected\n";
								my $Percentage_affected=100*$distance_affected/$TOTAL_distance_ENST;
								#~ $RESULT_hash{$string}{$Percentage_affected}=1;
								print OUTPUT "$string\t$Percentage_affected\n";
							}
							elsif($counter == scalar(@End_tmp)-1)
							{
								print "ERROR!:$string\tsplice_in_last_exon\n";
							}
						}
					}
					
					# Special conditions for splice positions outside coding region
					elsif($POS_tok > $End_tmp[0])
					{
						my $Percentage_affected="NaN";
						#~ $RESULT_hash{$string}{$Percentage_affected}=1;
						print OUTPUT "$string\t$Percentage_affected\n";
					}
					elsif($POS_tok < $End_tmp[scalar@End_tmp-1])
					{
						my $Percentage_affected="NaN";
						#~ $RESULT_hash{$string}{$Percentage_affected}=1;
						print OUTPUT "$string\t$Percentage_affected\n";
					}	
						
				}
			}
	}	
	}	
	}	
	}#	
	}	
	}
	}
}





	
sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
