##	Script to calculate the percentage and matched features for every protein domain or functional site. 2014.
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

my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $input3=$ARGV[2];
my $output1=$ARGV[3];

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
print "Start charging hash3 & PRINTING:$time\n";

if (open(INPUT3, $input3) && open(OUTPUT, '>'.$output1))
{
	# Input file= ALL_ISOFORMS_PROTEIN_table_full.txt
	
#	1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR013094**1    [12780885-12780948]     [12785189-12785455]     
#	1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR013094**2    [12785667-12785878]     
#	1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR029058**1    [12780885-12780948]     [12785189-12785617]     
#	1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR029058**2    [12785691-12785956]     
#	1       +       AADACL3 ENSG00000188984 ENST00000359318 IPR013094**1    [12779651-12779693]     [12780885-12780948]     [12785189-12785455]     
#	1       +       AADACL3 ENSG00000188984 ENST00000359318 IPR013094**2    [12785667-12785878]     
#	1       +       AADACL3 ENSG00000188984 ENST00000359318 IPR029058**1    [12779513-12779693]     [12780885-12780948]     [12785189-12785617]     
#	1       +       AADACL3 ENSG00000188984 ENST00000359318 IPR029058**2    [12785691-12785956]     
#	1       +       AADACL4 ENSG00000204518 ENST00000376221 IPR013094**1    [12711316-12711358]     [12721802-12721865]     [12725972-12726292]     
#	1       +       AADACL4 ENSG00000204518 ENST00000376221 IPR013094**2    [12726453-12726658]     
#	1       +       AADACL4 ENSG00000204518 ENST00000376221 IPR017157**1    [12704566-12704733]     [12711142-12711358]     [12721802-12721865]     [12725972-12726742]     
#	1       +       AADACL4 ENSG00000204518 ENST00000376221 IPR029058**1    [12711217-12711358]     [12721802-12721865]     [12725972-12726733]     

 #~ OR
 
 # Input file= ALL_ISOFORMS_DOMAIN_table_full.txt
 
#~ 1       +       AADACL3 ENSG00000188984 ENST00000332530 DOMAIN  [12780885-12780948]     [12785189-12785456]     [12780885-12780948]     [12785189-12785618]     [12785667-12785879]     [12785691-12785957]     
#~ 1       +       AADACL3 ENSG00000188984 ENST00000359318 DOMAIN  [12779513-12779693]     [12780885-12780948]     [12785189-12785618]     [12779651-12779693]     [12780885-12780948]     [12785189-12785456]     [12785667-12785879]     [12785691-12785957]     
#~ 1       +       AADACL4 ENSG00000204518 ENST00000376221 DOMAIN  [12704566-12704733]     [12711142-12711358]     [12721802-12721865]     [12725972-12726743]     [12711217-12711358]     [12721802-12721865]     [12725972-12726734]     [12711316-12711358]     [12721802-12721865]     [12725972-12726293]     [12726453-12726659]     
#~ 1       +       ABCD3   ENSG00000117528 ENST00000315713 DOMAIN  [94884053-94884144]     [94924163-94924199]     [94930331-94930429]     [94933475-94933563]     [94939322-94939391]     [94940699-94940796]     [94941170-94941293]     [94943815-94943871]     [94930358-94930429]     [94933475-94933563]     [94939322-94939391]     [94940699-94940796]     [94941170-94941293]     [94943815-94943871]     [94933487-94933563]     [94939322-94939391]     [94940699-94940796]     [94941170-94941293]     [94943815-94943871]     
#~ 1       +       ABCD3   ENSG00000117528 ENST00000315713 SITE    [94930364-94930366]     

	while(my $line=<INPUT3>)
	{
		chomp($line);
		#print "Hello_world_1:$line\n";
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)\t$/)
		{
			my %hash4=();
			
			my $CHROM=$1;
			my $strand=$2;
			my $SYMBOL=$3;
			my $ENSG=$4;
			my $ENST=$5;
			my $Feature=$6;
			my $fields=$7;
			if($SYMBOL eq 'A4GALT')
			{
				#~ print "**$CHROM**$strand**$SYMBOL**$ENST**$Feature**$fields\n";
			}
			
			my @coordinates_tmp=split(/\t/,$fields);
			foreach my $coordinates_tmp_tok(@coordinates_tmp)
			{
				if($coordinates_tmp_tok=~/\[(.+)\-(.+)\]/)
				{
					my $begin=$1;
					my $end=$2;
					my $distance=$end-$begin+1;
					if($SYMBOL eq 'ACTRT2')
					{
						#~ print "$CHROM\t$SYMBOL\t$ENST\t$strand\t$Feature\t$begin\t$end\t$distance\n";
					}
					# As there might be repeated distance values total distances should be kept in an array rahter than in a hash
					
					push(@{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'DISTANCE'}},$distance);
					$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'BEGIN'}{$begin}{$end}=1;
					#~ print "$CHROM\t$SYMBOL\t$ENST\t$Feature\t$strand\t'BEGIN'\t$begin\t$end\n";
					$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'END'}{$end}{$begin}=1;
					#~ print "$CHROM\t$SYMBOL\t$ENST\t$Feature\t$strand\t'END'\t$end\t$begin\n";
				}
			}
			
			foreach my $POS_tok(sort keys%{$hash1{$CHROM}{$SYMBOL}{$ENST}})
			{
			foreach my $REF_tok(sort keys%{$hash1{$CHROM}{$SYMBOL}{$ENST}{$POS_tok}})
			{
			foreach my $ALT_tok(sort keys%{$hash1{$CHROM}{$SYMBOL}{$ENST}{$POS_tok}{$REF_tok}})
			{
			foreach my $Effect_tok(sort keys%{$hash1{$CHROM}{$SYMBOL}{$ENST}{$POS_tok}{$REF_tok}{$ALT_tok}})
			{
				my $string=join("\t",$CHROM,$strand,$SYMBOL,$ENST,$POS_tok,$REF_tok,$ALT_tok,$Effect_tok,$Feature);
				#~ print "INICIO:$string\n";
				
				my @Begin_tmp=sort{$a<=>$b}keys%{$hash3{$CHROM}{$SYMBOL}{$ENST}{$strand}{'BEGIN'}};
				my @End_tmp=sort{$a<=>$b}keys%{$hash3{$CHROM}{$SYMBOL}{$ENST}{$strand}{'END'}};
				
				#~ print "El array begin del tránscrito es:@Begin_tmp\n";
				#~ print "El array end del tránscrito es:@End_tmp\n";
						
				if($POS_tok >= $Begin_tmp[0] && $POS_tok <= $End_tmp[scalar(@End_tmp)-1])
				{
					# Variant inside coding region
					
					my @distance_tmp=sort{$a<=>$b}@{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'DISTANCE'}};
					#~ print "El array es:$Feature\t@distance_tmp\n";
					my $TOTAL_distance_ENST = 0;
					for ( @distance_tmp )
					{
						$TOTAL_distance_ENST += $_;
					}
					#~ print "La distancia total es:$Feature\t$TOTAL_distance_ENST\n";
						if($strand eq '+')
						{
							my @positions=();
							my @downstream_distances=();
							
							my @Begin_Feature_tmp=sort{$a<=>$b}keys%{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'BEGIN'}};
							my @End_Feature_tmp=sort{$a<=>$b}keys%{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'END'}};
							
							#~ print "El array begin del feature es:@Begin_Feature_tmp\n";
							#~ print "El array end del feature es:@End_Feature_tmp\n";
							
							if($POS_tok >= $Begin_Feature_tmp[0] && $POS_tok <= $End_Feature_tmp[scalar(@End_Feature_tmp)-1])
							{
								my $contender=0;
								my $distance_affected_initial=0;
								my $counter=-1;
								my $matched=1;
								if(scalar(@Begin_Feature_tmp) == 1)
								{
									foreach my $end_tok(sort{$a<=>$b} keys%{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'BEGIN'}{$Begin_Feature_tmp[0]}})
									{
										my $distance_affected=$end_tok-$POS_tok+1;
										#~ print "La distancia afectada es:$end_tok-$POS_tok+1:$distance_affected\n";
										my $Percentage_affected=100*($distance_affected/$TOTAL_distance_ENST);
										print OUTPUT "$string\t$Percentage_affected\t$matched\n";
									}
								}
								elsif(scalar(@Begin_Feature_tmp) > 1)
								{
									my @Begin_shifted=@Begin_Feature_tmp;
									while(scalar(@Begin_shifted) !=0)
									{
										$contender=shift(@Begin_shifted);
										push(@positions,$contender);
										if($contender <= $POS_tok){$counter++;}
										#~ print "El array begin es:@Begin_shifted\t$counter\n";
										#~ print "El array contender es:@positions\t$counter\n";
									}
									foreach my $end_tok_affected(sort keys%{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'BEGIN'}{$Begin_Feature_tmp[$counter]}})
									{
										#~ print "***$Begin_Feature_tmp[$counter]\t$POS_tok\t$end_tok_affected****\n";
										
										# Condition for POS_toks arised from isoforms with coding exons absent in the displayed
										# isoform that interrupt a domain
										
										unless($POS_tok > $end_tok_affected)
										{
											$distance_affected_initial=$end_tok_affected-$POS_tok+1;
											#~ print "distance_affected_initial:$distance_affected_initial\t$end_tok_affected\t$POS_tok\n";
										}
									}
									if($counter < scalar(@Begin_Feature_tmp)-1)
									{
										for (my $i=$counter+1;$i < scalar(@Begin_Feature_tmp);$i++)
										{
											foreach my $end_downstream_tok(sort keys%{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'BEGIN'}{$Begin_Feature_tmp[$i]}})
											{
												my $distance=1+$end_downstream_tok-$Begin_Feature_tmp[$i];
												#~ print "$distance\t$end_downstream_tok\t$Begin_Feature_tmp[$i]\n";
												push(@downstream_distances,$distance);
											}
										}
										#~ print "Las distancias downstream son:@downstream_distances\n";
										my $Affected_downstream_distance = 0;
										for ( @downstream_distances )
										{
											$Affected_downstream_distance += $_;
										}
										#~ print "El total downstream afectado es:$Affected_downstream_distance\n";
										my $distance_affected=$distance_affected_initial+$Affected_downstream_distance;
										#~ print "La distancia afectada total es: $distance_affected_initial\t$Affected_downstream_distance\t$distance_affected\n";
										my $Percentage_affected=100*$distance_affected/$TOTAL_distance_ENST;
										print OUTPUT "$string\t$Percentage_affected\t$matched\n";
									}
									elsif($counter == scalar(@Begin_Feature_tmp)-1)
									{
										my $distance_affected=$distance_affected_initial;
										my $Percentage_affected=100*($distance_affected_initial/$TOTAL_distance_ENST);
										print OUTPUT "$string\t$Percentage_affected\t$matched\n";
									}
								}
							}
							elsif($POS_tok < $Begin_Feature_tmp[0])
							{
								# The domain is entirely downstream the variant
								my $Percentage_affected=100;
								my $matched=0;
								print OUTPUT "$string\t$Percentage_affected\t$matched\n";
								
							}
							elsif($POS_tok > $End_Feature_tmp[scalar(@End_Feature_tmp)-1])
							{
								# Do nothing, the domain is entirely downstream the variant
								
							}
						}
						elsif($strand eq '-')
						{
							my @positions=();
							my @upstream_distances=();
							my @Begin_Feature_tmp=sort{$a<=>$b}keys%{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'BEGIN'}};
							my @End_Feature_tmp=reverse sort{$a<=>$b}keys%{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'END'}};
							
							if($POS_tok >= $Begin_Feature_tmp[0] && $POS_tok <= $End_Feature_tmp[scalar(@End_Feature_tmp)-1])
							{
								my $contender=0;
								my $distance_affected_initial=0;
								my $counter=-1;
								my $matched=1;
								if(scalar(@End_Feature_tmp) == 1)
								{
									#~ print "El array es:$ENST_tok\t@distance_tmp\n";
									#~ print "$string\tEl array begin es:@End_Feature_tmp\n";
									foreach my $begin_tok(sort{$a<=>$b} keys%{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'END'}{$End_Feature_tmp[0]}})
									{
										my $distance_affected=$POS_tok-$begin_tok+1;
										#~ print "La distancia afectada es:$end_tok-$POS_tok+1:$distance_affected\n";
										my $Percentage_affected=100*($distance_affected/$TOTAL_distance_ENST);
										print OUTPUT "$string\t$Percentage_affected\t$matched\n";
									}
								}
								elsif(scalar(@End_Feature_tmp) > 1)
								{
									my @End_shifted=@End_Feature_tmp;
									#~ print "Array begin inicial:@End_Feature_tmp\n";
									while(scalar(@End_shifted) !=0)
									{
										$contender=shift(@End_shifted);
										push(@positions,$contender);
										if($contender >= $POS_tok){$counter++;}
										#~ print "El array end es:@End_shifted\t$counter\n";
										#~ print "El array contender es:@positions\t$counter\n";
									}
									foreach my $begin_tok(sort{$a<=>$b} keys%{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'END'}{$End_Feature_tmp[$counter]}})
									{
										#~ print "***$End_Feature_tmp[$counter]\t$POS_tok\t$end_tok_affected****\n";
										
										# Condition for POS_toks arised from isoforms with coding exons absent in the displayed
										# isoform that interrupt a domain (see issue in documentation)
										
										unless($POS_tok < $begin_tok)
										{
											$distance_affected_initial=$POS_tok-$begin_tok+1;
										}
										#~ print "distance_affected_initial:$distance_affected_initial\t$end_tok_affected\t$POS_tok\n";
									}
									if($counter < scalar(@End_Feature_tmp)-1)
									{
										for (my $i=$counter+1;$i < scalar(@End_Feature_tmp);$i++)
										{
											foreach my $begin_upstream_tok(sort keys%{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'END'}{$End_Feature_tmp[$i]}})
											{
												my $distance=1+$End_Feature_tmp[$i]-$begin_upstream_tok;
												#~ print "$distance\t$end_upstream_tok\t$End_Feature_tmp[$i]\n";
												push(@upstream_distances,$distance);
											}
										}
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
										print OUTPUT "$string\t$Percentage_affected\t$matched\n";
									}
									elsif($counter == scalar(@End_Feature_tmp)-1)
									{
										my $distance_affected=$distance_affected_initial;
										my $Percentage_affected=100*($distance_affected_initial/$TOTAL_distance_ENST);
										print OUTPUT "$string\t$Percentage_affected\t$matched\n";
									}
								}
							}
							elsif($POS_tok < $Begin_Feature_tmp[0])
							{
								# Do nothing, The domain is entirely downstream the variant
							}
							elsif($POS_tok > $End_Feature_tmp[scalar(@End_Feature_tmp)-1])
							{
								# The domain is entirely upstream the variant
								my $Percentage_affected=100;
								my $matched=0;
								print OUTPUT "$string\t$Percentage_affected\t$matched\n";
							}
						}
				}
				elsif($POS_tok < $Begin_tmp[0] || $POS_tok > $End_tmp[scalar(@End_tmp)-1])
				{
					# Variant out of coding region
					
					my $Percentage_affected="NaN";
					my $matched="NaN";
					print OUTPUT "$string\t$Percentage_affected\t$matched\n";
				}
			}
			}
			}
			}
			foreach my $POS_tok(sort keys%{$hash_splice{$CHROM}{$SYMBOL}{$ENST}})
			{
			foreach my $REF_tok(sort keys%{$hash_splice{$CHROM}{$SYMBOL}{$ENST}{$POS_tok}})
			{
			foreach my $ALT_tok(sort keys%{$hash_splice{$CHROM}{$SYMBOL}{$ENST}{$POS_tok}{$REF_tok}})
			{
			foreach my $Effect_tok(sort keys%{$hash_splice{$CHROM}{$SYMBOL}{$ENST}{$POS_tok}{$REF_tok}{$ALT_tok}})
			{
				my $string=join("\t",$CHROM,$strand,$SYMBOL,$ENST,$POS_tok,$REF_tok,$ALT_tok,$Effect_tok,$Feature);

				my @Begin_tmp=sort{$a<=>$b}keys%{$hash3{$CHROM}{$SYMBOL}{$ENST}{$strand}{'BEGIN'}};
				my @End_tmp=sort{$a<=>$b}keys%{$hash3{$CHROM}{$SYMBOL}{$ENST}{$strand}{'END'}};
						
				if($POS_tok >= $Begin_tmp[0] && $POS_tok <= $End_tmp[scalar(@End_tmp)-1])
				{
					# Variant inside coding region
					
					my @distance_tmp=sort{$a<=>$b}@{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'DISTANCE'}};
					#~ print "El array es:$ENST\t@distance_tmp\n";
					my $TOTAL_distance_ENST = 0;
					for ( @distance_tmp )
					{
						$TOTAL_distance_ENST += $_;
					}
					#~ print "La distancia total es:$ENST_tok\t$TOTAL_distance_ENST\n";
						if($strand eq '+')
						{
							my @positions=();
							my @downstream_distances=();
							
							my @Begin_Feature_tmp=sort{$a<=>$b}keys%{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'BEGIN'}};
							my @End_Feature_tmp=sort{$a<=>$b}keys%{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'END'}};
							
							if($POS_tok >= $Begin_Feature_tmp[0] && $POS_tok <= $End_Feature_tmp[scalar(@End_Feature_tmp)-1])
							{
								my $contender=0;
								my $distance_affected_initial=0;
								my $counter=-1;
								my $matched=1;
								if(scalar(@Begin_Feature_tmp) == 1)
								{
									print "ERROR!:$string\tsplice_in_mono_component\n";
								}
								elsif(scalar(@Begin_Feature_tmp) > 1)
								{
									my @Begin_shifted=@Begin_Feature_tmp;
									#~ print "Array begin inicial:@Begin_Feature_tmp\n";
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
										for (my $i=$counter+1;$i < scalar(@Begin_Feature_tmp);$i++)
										{
											foreach my $end_downstream_tok(sort keys%{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'BEGIN'}{$Begin_Feature_tmp[$i]}})
											{
												my $distance=1+$end_downstream_tok-$Begin_Feature_tmp[$i];
												#~ print "$distance\t$end_downstream_tok\t$Begin_Feature_tmp[$i]\n";
												push(@downstream_distances,$distance);
											}
										}
										#~ print "Las distancias downstream son:@downstream_distances\n";
										my $Affected_downstream_distance = 0;
										for ( @downstream_distances )
										{
											$Affected_downstream_distance += $_;
										}
										#~ print "El total downstream afectado es:$Affected_downstream_distance\n";
										my $distance_affected=$distance_affected_initial+$Affected_downstream_distance;
										#~ print "La distancia afectada total es: $distance_affected_initial\t$Affected_downstream_distance\t$distance_affected\n";
										my $Percentage_affected=100*$distance_affected/$TOTAL_distance_ENST;
										print OUTPUT "$string\t$Percentage_affected\t$matched\n";
									}
									elsif($counter == scalar(@Begin_Feature_tmp)-1)
									{
										print "ERROR!:$string\tsplice_in_last_component\n";
									}
								}
							}
							elsif($POS_tok < $Begin_Feature_tmp[0])
							{
								# The domain is entirely downstream the variant
								my $Percentage_affected=100;
								my $matched=0;
								print OUTPUT "$string\t$Percentage_affected\t$matched\n";
								
							}
							elsif($POS_tok > $End_Feature_tmp[scalar(@End_Feature_tmp)-1])
							{
								# Do nothing, the domain is entirely downstream the variant
								
							}
						}
						elsif($strand eq '-')
						{
							my @positions=();
							my @upstream_distances=();
							my @Begin_Feature_tmp=sort{$a<=>$b}keys%{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'BEGIN'}};
							my @End_Feature_tmp=reverse sort{$a<=>$b}keys%{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'END'}};
							
							if($POS_tok >= $Begin_Feature_tmp[0] && $POS_tok <= $End_Feature_tmp[scalar(@End_Feature_tmp)-1])
							{
								my $contender=0;
								my $distance_affected_initial=0;
								my $counter=-1;
								my $matched=1;
								if(scalar(@End_Feature_tmp) == 1)
								{
									print "ERROR!:$string\tsplice_in_mono_component\n";
								}
								elsif(scalar(@End_Feature_tmp) > 1)
								{
									my @End_shifted=@End_Feature_tmp;
									#~ print "Array begin inicial:@End_Feature_tmp\n";
									while(scalar(@End_shifted) !=0)
									{
										$contender=shift(@End_shifted);
										push(@positions,$contender);
										if($contender >= $POS_tok){$counter++;}
										#~ print "El array end es:@End_shifted\t$counter\n";
										#~ print "El array contender es:@positions\t$counter\n";
									}
									
									if($counter < scalar(@End_Feature_tmp)-1)
									{
										for (my $i=$counter+1;$i < scalar(@End_Feature_tmp);$i++)
										{
											foreach my $begin_upstream_tok(sort keys%{$hash4{$CHROM}{$SYMBOL}{$ENST}{$strand}{$Feature}{'END'}{$End_Feature_tmp[$i]}})
											{
												my $distance=1+$End_Feature_tmp[$i]-$begin_upstream_tok;
												#~ print "$distance\t$end_upstream_tok\t$End_Feature_tmp[$i]\n";
												push(@upstream_distances,$distance);
											}
										}
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
										print OUTPUT "$string\t$Percentage_affected\t$matched\n";
									}
									elsif($counter == scalar(@End_Feature_tmp)-1)
									{
										print "ERROR!:$string\tsplice_in_last_component\n";
									}
								}
							}
							elsif($POS_tok < $Begin_Feature_tmp[0])
							{
								# Do nothing, The domain is entirely downstream the variant
							}
							elsif($POS_tok > $End_Feature_tmp[scalar(@End_Feature_tmp)-1])
							{
								# The domain is entirely upstream the variant
								my $Percentage_affected=100;
								my $matched=0;
								print OUTPUT "$string\t$Percentage_affected\t$matched\n";
							}
						}
				}
				elsif($POS_tok < $Begin_tmp[0] || $POS_tok > $End_tmp[scalar(@End_tmp)-1])
				{
					# Variant out of coding region
					
					my $Percentage_affected="NaN";
					my $matched="NaN";
					print OUTPUT "$string\t$Percentage_affected\t$matched\n";
				}
			}
			}
			}
			}
		}
	}
}else {print "impossible to open INPUT3\n";die;}

sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
