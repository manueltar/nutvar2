##	Script to claculate the total domain and site positions affected by the variant. 2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne) 


use strict;
use warnings;
use Time::localtime;

my %hash1=();
my %Feature_hash=();
my %RESULT_hash=();
my %RESULT_hash2=();


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
			
			$hash1{$CHROM}{$SYMBOL}{$ENST}{$Effect}{$POS}{$REF}{$ALT}=1;
			#print "hello_world:$CHROM\t$POS\t$REF\t$ALT\t$Effect\t$SYMBOL\t$ENST\n";
		}
	}
}else{print "Unable to open INPUT1\n";}

$time='['. timestamp(). ']'."\n";
print "Start charging hash2:$time\n";

if (open(INPUT2, $input2))
{
	# Input file= ALL_ISOFORMS_DOMAIN_table_midC.txt
	
#	1       +       AADACL3 ENSG00000188984 ENST00000332530 DOMAIN  [12780885-12780948]     [12785189-12785617]     [12785667-12785956]
#	1       +       AADACL3 ENSG00000188984 ENST00000359318 DOMAIN  [12779513-12779693]     [12780885-12780948]     [12785189-12785617]
#	1       +       AADACL4 ENSG00000204518 ENST00000376221 DOMAIN  [12704566-12704733]     [12711142-12711358]     [12721802-12721865]
#	1       +       ABCD3   ENSG00000117528 ENST00000315713 DOMAIN  [94884053-94884144]     [94924163-94924199]     [94930331-94930429]
#	1       +       ABCD3   ENSG00000117528 ENST00000315713 SITE    [94884038-94884144]     [94924163-94924199]     [94930331-94930429]
#	1       +       ABCD3   ENSG00000117528 ENST00000370214 DOMAIN  [94884053-94884144]     [94924163-94924199]     [94930331-94930429]
#	1       +       ABCD3   ENSG00000117528 ENST00000370214 SITE    [94884038-94884144]     [94924163-94924199]     [94930331-94930429]
#	1       +       ABCD3   ENSG00000117528 ENST00000394233 DOMAIN  [94884053-94884144]     [94924163-94924199]     [94930331-94930429]

	while(my $line=<INPUT2>)
	{
		chomp($line);
		#print "Hello_world_1:$line\n";
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)/)
		{
			my $CHROM=$1;
			my $strand=$2;
			my $SYMBOL=$3;
			my $ENSG=$4;
			my $ENST=$5;
			my $Feature=$6;
			my $fields=$7;
			# print "**$CHROM**$strand**$SYMBOL**$ENST**$Feature**$fields\n";
				if(exists($hash1{$CHROM}{$SYMBOL}{$ENST}))
				{
					# print "Hello_world_1:$CHROM\t$SYMBOL\n";
					my %ordered_coordinates=();
					my @coordinates_tmp=split(/\t/,$fields);
					foreach my $coordinates_tmp_tok(@coordinates_tmp)
					{
						if($coordinates_tmp_tok=~/\[(.+)\-(.+)\]/)
						{
							my $begin=$1;
							my $end=$2;
							$ordered_coordinates{$begin}{$end}=1;
						}
					}
					if($strand eq '+')
						{
							my @begin_tmp=sort{$a<=>$b}keys%ordered_coordinates;
							## print "Hello_world_III:$begin_tmp[$i]\t$Pervasive\n";
							for(my$i=0;$i<scalar(@begin_tmp);$i++)
							{
								foreach my $end_tok(sort{$a<=>$b}keys%{$ordered_coordinates{$begin_tmp[$i]}})
								{
									$Feature_hash{$CHROM}{$strand}{$SYMBOL}{$ENST}{$Feature}{$begin_tmp[$i]}{$end_tok}=1;
									if($SYMBOL eq 'A4GALT')
									{
										print "Feature_hash:$CHROM\t$strand\t$SYMBOL\t$ENST\t$Feature\t$begin_tmp[$i]\t$end_tok\n";
									}
									# print "********************************$CHROM\n";
								}
							}
						}
					elsif($strand eq '-')
						{
							my @begin_tmp=reverse sort{$a<=>$b}keys%ordered_coordinates;
							## print "Hello_world_III:$begin_tmp[$i]\t$Pervasive\n";
							for(my$i=0;$i<scalar(@begin_tmp);$i++)
							{
								foreach my $end_tok(sort{$a<=>$b}keys%{$ordered_coordinates{$begin_tmp[$i]}})
								{
									$Feature_hash{$CHROM}{$strand}{$SYMBOL}{$ENST}{$Feature}{$begin_tmp[$i]}{$end_tok}=1;
									if($SYMBOL eq 'A4GALT')
									{
										print "Feature_hash:$CHROM\t$strand\t$SYMBOL\t$ENST\t$Feature\t$begin_tmp[$i]\t$end_tok\n";
									}
								}
							}
						}
				}
		}
	}
}else {print "impossible to open INPUT2\n";die;}

$time='['. timestamp(). ']'."\n";
print "Start charging PROCESSING:$time\n";

if (open(OUTPUT, '>'.$output))
{
foreach my $CHROM_tok(sort keys %hash1)
{
foreach my $SYMBOL_tok(sort keys %{$hash1{$CHROM_tok}})
{
foreach my $strand_tok(sort keys %{$Feature_hash{$CHROM_tok}})
{
foreach my $ENST_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}})
{					
						foreach my $Effect_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}})
						{
						foreach my $POS_tok(sort {$a<=>$b}keys%{$hash1{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$Effect_tok}})
						{
						foreach my $REF_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$Effect_tok}{$POS_tok}})
						{
						foreach my $ALT_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$Effect_tok}{$POS_tok}{$REF_tok}})
						{
							foreach my $Feature_tok(sort keys %{$Feature_hash{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENST_tok}})
							{
								my %affected_feature_intervals=();
								my %feature_distances=();
								my @POS_affected_tmp=();
								my @POS_feature=();
								my $counter=0;
								my $Affected_percentage=0;
								
								foreach my $begin_tok(sort{$a<=>$b} keys %{$Feature_hash{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENST_tok}{$Feature_tok}})
								{
								foreach my $end_tok(sort{$a<=>$b}  keys %{$Feature_hash{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENST_tok}{$Feature_tok}{$begin_tok}})
								{
									
									my $distance_feature=$end_tok-$begin_tok+1;
									$affected_feature_intervals{'TOTAL'}{$distance_feature}=1;
									
									###################################################################################################3
									
									if($strand_tok eq '+')
									{
										if ($POS_tok > $end_tok)
										{
											# Do nothing
										}
										
										else
										{
											if($POS_tok < $begin_tok)
											{
												my $location_tag=0;
												my $distance_affected=$end_tok-$begin_tok+1;
												$affected_feature_intervals{'AFFECTED'}{$distance_affected}=1;
												$affected_feature_intervals{'TAG'}{$location_tag}=1;
											}
											
											elsif($POS_tok >= $begin_tok)
											{
												my $location_tag=1;
												my $distance_affected=$end_tok-$POS_tok+1;
												$affected_feature_intervals{'AFFECTED'}{$distance_affected}=1;
												$affected_feature_intervals{'TAG'}{$location_tag}=1;
											}
										}
									}
									
									if($strand_tok eq '-')
									{
										if ($POS_tok < $begin_tok)
										{
											# Do nothing
										}
										
										else
										{
											if($POS_tok > $end_tok)
											{
												my $location_tag=0;
												my $distance_affected=$end_tok-$begin_tok+1;
												$affected_feature_intervals{'AFFECTED'}{$distance_affected}=1;
												$affected_feature_intervals{'TAG'}{$location_tag}=1;
											}
											
											elsif($POS_tok <= $end_tok)
											{
												my $location_tag=1;
												my $distance_affected=$POS_tok-$begin_tok+1;
												$affected_feature_intervals{'AFFECTED'}{$distance_affected}=1;
												$affected_feature_intervals{'TAG'}{$location_tag}=1;
											}
										}
									}
								}
								}
								#######################################################################################################################
								my @distance_feature_tmp=sort keys%{$affected_feature_intervals{'TOTAL'}};
								my @distance_feature_affected_tmp=sort keys%{$affected_feature_intervals{'AFFECTED'}};
								unless(scalar(@distance_feature_affected_tmp) == 0)
								{	
									my $TOTAL_feature_POS = 0;
									for ( @distance_feature_tmp )
									{
										$TOTAL_feature_POS += $_;
									}
									
									my $AFFECTED_feature_POS = 0;
									for ( @distance_feature_affected_tmp )
									{
										$AFFECTED_feature_POS += $_;
									}
																								
									my $percentage_of_feature_affected=0;
									unless($AFFECTED_feature_POS==0)
									{
										$percentage_of_feature_affected=100*($AFFECTED_feature_POS/$TOTAL_feature_POS)
									}
									
									#~ if($SYMBOL_tok eq 'A4GALT')
									#~ {
										#~ print "AAA:$CHROM_tok\t$strand_tok\t$SYMBOL_tok\t$ENST_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$Effect_tok\t$Feature_tok\t$AFFECTED_feature_POS***$TOTAL_feature_POS\t$percentage_of_feature_affected:AA\n";
									#~ }
									my $presence=1;
									my $domain_tag=0;
									my @location_tag_tmp=sort keys%{$affected_feature_intervals{'TAG'}};
									# print "TAGTAGTAG: @location_tag_tmp\n";
									if (grep { $presence eq $_ } @location_tag_tmp) 
									{
									   $domain_tag++;
									}
									# print "***$ENST_tok\t$Feature_tok\tArray_total=@distance_feature_tmp\t$TOTAL_feature_POS\t$percentage_of_feature_affected\t$domain_tag\n";
									# print "***$ENST_tok\t$Feature_tok\tArray_afectado=@distance_feature_affected_tmp\t$AFFECTED_feature_POS\t$percentage_of_feature_affected\t$domain_tag\n";
									
									#~ $RESULT_hash{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$Feature_tok}{$percentage_of_feature_affected}{$domain_tag}=1;
									#~ if($SYMBOL_tok eq 'A4GALT')
									#~ {
		print OUTPUT "$CHROM_tok\t$strand_tok\t$SYMBOL_tok\t$ENST_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$Effect_tok\t$Feature_tok\t$percentage_of_feature_affected\t$domain_tag\n";
									#~ }
								}
							}#
						}
						}
						}
						}
}
}
}
}
}

#~ $time='['. timestamp(). ']'."\n";
#~ print "Start printing:$time\n";
#~ 
#~ if (open(OUTPUT, '>'.$output))
#~ {
	#~ foreach my $CHROM_tok(sort{$a<=>$b}keys%RESULT_hash)
	#~ {
	#~ foreach my $strand_tok(sort keys%{$RESULT_hash{$CHROM_tok}})
	#~ {
	#~ foreach my $SYMBOL_tok(sort keys%{$RESULT_hash{$CHROM_tok}{$strand_tok}})
	#~ {
	#~ foreach my $ENST_tok(sort keys%{$RESULT_hash{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}})
	#~ {
	#~ foreach my $POS_tok(sort keys%{$RESULT_hash{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENST_tok}})
	#~ {
	#~ foreach my $REF_tok(sort keys%{$RESULT_hash{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}})
	#~ {
	#~ foreach my $ALT_tok(sort keys%{$RESULT_hash{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}{$REF_tok}})
	#~ {
	#~ foreach my $Effect_tok(sort keys%{$RESULT_hash{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}{$REF_tok}{$ALT_tok}})
	#~ {
		#~ #foreach my $Percentage_ENST_affected_feature_positions(sort keys %{$RESULT_hash2{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}})
		#~ #{
				#~ #print OUTPUT "$CHROM_tok\t$strand_tok\t$SYMBOL_tok\t$ENST_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$Effect_tok\t$Percentage_ENST_affected_feature_positions\t";
				#~ print OUTPUT "$CHROM_tok\t$strand_tok\t$SYMBOL_tok\t$ENST_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$Effect_tok\t";
#~ 
				#~ my @tmp_print=();
				#~ foreach my $Feature_tok(sort keys%{$RESULT_hash{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}})
				#~ {
				#~ foreach my $Affected_percentage_tok(sort keys%{$RESULT_hash{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$Feature_tok}})
				#~ {
				#~ foreach my $location_tag_tok(sort keys%{$RESULT_hash{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$Feature_tok}{$Affected_percentage_tok}})
				#~ {
					#~ my $string1=join("|",$Feature_tok,$Affected_percentage_tok,$location_tag_tok);
					#~ push(@tmp_print,$string1);
					#~ #print OUTPUT ";$Feature_tok|$Affected_percentage_tok";
				#~ }
				#~ }
				#~ }
				#~ my $string2=join(";",@tmp_print);
				#~ print OUTPUT "$string2\n";
		#~ #}
	#~ }
	#~ }
	#~ }
	#~ }		
	#~ }
	#~ }
	#~ }
	#~ }
#~ } else{print "Impossible to open OUTPUT\n";}				
	
sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
