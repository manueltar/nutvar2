##	Script to calculate the new reading frame after a framshift variant and ascertain if there is derived NMD. 2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne) 


use strict;
use warnings;
use Time::localtime;

my $time='['. timestamp(). ']'."\n";
print "Start charging hash chromosomes:$time\n";

my %AcceptedChromosomes=();
my @AcceptChro=('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y');
foreach my $AcceptedChro_tok(@AcceptChro)
{
	$AcceptedChromosomes{$AcceptedChro_tok}=1;
}

my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $output=$ARGV[2];
my %hash1=();
my %RESULTS_hash=();
my %RESULTS_hash_previous_stop_codon=();
my %RESULTS_hash_stop_lost=();

# Here we parse the variant fields

$time='['. timestamp(). ']'."\n";
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
		#~ print "$line\n";
		if($line=~/^([^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)\t[^\t]+\t[^\t]+\t(.+)/)
		{
			my $CHROM=$1;
			my $POS=$2;
			my $REF=$3;
			my $ALT=$4;
			my $fields=$5;
			#~ print "hello_world:$CHROM\t$POS\t$fields\t$REF\t$ALT\n";
			my @fields_tmp=split(";",$fields);
			my $Effect=$fields_tmp[0];
			my $SYMBOL=$fields_tmp[1];
			my $ENST=$fields_tmp[2];
			
			# We restrict to frameshift variants to calculate the new reading frame
			
			if ($Effect =~/frameshift_variant/)
			{
				#~ print "hello_world:$CHROM\t$POS\t$fields\t$REF\t$ALT\n";
				my $length_REF=length($REF);
				my $length_ALT=length($ALT);
				
				# If the alternative allele is shorter than the reference allele -> DELETION
				
				if($length_REF > $length_ALT)
				{
					#~  ATG	->	AC ; 3-2; $length_ALT +1
					my $del_expand=$length_REF - $length_ALT;
					$hash1{$CHROM}{$SYMBOL}{$ENST}{$POS}{$REF}{$ALT}{$Effect}{'DEL'}{$length_ALT}{$del_expand}=1;
					#~ print "$CHROM\t$SYMBOL\t$ENST\t$POS\t$REF\t$ALT\t$Effect\t'DEL'\t$length_ALT\t$del_expand\n";
				}
				
				# If the alternative allele is longer than the reference allele -> INSERTION
				
				elsif($length_REF < $length_ALT)
				{
					my $ins_expand=$length_ALT - $length_REF;
					$hash1{$CHROM}{$SYMBOL}{$ENST}{$POS}{$REF}{$ALT}{$Effect}{'INS'}{$length_REF}{$ins_expand}=1;
					#~ print "$CHROM\t$SYMBOL\t$ENST\t$POS\t$REF\t$ALT\t$Effect\t'INS'\t$length_REF\t$ins_expand\n";
				}
				
			}
		}
	}
}else{print "Unable to open INPUT1\n";}

my $CHROM="NaN";
my $ENSG="NaN";
my $ENST="NaN";
my $SYMBOL="NaN";
my $strand="NaN";
my @seq_tmp=();
my $seq="NaN";
my %position_hash=();

# Here we parse the file with all the positions of the coding sequence and the nucleotide sequence for every transcript

$time='['. timestamp(). ']'."\n";
print "Start charging hash2:$time\n";


if(open (INPUT2, $input2) && open (OUTPUT,'>'.$output))
{
	
	#~ Input:CDS_genomic_coordinates_full_compresed.txt
	#~ 
	#~ >X      ENSG00000071553 ATP6AP1 ENST00000449556 +
	#~ seq>ATGATGGCGGCCATGGCGACGGCTCGAGTGCGGATGGGGCCGCGATGCGCCCAGGCGCTCTGGCGCATGCCGTGGCTGCCGGTGTTTTTGTCGTTGGCGGCGGCGGCGGCGGCGGCAGCGGCGGAGCAGCAGGTCCCGCTGGTGCTGTGGTCGAGTGACCGGGACTTGTGGGCTCCTGCGGCCGACACTCATGAAGGCCAC
	#~ POS>[153657039-153657199]       [153657394-153657520]   [153660176-153660250]   [153660612-153660805]   [153661277-153661311]   [153661981-153662066]   [153662554-153662565]
	#~ //


	while(my $line=<INPUT2>)
	{
		chomp $line;
		#~ print "$line\n";
		
		# If we are in the last line
		
		if ($line=~/^\/\//)
		{
			if(exists($hash1{$CHROM}{$SYMBOL}{$ENST}))
			{
				# We obtain all the coding position
				
				my @base_positions_tmp=sort{$a<=>$b}keys%position_hash;
				my @base_positions_tmp_shifted=@base_positions_tmp;
				
				foreach my $POS_tok(sort{$a<=>$b} keys%{$hash1{$CHROM}{$SYMBOL}{$ENST}})
				{
				foreach my $REF_tok(sort keys%{$hash1{$CHROM}{$SYMBOL}{$ENST}{$POS_tok}})
				{
				foreach my $ALT_tok(sort keys%{$hash1{$CHROM}{$SYMBOL}{$ENST}{$POS_tok}{$REF_tok}})
				{
				foreach my $effect_tok(sort keys%{$hash1{$CHROM}{$SYMBOL}{$ENST}{$POS_tok}{$REF_tok}{$ALT_tok}})
				{
					#~ print ">>$CHROM\t$SYMBOL\t$ENST\t$strand\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\n";
					
				# IN CASE POS is inside coding but is not coding; splice positions; or outside coding positions.
				# We do not calculate downstream PTC for non coding position whatsoever
					
				if(grep /$POS_tok/, @base_positions_tmp)
				{
						#~ print "$POS_tok\n";
					
					# If the variant is a DELETION
					
					foreach my $length_ALT_tok(sort keys%{$hash1{$CHROM}{$SYMBOL}{$ENST}{$POS_tok}{$REF_tok}{$ALT_tok}{$effect_tok}{'DEL'}})
					{
					foreach my $del_expand_tok(sort keys%{$hash1{$CHROM}{$SYMBOL}{$ENST}{$POS_tok}{$REF_tok}{$ALT_tok}{$effect_tok}{'DEL'}{$length_ALT_tok}})
					{
						my @affected_positions=();
						my %modified_sequence=();
						my %codons=();
						
						# We locate the affected positions and its extension and put them in an array
								
						my $begin_deletion=$POS_tok+$length_ALT_tok;
						for(my $i=$begin_deletion;$i<=$begin_deletion+$del_expand_tok-1;$i++)
						{
							push(@affected_positions,$i);
						}
						#~ print "Posiones de partida:@base_positions_tmp\n";
						#~ print "Las posiciones afectadas son:@affected_positions\n";
						#~ print "Las posiciones no afectadas son:";
									
						foreach my $base_positions_tmp_tok(@base_positions_tmp)
						{
							#~ print "$base_positions_tmp_tok\t";
							
							# We discard all the bases and nucleotides afected for the deletion
							
							unless(grep /$base_positions_tmp_tok/, @affected_positions)
							{
								foreach my $nucleotide_tok(sort keys%{$position_hash{$base_positions_tmp_tok}})
								{
									$modified_sequence{'POS'}{$base_positions_tmp_tok}=1;
									push(@{$modified_sequence{'NUCLEOTIDE'}},$nucleotide_tok);
								}
							}
						}
						
						# Now we obtain the new arrays of nucleotide and positions without the affected positions and nucleotides
						
						my @tmp_pos_modified_sequence=sort{$a<=>$b}keys%{$modified_sequence{'POS'}};
						#~ print "*********************@tmp_pos_modified_sequence\n";
						my @tmp_nucleotide_modified_sequence=@{$modified_sequence{'NUCLEOTIDE'}};
						#~ print "Modified_sequence>".join("",@tmp_nucleotide_modified_sequence)."\n";
						
						# Here we calculate if the reading frame is 1-based displaced, 2-based displaced or the original frame
									
						my $remainder=scalar(@tmp_pos_modified_sequence) % 3;
						my $FLAG=0;
						
						# Original frame
						
						if($remainder == 0)
						{ 
							# Do nothing
						}
						
						# 1-based displaced
						
						elsif($remainder == 1){$FLAG=1;}
						
						# 2-based displaced
						
						elsif($remainder == 2){$FLAG=2;}
						#~ print "$remainder\t$FLAG\n";
						
						# The calculation of the new reading frame is strand dependent
						
						if($strand eq '+')
						{
							# We join the nucleotides in triplets
							
							# We advance every three positions
							
							for(my $h=0;$h<scalar(@tmp_pos_modified_sequence);$h+=3)
							{
								if($FLAG == 0)
								{
									# Join the three previous nucleotides
									
									#~ print "$h\t".scalar(@tmp_pos_modified_sequence)."\n";
									my $codon=join("",$tmp_nucleotide_modified_sequence[$h],$tmp_nucleotide_modified_sequence[$h+1],$tmp_nucleotide_modified_sequence[$h+2]);
									#~ print "CODON:$codon\t$tmp_pos_modified_sequence[$h]\n";
									
									# We locate the stop codons
									
									if($codon eq 'TAG'||$codon eq 'TAA'||$codon eq 'TGA')
									{
										$codons{'STOP'}{$tmp_pos_modified_sequence[$h]}{$codon}=1;
										#~ print "*********************STOP:$codon\n";
									}			
								}
								elsif($FLAG == 1)
								{
									# Join the three previous nucleotides
									
									if($h < scalar(@tmp_pos_modified_sequence) -1)
									{
										#~ print "$h\t".scalar(@tmp_pos_modified_sequence)."\n";
										my $codon=join("",$tmp_nucleotide_modified_sequence[$h],$tmp_nucleotide_modified_sequence[$h+1],$tmp_nucleotide_modified_sequence[$h+2]);
										#~ print "CODON:$codon\t$tmp_pos_modified_sequence[$h]\n";
										
										# We locate the stop codons
										
										if($codon eq 'TAG'||$codon eq 'TAA'||$codon eq 'TGA')
										{
											$codons{'STOP'}{$tmp_pos_modified_sequence[$h]}{$codon}=1;
											#~ print "*********************STOP:$codon\n";
										}
									}
									elsif ($h >= scalar(@tmp_pos_modified_sequence) -1)
									{
										# Do nothing
									}
								}
								elsif($FLAG == 2)
								{
									# Join the three previous nucleotides
									
									if($h < scalar(@tmp_pos_modified_sequence) -2)
									{
										my $codon=join("",$tmp_nucleotide_modified_sequence[$h],$tmp_nucleotide_modified_sequence[$h+1],$tmp_nucleotide_modified_sequence[$h+2]);
										#~ print "CODON:$codon\t$tmp_pos_modified_sequence[$h]\n";
										
										# We locate the stop codons
										
										if($codon eq 'TAG'||$codon eq 'TAA'||$codon eq 'TGA')
										{
											$codons{'STOP'}{$tmp_pos_modified_sequence[$h]}{$codon}=1;
											#~ print "*********************STOP:$codon\n";
										}
									}
									elsif ($h >= scalar(@tmp_pos_modified_sequence) -2)
									{
										# Do nothing
									}			
								}	
							}
							my @STOP_tmp=sort{$a<=>$b} keys %{$codons{'STOP'}};
							
							# If there is a derived STOP, print the first element of the array
										
							if(scalar(@STOP_tmp) > 0)
							{
								#~ print "$CHROM\t$SYMBOL\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\t$ENST\t->DERIVED_PTC:$STOP_tmp[0]\n";
								print OUTPUT "$CHROM\t$SYMBOL\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\t$ENST\t->DERIVED_PTC:$STOP_tmp[0]\n";	
							}
							
							# If there is no derived STOP print stop lost
							
							elsif(scalar(@STOP_tmp) == 0)
							{
								#~ print "$CHROM\t$SYMBOL\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\t$ENST\t->STOP_LOST\n";
								print OUTPUT "$CHROM\t$SYMBOL\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\t$ENST\t->STOP_LOST\n";
							}
						}
						
						# We repeat the same scheme as in the positive strand, BUT this time we use the reverse the arrays of positions and nucleotides
						
						elsif($strand eq '-')
						{
							my @tmp_pos_reversed=reverse @tmp_pos_modified_sequence;
							my @tmp_nucleotide_reversed=reverse @tmp_nucleotide_modified_sequence;			
							for(my $h=0;$h<scalar(@tmp_pos_reversed);$h+=3)
							{
								if($FLAG == 0)
								{
									#~ print "$h\t".scalar(@tmp_pos_reversed)."\n";
									my $codon=join("",$tmp_nucleotide_reversed[$h],$tmp_nucleotide_reversed[$h+1],$tmp_nucleotide_reversed[$h+2]);
									#~ print "CODON:$codon\t$tmp_pos_reversed[$h]\n";
									if($codon eq 'TAG'||$codon eq 'TAA'||$codon eq 'TGA')
									{
										$codons{'STOP'}{$tmp_pos_reversed[$h]}{$codon}=1;
										#~ print "*********************STOP:$codon\n";
									}
								}
								elsif($FLAG == 1)
								{
									if($h < scalar(@tmp_pos_reversed) -1)
									{
										#~ print "$h\t".scalar(@tmp_pos_reversed)."\n";
										my $codon=join("",$tmp_nucleotide_reversed[$h],$tmp_nucleotide_reversed[$h+1],$tmp_nucleotide_reversed[$h+2]);
										#~ print "CODON:$codon\t$tmp_pos_reversed[$h]\n";
										if($codon eq 'TAG'||$codon eq 'TAA'||$codon eq 'TGA')
										{
											$codons{'STOP'}{$tmp_pos_reversed[$h]}{$codon}=1;
											#~ print "*********************STOP:$codon\n";
										}
									}
									elsif ($h >= scalar(@tmp_pos_reversed) -1)
									{
										# Do nothing
									}
								}
								elsif($FLAG == 2)
								{
									if($h < scalar(@tmp_pos_reversed) -2)
									{
										my $codon=join("",$tmp_nucleotide_reversed[$h],$tmp_nucleotide_reversed[$h+1],$tmp_nucleotide_reversed[$h+2]);
										#~ print "CODON:$codon\t$tmp_pos_reversed[$h]\n";
										if($codon eq 'TAG'||$codon eq 'TAA'||$codon eq 'TGA')
										{
											$codons{'STOP'}{$tmp_pos_reversed[$h]}{$codon}=1;
											#~ print "*********************STOP:$codon\n";
										}
									}
									elsif ($h >= scalar(@tmp_pos_reversed) -2)
									{
										# Do nothing
									}
								}	
						}
						my @STOP_tmp=sort{$a<=>$b} keys %{$codons{'STOP'}};				
						if(scalar(@STOP_tmp) > 0)
						{
							#~ print "$CHROM\t$SYMBOL\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\t$ENST\t->DERIVED_PTC:$STOP_tmp[0]\n";
							print OUTPUT "$CHROM\t$SYMBOL\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\t$ENST\t->DERIVED_PTC:$STOP_tmp[0]\n";
						}
						elsif(scalar(@STOP_tmp) == 0)
						{
							#~ print "$CHROM\t$SYMBOL\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\t$ENST\t->STOP_LOST\n";
							print OUTPUT "$CHROM\t$SYMBOL\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\t$ENST\t->STOP_LOST\n";
						}
					}
				}
				} # DELETION
				
				# Now in the case of INSERTIONS we repeat the same scheme
				
				foreach my $length_REF_tok(sort keys%{$hash1{$CHROM}{$SYMBOL}{$ENST}{$POS_tok}{$REF_tok}{$ALT_tok}{$effect_tok}{'INS'}})
				{
				foreach my $ins_expand_tok(sort keys%{$hash1{$CHROM}{$SYMBOL}{$ENST}{$POS_tok}{$REF_tok}{$ALT_tok}{$effect_tok}{'INS'}{$length_REF_tok}})
				{
						
						my @affected_positions=();
						my %modified_sequence=();
						my %codons=();
								
						my $begin_insertion=$POS_tok+$length_REF_tok-1;
						push(@affected_positions,$begin_insertion);

						#~ print "Posiones de partida:@base_positions_tmp\n";
						#~ print "Las posiciones afectadas son:@affected_positions\n";
						#~ print "Las posiciones no afectadas son:";
									
						foreach my $base_positions_tmp_tok(@base_positions_tmp)
						{
							#~ print "$base_positions_tmp_tok\t";
							unless(grep /$base_positions_tmp_tok/, @affected_positions)
							{
								#~ print "********************$base_positions_tmp_tok\t";
								foreach my $nucleotide_tok(sort keys%{$position_hash{$base_positions_tmp_tok}})
								{
									push(@{$modified_sequence{'NUCLEOTIDE'}},$nucleotide_tok);
									push(@{$modified_sequence{'POS'}},$base_positions_tmp_tok);
								}
							}
							
							# As it is an INSERTION we add the inderted nucleotides with the genomic coordinate corresponding to the point of insertion 
							
							elsif(grep /$base_positions_tmp_tok/, @affected_positions)
							{
								#~ print "********************$base_positions_tmp_tok\t";
								my @ALT_split=split("",$ALT_tok);
								#~ print "La REF es:@ALT_split\n";
								for(my$i=$length_REF_tok-1;$i<scalar(@ALT_split);$i++)
								{
									push(@{$modified_sequence{'POS'}},$base_positions_tmp_tok);
									push(@{$modified_sequence{'NUCLEOTIDE'}},$ALT_split[$i]);
									#~ print "Base_inserted:$ALT_split[$i]\n";
								}
							}
						}
						
						# From now onwards we repeat the same scheme as for deletions, if the stop codon is found within the position
						# of the insertion (insertion carrying within the STOP codon) it will be reported to that position
						
						my @tmp_pos_modified_sequence=sort @{$modified_sequence{'POS'}};
						#~ print "****POS*@tmp_pos_modified_sequence\n";
						#~ print "****POS*".scalar(@tmp_pos_modified_sequence)."\n";
						my @tmp_nucleotide_modified_sequence=@{$modified_sequence{'NUCLEOTIDE'}};
						#~ print"****NUCLEOTIDE_modified>".join("",@tmp_nucleotide_modified_sequence)."\n";
						#~ print"****NUCLEOTIDE*". scalar(@tmp_nucleotide_modified_sequence)."\n";
									
						my $remainder=scalar(@tmp_nucleotide_modified_sequence) % 3;
						my $FLAG=0;
						if($remainder == 0)
						{ 
							# Do nothing
						}
						elsif($remainder == 1){$FLAG=1;}
						elsif($remainder == 2){$FLAG=2;}
						#~ print "$remainder\t$FLAG\n";
						if($strand eq '+')
						{
							for(my $h=0;$h<scalar(@tmp_pos_modified_sequence);$h+=3)
							{
								if($FLAG == 0)
								{
									#~ print "$h\t".scalar(@tmp_pos_modified_sequence)."\n";
									my $codon=join("",$tmp_nucleotide_modified_sequence[$h],$tmp_nucleotide_modified_sequence[$h+1],$tmp_nucleotide_modified_sequence[$h+2]);
									#~ print "CODON:$codon\t$tmp_pos_modified_sequence[$h]\n";
									if($codon eq 'TAG'||$codon eq 'TAA'||$codon eq 'TGA')
									{
										$codons{'STOP'}{$tmp_pos_modified_sequence[$h]}{$codon}=1;
										#~ print "*********************STOP:$codon\n";
									}			
								}
								elsif($FLAG == 1)
								{
									if($h < scalar(@tmp_pos_modified_sequence) -1)
									{
										#~ print "$h\t".scalar(@tmp_pos_modified_sequence)."\n";
										my $codon=join("",$tmp_nucleotide_modified_sequence[$h],$tmp_nucleotide_modified_sequence[$h+1],$tmp_nucleotide_modified_sequence[$h+2]);
										#~ print "CODON:$codon\t$tmp_pos_modified_sequence[$h]\n";
										if($codon eq 'TAG'||$codon eq 'TAA'||$codon eq 'TGA')
										{
											$codons{'STOP'}{$tmp_pos_modified_sequence[$h]}{$codon}=1;
											#~ print "*********************STOP:$codon\n";
										}
									}
									elsif ($h >= scalar(@tmp_pos_modified_sequence) -1)
									{
										# Do nothing
									}
								}
								elsif($FLAG == 2)
								{
									if($h < scalar(@tmp_pos_modified_sequence) -2)
									{
										my $codon=join("",$tmp_nucleotide_modified_sequence[$h],$tmp_nucleotide_modified_sequence[$h+1],$tmp_nucleotide_modified_sequence[$h+2]);
										#~ print "CODON:$codon\t$tmp_pos_modified_sequence[$h]\n";
										if($codon eq 'TAG'||$codon eq 'TAA'||$codon eq 'TGA')
										{
											$codons{'STOP'}{$tmp_pos_modified_sequence[$h]}{$codon}=1;
											#~ print "*********************STOP:$codon\n";
										}
									}
									elsif ($h >= scalar(@tmp_pos_modified_sequence) -2)
									{
										# Do nothing
									}			
								}	
							}
							my @STOP_tmp=sort{$a<=>$b} keys %{$codons{'STOP'}};
										
							if(scalar(@STOP_tmp) > 0)
							{
								#~ print "$CHROM\t$SYMBOL\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\t$ENST\t->DERIVED_PTC:$STOP_tmp[0]\n";
								print OUTPUT "$CHROM\t$SYMBOL\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\t$ENST\t->DERIVED_PTC:$STOP_tmp[0]\n";	
							}
							elsif(scalar(@STOP_tmp) == 0)
							{
								#~ print "$CHROM\t$SYMBOL\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\t$ENST\t->STOP_LOST\n";
								print OUTPUT "$CHROM\t$SYMBOL\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\t$ENST\t->STOP_LOST\n";
							}
						}
						elsif($strand eq '-')
						{
							my @tmp_pos_reversed=reverse @tmp_pos_modified_sequence;
							my @tmp_nucleotide_reversed=reverse @tmp_nucleotide_modified_sequence;			
							for(my $h=0;$h<scalar(@tmp_pos_reversed);$h+=3)
							{
								if($FLAG == 0)
								{
									#~ print "$h\t".scalar(@tmp_pos_reversed)."\n";
									my $codon=join("",$tmp_nucleotide_reversed[$h],$tmp_nucleotide_reversed[$h+1],$tmp_nucleotide_reversed[$h+2]);
									#~ print "CODON:$codon\t$tmp_pos_reversed[$h]\n";
									if($codon eq 'TAG'||$codon eq 'TAA'||$codon eq 'TGA')
									{
										$codons{'STOP'}{$tmp_pos_reversed[$h]}{$codon}=1;
										#~ print "*********************STOP:$codon\n";
									}
								}
								elsif($FLAG == 1)
								{
									if($h < scalar(@tmp_pos_reversed) -1)
									{
										#~ print "$h\t".scalar(@tmp_pos_reversed)."\n";
										my $codon=join("",$tmp_nucleotide_reversed[$h],$tmp_nucleotide_reversed[$h+1],$tmp_nucleotide_reversed[$h+2]);
										#~ print "CODON:$codon\t$tmp_pos_reversed[$h]\n";
										if($codon eq 'TAG'||$codon eq 'TAA'||$codon eq 'TGA')
										{
											$codons{'STOP'}{$tmp_pos_reversed[$h]}{$codon}=1;
											#~ print "*********************STOP:$codon\n";
										}
									}
									elsif ($h >= scalar(@tmp_pos_reversed) -1)
									{
										# Do nothing
									}
								}
								elsif($FLAG == 2)
								{
									if($h < scalar(@tmp_pos_reversed) -2)
									{
										my $codon=join("",$tmp_nucleotide_reversed[$h],$tmp_nucleotide_reversed[$h+1],$tmp_nucleotide_reversed[$h+2]);
										#~ print "CODON:$codon\t$tmp_pos_reversed[$h]\n";
										if($codon eq 'TAG'||$codon eq 'TAA'||$codon eq 'TGA')
										{
											$codons{'STOP'}{$tmp_pos_reversed[$h]}{$codon}=1;
											#~ print "*********************STOP:$codon\n";
										}
									}
									elsif ($h >= scalar(@tmp_pos_reversed) -2)
									{
										# Do nothing
									}
								}	
						}
						my @STOP_tmp=sort{$a<=>$b} keys %{$codons{'STOP'}};				
						if(scalar(@STOP_tmp) > 0)
						{
							#~ print "$CHROM\t$SYMBOL\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\t$ENST\t->DERIVED_PTC:$STOP_tmp[0]\n";
							print OUTPUT "$CHROM\t$SYMBOL\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\t$ENST\t->DERIVED_PTC:$STOP_tmp[0]\n";
						}
						elsif(scalar(@STOP_tmp) == 0)
						{
							#~ print "$CHROM\t$SYMBOL\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\t$ENST\t->STOP_LOST\n";
							print OUTPUT "$CHROM\t$SYMBOL\t$POS_tok\t$REF_tok\t$ALT_tok\t$effect_tok\t$ENST\t->STOP_LOST\n";
						}
					}
				}
				}
			}
			}
			}
			}	
			}
				
			}##exist
			%position_hash=();
			@seq_tmp=();
		}
		#~ Input:CDS_genomic_coordinates_full_compresed.txt
	#~ 
	#~ >X      ENSG00000071553 ATP6AP1 ENST00000449556 +
	#~ seq>ATGATGGCGGCCATGGCGACGGCTCGAGTGCGGATGGGGCCGCGATGCGCCCAGGCGCTCTGGCGCATGCCGTGGCTGCCGGTGTTTTTGTCGTTGGCGGCGGCGGCGGCGGCGGCAGCGGCGGAGCAGCAGGTCCCGCTGGTGCTGTGGTCGAGTGACCGGGACTTGTGGGCTCCTGCGGCCGACACTCATGAAGGCCAC
	#~ POS>[153657039-153657199]       [153657394-153657520]   [153660176-153660250]   [153660612-153660805]   [153661277-153661311]   [153661981-153662066]   [153662554-153662565]
	#~ //
		
		elsif ($line=~/^>([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			$CHROM=$1;
			$ENSG=$2;
			$SYMBOL=$3;
			$ENST=$4;
			$strand=$5;
			#~ print ">$CHROM\t$ENSG\t$SYMBOL\t$ENST\t$strand\n";
			#~ if(exists($hash1{$CHROM}{$SYMBOL}{$ENST})){print "Hello_world:>$CHROM\t$ENSG\t$SYMBOL\t$ENST\t$strand\n";}
		}
		
		# Here we parse the sequence
		
		if(exists($hash1{$CHROM}{$SYMBOL}{$ENST}))
		{
			if($line =~ /^seq>(.+)/)
			{
				$seq=$1;
				@seq_tmp=split("",$seq);
				if($strand eq '-')
				{
					@seq_tmp=reverse(@seq_tmp);
				}
				
				#~ print "El array es:@seq_tmp\n";
			}
			
			# Here we parse the positions
			
			elsif ($line=~/^POS>(.+)/)
			{
				my @POS_def=();
				
				my $POS=$1;
				#~ print "****************************$POS\n";
				my @POS_tmp=split(/\t/,$POS);
				#~ print "****************************".join("**",@POS_tmp)."\n";
				foreach my $POS_tmp_tok(@POS_tmp)
				{
					if($POS_tmp_tok=~/\[(.+)\-(.+)\]/)
					{
						my $begin=$1;
						my $end=$2;
						my $distance=$end-$begin+1;
						for(my $h=$begin;$h<$begin+$distance;$h++)
						{
							push(@POS_def,$h)
						}	
					}
				}
				for(my$i=0;$i<scalar(@POS_def);$i++)
				{
					$position_hash{$POS_def[$i]}{$seq_tmp[$i]}=1;
					#~ print "$POS_def[$i]\t$seq_tmp[$i]\n";
				}
			}
		}
	}
}else{print "Unable to open INPUT2\n";}

$time='['. timestamp(). ']'."\n";
print "Print FIN:$time\n";


sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
