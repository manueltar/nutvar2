##	Script to build a Pre-C with DNA sequence and genomic positions of all the CDS positions in the transcript in ENSEMBL.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

use strict;
use warnings;
use Time::localtime;
use Memory::Usage;

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

# Here we obtain the CDS sequence

$time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";
my $CHROM_2="NaN";
my $ENSG_2="NaN";
my $ENST_2="NaN";
my $ENSG_2_biotype="NaN";
my $ENST_2_biotype="NaN";
my %hash_CDS=();
if(open (INPUT, $input1))
{	
	#~ Input:Homo_sapiens.GRCh37.75.cds.all.fa
	#~ 
	#~ >ENST_200000324659 cds:known chromosome:GRCh37:16:3612901:3627401:-1 gene:ENSG_200000167984 gene_biotype:polymorphic_pseudogene transcript_biotype:protein_coding
	#~ CTCAGCCTGCCCTCAGTCACCTATCTGCTCCTGGAGGTGATCCCCGACTCCATGAGGAAG
	#~ CAAGAGGTGCGGACGGGCAGGGAGGCCGGCCAGGGCCACGGTACGGGCTCCCCAGCCGAG
	#~ CAGGTGAAAGCCCTCATGGATCTGCTGGCTGGGAAGGGCAGTCAAGGCTCCCAGGCCCCG
	#~ CAGGCCCTGGATAGGACACCGGATGCCCCGCTGGGGCCCTGCAGCAATGACTCAAGGATA
	#~ CAGAGGCACCGCAAGGCCCTGCTGAGCAAGGTGGAGGCCACCCGCGGGGGCGGGCACCCC
	#~ GCCAGGACCGTCGCCCTGGACCGGCTCTTCCTGCCTCTCTCCCGGGTGTCTGTCCCACCC
	#~ CGGGTCTCCATCACTATCGGGGTGGCCGGCATGGGCAAGACCACCCTGGTGAGGCACTTC
	#~ GTCCGCCTCTGGGCCCATGGGCAGGTCGGCAAGGACTTCTCGCTGGTGCTGCCTCTGACC
	
	while(my $line=<INPUT>)
	{
		chomp $line;
		if ($line=~/^>([^\s]+)\s+(.+)/)
		{
			$ENST_2=$1;
			#~ print "El  ENST_2 es:$ENST_2\n";
			my $fields=$2;
			my @tmp=split(/\s+/,$fields);
			foreach my $tmp_tok(@tmp)
			{
				if($tmp_tok =~ /chromosome:(.+)/)
				{
					my $sub_fields=$1;
					my @sub_fields_tmp=split(":",$sub_fields);
					$CHROM_2=$sub_fields_tmp[1];
					#~ print "El CHROM_2 es:$CHROM_2\n";
					#~ print "El CHROM_2 es:$CHROM_2\t$begin_CDS\t$end_CDS\t$orientation\n";
				}
				
				elsif($tmp_tok =~ /gene:(.+)/)
				{
					$ENSG_2=$1;
					#~ print "El ENSG_2 es:$ENSG_2\n";
				}
				elsif($tmp_tok =~ /gene_biotype:(.+)/)
				{
					$ENSG_2_biotype=$1;
					#~ print "El ENSG_2:_biotype es:$ENSG_2_biotype\n";
				}
				elsif($tmp_tok =~ /transcript_biotype:(.+)/)
				{
					$ENST_2_biotype=$1;
					#~ print "El ENST_2_biotype es:$ENST_2_biotype\n";
				}
			}
		}
		else
		{
			#~ print "**************$CHROM_2\t$ENSG_2\t$ENST_2\n";
			if(exists($AcceptedChromosomes{$CHROM_2}))
			{
				#~ if($ENSG_2_biotype eq 'protein_coding' && $ENST_2_biotype eq 'protein_coding')
				#~ {
						#~ print "$CHROM_2\t$ENSG_2\t$ENST_2\t'SEQ'\t$line\n";
						push(@{$hash_CDS{$CHROM_2}{$ENSG_2}{$ENST_2}{'SEQ'}},$line);
				#~ }
			}
		}
	}
}
my $CHROM="NaN";
my $ENSG="NaN";
my $ENST="NaN";
my $SYMBOL="NaN";
my $strand="NaN";
my $seq="NaN";
my %position_hash=();

$time='['. timestamp(). ']'."\n";
print "Start charging hash2:$time\n";

# From this file we obtain the cDNA sequence and the genomic coordinates occupied by it.

if(open (INPUT2, $input2) && open(OUTPUT,'>'.$output))
{
	
	#~ Input:cDNA_genomic_coordinates.txt
	#~ 
	#~ >X      ENSG00000000003 TSPAN6  ENST00000373020 -
	#~ seq>AGTTGTGGACGCTCGTAAGTTTTCGGCAGTTTCCGGGGAGACTCGGGGACTCCGCGTCTCGCTCTCTGTGTTCCAATCGCCCGGTGCGGTGGTGCAGGGTCTCGGGCTAGTCATGGCGTCCCCGTCTCGGAGACTGCAGACTAAACCAGTCATTACTTGTTTCAAGAGCGTTCTGCTAATCTACACTTTTATTTTCTGGATCACTGGCGTTATCCTTCTTGCAGTTGGCATTTGGGGCAAGGTGAGCCTGGAGAATTACTTTTCTCTTTTAAATGAGAAGGCCACCAATGTCCCCTTCGTGCTCATTGCTACTGGTACCGTCATTATTCTTTTGGGCACCTTTGGTTGTTTTGCTACCTGCCGAGCTTCTGCATGGATGCTAAAACTGTATGCAATGTTTCTGACTCTCGTTTTTTTGGTCGAACTGGTCGCTGCCATCGTAGGATTTGTTTTCAGACATGAGATTAAGAACAGCTTTAAGAATAATTATGAGAAGGCTTTGAAGCAGTATAACTCTACAGGAGATTATAGAAGCCATGCAGTAGACAAGATCCAAAATACGTTGCATTGTTGTGGTGTCACCGATTATAGAGATTGGACAGATACTAATTATTACTCAGAAAAAGGATTTCCTAAGAGTTGCTGTAAACTTGAAGATTGTACTCCACAGAGAGATGCAGACAAAGTAAACAATGAAGGTTGTTTTATAAAGGTGATGACCATTATAGAGTCAGAAATGGGAGTCGTTGCAGGAATTTCCTTTGGAGTTGCTTGCTTCCAACTGATTGGAATCTTTCTCGCCTACTGCCTCTCTCGTGCCATAACAAATAACCAGTATGAGATAGTGTAACCCAATGTATCTGTGGGCCTATTCCTCTCTACCTTTAAGGACATTTAGGGTCCCCCCTGTGAATTAGAAAGTTGCTTGGCTGGAGAACTGACAACACTACTTACTGATAGACCAAAAAACTACACCAGTAGGTTGATTCAATCAAGATGTATGTAGACCTAAAACTACACCAATAGGCTGATTCAATCAAGATCCGTGCTCGCAGTGGGCTGATTCAATCAAGATGTATGTTTGCTATGTTCTAAGTCCACCTTCTATCCCATTCATGTTAGATCGTTGAAACCCTGTATCCCTCTGAAACACTGGAAGAGCTAGTAAATTGTAAATGAAGTAATACTGTGTTCCTCTTGACTGTTATTTTTCTTAGTAGGGGGCCTTTGGAAGGCACTGTGAATTTGCTATTTTGATGTAGTGTTACAAGATGGAAAATTGATTCCTCTGACTTTGCTATTGATGTAGTGTGATAGAAAATTCACCCCTCTGAACTGGCTCCTTCCCAGTCAAGGTTATCTGGTTTGATTGTATAATTTGCACCAAGAAGTTAAAATGTTTTATGACTCTCTGTTCTGCTGACAGGCAGAGAGTCACATTGTGTAATTTAATTTCAGTCAGTCAATAGATGGCATCCCTCATCAGGGTTGCCAGATGGTGATAACAGTGTAAGGCCTTGGGTCTAAGGCATCCACGACTGGAAGGGACTACTGATGTTCTGTGATACATCAGGTTTCAGCACACAACTTACATTTCTTTGCCTCCAAATTGAGGCATTTATTATGATGTTCATACTTTCCCTCTTGTTTGAAAGTTTCTAATTATTAAATGGTGTCGGAATTGTTGTATTTTCCTTAGGAATTCAGTGGAACTTATCTTCATTAAATTTAGCTGGTACCAGGTTGATATGACTTGTCAATATTATGGTCAACTTTAAGTCTTAGTTTTCGTTTGTGCCTTTGATTAATAAGTATAACTCTTATACAATAAATACTGCTTTCCTCTAAAAAGATCGTGTTTAAATTAACTTGTAGAAAATCTGCTGGAATGGTTGTTGTTTTCCACTGAGAAAGCTAAGCCCTACATTTCTATTCAGAGTACTGTTTTTAGATGTGAAATATAAGCCTGCGGCCTTAACTCTGTATTAAAAAAAATGTTTTTGTTTAAAAAAAACTGTTCCCATAGGTGCAGCAAACCACCATGGCACATGTATACCTATGTAACAAACCTGCACATTCTGCACATGTATCCCAGAACTTAATGTAAACAAAAAAATCTTAAAGTGCAAATATTAAAAAAAACTGTTCTCTGTGAAAAAAATTATATTCCATGTTATAAAGTAGCATATGACTAGTGTTCTCCTAG
	#~ seq_length>2206
	#~ POS>    99883667__1317  99885756__108   99887482__84    99888402__135   99888928__99    99890175__75    99890555__189   99891605__199   
	#~ total_distance>2206
	#~ //
	while(my $line=<INPUT2>)
	{
		chomp $line;
		#~ print "$line\n";
		# If we are in the last line of the register
		
		if ($line=~/^\/\//)
		{
			my %base_position_hash=();
				print OUTPUT ">$CHROM\t$ENSG\t$SYMBOL\t$ENST\t$strand\n";
				my @CDS_seq=@{$hash_CDS{$CHROM}{$ENSG}{$ENST}{'SEQ'}};
				#~ print join("**",@CDS_seq)."\n";
				my $seq_CDS=join("",@CDS_seq);
				#~ print "CDS>$seq_CDS\n";
				#~ print "cDNA>$seq\n";
				my $upstream=0;
				my $downstream=0;
				my $CDS="NaN";
				
				# The cDNA always includes the CDS but ENSEMBL adjust the reading frame adding N or NN to the beggining of the CDS nucleotide seq.
				# Therefore in these cases the CDS containing N's does not match the cDNA. Hence the else in this conditional IF
				
				if($seq =~/$seq_CDS/)
				{
					#~ print "Hello_world\n";
					
					# In case there are no UTRs
					if($seq eq $seq_CDS)
					{
						# Do nothing
					}
					
					# In case there are both 5 and 3 prime UTRS
					
					elsif($seq =~/(.+)$seq_CDS(.+)/)
					{
						my $upstream_seq=$1;
						$upstream=length($upstream_seq);
						my $downstream_seq=$2;
						$downstream=length($downstream_seq);
					}
					# Only 5 prime UTR
					elsif($seq =~/(.+)$seq_CDS/)
					{
						my $upstream_seq=$1;
						$upstream=length($upstream_seq);
					}
					# Only 3 prime UTR
					elsif($seq =~/$seq_CDS(.+)/)
					{
						my $downstream_seq=$1;
						$downstream=length($downstream_seq);
					}
					else
					{
						print "************ERROR2**************$CHROM\t$ENSG\t$SYMBOL\t$ENST\t$strand\n";
					}
					$CDS=$seq_CDS;
				}
				
				# This is for the N cases mentioned in the IF
				
				else
				{
					print ">************ERROR1**************$CHROM\t$ENSG\t$SYMBOL\t$ENST\t$strand\tcDNA\n";
					print "$seq\n";
					print ">************ERROR1**************$CHROM\t$ENSG\t$SYMBOL\t$ENST\t$strand\tCDS\n";
					print "$seq_CDS\n";
					
					my @CDS_split_tmp=split("",$seq_CDS);
					my @CDS_split_tmp_shifted=@CDS_split_tmp;
					my $counter=0;
					my @corrected_sequence=();
					while(scalar(@CDS_split_tmp_shifted) != 0)
					{
						my $contender=shift(@CDS_split_tmp_shifted);
						$counter++;
						# This is to eliminate 'N' letters introduced by ENSEMBL to adjust the reading frame
						if($contender eq 'N')
						{
							print "$contender\t$counter\n";
						}
						else{push(@corrected_sequence,$contender);}
					}
					
					my $seq_CDS_corrected=join("",@corrected_sequence);
					print ">************ERROR1**************$CHROM\t$ENSG\t$SYMBOL\t$ENST\t$strand\tseq_corrected\n";
					print "$seq_CDS_corrected\n";

					if($seq =~/$seq_CDS_corrected/)
					{
						# Here as mentioned above we ascertain the length of UTRs
						
						#~ print "Hello_world\n";
						if($seq eq $seq_CDS_corrected)
						{
							# Do nothing
							#~ print "Hello_world_1\n";
						}
						elsif($seq =~/(.+)$seq_CDS_corrected(.+)/)
						{
							my $upstream_seq=$1;
							$upstream=length($upstream_seq);
							my $downstream_seq=$2;
							$downstream=length($downstream_seq);
							#~ print "Hello_world_2:$upstream_seq\t$upstream\t$downstream_seq\t$downstream\n";
						}
						elsif($seq =~/(.+)$seq_CDS_corrected/)
						{
							my $upstream_seq=$1;
							$upstream=length($upstream_seq);
							#~ print "Hello_world_3:$upstream_seq\t$upstream\n";
						}
						elsif($seq =~/$seq_CDS_corrected(.+)/)
						{
							my $downstream_seq=$1;
							$downstream=length($downstream_seq);
							#~ print "Hello_world_4:$downstream_seq\t$downstream\n";
						}
						else
						{
							print "************ERROR2**************$CHROM\t$ENSG\t$SYMBOL\t$ENST\t$strand\n";
						}
						$CDS=$seq_CDS_corrected;
					}
					else
					{
						print "Still_error1\n";
					}
				}
				
				# Hash of positions
				
				my @positions=sort{$a<=>$b}keys%position_hash;
				
				
				# Eliminate positions upstream the CDS (as many positions as the length of the UTR)
				
				for(my $i=1;$i<=$upstream;$i++)
				{
					shift(@positions);
				}
				for(my $j=1;$j<=$downstream;$j++)
				{
					pop(@positions);
				}
				
				my $check=scalar(@positions);
				my $CDS_length=length($CDS);
				if($CDS_length ne $check){print "************ERROR3**************$CHROM\t$ENSG\t$SYMBOL\t$ENST\t$strand\n";}
				#~ print "El array es:@positions\n";
				#~ print "La secuencia es:$CDS\n";
				
				# Now we print the sequence; this is strand dependent.
				
				if($strand eq '+')
				{
					my @sequence_split=split("",$CDS);
					for(my$i=0;$i<scalar(@sequence_split);$i++)
					{
						print OUTPUT "$positions[$i]"."__"."$sequence_split[$i]\t";
					}
				}
				elsif($strand eq '-')
				{
					my @sequence_split=split("",$CDS);
					my @sequence_split_reversed=reverse @sequence_split;
					for(my$i=0;$i<scalar(@sequence_split_reversed);$i++)
					{
						print OUTPUT "$positions[$i]"."__"."$sequence_split_reversed[$i]\t";
					}
				}
				print OUTPUT "\n";	
				print OUTPUT "//\n";	
		%position_hash=();
		}
		
		# Here we obtain the IDs
		
		elsif ($line=~/^>([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			$CHROM=$1;
			$ENSG=$2;
			$SYMBOL=$3;
			$ENST=$4;
			$strand=$5;
			#~ print "$CHROM\t$ENSG\t$SYMBOL\t$ENST\t$strand\n";
		}
		
		# Here we parse the sequence
		
		elsif ($line=~/^seq>(.+)/)
		{
			$seq=$1;
			#~ print "$seq\n";
		}
		
		# Here we parse the begin positions and the distance and converte them to a hash of positions
		
		elsif ($line=~/^POS>(.+)\t$/)
		{
			my $POS=$1;
			#~ print "****************************$POS\n";
			my @POS_tmp=split(/\t/,$POS);
			#~ print "****************************".join("**",@POS_tmp)."\n";
			foreach my $POS_tmp_tok(@POS_tmp)
			{
				my @tmp=split("__",$POS_tmp_tok);
				#~ print "****************************".join("**",@tmp)."\n";
				my $begin=$tmp[0];
				my $distance=$tmp[1];
				for(my $i=$begin;$i<$begin+$distance;$i++)
				{
					$position_hash{$i}=1;
					#~ print "$i\n";
				}
			}
		}
	}
}else{print "Unable to open INPUT2\n";}

$time='['. timestamp(). ']'."\n";
print "THE END:$time\n";

sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
