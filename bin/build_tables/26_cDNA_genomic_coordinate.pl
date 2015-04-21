##	Script to build a Pre-C with DNA sequence and genomic positions of all the cDNA in ENSEMBL.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

use strict;
use warnings;
use Time::localtime;



my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $output=$ARGV[2];
my %hash1=();
my %RESULTS_hash=();
my %RESULTS_hash_previous_stop_codon=();
my %RESULTS_hash_stop_lost=();


my $time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

my %hash_begin=();

# Here we extract genomic coordinates of all the exons in ENSEMBL protein coding transcripts

if (open(INPUT1, $input1))
{
	## Input file= NMD_table.txt

# X       +       ABCD1   ENSG00000101986 ENST00000218104 10      [153008943]-[153010216] [152990323]-[152991621] [152994687]-[152994867] [153001566]-[153001708] [153001799]-[153001967] [153002611]-[153002705] [153005546]-[153005691] [153006028]-[153006173] [153008441]-[153008525] [153008675]-[153008800] 
# X       +       ABCD1   ENSG00000101986 ENST00000370129 2       [152991272]-[152991621] [152994687]-[152995352] 
# X       +       ABCD1   ENSG00000101986 ENST00000443684 6       [153000798]-[153000881] [153001566]-[153001708] [153001799]-[153001967] [153002611]-[153002705] [153005546]-[153005691] [153006028]-[153006058] 
# X       +       AC002365.1      ENSG00000234469 ENST00000445307 2       [9935392]-[9936005]     [9936107]-[9936134]     
# X       +       ACRC    ENSG00000147174 ENST00000373695 12      [70830531]-[70830669]   [70832205]-[70832409]   [70832712]-[70833432]
#   [70800135]-[70800730]   [70811972]-[70812018]   [70812387]-[70812420]   [70814181]-[70814227]   [70814596]-[70814629]   [70817800]-[70817888]   [70823438]-[70824526]   [70825513]-[70825579]   [70828823]-[70828967]   
# X       +       ACRC    ENSG00000147174 ENST00000373696 13      [70828823]-[70828967]   [70830531]-[70830669]   [70832205]-[70832409]
#   [70832712]-[70833433]   [70798261]-[70798372]   [70800670]-[70800730]   [70811972]-[70812018]   [70812387]-[70812420]   [70814181]-[70814227]   [70814596]-[70814629]   [70817800]-[70817888]   [70823438]-[70824526]   [70825513]-[70825579]   
# X       +       AFF2    ENSG00000155966 ENST00000286437 18      [148044245]-[148044467] [148048320]-[148048609] [148049159]-[148049222] [148055001]-[148055137] [148059463]-[148059534] [148059892]-[148059985] [148062268]-[148062320] [148068897]-[148069087] [148072741]-[148082193] [147800615]-[147800755] [147891400]-[147891444] [147924490]-[147924526] [147924906]-[147924957] [147967419]-[147967515] [147985751]-[147985788] [148035110]-[148035269] [148037133]-[148038143] [148039867]-[148039988] 
# X       +       AFF2    ENSG00000155966 ENST00000342251 20      [148037133]-[148038143] [148039867]-[148039988] [148044245]-[148044467] [148048320]-[148048609] [148049159]-[148049222] [148055001]-[148055137] [148059463]-[148059534] [148059892]-[148059985] [148062268]-[148062320] [148068897]-[148069087] [147582228]-[147582664] [148072741]-[148075954] [147733520]-[147733640] [147743429]-[147744289] [147891400]-[147891444] [147924490]-[147924526] [147924906]-[147924957] [147967419]-[147967515] [147985751]-[147985788] [148035110]-[148035269] 
# X       +       AFF2    ENSG00000155966 ENST00000370457 20      [148037133]-[148038143] [148039867]-[148039988] [148044245]-[148044461] [148048320]-[148048609] [148049159]-[148049222] [148055001]-[148055137] [148059463]-[148059534] [148059892]-[148059985] [148062268]-[148062320] [148068897]-[148069087] [147582139]-[147582664] [148072741]-[148082193] [147733520]-[147733640] [147743429]-[147744289] [147891400]-[147891444] [147924490]-[147924526] [147924906]-[147924957] [147967419]-[147967515] [147985751]-[147985788] [148035110]-[148035269] 

	while(my $line=<INPUT1>)
	{
		chomp($line);
		#~ print "$line\n";
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)/)
		{
			#~ print "HELLO_WORLD:$line\n";
			my $CHROM=$1;
			my $strand=$2;
			my $SYMBOL=$3;
			my $ENSG=$4;
			my $ENST=$5;
			my $exon_number=$6;
			my $fields=$7;
				
			my @fields_tmp=split(/\t/,$fields);
			foreach my $fields_tmp_tok(@fields_tmp)
			{
				if($fields_tmp_tok=~/\[(.+)\]\-\[(.+)\]/)
				{
					my $begin=$1;
					my $end=$2;
					$hash_begin{$CHROM}{$ENSG}{$ENST}{'BEGIN'}{$SYMBOL}{$strand}{$begin}{$end}=1;
					#~ print "$CHROM\t$ENSG\t$ENST\n";
				}
			}
		}
	}
}

#~ exit;

$time='['. timestamp(). ']'."\n";
print "Start charging hash2:$time\n";

my $CHROM="NaN";
my $ENSG="NaN";
my $ENST="NaN";
my $orientation="NaN";
my $begin_CDS="NaN";
my $end_CDS="NaN";
my $ENSG_biotype="NaN";
my $ENST_biotype="NaN";
my %sequence_hash=();

# Here we open a file from ENSEMBL containing all the cDNA nucleotide sequences

if(open (INPUT2, $input2))
{
	#~ Input:Homo_sapiens.GRCh37.75.cdna.all.fa
	#~ 
	#~ >ENST00000567516 cdna:novel chromosome:GRCh37:HG706_PATCH:75612096:75635873:-1 gene:ENSG00000261530 gene_biotype:protein_coding transcript_biotype:protein_coding
	#~ GTCTAGTGATCCTTCACTGTGTGGTGGCAGATGGGAATTCCACCAGAAGTCCTGAAACTA
	#~ ATGGCCTCCTCTGTGGAGACCCTGAGGAAAACTGTGCAGCTACCACCACACAATCAAAGC
	#~ GGAAAGGCCACTTCTCTAGGTGCCCCAAGCAATACAAGCATTACTGCATCAAAGGGAGAT
	#~ GCCGCTTCGTGGTGGCCGAGCAGACGCCCTCCTGTGTCCCTCTTCGGAAACGTCGTAAAA
	#~ GAAAGAAGAAAGAAGAAGAAATGGAAACTCTGGGTAAAGATATAACTCCTATCAATGAAG
	#~ ATATTGAAGAGACAAATATTGCTTAAAAGGCTATGAAGTTACCTCCAGGTTGGTGGCAAG
	#~ CTGCAAAGTGCCTTGCTCATTTGAAAATGGACAGAATGTGTCTCAGGAAAACAGCTAGTA
	#~ GACATGAATTTTAAATAATGTATTTACTTTTTATTTGCAACTTTAGTTTGTGTTATTATT
	#~ TTTTAATAAGAACATTAATTATATGTATATTGTCTAGTAATTGGGAAAAAAGCAACTGGT
	#~ TAGGTAGCAACAACAGAAGGGAAATTTCAATAACCTTTCACTTAAGTATTGTCACCAGGA
	#~ TTACTAGTCAAACAAAAAAGAAAA

	
	while(my $line=<INPUT2>)
	{
		chomp $line;
		if ($line=~/^>([^\s]+)\s+(.+)/)
		{
			$ENST=$1;
			#~ print "El  ENST es:$ENST\n";
			my $fields=$2;
			my @tmp=split(/\s+/,$fields);
			foreach my $tmp_tok(@tmp)
			{
				if($tmp_tok =~ /chromosome:(.+)/)
				{
					my $sub_fields=$1;
					my @sub_fields_tmp=split(":",$sub_fields);
					$CHROM=$sub_fields_tmp[1];
					$begin_CDS=$sub_fields_tmp[2];
					$end_CDS=$sub_fields_tmp[3];
					$orientation=$sub_fields_tmp[4];
					#~ print "El CHROM es:$CHROM\n";
					#~ print "El CHROM es:$CHROM\t$begin_CDS\t$end_CDS\t$orientation\n";
				}
				
				elsif($tmp_tok =~ /gene:(.+)/)
				{
					$ENSG=$1;
					#~ print "El ENSG es:$ENSG\n";
				}
				elsif($tmp_tok =~ /gene_biotype:(.+)/)
				{
					$ENSG_biotype=$1;
					#~ print "El ENSG:_biotype es:$ENSG_biotype\n";
				}
				elsif($tmp_tok =~ /transcript_biotype:(.+)/)
				{
					$ENST_biotype=$1;
					#~ print "El ENST_biotype es:$ENST_biotype\n";
				}
			}
		}
		else
		{
			#~ print "**************$CHROM\t$ENSG\t$ENST\n";
			
			# Here we select only the ones of protein coding transcripts (hash_begin)
			
			if(exists($hash_begin{$CHROM}{$ENSG}{$ENST}))
			{
				#~ if($ENSG_biotype eq 'protein_coding' && $ENST_biotype eq 'protein_coding')
				#~ {
						#~ print "$CHROM\t$ENSG\t$ENST\t'SEQ'\t$begin_CDS\t$end_CDS\t$orientation\t$line\n";
						push(@{$hash_begin{$CHROM}{$ENSG}{$ENST}{'SEQ'}{$begin_CDS}{$end_CDS}{$orientation}},$line);
				#~ }
			}
		}
	}
}

# Here we print a file with the cDNA nucleotide sequence and the positions it occupies in terms of beginEXON_distanceEXON


$time='['. timestamp(). ']'."\n";
print "Start PROCESSING:$time\n";

if(open(OUTPUT,'>'.$output))
{
foreach my $CHROM_tok(sort {$a<=>$b} keys %hash_begin)
{
foreach my $ENSG_tok(sort keys %{$hash_begin{$CHROM_tok}})
{	
foreach my $ENST_tok(sort keys %{$hash_begin{$CHROM_tok}{$ENSG_tok}})
{	
foreach my $begin_cDNA_tok(sort keys %{$hash_begin{$CHROM_tok}{$ENSG_tok}{$ENST_tok}{'SEQ'}})
{	
foreach my $end_cDNA_tok(sort keys %{$hash_begin{$CHROM_tok}{$ENSG_tok}{$ENST_tok}{'SEQ'}{$begin_cDNA_tok}})
{	
foreach my $orientation_cDNA_tok(sort keys %{$hash_begin{$CHROM_tok}{$ENSG_tok}{$ENST_tok}{'SEQ'}{$begin_cDNA_tok}{$end_cDNA_tok}})
{	
	my @seq=@{$hash_begin{$CHROM_tok}{$ENSG_tok}{$ENST_tok}{'SEQ'}{$begin_cDNA_tok}{$end_cDNA_tok}{$orientation_cDNA_tok}};
	my $seq=join("",@seq);
	my $length_seq=length($seq);
	#~ print "$seq\n";
	#~ print "//\n";
	
	foreach my $SYMBOL_tok(sort keys %{$hash_begin{$CHROM_tok}{$ENSG_tok}{$ENST_tok}{'BEGIN'}})
	{
	foreach my $strand_tok(sort keys %{$hash_begin{$CHROM_tok}{$ENSG_tok}{$ENST_tok}{'BEGIN'}{$SYMBOL_tok}})
	{
		print OUTPUT ">$CHROM_tok\t$ENSG_tok\t$SYMBOL_tok\t$ENST_tok\t$strand_tok\n";
		print OUTPUT "seq>$seq\n";
		print OUTPUT "seq_length>$length_seq\n";
		
		# Now; first exclude non_coding positions
		
		my @positions=();
		my @positions_denied=();
		my @distance=();
	
		my @begin_tmp=sort{$a<=>$b}keys%{$hash_begin{$CHROM_tok}{$ENSG_tok}{$ENST_tok}{'BEGIN'}{$SYMBOL_tok}{$strand_tok}};
		my @begin_tmp_shifted=@begin_tmp;
		print OUTPUT "POS>";
		foreach my $begin_tok(@begin_tmp)
		{
			foreach my $end_tok(sort{$a<=>$b}keys%{$hash_begin{$CHROM_tok}{$ENSG_tok}{$ENST_tok}{'BEGIN'}{$SYMBOL_tok}{$strand_tok}{$begin_tok}})
			{
				my $distance=$end_tok-$begin_tok+1;
				print OUTPUT "$begin_tok"."__"."$distance\t";
				push(@distance,$distance);		
			}
		}
		print OUTPUT  "\n";
		my $TOTAL_distance= 0;
			for ( @distance )
			{
				$TOTAL_distance += $_;
			}
		print OUTPUT "total_distance>$TOTAL_distance\n";
		print OUTPUT "//\n";
		
		# Here we check that the length of the nucleotide sequence and the total amount of positions are equal.
		
		if ($TOTAL_distance ne $length_seq)
		{
			print "ERROR in correspondence:$CHROM_tok\t$ENSG_tok\t$SYMBOL_tok\t$ENST_tok\t$strand_tok\n";
		}
	}
	}
}		
}		
}	
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
