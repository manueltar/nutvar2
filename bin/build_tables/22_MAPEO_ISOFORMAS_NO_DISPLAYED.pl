##	Script to transfer the protein features condensed in genomic intervals mapped to the displayed isoform to
## other isoforms non-displayed.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

use strict;
use warnings;
use Time::localtime;

my %hash1=();
my %hash2=();
my %hash3=();
my %hash4=();

my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $input3=$ARGV[2];
my $output=$ARGV[3];

my $time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

# Here we obtain the equivalence SYMBOL-ENSG-strand

if (open (INPUT1, $input1))
{
	#INPUT1=gtf_output_ENSG.txt
	
	#~ ENSG00000186092 OR4F5   +       69091   70008
	#~ ENSG00000237683 AL627309.1      -       134901  139379
	#~ ENSG00000235249 OR4F29  +       367640  368634
	#~ ENSG00000185097 OR4F16  -       621059  622053
	#~ ENSG00000269831 AL669831.1      -       738532  739137
	#~ ENSG00000269308 AL645608.2      +       818043  819983
	
while(my $line=<INPUT1>)
	{
		chomp $line;
		#print "$line\n";
		if($line=~/^(ENSG[^\t]+)\t([^\t]+)\t([^\t]+)\t.+/)
		{
			my $ENSG=$1;
			my $SYMBOL=$2;
			my $strand=$3;
			$hash1{$SYMBOL}{$strand}{$ENSG}=1;
			#~ print "$SYMBOL\t$strand\t$ENSG**\n";
		}
	}
}else{print "Unable to open INPUT1\n";}

# Here we obtain the genomic coordinates intervals for the transcripts in the trasncript corresponding to the displayed isoform of ENSEMBL.

$time='['. timestamp(). ']'."\n";
print "Start charging hash2:$time\n";

if (open(INPUT2, $input2))
{
	## Input file= PROTEIN_full_condensed_feature.txt

#~ 1       +       ENST00000009105 ENSG00000008118 ACT_SITE**1     [209779011-209779013]   
#~ 1       +       ENST00000009105 ENSG00000008118 BINDING**1      [209773389-209773391]   
#~ 1       +       ENST00000009105 ENSG00000008118 IPR000719**1    [209768395-209768412]   [209768413-209768420]   [209773328-209773346]   [209773347-209773388]   [209773389-209773391]   [209773392-209773403]   [209773404-209773456]   [209776559-209776633]   [209778881-209778998]   [209778999-209779010]   [209779011-209779013]   [209779014-209779019]   [209779665-209779682]   [209779683-209779788]   [209781203-209781278]   [209782325-209782437]   [209783196-209783274]   [209784810-209784813] 

	while(my $line=<INPUT2>)
	{
		chomp($line);
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)\t$/)
		{
			my $CHROM=$1;
			my $strand=$2;
			my $ENST=$3;
			my $ENSG=$4;
			my $FEATURE=$5;
			my $fields=$6;
			$hash2{$CHROM}{$strand}{$ENSG}{$ENST}{$FEATURE}{$fields}=1;
			#~ print "$CHROM\t$strand\t$ENSG\t$ENST\t$FEATURE\t$fields**\n";
		}
	}
}else {print "impossible to open INPUT2\n";die;}

# Here we parse all the coding intervals of all the protein coding transcripts in the genome

$time='['. timestamp(). ']'."\n";
print "Start charging hash3:$time\n";

if (open(INPUT3, $input3))
{
	while(my $line=<INPUT3>)
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
				
				$hash3{$CHROM}{$strand}{$SYMBOL}{$ENST_tok}{$interval_tok}=1;
				#~ print "$CHROM\t$strand\t$SYMBOL\t$ENST_tok\t$interval_tok\n";
			}	
			}
		}	
	}
}else {print "impossible to open INPUT3\n";die;}

# Now we start to transfer the features from the transcript correspondent to the displayed isoforms to the rest of non-displayed isoforms

$time='['. timestamp(). ']'."\n";
print "Start printing:$time\n";

if (open(OUTPUT, '>'.$output))
{
foreach my $CHROM_tok(sort keys %hash3)
{
foreach my $strand_tok(sort keys%{$hash3{$CHROM_tok}})
{
foreach my $SYMBOL_tok(sort keys%{$hash3{$CHROM_tok}{$strand_tok}})
{
	foreach my $ENSG_tok(sort keys%{$hash2{$CHROM_tok}{$strand_tok}})
	{
		# If the gene is in autosomal or sexual chromosomes continue. This is a check point.
		
		if(exists($hash1{$SYMBOL_tok}{$strand_tok}{$ENSG_tok}))
		{
			# On a per gene basis we create position hashes with all the coding positions in a transcript
			
			my $string=join("\t",$CHROM_tok,$strand_tok,$SYMBOL_tok,$ENSG_tok);
			my %hash4=();
			my %hash5=();
			
			foreach my $ENST_tok(sort keys%{$hash3{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}})
			{
			foreach my $fields_ENST_tok(sort keys%{$hash3{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENST_tok}})
			{
				my @fields_ENST_tok_tmp=split(/\t/,$fields_ENST_tok);
				foreach my $fields_ENST_tok_tmp_tok(@fields_ENST_tok_tmp)
				{
					if($fields_ENST_tok_tmp_tok=~/\[(.+)\-(.+)\]/)
					{
						my $begin=$1;
						my $end=$2;
						for(my $i=$begin;$i<=$end;$i++)
						{
							# This hash includes all the protein coding transcripts in a gene, the one corresponding to the displayed isoform as well.
							$hash4{$string}{$ENST_tok}{$i}=1;
						}
					}
				}
				
				# On a per gene basis we create position hashes with all the coding positions occupied by the feature
				
				foreach my $FEATURE_tok(sort keys%{$hash2{$CHROM_tok}{$strand_tok}{$ENSG_tok}{$ENST_tok}})
				{
				foreach my $fields_tok(sort keys%{$hash2{$CHROM_tok}{$strand_tok}{$ENSG_tok}{$ENST_tok}{$FEATURE_tok}})
				{
					my @fields_tok_tmp=split(/\t/,$fields_tok);
					foreach my $fields_tok_tmp_tok(@fields_tok_tmp)
					{
						if($fields_tok_tmp_tok=~/\[(.+)\-(.+)\]/)
						{
							my $begin=$1;
							my $end=$2;
							for(my $i=$begin;$i<=$end;$i++)
							{
								# Note here the coding positions are associated with the gene not the transcript
								
								$hash5{$string}{$FEATURE_tok}{$i}=1;
							}
						}
					}
				}	
				}		
			}	
			}
			foreach my $string_tok(sort keys%hash4)
			{
				#print "$string_tok\t";
				# This includes all the protein coding transcripts in the gene, (the one corresponding to the displayed isoform as well)
				foreach my $ENST_tok(sort keys%{$hash4{$string_tok}})
				{
				 foreach my $POS_ENST_tok(sort {$a<=>$b}keys%{$hash4{$string_tok}{$ENST_tok}})
				 {
					foreach my $FEATURE_tok(sort keys %{$hash5{$string_tok}})
					{
						# If the coding position of the transcript exists in the array of positions occupied by the feature
						# we map the feature to the transcript
						
						if(exists($hash5{$string_tok}{$FEATURE_tok}{$POS_ENST_tok}))
						{
							print OUTPUT "$POS_ENST_tok\t$string_tok\t$ENST_tok\t$FEATURE_tok\n";
							
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
}

$time='['. timestamp(). ']'."\n";
print "END:$time\n";
	
sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
