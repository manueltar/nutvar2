##	Script to transfer domain and site ENSEMBL codon coordinates to ENSEMBL genomic coordinates.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;


my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $output=$ARGV[2];

## Open input1 carrying protein information: equivalence ENST-AC_isoform, Domain_IDs, Domain_Coordinates, Postranslational sites and 'Active'sites and
## their coordinates. Position and extension of the features are already transformed to 3-pb codons. 


# Print Time log

my $time='['. timestamp(). ']'."\n";
print "Start charging hash 1:$time\n";
my %hash1=();
my %hash2=();
my %hash3=();
my %hash4=();

if(open (INPUT1, $input1))
{
	while(my $line=<INPUT1>)
	{
		## Input1= Re_mapping/UNIPROT_PLUS_GLOBAL.txt

	#~ A0A183  ENSG00000235942 ENST00000431011 
	#~ A0AUZ9  ENSG00000144445 ENST00000281772 FEATURE:MOD_RES**1__2574__3     FEATURE:IPR026180**1__0__2961   FEATURE:IPR029332**1__2382__363 
	#~ A0AV02  ENSG00000221955 ENST00000393469 FEATURE:IPR004841**1__129__1104 
	#~ A0AV96  ENSG00000163694 ENST00000295971 FEATURE:IPR000504**1__210__225  FEATURE:IPR000504**2__450__249  FEATURE:IPR000504**3__735__219  FEATURE:IPR006535**1__42__1737  FEATURE:IPR012677**1__219__480  FEATU
	#~ A0AVF1  ENSG00000105948 ENST00000464848 FEATURE:REPEAT**1__168__102     FEATURE:REPEAT**2__273__102     FEATURE:REPEAT**3__450__102     FEATURE:REPEAT**4__1401__102    FEATURE:IPR011990**1__99__264   FEATU
	#~ A0AVI4  ENSG00000168936 ENST00000382936 FEATURE:IPR018801**1__48__1035

		my @tmp=split("\t",$line);
		my $ENST=$tmp[2];
		my $ENSG=$tmp[1];
		foreach my $tmp_tok(@tmp)
		{
			if ($tmp_tok=~/^FEATURE:(.+)/)
			{
				my $FEATURE_tok=$1;
				my @FEATURE_tmp=split("__",$FEATURE_tok);
				my $distance_begin=$FEATURE_tmp[1];
				my $distance_feature=$FEATURE_tmp[2];
				my $FEATURE_ID=$FEATURE_tmp[0];
				
				# Here we sort out some features (CHAIN, REGION, REPEAT, DISULFID and CROSSLNK) that are very big, overlap with protein domains
				# and do not fulfil our needs
				
				unless($FEATURE_ID=~/CHAIN/ || $FEATURE_ID=~/REGION/ || $FEATURE_ID=~/REPEAT/ || $FEATURE_ID=~/DISULFID/ || $FEATURE_ID=~/CROSSLNK/)
				{
					$hash1{$ENST}{$ENSG}{$FEATURE_ID}{$distance_begin}{$distance_feature}=1;
					#~ print "$ENST\t$ENSG\t$FEATURE_ID\t$distance_begin\t$distance_feature\n";
				}
			}
		}
	}
}else {print "Unable to open INPUT1\n";}

# Here we parse the file containing the protein cofding transcripts mapped to genomic coordinates on aper base basis.

$time='['. timestamp(). ']'."\n";
print "Tiempo de carga hash_2:$time\n";

if (open(INPUT2, $input2))
{
	#~ Input=Escritorio/Proyecto_clasificador/Raw_Data/TRANSCRIPTS_table.txt
	#~ 
#~ 33772370        1       A3GALT2 -       ENST00000330379 ENST00000442999 
#~ 33772371        1       A3GALT2 -       ENST00000330379 ENST00000442999 
#~ 33772372        1       A3GALT2 -       ENST00000330379 ENST00000442999 
#~ 33772373        1       A3GALT2 -       ENST00000330379 ENST00000442999 
#~ 33772374        1       A3GALT2 -       ENST00000330379 ENST00000442999 
#~ 33772375        1       A3GALT2 -       ENST00000330379 ENST00000442999 
#~ 33772376        1       A3GALT2 -       ENST00000330379 ENST00000442999 
#~ 33772377        1       A3GALT2 -       ENST00000330379 ENST00000442999 

	
	#my @tmp_START_positions=();
	while(my $line=<INPUT2>)
	{
		chomp $line;
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)\t$/)
		{
			my $POS=$1;
			my $CHROM=$2;
			my $SYMBOL=$3;
			my $strand=$4;
			my $ENST_fields=$5;
			
			my @ENST_tmp=split(/\t/,$ENST_fields);
			
			foreach my $ENST_tmp_tok(@ENST_tmp)
			{
				## print "$POS\t$CHROM\t$strand\t$tmp_tok\n";
				
				# Note that we select the trasncripts we have ascertained correspond to the displayed UniProt isoform. 
				
				if(exists($hash1{$ENST_tmp_tok}))
				{
					# We create and array with all the CDS positions for the chosen transcripts of every gen
					$hash2{$ENST_tmp_tok}{$CHROM}{$strand}{$POS}=1;
					#~ print "$ENST_tmp_tok\t$CHROM\t$strand\t$POS\n";	
				}
			}
		}
	}
}else{print "Unable to open INPUT2\n";}

#print "La posición es:$tmp_START_positions[0] \n";

# Now we start the transfer

$time='['. timestamp(). ']'."\n";
print "Tiempo de creación hash3:$time\n";

if(open(OUTPUT, '>'.$output))
{
foreach my $ENST_tok(sort keys%hash1)
{
	my %order_hash=();
	foreach my $ENSG_tok(sort keys%{$hash1{$ENST_tok}})
	{
	foreach my $IPR_OR_SITE_tok(sort keys%{$hash1{$ENST_tok}{$ENSG_tok}})
	{
	foreach my $begin_tok(sort keys%{$hash1{$ENST_tok}{$ENSG_tok}{$IPR_OR_SITE_tok}})
	{
	foreach my $distance_tok(sort keys%{$hash1{$ENST_tok}{$ENSG_tok}{$IPR_OR_SITE_tok}{$begin_tok}})
	{
		foreach my $CHROM_tok(sort keys%{$hash2{$ENST_tok}})
		{
			# The trasnference is strand dependent
			
				foreach my $strand_tok(sort keys%{$hash2{$ENST_tok}{$CHROM_tok}})
				{
					#~ print "$CHROM_tok\t$strand_tok\t$ENST_tok\t$ENSG_tok\tFEATURE:$IPR_OR_SITE_tok\t$begin_tok\t$distance_tok\n";
					
					# Here we retrieve all the coding positions of the transcript and we calculate the total amount of them (max)
					
					my @POS_tmp=sort{ $a <=> $b }(keys%{$hash2{$ENST_tok}{$CHROM_tok}{$strand_tok}});
					#~ print "El array es:@POS_tmp\n";
					#~ print "El scalar es:".scalar(@POS_tmp)."\n";
					my $max=scalar(@POS_tmp);
					 if(scalar(@POS_tmp)==0){print "*******************VOID_ARRAY_IN:$CHROM_tok\t$strand_tok\t$ENST_tok\n";}
					 else
					 {
						if ($strand_tok eq '+')
						{
							# Now we select from the array of positions those covered by the features. We give as the initial index for the array
							# the begin of the feature and the extension is all the indexes smaller than the addition of begin plus distance
							
							for (my $i=$begin_tok; $i<=($begin_tok+$distance_tok-1);$i++)
							{
								# This is just to secure that we won't extend features beyond the coding positions. It is a check point to avoid errors.
								
								unless ($i >= $max)
								{
									#~ print "El índice es:$i\n";
									$order_hash{$POS_tmp[$i]}{$CHROM_tok}{$strand_tok}{$ENST_tok}{$ENSG_tok}{$IPR_OR_SITE_tok}=1;
									#~ print OUTPUT "$POS_tmp[$i]\tchr:$CHROM_tok\t$strand_tok\t$ENST_tok\t$ENSG_tok\t$IPR_OR_SITE_tok\n";
								}
							}
						}
						elsif ($strand_tok eq '-')
						{
							# Note that when the strand is negative we reverse the position array. In this case protein features in the N-ter
							# map to the highest genomic coordinates of the coding positions of the transcript.
							
							my @POS_tmp_reversed=reverse(@POS_tmp);
							
							# Now we select from the array of positions those covered by the features. We give as the initial index for the array
							# the begin of the feature and the extension is all the indexes smaller than the addition of begin plus distance
							
							for (my $i=$begin_tok; $i<=($begin_tok+$distance_tok-1);$i++)
							{
								# This is just to secure that we won't extend features beyond the coding positions. It is a check point to avoid errors
								
								unless ($i >= $max)
								{
									#print "$POS_tmp[$i]\t$CHROM_tok\t$strand_tok\t$ENST_tok\t$AC_tok\t$IPR_ID_OR_SITE_tok\n";
									$order_hash{$POS_tmp_reversed[$i]}{$CHROM_tok}{$strand_tok}{$ENST_tok}{$ENSG_tok}{$IPR_OR_SITE_tok}=1;
									#~ print OUTPUT "$POS_tmp_reversed[$i]\tchr:$CHROM_tok\t$strand_tok\t$ENST_tok\t$ENSG_tok\t$IPR_OR_SITE_tok\n";
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
	
	# Here we print on a per base basis the features mapped to genomic coordinates
	
	foreach my $POS_tok(sort{$a<=>$b}keys%order_hash)
	{
	foreach my $CHROM_tok(sort keys%{$order_hash{$POS_tok}})
	{
	foreach my $strand_tok(sort keys%{$order_hash{$POS_tok}{$CHROM_tok}})
	{
	foreach my $ENST_tok(sort keys%{$order_hash{$POS_tok}{$CHROM_tok}{$strand_tok}})
	{
	foreach my $ENSG_tok(sort keys%{$order_hash{$POS_tok}{$CHROM_tok}{$strand_tok}{$ENST_tok}})
	{
		print OUTPUT "$POS_tok\t$CHROM_tok\t$strand_tok\t$ENST_tok\t$ENSG_tok\t";
		foreach my $Feature_tok(sort keys%{$order_hash{$POS_tok}{$CHROM_tok}{$strand_tok}{$ENST_tok}{$ENSG_tok}})
		{
			print OUTPUT "$Feature_tok\t";
			
		}
		print OUTPUT "\n";
	}
	}	
	}	
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
