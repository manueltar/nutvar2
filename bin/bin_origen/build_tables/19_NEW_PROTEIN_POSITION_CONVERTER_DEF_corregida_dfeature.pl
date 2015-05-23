###############
###############

#!/usr/bin/perl

use strict;
# use warnings;
use Time::localtime;


my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $output=$ARGV[2];

## Open input1 carrying protein information: equivalence ENST-AC_isoform, Domain_IDs, Domain_Coordinates, Postranslational sites and 'Active'sites and
## their coordinates.


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
		## Input1= 

		#	ENSG00000000003 ENST00000373020 FEATURE:IPR000301__168__81      FEATURE:IPR000301__21__714      FEATURE:IPR000301__249__87      FEATUR
		#	ENSG00000000005 ENST00000373031 FEATURE:DISULFID__357__177      FEATURE:IPR007084__276__282     FEATURE:IPR007084__279__279     
		#	ENSG00000000419 ENST00000371588 FEATURE:IPR001173__81__513      FEATURE:IPR029044__27__702      FEATURE:IPR029044__72__603      FEATUR
		#	ENSG00000000457 ENST00000367770 FEATURE:IPR000719__0__735       FEATURE:IPR000719__105__603     FEATURE:IPR011009__99__648      FEATUR
		#	ENSG00000000460 ENST00000286031 FEATURE:IPR027902__525__1659    FEATURE:MOD_RES__2373__3        
		#	ENSG00000000938 ENST00000374003 FEATURE:ACT_SITE__1143__3       FEATURE:BINDING__870__3 FEATURE:IPR000719__786__762     FEATURE:IPR000
		#	ENSG00000001036 ENST00000002165 FEATURE:CARBOHYD__714__3        FEATURE:IPR000933__0__1398      FEATURE:IPR000933__60__1047     FEATUR
		
		#~ A0A183  ENSG00000235942 ENST00000431011 ENSP00000411070 CCDS:CCDS44227.1        
		#~ A0AUZ9  ENSG00000144445 ENST00000281772 ENSP00000281772 FEATURE:N6-acetyllysine.__2574__3__MOD_RES**1   IPR:IPR026180__NSL1__1__0__2961 IPR:IPR029332__PEHE_dom__1__2382__363   CCDS:CCDS33370.1  

		my @tmp=split("\t",$line);
		my $ENST=$tmp[1];
		my $ENSG=$tmp[0];
		foreach my $tmp_tok(@tmp)
		{
			if ($tmp_tok=~/^FEATURE:(.+)/)
			{
				my $FEATURE_tok=$1;
				my @FEATURE_tmp=split("__",$FEATURE_tok);
				my $distance_begin=$FEATURE_tmp[1];
				my $distance_feature=$FEATURE_tmp[2];
				my $FEATURE_ID=$FEATURE_tmp[0];
				unless($FEATURE_ID=~/CHAIN/ || $FEATURE_ID=~/REGION/ || $FEATURE_ID=~/REPEAT/ || $FEATURE_ID=~/DISULFID/ || $FEATURE_ID=~/CROSSLNK/)
				{
					$hash1{$ENST}{$ENSG}{$FEATURE_ID}{$distance_begin}{$distance_feature}=1;
					#~ print "$ENST\t$ENSG\t$FEATURE_ID\t$distance_begin\t$distance_feature\n";
				}
			}
		}
	}
}else {print "Unable to open INPUT1\n";}

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
				foreach my $strand_tok(sort keys%{$hash2{$ENST_tok}{$CHROM_tok}})
				{
				 
				 my @POS_tmp=sort{ $a <=> $b }(keys%{$hash2{$ENST_tok}{$CHROM_tok}{$strand_tok}});
					 if(scalar(@POS_tmp)==0){print "*******************VOID_ARRAY_IN:$CHROM_tok\t$strand_tok\t$ENST_tok\n";}
					 else
					 {
						if ($strand_tok eq '+')
						{
							################	CORRECCIÓN DE LA EXTENSIÓN DEL FEATURE 1 base menos	
							for (my $i=$begin_tok; $i<=($begin_tok+$distance_tok-1);$i++)
							{
								#print "$POS_tmp[$i]\t$CHROM_tok\t$strand_tok\t$ENST_tok\t$AC_tok\t$IPR_ID_OR_SITE_tok\n";
								$order_hash{$POS_tmp[$i]}{$CHROM_tok}{$strand_tok}{$ENST_tok}{$ENSG_tok}{$IPR_OR_SITE_tok}=1;
								#~ print OUTPUT "$POS_tmp[$i]\tchr:$CHROM_tok\t$strand_tok\t$ENST_tok\t$ENSG_tok\t$IPR_OR_SITE_tok\n";
							}
						}
						elsif ($strand_tok eq '-')
						{
							my @POS_tmp_reversed=reverse(@POS_tmp);
							for (my $i=$begin_tok; $i<=($begin_tok+$distance_tok-1);$i++)
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
