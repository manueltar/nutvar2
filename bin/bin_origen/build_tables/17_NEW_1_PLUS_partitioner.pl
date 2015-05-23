###############
###############

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;


my $input1=$ARGV[0];
my $output1=$ARGV[1];

## Open input1 carrying protein information: equivalence ENST-AC_isoform, Domain_IDs, Domain_Coordinates, Postranslational sites and 'Active'sites and
## their coordinates.


# Print Time log

my $time='['. timestamp(). ']'."\n";
print "Start charging hash 1:$time\n";
my %hash1=();

if(open (INPUT1, $input1))
{
	while(my $line=<INPUT1>)
	{
		## Input1= UNIPROT_GLOBAL_POST_UNIQUE_PLUS_COORDINATES_rereferenced_2.txt

		# P00519	ENSG00000097007	ENST00000318560	ENSP00000323315	FEATURE:ATP.__741__27__NP_BIND**1	FEATURE:ATP.__810__3__BINDING**1	FEATURE:ATP.__945__21__NP_BIND**2	FEATURE:CAP.__0__180__REGION**1	FEATURE:DNA-binding (By similarity).__2604__300__REGION**2	FEATURE:F-actin-binding.__2856__534__REGION**3	FEATURE:N6-acetyllysine; by EP300.__2130__3__MOD_RES**25	FEATURE:Phosphoserine (By similarity).__1335__3__MOD_RES**16	FEATURE:Phosphoserine.__147__3__MOD_RES**1	FEATURE:Phosphoserine.__1704__3__MOD_RES**18	FEATURE:Phosphoserine.__1857__3__MOD_RES**22	FEATURE:Phosphoserine.__1974__3__MOD_RES**23	FEATURE:Phosphoserine.__2046__3__MOD_RES**24	FEATURE:Phosphoserine.__2151__3__MOD_RES**26	FEATURE:Phosphoserine.__2412__3__MOD_RES**29	FEATURE:Phosphoserine.__2424__3__MOD_RES**30	FEATURE:Phosphoserine.__2562__3__MOD_RES**34	FEATURE:Phosphoserine.__2748__3__MOD_RES**35	FEATURE:Phosphoserine.__2754__3__MOD_RES**36	FEATURE:Phosphoserine.__2805__3__MOD_RES**37	FEATURE:Phosphoserine.__2844__3__MOD_RES**38	FEATURE:Phosphoserine.__2928__3__MOD_RES**39	FEATURE:Phosphoserine; by PAK2.__1851__3__MOD_RES**20	FEATURE:Phosphoserine; by PAK2.__1854__3__MOD_RES**21	FEATURE:Phosphothreonine.__1173__3__MOD_RES**13	FEATURE:Phosphothreonine.__1179__3__MOD_RES**15	FEATURE:Phosphothreonine.__1836__3__MOD_RES**19	FEATURE:Phosphothreonine.__2202__3__MOD_RES**27	FEATURE:Phosphothreonine.__2340__3__MOD_RES**28	FEATURE:Phosphothreonine.__2439__3__MOD_RES**31	FEATURE:Phosphothreonine.__2529__3__MOD_RES**32	FEATURE:Phosphothreonine.__2553__3__MOD_RES**33	FEATURE:Phosphotyrosine.__1404__3__MOD_RES**17	FEATURE:Phosphotyrosine.__342__3__MOD_RES**3	FEATURE:Phosphotyrosine.__381__3__MOD_RES**4	FEATURE:Phosphotyrosine.__414__3__MOD_RES**5	FEATURE:Phosphotyrosine.__513__3__MOD_RES**6	FEATURE:Phosphotyrosine.__552__3__MOD_RES**7	FEATURE:Phosphotyrosine.__642__3__MOD_RES**8	FEATURE:Phosphotyrosine.__756__3__MOD_RES**10	FEATURE:Phosphotyrosine.__768__3__MOD_RES**11	FEATURE:Phosphotyrosine.__789__3__MOD_RES**12	FEATURE:Phosphotyrosine; by autocatalysis and__1176__3__MOD_RES**14	FEATURE:Phosphotyrosine; by autocatalysis.__207__3__MOD_RES**2	FEATURE:Phosphotyrosine; by autocatalysis.__675__3__MOD_RES**9	FEATURE:Proton acceptor (By similarity).__1086__3__ACT_SITE**1	IPR:IPR000719__Prot_kinase_dom__1__723__756	IPR:IPR000980__SH2__1__354__432	IPR:IPR001245__Ser-Thr/Tyr_kinase_cat_dom__1__723__753	IPR:IPR001452__SH3_domain__1__180__261	IPR:IPR008266__Tyr_kinase_AS__1__1074__39	IPR:IPR011009__Kinase-like_dom__1__690__804	IPR:IPR015015__F-actin_binding__1__3009__381	IPR:IPR017441__Protein_kinase_ATP_BS__1__741__72	IPR:IPR020635__Tyr_kinase_cat_dom__1__723__756	CCDS:CCDS35165.1	CCDS:CCDS35166.1	
		
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
				my $FEATURE_ID=$FEATURE_tmp[3];
				unless($FEATURE_ID eq 'CHAIN')
				{
					$hash1{$ENSG}{$ENST}{$FEATURE_ID}{$distance_begin}{$distance_feature}=1;
				}
			}
			elsif ($tmp_tok=~/^IPR:(.+)/)
			{
				my $IPR_tok=$1;
				my @IPR_tmp=split("__",$IPR_tok);
				my $IPR_ID=join("**",$IPR_tmp[0],$IPR_tmp[2]);
				my $distance_begin=$IPR_tmp[3];
				my $distance_feature=$IPR_tmp[4];
				
				$hash1{$ENSG}{$ENST}{$IPR_ID}{$distance_begin}{$distance_feature}=1;
			}
			#elsif ($tmp_tok=~/^CCDS:/)
			#{
			#	my $CCDS=$tmp_tok;
			#	$hash1{$ENST}{$ENSG}{'CCDS'}{$CCDS}=1;
			#}
		}
	}
}else {print "Unable to open INPUT1\n";}

$time='['. timestamp(). ']'."\n";
print "START PRINTING:$time\n";

if(open(OUTPUT1, '>'.$output1))
{

	my @ENSG_tmp=sort keys %hash1;
	print "El scalar es:".scalar(@ENSG_tmp)."\n";
	for (my $i=0; $i<scalar(@ENSG_tmp); $i++)
	{
		
	foreach my $ENST_tok(sort keys %{$hash1{$ENSG_tmp[$i]}})
	{
			print OUTPUT1 "$ENSG_tmp[$i]\t$ENST_tok\t";
			
			foreach my $IPR_ID_tok(sort keys %{$hash1{$ENSG_tmp[$i]}{$ENST_tok}})
			{
			foreach my $distance_begin(sort keys %{$hash1{$ENSG_tmp[$i]}{$ENST_tok}{$IPR_ID_tok}})
			{
			foreach my $distance_feature(sort keys %{$hash1{$ENSG_tmp[$i]}{$ENST_tok}{$IPR_ID_tok}{$distance_begin}})
			{
				print OUTPUT1 "FEATURE:$IPR_ID_tok"."__"."$distance_begin"."__"."$distance_feature\t";
			}
			}
			}
			
			print OUTPUT1 "\n";
	}	
	}		
}else{print "Impossible to open OUTPUTS\n"; die;}

sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
