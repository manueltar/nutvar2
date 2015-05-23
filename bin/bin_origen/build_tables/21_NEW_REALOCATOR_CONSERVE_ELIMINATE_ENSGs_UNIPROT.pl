###############
###############

#!/usr/bin/perl

use strict;
# use warnings;
use Time::localtime;
use Bio::EnsEMBL::Registry;



my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $input3=$ARGV[2];
my $input4=$ARGV[3];
my $output=$ARGV[4];

## Open input1 carrying protein information: equivalence ENST-AC_isoform, Domain_IDs, Domain_Coordinates, Postranslational sites and 'Active'sites and
## their coordinates.


# Print Time log

my $time='['. timestamp(). ']'."\n";
print "Start charging hash 1:$time\n";
my %hash1=();
my %hash2=();
my %hash3=();
my %hash4=();
my %hash5=();

if (open(INPUT1, $input1))
{
	my $ENSG="NaN";
	my $ENST="NaN";
	while(my $line=<INPUT1>)
	{
		chomp($line);
		
		#Input1 = UNIPROTseq_eq_ENSEMBL_seq.txt
		
		#	>ENSG00000000003        ENST00000373020
		#	MASPSRRLQTKPVITCFKSVLLIYTFIFWITGVILLAVGIWGKVSLENYFSLLNEKATNVPFVLIATGTVIILLGTFGCFATCRASAWMLKLYAMFLTLVFLVELVAAIVGFVFRHEIKNSFKNNYEKALKQYNSTGDYRSHAVDKIQNTLHCCGVTDYRDWTDTNYYSEKGFPKSCCKLEDCTPQRDADKVNNEGCFIKVMTIIESEMGVVAGISFGVACFQLIGIFLAYCLSRAITNNQYEIV
		#	MASPSRRLQTKPVITCFKSVLLIYTFIFWITGVILLAVGIWGKVSLENYFSLLNEKATNVPFVLIATGTVIILLGTFGCFATCRASAWMLKLYAMFLTLVFLVELVAAIVGFVFRHEIKNSFKNNYEKALKQYNSTGDYRSHAVDKIQNTLHCCGVTDYRDWTDTNYYSEKGFPKSCCKLEDCTPQRDADKVNNEGCFIKVMTIIESEMGVVAGISFGVACFQLIGIFLAYCLSRAITNNQYEIV
		#	>ENSG00000000005        ENST00000373031
		#	MAKNPPENCEDCHILNAEAFKSKKICKSLKICGLVFGILALTLIVLFWGSKHFWPEVPKKAYDMEHTFYSNGEKKKIYMEIDPVTRTEIFRSGNGTDETLEVHDFKNGYTGIYFVGLQKCFIKTQIKVIPEFSEPEEEIDENEEITTTFFEQSVIWVPAEKPIENRDFLKNSKILEICDNVTMYWINPTLISVSELQDFEEEGEDLHFPANEKKGIEQNEQWVVPQVKVEKTRHARQASEEELPINDYTENGIEFDPMLDERGYCCIYCRRGNRYCRRVCEPLLGYYPYPYCYQGGRVICRVIMPCNWWVARMLGRV
		#	MAKNPPENCEDCHILNAEAFKSKKICKSLKICGLVFGILALTLIVLFWGSKHFWPEVPKKAYDMEHTFYSNGEKKKIYMEIDPVTRTEIFRSGNGTDETLEVHDFKNGYTGIYFVGLQKCFIKTQIKVIPEFSEPEEEIDENEEITTTFFEQSVIWVPAEKPIENRDFLKNSKILEICDNVTMYWINPTLISVSELQDFEEEGEDLHFPANEKKGIEQNEQWVVPQVKVEKTRHARQASEEELPINDYTENGIEFDPMLDERGYCCIYCRRGNRYCRRVCEPLLGYYPYPYCYQGGRVICRVIMPCNWWVARMLGRV
		if($line=~/^>(ENSG[^\t]+)\t(ENST[^\t]+)/)
		{
			$ENSG=$1;
			$ENST=$2;
			# print "$ENSG\t$ENST\t$FLAG\n";	
			$hash1{$ENSG}{$ENST}=1;
		}
	}
}
$time='['. timestamp(). ']'."\n";
print "Tiempo de carga hash_2:$time\n";

if (open(INPUT2, $input2))
{
	my $ENSG="NaN";
	my $ENST="NaN";
	my $FLAG=0;
	
	while(my $line=<INPUT2>)
	{
		chomp($line);
		
		#Input2 = UNIPROTseq_ne_ENSEMBL_seq_contenidos.txt
		
		#	>ENSG00000084628        ENST00000263693
		#	MGKCSGRCTLVAFCCLQLVAALERQIFDFLGYQWAPILANFLHIMAVILGIFGTVQYRSRYLILYAAWLVLWVGWNAFIICFYLEVGQLSQDRDFIMTFNTSLHRSWWMENGPGCLVTPVLNSRLALEDHHVISVTGCLLDYPYIEALSSALQIFLALFGFVFACYVSKVFLEEEDSFDFIGGFDSYGYQAPQKTSHLQLQPLYTSG
		#	MAVILGIFGTVQYRSRYLILYAAWLVLWVGWNAFIICFYLEVGQLSQDRDFIMTFNTSLHRSWWMENGPGCLVTPVLNSRLALEDHHVISVTGCLLDYPYIEALSSALQIFLALFGFVFACYVSKVFLEEEDSFDFIGGFDSYGYQAPQKTSHLQLQPLYTSG

		if($line=~/^>(ENSG[^\t]+)\t(ENST[^\t]+)/)
		{
			$FLAG=1;
			$ENSG=$1;
			$ENST=$2;
			# print "$ENSG\t$ENST\t$FLAG\n";	
		}
		elsif($line=~/^\w+/ && $FLAG==1)
		{
			## print "LA LÍNEA ES: $line\n";
			my $UNIPROTseq=$line;
			# print "HELLO_WORLD_I:$UNIPROTseq\n";
			$hash2{$ENSG}{$ENST}{'UNIPROT'}{$UNIPROTseq}=1;
			$FLAG++;
		}
		elsif($line=~/^\w+/ && $FLAG==2)
		{
			## print "LA LÍNEA ES: $line\n";
			my $ENSEMBLseq=$line;
			# print "HELLO_WORLD_II:$ENSEMBLseq\n";
			$hash2{$ENSG}{$ENST}{'ENSEMBL'}{$ENSEMBLseq}=1;
			#exit;
			$FLAG=0;
		}
		
	}

$time='['. timestamp(). ']'."\n";
print "Tiempo de generación del offset:$time\n";

foreach my $ENSG_tok(sort keys%hash2)
{
foreach my $ENST_tok(sort keys%{$hash2{$ENSG_tok}})
{
	my @UNIPROTseq_non_overlap=();
	
	foreach my $UNIPROTseq_tok(sort keys%{$hash2{$ENSG_tok}{$ENST_tok}{'UNIPROT'}})
	{
	foreach my $ENSEMBLseq_tok(sort keys%{$hash2{$ENSG_tok}{$ENST_tok}{'ENSEMBL'}})
	{
		if($UNIPROTseq_tok=~/^(.+)$ENSEMBLseq_tok.*$/)
		{
			my $UNIPROTseq_non_overlap=$1;
			#print "$UNIPROTseq_non_overlap\n";
			@UNIPROTseq_non_overlap=split("",$UNIPROTseq_non_overlap);
			my $offset=scalar(@UNIPROTseq_non_overlap);
			#my $string=join("***",@UNIPROTseq_non_overlap);
			#print "$ENSG_tok\t$ENST_tok\t$string\n";
			$hash3{$ENSG_tok}{$ENST_tok}{$offset}=1;
		}
	}
	}		
	
}

}

}

$time='['. timestamp(). ']'."\n";
print "Tiempo de carga hash_4:$time\n";

if (open (INPUT3, $input3))
{
	#Input3=UNIPROTseq_ne_ENSEMBL_seq_no_contenidos.txt
	
	#	>ENSG00000005243        ENST00000006101
	#	MQRPEAWPRPHPGEGAAAAQAGGPAPPARAGEPSGLRLQEPSLYTIKAVFILDNDGRRLLAKYYDDTFPSMKEQMVFEKNVFNKTSRTESEIAFFGGMTIVYKNSIDLFLYVVGSSYENELMLMSVLTCLFESLNHMLRKNVEKRWLLENMDGAFLVLDEIVDGGVILESDPQQVIQKVNFRADDGGLTEQSVAQVLQSAKEQIKWSLLK
	#	MQRPEWPRPHPGEGAAAAGRGPAPPARAGEPSGLRLQEPSLYTIKAVFILDNDGRRLLAKYYDDTFPSMKEQMVFEKNVFNKTSRTESEIAFFGGMTIVYKNSIDLFLYVVGSSYENELMLMSVLTCLFESLNHMLRKNVEKRWLLENMDGAFLVLDEIVDGGVILESDPQQVIQKVNFRADDGGLTEQSVAQVLQSAKEQIKWSLLK

	while (my $line=<INPUT3>)
	{
		chomp $line;
		if($line=~/^>(ENSG[^\t]+)\t(ENST[^\t]+)/)
		{
			my $ENSG=$1;
			my $ENST=$2;
			# print "$ENSG\t$ENST\t$FLAG\n";	
			$hash4{$ENSG}{$ENST}=1;
		}
		
	}
}

$time='['. timestamp(). ']'."\n";
print "Tiempo de carga impresión líneas no afectadas y creación de hash_rereferenced:$time\n";

if(open (OUTPUT, '>'.$output))
{
	if(open (INPUT4, $input4))
	{
		while(my $line=<INPUT4>)
		{
			chomp($line);
			## Input4= UNIPROT_GLOBAL_POST_UNIQUE_PLUS_COORDINATES.txt

			#	A0A183  ENSG00000235942 ENST00000431011 ENSP00000411070 FEATURE:Late cornified envelope protein 6A.__0__240__CHAIN      CCDS:CCDS44227
			#	A0AUZ9  ENSG00000144445 ENST00000281772 ENSP00000281772 FEATURE:KAT8 regulatory NSL complex subunit 1-__0__2961__CHAIN  FEATURE:N6-ace
			#	A0AV02  ENSG00000221955 ENST00000393469 ENSP00000377112 FEATURE:Solute carrier family 12 member 8.__0__2142__CHAIN      IPR:IPR004841_
			
			my @tmp=split("\t",$line);
			my $AC=$tmp[0];
			my $ENSG=$tmp[1];
			my $ENST=$tmp[2];
			my $ENSP=$tmp[3];
			
			if(exists($hash1{$ENSG}{$ENST}))
			{
				print OUTPUT "$line\n";
			}
			elsif (exists($hash4{$ENSG}{$ENST}))
			{
				#do nothing eliminate ENSG with no sequence correspondence ENSEMBL UNIPROT
			}
			elsif (exists($hash3{$ENSG}{$ENST}))
			{
				foreach my $offset_tok(sort keys %{$hash3{$ENSG}{$ENST}})
				{
					my $offset_tok_base=$offset_tok*3;
					foreach my $tmp_tok(@tmp)
					{		
						if ($tmp_tok=~/^FEATURE:/)
						{
							my @FEATURE_tmp=split("__",$tmp_tok);
							my $FEATURE_ID=$FEATURE_tmp[0];
							my $distance_begin=$FEATURE_tmp[1];
							my $distance_feature=$FEATURE_tmp[2];
							my $FEATURE_name=$FEATURE_tmp[3];
							if($distance_begin < $offset_tok_base)
							{
								#do nothing, eliminate features in the offset zone
							}
							elsif($distance_begin => $offset_tok_base)
							{
								#Re-reference everything
								my $rereferenced_begin=$distance_begin-$offset_tok_base;
								my $string=join("__",$FEATURE_ID,$rereferenced_begin,$distance_feature,$FEATURE_name);
								$hash5{$AC}{$ENSG}{$ENST}{$ENSP}{$string}=1;
							}	
						}
						elsif ($tmp_tok=~/^IPR:/)
						{
							my @IPR_tmp=split("__",$tmp_tok);
							my $IPR_ID=$IPR_tmp[0];
							my $description=$IPR_tmp[1];
							my $counter=$IPR_tmp[2];
							my $distance_begin=$IPR_tmp[3];
							my $distance_IPR=$IPR_tmp[4];
							if($distance_begin < $offset_tok_base)
							{
								#do nothing, eliminate IPRs in the offset zone
							}
							elsif($distance_begin => $offset_tok_base)
							{
								#Re-reference everything
								my $rereferenced_begin=$distance_begin-$offset_tok_base;
								my $string=join("__",$IPR_ID,$description,$counter,$rereferenced_begin,$distance_IPR);
								$hash5{$AC}{$ENSG}{$ENST}{$ENSP}{$string}=1;
							}	
						}
						elsif ($tmp_tok=~/^CCDS:/)
						{
							$hash5{$AC}{$ENSG}{$ENST}{$ENSP}{$tmp_tok}=1;
						}
					}
				}
			}
		}
	}else {print "Unable to open INPUT or OUTPUT\n";}

	foreach my $AC_tok(sort keys %hash5)
	{
	foreach my $ENSG_tok(sort keys %{$hash5{$AC_tok}})
	{
	foreach my $ENST_tok(sort keys %{$hash5{$AC_tok}{$ENSG_tok}})
	{
	foreach my $ENSP_tok(sort keys %{$hash5{$AC_tok}{$ENSG_tok}{$ENST_tok}})
	{
		print OUTPUT "$AC_tok\t$ENSG_tok\t$ENST_tok\t$ENSP_tok\t";
		foreach my $item_tok(sort keys %{$hash5{$AC_tok}{$ENSG_tok}{$ENST_tok}{$ENSP_tok}})
		{
			print OUTPUT "$item_tok\t";
		}
		print OUTPUT "\n";
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
