##	Script to put together the Percentage of sequence damged, the protein and the site features. 2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

use strict;
use warnings;
use Time::localtime;

my %hash1=();
my %primer_hash=();


my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $input3=$ARGV[2];
my $input4=$ARGV[3];
my $input5=$ARGV[4];
my $output=$ARGV[5];



my $time='['. timestamp(). ']'."\n";
print "Start charging hash3:$time\n";

if (open(INPUT3, $input3))
{
	while(my $line=<INPUT3>)
	{
		chomp($line);
		
		## Input file= 1GK_PROTEINS_GLOBAL.txt 
		
		#  X       +       ACRC    ENST00000373696 70832274        A       G       missense_variant        IPR006640**1|44.327731092437|1
		#  X       +       ACRC    ENST00000373696 70832724        C       T       synonymous_variant      IPR006640**1|13.2352941176471|1
		#  X       +       ACRC    ENST00000373696 70832741        C       T       missense_variant        IPR006640**1|9.66386554621849|1
		#  X       +       ACRC    ENST00000373696 70832746        C       T       missense_variant        IPR006640**1|8.61344537815126|1
		#  X       +       AFF2    ENST00000286437 147800737       G       A       synonymous_variant      IPR007797**1|100|0;MOD_RES**1|100|0;MO
		#  X       +       AFF2    ENST00000286437 147800748       A       G       missense_variant        IPR007797**1|100|0;MOD_RES**1|100|0;MO
		#  X       +       AFF2    ENST00000286437 147800752       G       C       missense_variant        IPR007797**1|100|0;MOD_RES**1|100|0;MO
		
		#~ print "$line\n";
		if($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)/)
		{
			#~ print "Hello_world$line\n";
			my $CHROM=$1;
			my $strand=$2;
			my $SYMBOL=$3;
			my $ENST=$4;
			my $POS=$5;
			my $REF=$6;
			my $ALT=$7;
			my $Effect=$8;
			my $fields=$9;
			#~ if(exists($primer_hash{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}))
			#~ {
				my @tmp=split(";",$fields);
				#print "hello_world:$CHROM\t$strand\t$SYMBOL\t$ENST\t$POS\t$REF\t$ALT\t$Effect\t$Percentage_affected\t$fields\n";
				foreach my $tmp_tok(@tmp)
				{
					$hash1{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{'TRANSCRIPT'}{$ENST}{'DETAIL_DOMAIN'}{$tmp_tok}=1;
					#~ print"$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t$Effect\t'TRANSCRIPT'\t$ENST\t'DETAIL_DOMAIN'\t$tmp_tok\n";
				}
			#~ }
		}
	}
}else {print "impossible to open INPUT3\n";die;}


$time='['. timestamp(). ']'."\n";
print "Start charging hash4:$time\n";

if (open(INPUT4, $input4))
{
	while(my $line=<INPUT4>)
	{
		chomp($line);
		
		## Input file= 1GK_DOMAINS.txt 
		
		#  X       +       ABCD1   ENST00000218104 152990759       A       C       missense_variant        DOMAIN|98.1141692150866|1;SITE|100|0
		#  X       +       ABCD1   ENST00000218104 152990979       C       T       synonymous_variant      DOMAIN|86.901121304791|1;SITE|84.63541
		#  X       +       ABCD1   ENST00000218104 152991192       A       G       synonymous_variant      DOMAIN|76.0448521916412|1;SITE|29.1666
		#  X       +       ABCD1   ENST00000218104 152991322       G       A       missense_variant        DOMAIN|69.4189602446483|1;SITE|6.51041
		#  X       +       ABCD1   ENST00000218104 152991417       G       T       synonymous_variant      DOMAIN|64.5769622833843|1;SITE|6.51041
		#  X       +       ABCD1   ENST00000218104 152991477       C       A       missense_variant        DOMAIN|61.5188583078491|1;SITE|6.51041
		#  X       +       ABCD1   ENST00000218104 152991478       C       G       missense_variant        DOMAIN|61.4678899082569|1;SITE|6.51041
		#  X       +       ABCD1   ENST00000218104 152991616       C       T       missense_variant        DOMAIN|54.434250764526|1;SITE|6.510416
		#  X       +       ABCD1   ENST00000218104 152994731       G       A       synonymous_variant      DOMAIN|51.8858307849134|1;SITE|6.51041
		#  X       +       ABCD1   ENST00000218104 152994833       C       A       synonymous_variant      DOMAIN|46.6870540265036|1;SITE|6.51041
		#  X       +       ABCD1   ENST00000218104 152994857       C       T       synonymous_variant      DOMAIN|45.4638124362895|1;SITE|6.51041
		#  X       +       ABCD1   ENST00000218104 153001567       T       A       missense_variant&splice_region_variant  DOMAIN|44.852191641182
		#  X       +       ABCD1   ENST00000218104 153001699       G       A       synonymous_variant      DOMAIN|38.1243628950051|1;SITE|6.51041

		#print "$line\n";
		if($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)/)
		{
			my $CHROM=$1;
			my $strand=$2;
			my $SYMBOL=$3;
			my $ENST=$4;
			my $POS=$5;
			my $REF=$6;
			my $ALT=$7;
			my $Effect=$8;
			my $fields=$9;
			my @tmp=split(";",$fields);
			#~ if(exists($primer_hash{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}))
			#~ {
				#print "hello_world:$CHROM\t$strand\t$SYMBOL\t$ENST\t$POS\t$REF\t$ALT\t$Effect\t$Percentage_affected\t$fields\n";
				foreach my $tmp_tok(@tmp)
				{
					$hash1{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{'TRANSCRIPT'}{$ENST}{'GLOBAL_DOMAIN'}{$tmp_tok}=1;
					#print"$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t$Effect\t'TRANSCRIPT'\t$ENST\t'GLOBAL_DOMAIN'\t$tmp_tok\n";
				}
			#~ }
		}
	}
}else {print "impossible to open INPUT4\n";die;}

$time='['. timestamp(). ']'."\n";
print "Start charging hash5:$time\n";

if (open(INPUT5, $input5))
{
	while(my $line=<INPUT5>)
	{
		chomp($line);
		
		## Input file= 1GK_sequence_percentage_def.txt
		
		#	1       A3GALT2 33772799        G       T       stop_gained     ENST00000330379 50.2923976608187
		#	1       A3GALT2 33772799        G       T       stop_gained     ENST00000442999 42.156862745098
		#	1       A3GALT2 33772876        C       G       missense_variant        ENST00000330379 59.2982456140351
		#	1       A3GALT2 33772876        C       G       missense_variant        ENST00000442999 49.7058823529412
		#	1       A3GALT2 33773016        G       A       missense_variant        ENST00000330379 75.672514619883
		#	1       A3GALT2 33773016        G       A       missense_variant        ENST00000442999 63.4313725490196
		#	1       A3GALT2 33777669        TAGTC   T       frameshift_variant&feature_truncation   ENST00000330379 82.1052631578947
		#	1       A3GALT2 33777669        TAGTC   T       frameshift_variant&feature_truncation   ENST00000442999 68.8235294117647


		#print "$line\n";
		if($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			my $CHROM=$1;
			my $SYMBOL=$2;
			my $POS=$3;
			my $REF=$4;
			my $ALT=$5;
			my $Effect=$6;
			my $ENST=$7;
			my $Percentage=$8;
			#~ if(exists($primer_hash{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}))
			#~ {
				#print "hello_world:$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t$Effect\t$ENST\t$Percentage\n";
				$hash1{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{'TRANSCRIPT'}{$ENST}{'PERCENTAGE_AFFECTED'}{$Percentage}=1;
					#print"$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t$Effect\t'TRANSCRIPT'\t$ENST\t'GLOBAL_DOMAIN'\t$tmp_tok\n";
			#~ }
		}
	}
}else {print "impossible to open INPUT5\n";die;}

$time='['. timestamp(). ']'."\n";
print "Start printing:$time\n";

if (open(OUTPUT, '>'.$output))
{
	foreach my $CHROM_tok(sort keys %hash1)
	{
	foreach my $SYMBOL_tok(sort keys %{$hash1{$CHROM_tok}})
	{
	foreach my $POS_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}})
	{
	foreach my $REF_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}})
	{
	foreach my $ALT_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}})
	{
	foreach my $Effect_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}})
	{
		
		my @tmp_impresion=("NaN","NaN","NaN");
		
		my $string_impresion=join(";",@tmp_impresion);
		# print "$string_impresion\t";
		# print OUTPUT "$string_impresion\t";
		
		foreach my $ENST_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'TRANSCRIPT'}})
		{
			my @tmp_impresion2=();
			my @tmp_impresion3=();
			
			#print "$CHROM_tok\t$SYMBOL_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$Effect_tok\t$string_impresion\t$ENST_tok\t";
			print OUTPUT "$CHROM_tok\t$SYMBOL_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$Effect_tok\t$string_impresion\t$ENST_tok\t";
			
			my $derived_NMD="NaN";
			my $NMD_tok="NaN";
			
			foreach my $Percentage_affected_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'TRANSCRIPT'}{$ENST_tok}{'PERCENTAGE_AFFECTED'}})
				{
					push(@tmp_impresion3, $Percentage_affected_tok);
				}
			
			#~ foreach my $NMD_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'TRANSCRIPT'}{$ENST_tok}{'NMD'}})
				#~ {
					push(@tmp_impresion2, $NMD_tok);
				#~ }
			foreach my $domains_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'TRANSCRIPT'}{$ENST_tok}{'DETAIL_DOMAIN'}})
				{
					#print "DETAIL:$domains_tok\n";
					push(@tmp_impresion2, $domains_tok);
				}
			foreach my $global_domains_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'TRANSCRIPT'}{$ENST_tok}{'GLOBAL_DOMAIN'}})
				{
					#print "GLOBAL:$global_domains_tok\n";
					push(@tmp_impresion3, $global_domains_tok);
				}
			#~ unless (exists($hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'TRANSCRIPT'}{$ENST_tok}{'DERIVED_PTC'}))
				#~ {
					push(@tmp_impresion2, $derived_NMD);
				#~ }
			
			#~ else
			#~ {
				#~ foreach my $PTC_derived_NMD_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'TRANSCRIPT'}{$ENST_tok}{'DERIVED_PTC'}})
				#~ {
					#~ push(@tmp_impresion2, $PTC_derived_NMD_tok);
				#~ }
			#~ }
			
			my $string_impresion2=join(";",@tmp_impresion2);
			my $string_impresion3=join(";",@tmp_impresion3);
			#print "$string_impresion2\n";
			print OUTPUT "$string_impresion3\t$string_impresion2\n";
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
