###############
###############

#!/usr/bin/perl

use strict;
# use warnings;
use Time::localtime;
use Bio::EnsEMBL::Registry;



my $input1=$ARGV[0];
my $output1=$ARGV[1];
my $output2=$ARGV[2];

## Open input1 carrying protein information: equivalence ENST-AC_isoform, Domain_IDs, Domain_Coordinates, Postranslational sites and 'Active'sites and
## their coordinates.


# Print Time log

my $time='['. timestamp(). ']'."\n";
print "Start charging hash 1:$time\n";
my %hash1=();
my %hash2=();
my $FLAG=0;


if(open (INPUT1, $input1))
{
	my $ENSG="NaN";
	my $ENST="NaN";
	
	while(my $line=<INPUT1>)
	{
		chomp($line);
		## Input1= UNIPROTseq_ne_ENSEMBL_seq.txt

		#>ENSG00000005243        ENST00000006101
#UNIPROT#MQRPEAWPRPHPGEGAAAAQAGGPAPPARAGEPSGLRLQEPSLYTIKAVFILDNDGRRLLAKYYDDTFPSMKEQMVFEKNVFNKTSRTESEIAFFGGMTIVYKNSIDLFLYVVGSSYENELMLMSVLTCLFESL
#ENSEMBL#MQRPEWPRPHPGEGAAAAGRGPAPPARAGEPSGLRLQEPSLYTIKAVFILDNDGRRLLAKYYDDTFPSMKEQMVFEKNVFNKTSRTESEIAFFGGMTIVYKNSIDLFLYVVGSSYENELMLMSVLTCLFESLNH
		#>ENSG00000007545        ENST00000262317
#UNIPROT#MTVKLGDGGSGEDGLKKLGKRAADEESLEGEGAGGADAAEESSGTKRDEKTPRAGADGPPAPPGAPQAPSPPQGSPQDQHHFLRSSVRPQSKRPRKDPPSAVGSGNAGGSGPRGKGAEGGGSSSGNVSGVAPAA
#ENSEMBL#MLRTDAGEDLADAPAEELQEKGSPAGPPPSQGQPAARPPKEVPASRLAQQLREEGWNLQTSESLTLAEVYLMMGKPSKLQLEYDWLGPGRQDPRPGSLPTALHKQRLLSCLLKLISTEVNPKLALEANTISTAS
		#>ENSG00000007866        ENST00000338863
#UNIPROT#MASNSWNASSSPGEAREDGPEGLDKGLDNDAEGVWSPDIEQSFQEALAIYPPCGRRKIILSDEGKMYGRNELIARYIKLRTGKTRTRKQVSSHIQVLARKKVREYQVGIKAMNLDQVSKDKALQSMASMSSAQI
#ENSEMBL#IASNSWNASSSPGEAREDGPEGLDKGLDNDAEGVWSPDIEQSFQEALAIYPPCGRRKIILSDEGKMYGRNELIARYIKLRTGKTRTRKQVSSHIQVLARKKVREYQVGIKAMNLDQVSKDKALQSMASMSSAQI
		
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
			$hash1{$ENSG}{$ENST}{'UNIPROT'}{$UNIPROTseq}=1;
			$FLAG++;
		}
		elsif($line=~/^\w+/ && $FLAG==2)
		{
			## print "LA LÍNEA ES: $line\n";
			my $ENSEMBLseq=$line;
			# print "HELLO_WORLD_II:$ENSEMBLseq\n";
			$hash1{$ENSG}{$ENST}{'ENSEMBL'}{$ENSEMBLseq}=1;
			#exit;
			$FLAG=0;
		}
	}

}else {print "Unable to open INPUT1\n";}

$time='['. timestamp(). ']'."\n";
print "Tiempo de impresión:$time\n";

if (open(OUTPUT,'>'.$output1) && open(OUTPUT2, '>'.$output2))
{
	foreach my $ENSG_tok(sort keys %hash1)
	{
	foreach my $ENST_tok(sort keys %{$hash1{$ENSG_tok}})
	{
		foreach my $UNIPROTseq_tok(sort keys %{$hash1{$ENSG_tok}{$ENST_tok}{'UNIPROT'}})
		{
			#print "AAAAA:$UNIPROTseq_tok\n";	
		foreach my $ENSEMBLseq_tok(sort keys %{$hash1{$ENSG_tok}{$ENST_tok}{'ENSEMBL'}})
		{
			#print "HHH:$ENSEMBLseq_tok\n";	
			if($UNIPROTseq_tok=~/$ENSEMBLseq_tok/)
			{
				print OUTPUT  ">$ENSG_tok\t$ENST_tok\n";
				print OUTPUT  "$UNIPROTseq_tok\n";
				print OUTPUT  "$ENSEMBLseq_tok\n";
			}
			else
			{
				print OUTPUT2  ">$ENSG_tok\t$ENST_tok\n";
				print OUTPUT2  "$UNIPROTseq_tok\n";
				print OUTPUT2  "$ENSEMBLseq_tok\n";
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
