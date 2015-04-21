##	Script to condensate the protein features on a per variant basis. 2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne) 

use strict;
use warnings;
use Time::localtime;

my %hash1=();
my %hash2=();



my $input1=$ARGV[0];
my $output1=$ARGV[1];

my $time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

my $counter=1;
my $Var="NaN";
my $Char="NaN";

if (open (INPUT1, $input1) && open (OUTPUT1, '>'.$output1))
{
	
#~ 1       +       ACTRT2  ENST00000378404 2938261 C       T       missense_variant        IPR004000**1    99.1158267020336        1
#~ 1       +       ACTRT2  ENST00000378404 2938265 C       T       synonymous_variant      IPR004000**1    98.762157382847 1
#~ 1       +       ACTRT2  ENST00000378404 2938266 G       A       missense_variant        IPR004000**1    98.6737400530504        1
#~ 1       +       ACTRT2  ENST00000378404 2938276 C       A       missense_variant        IPR004000**1    97.789566755084 1
#~ 1       +       ACTRT2  ENST00000378404 2938318 C       T       missense_variant        IPR004000**1    94.0760389036251        1
#~ 1       +       ACTRT2  ENST00000378404 2938319 G       A       synonymous_variant      IPR004000**1    93.9876215738285        1
#~ 1       +       ACTRT2  ENST00000378404 2938365 G       A       missense_variant        IPR004000**1    89.920424403183 1
#~ 1       +       ACTRT2  ENST00000378404 2938395 G       A       missense_variant        IPR004000**1    87.2679045092838        1
#~ 1       +       ACTRT2  ENST00000378404 2938424 G       T       synonymous_variant      IPR004000**1    84.7038019451813        1
#~ 1       +       ACTRT2  ENST00000378404 2938480 G       A       missense_variant        IPR004000**1    79.7524314765694        1
#~ 1       +       ACTRT2  ENST00000378404 2938505 C       T       synonymous_variant      IPR004000**1    77.5419982316534        1

while(my $line=<INPUT1>)
	{
		chomp $line;
		#~ print "INITIO\t$counter:$line\n";
		if($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/ && $counter == 1)
		{
			
			my $CHROM=$1;
			my $strand=$2;
			my $SYMBOL=$3;
			my $ENST=$4;
			my $POS=$5;
			my $REF=$6;
			my $ALT=$7;
			my $Effect=$8;
			
			$Var=join("\t",$CHROM,$strand,$SYMBOL,$ENST,$POS,$REF,$ALT,$Effect);
			
			my $Feature=$9;
			my $Percentage=$10;
			my $matched=$11;
			
			#~ print "hello_world:**$Var\t**$Feature\t**$Percentage\t**$matched**\n";
			
			$Char=join("|",$Feature,$Percentage,$matched);
			$hash1{$Var}{$Char}=1;
			#~ print "$Var\t$Char\t####$counter\n";
			$counter++;
		}
		if($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/ && $counter > 1)
		{
			my $CHROM=$1;
			my $strand=$2;
			my $SYMBOL=$3;
			my $ENST=$4;
			my $POS=$5;
			my $REF=$6;
			my $ALT=$7;
			my $Effect=$8;
			
			my $NewVar=join("\t",$CHROM,$strand,$SYMBOL,$ENST,$POS,$REF,$ALT,$Effect);
			
			my $Feature=$9;
			my $Percentage=$10;
			my $matched=$11;
			
			my $NewChar=join("|",$Feature,$Percentage,$matched);
			
			if($NewVar eq $Var)
			{
				$hash1{$NewVar}{$NewChar}=1;
				#~ print "$NewVar\t$NewChar\t####$counter\n";
				$counter++;
			}
			else
			{
				foreach my $Var_tok(sort keys %hash1)
				{
					my @Char_tmp=sort keys %{$hash1{$Var_tok}};
					my $string=join(";",@Char_tmp);
					print OUTPUT1 "$Var_tok\t$string\n";
				}
				
				%hash1=();
				$Var=$NewVar;
				$Char=$NewChar;
				$hash1{$Var}{$Char}=1;
				$counter=1;
				#~ print "$Var\t$Char\t####$counter\n";
			}
		}
	}
}else{print "Unable to open INPUT1\n";}
	
sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
