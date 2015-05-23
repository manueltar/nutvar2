###############
###############

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;


my $input=$ARGV[0];
my $input2=$ARGV[1];
my $output1=$ARGV[2];
my %hash1=();
my %hash2=();
my %hash3=();


my $time='['. timestamp(). ']'."\n";
print "Tiempo de incio de carga_hash1:$time\n";

if(open(INPUT,$input))# && open (OUTPUT2, '>'.$output2))
{
	# Input= HUMAN.fa (fasta de UNIPROT para human del 28/07/2014) o paralelizado FASTA_chunk_1.fa
	
	#>tr|A0A024R3B9|A0A024R3B9_HUMAN Crystallin, alpha B, isoform CRA_c OS=Homo sapiens GN=CRYAB PE=4 SV=1
	#MRLEKDRFSVNLDVKHFSPEELKVKVLGDVIEVHGKHEERQDEHGFISREFHRKYRIPAD
	#VDPLTITSSLSSDGVLTVNGPRKQVSGPERTIPITREEKPAVTAAPKK
	#>tr|A0A024R7E8|A0A024R7E8_HUMAN Elongation factor 1 homolog (S. cerevisiae), isoform CRA_a OS=Homo sapiens GN=ELOF1 PE=4 SV=1
	#MVRSRLTAVSASWVQAHPPADMGRRKSKRKPPPKKKMTGTLETQFTCPFCNHEKSCDVKM
	
	
	while(my $line=<INPUT>)
	{
		chomp $line;
		if ($line=~/^>/)
		{
			$line=~s/>//ig;
			## print "$line\n";
			my @tmp=split('\|',$line);
			my $AC=$tmp[1];
			my $description=$tmp[2];
			my $SYMBOL="NaN";
			
			my @description_tmp=split(/\s+/,$description);
			foreach my $description_tmp_tok(@description_tmp)
			{
				if($description_tmp_tok =~ /GN=(.+)/)
				{
					$SYMBOL=$1;
				}
			}
			$hash1{$SYMBOL}{$AC}=1;
		}
	
	}	
}else{print "Unable to open $input\n";}

$time='['. timestamp(). ']'."\n";
print "Tiempo de incio de carga_hash2:$time\n";

if(open(INPUT2,$input2))
{
	# Input=  gtf_output_ENSG_added_chrom.txt
	
	#~ 1       ENSG00000186092 OR4F5   +       69091   70008
	#~ 1       ENSG00000237683 AL627309.1      -       134901  139379
	#~ 1       ENSG00000235249 OR4F29  +       367640  368634
	#~ 1       ENSG00000185097 OR4F16  -       621059  622053
	#~ 1       ENSG00000269831 AL669831.1      -       738532  739137
	#~ 1       ENSG00000269308 AL645608.2      +       818043  819983
	#~ 1       ENSG00000187634 SAMD11  +       860260  879955
	#~ 1       ENSG00000268179 AL645608.1      -       861264  866445

	
	while(my $line=<INPUT2>)
	{
		chomp $line;
		#~ print "$line\n";
		my @tmp=split(/\t/,$line);
		#~ print "El array es:@tmp\n";
		my $ENSG=$tmp[1];
		#~ print "************$ENSG\n";
		my $SYMBOL=$tmp[2];
		
		if(exists($hash1{$SYMBOL}))
		{
			foreach my $AC_tok(sort keys %{$hash1{}})
			$hash2{$SYMBOL}
		}
	}	
}else{print "Unable to open $input2\n";}

if(open(OUTPUT,'>'.$output1))
{
foreach my $ENSG_tok(sort keys %hash1)
{
	unless (exists($hash2{$ENSG_tok}) || exists($hash3{$ENSG_tok}))
	{
		foreach my $SYMBOL_tok(sort keys %{$hash1{$ENSG_tok}})
		{
			print OUTPUT "$ENSG_tok\t$SYMBOL_tok\n";
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
