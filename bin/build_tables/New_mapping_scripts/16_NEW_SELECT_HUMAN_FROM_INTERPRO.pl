##	Script to select from the INTERPRO file only those domains belonging to human proteins.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;


my $input=$ARGV[0];
my $input2=$ARGV[1];
my $output1=$ARGV[2];

my %initial_hash=();
my %ID_hash=();
my $time="NaN";

$time='['. timestamp(). ']'."\n";
	print "Tiempo de incio de carga_hash1:$time\n";

# We open the file with all the human identifiers of UniProt

if(open(INPUT,$input))#&& open (OUTPUT2, '>'.$output2))
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
			my $db=$tmp[0];
			## print "$AC\t$db\t$description\n";
			$initial_hash{$AC}{$db}{$description}=1;
			#print OUTPUT2 "$AC\n";
		}
	
	}	
}else{print "Unable to open $input\n";}




$time='['. timestamp(). ']'."\n";
print "Carga del hash2:$time\n";

my %hash_IPR_coordinates=();
my $counter_INTERPRO=0;

# We open the file with hole INTERPRO (27GB)

if(open(INPUT2,$input2) && open(OUTPUT,'>'.$output1))
{	
	#~ Input=protein2ipr.dat
#~ 
	#~ P24821  IPR000742       Epidermal growth factor-like domain     PS50026 186     217
	#~ P24821  IPR000742       Epidermal growth factor-like domain     PS50026 280     311
	#~ P24821  IPR000742       Epidermal growth factor-like domain     PS50026 373     404
	#~ P24821  IPR000742       Epidermal growth factor-like domain     PS50026 466     497
	#~ P24821  IPR000742       Epidermal growth factor-like domain     PS50026 559     590
	#~ P24821  IPR000742       Epidermal growth factor-like domain     SM00181 189     217
	#~ P24821  IPR000742       Epidermal growth factor-like domain     SM00181 220     248
	#~ P24821  IPR000742       Epidermal growth factor-like domain     SM00181 251     280
	#~ P24821  IPR000742       Epidermal growth factor-like domain     SM00181 283     311
	#~ P24821  IPR000742       Epidermal growth factor-like domain     SM00181 376     404
	#~ P24821  IPR000742       Epidermal growth factor-like domain     SM00181 407     435
	#~ P24821  IPR000742       Epidermal growth factor-like domain     SM00181 438     466
	#~ P24821  IPR000742       Epidermal growth factor-like domain     SM00181 469     497
	#~ P24821  IPR000742       Epidermal growth factor-like domain     SM00181 500     528
	#~ P24821  IPR000742       Epidermal growth factor-like domain     SM00181 531     559
	#~ P24821  IPR000742       Epidermal growth factor-like domain     SM00181 562     590
	#~ P24821  IPR000742       Epidermal growth factor-like domain     SM00181 593     621
	#~ P24821  IPR002181       Fibrinogen, alpha/beta/gamma chain, C-terminal globular domain  PF00147 1981    2188
	#~ P24821  IPR002181       Fibrinogen, alpha/beta/gamma chain, C-terminal globular domain  PS51406 1975    2190
	#~ P24821  IPR002181       Fibrinogen, alpha/beta/gamma chain, C-terminal globular domain  SM00186 1979    2189
	#~ P24821  IPR002181       Fibrinogen, alpha/beta/gamma chain, C-terminal globular domain  SSF56496        1971    2194

while(my $line=<INPUT2>)
	{
		#print "Hello_world1\n";
		chomp $line;
		#~ print "**$line**\n";
		if ($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			#~ print "Hello_world:$line**\n";
			my $AC_INTERPRO=$1;
			my $IPR_INTERPRO=$2;
			my $description=$3;
			my $family=$4;
			my $begin=$5;
			my $end=$6;
			
			# We restrict the file of INTERPRO to domains with a human UniProt identifier

			if(exists($initial_hash{$AC_INTERPRO}))
			{
				print OUTPUT "$AC_INTERPRO\t$IPR_INTERPRO\t$family\t$begin\t$end\n";
				print "**********************$AC_INTERPRO\t$IPR_INTERPRO\t$family\t$begin\t$end\n";
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
