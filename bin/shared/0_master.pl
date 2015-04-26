
# Manuel Tardaguila, PhD. University of Florida.
# The aim of this script is creating a hash with all the short reads for a certain given sample

use strict;
use warnings;
use Time::localtime;

my %hash1=();
my %hash2=();

my $time='['. timestamp(). ']'."\n";
print "Loading_hash1:$time\n";

print "Type the name of the vcf file you want to analyse (p.e. test.vcf)\n";
my $name=<STDIN>;

my $input="./input/$name";
my $output1="./data/intermediate/user_vcf_eff.vcf";

my $counter=0;

# Here I can do a checking of the format of the user file

if (open(INPUT1, $input))
{
	while(my $line=<INPUT1>)
	{
		chomp($line);
		$counter++;
		if($counter<=10)
		{
			print "Hello_world:$counter\t$line\n";
		}
	}
}



system ("perl name.pl");
exit;
system("java -Xmx4g -jar SOFTWARE/snpEff/snpEff.jar eff -c SOFTWARE/snpEff/snpEff.config -v GRCh37.75 -lof -csvStats -nextProt -sequenceOntology $input > here.txt");

exit;



## Input file= /home/manueltar/Desktop/Proyect_PacBio_LHMA/Raw_data/SR_correction/prueba.fq

#~ @NS500169:5:H0JV9AGXX:1:11101:23399:1130 1:N:0:1
#~ TGCCAGAGATTGCCGTGTACCCTGCCTTTGAAACACCTACACAGTATGTTTTGCCA
#~ +
#~ <A)A<FFF//F/FFFF/FFF<FA/FFA)A7FFFFFAFA)FFFF)7/7A)FA)FFF<
#~ @NS500169:5:H0JV9AGXX:1:11101:8262:1130 1:N:0:1
#~ GCAGTGAATGATGAGGGGTTGGGCAGTGGATGAAGGCAGAGGGGGGCAGGACAGAA
#~ +
#~ <//7AFF<FFFFFFFFFFFFFFFF<FFFFFAFFAFFF<FFFFFFFFFAFFAFFF7F


# We set a  counter variable

my $counter=0;

my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $output=$ARGV[2];

if (open(INPUT1, $input1))
{
	while(my $line=<INPUT1>)
	{
		chomp($line);
	#print "hello_world_1\n";
		if ($line=~/^\@NS5/)
		{
			#~ print "$line\n";
			$counter=1;
		}
		else{$counter++;}
		if($counter==2)
		{
			#~ print "Hello_world:$line\n";
			$hash1{$line}=1;
		}
	}
}

#~ exit;

# Reseting the variable to 0

$counter=0;

my $ID=0;
my $SR=0;
my $strand=0;
my $quality=0;

$time='['. timestamp(). ']'."\n";
print "Matching_hash1_against_hash2:$time\n";

if (open(INPUT2, $input2) && open(OUTPUT, '>'.$output))
{
	## Input file= /home/manueltar/Desktop/Proyect_PacBio_LHMA/Raw_data/SR_correction/NSC1_S1_L002_R1_001_trimmed3.fq
	
	#~ @NS500169:5:H0JV9AGXX:2:11101:16567:1037 1:N:0:1
	#~ GAGTGAAGCTCTTCCTGGGGACAATGTGGGCTTCAATGTAAAG
	#~ +
	#~ A7AA<FFFFFAFFFFFFAFFFFFFAFFFFFFFFFFAF<FFAFF
	while(my $line=<INPUT2>)
		{
			chomp($line);
		#print "hello_world_1\n";
			if ($line=~/^\@NS5/)
			{
				#~ print "$line\n";
				$counter=1;
				$ID=$line;
			}
			else{$counter++;}
			if($counter==2)
			{
				$SR=$line;
			}
			elsif($counter==3)
			{
				$strand=$line;
			}
			elsif($counter==4)
			{
				$quality=$line;
				unless(exists($hash1{$SR}))
				{
					print OUTPUT "$ID\n";
					print OUTPUT "$SR\n";
					print OUTPUT "$strand\n";
					print OUTPUT "$quality\n";
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
