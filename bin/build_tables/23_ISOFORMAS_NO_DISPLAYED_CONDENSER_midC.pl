##	Script to condensed the features mapped to all the isoforms on a per line basis.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

use strict;
use warnings;
use Time::localtime;

my %hash1=();
my %hash2=();
my %hash3=();
my %hash4=();

my $input1=$ARGV[0];
my $output1=$ARGV[1];
my $output2=$ARGV[2];
my $output3=$ARGV[3];

# Here we extract all the genomic coordinates of all the transcipts in which protein or site features map on a per base basis

my $time='['. timestamp(). ']'."\n";
print "Start charging hash1 & PRINTING:$time\n";
my $counter2=0;

if (open (INPUT1, $input1) && open(OUTPUT, '>'.$output1))
{
	#INPUT1=ALL_ISOFORMS_PROTEIN_table.txt
	
#~ 12780885        1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR013094**1
#~ 12780885        1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR029058**1
#~ 12780886        1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR013094**1
#~ 12780886        1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR029058**1
#~ 12780887        1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR013094**1
#~ 12780887        1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR029058**1
#~ 12780888        1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR013094**1


	while(my $line=<INPUT1>)
	{
		chomp($line);
		#~ print "AAAAA:$line:HHHH\n";
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			$counter2++;
			#~ print "Hello_world:$line****************\n";
			my $POS=$1;
			my $CHROM=$2;
			my $strand=$3;
			my $SYMBOL=$4;
			my $ENSG=$5;
			my $ENST=$6;
			my $FEATURE=$7;
			my $string=join("\t",$CHROM,$strand,$SYMBOL,$ENSG,$ENST);
			if($counter2 ==1)
			{
				$hash1{$string}=1;
				#~ print "HELLO_WORLD_I:$counter2\t$string\t$FEATURE\t$POS\n"; 
			}
			foreach my $string_tok (sort keys %hash1)
			{
				if($string eq $string_tok)
				{
					#~ $hash2{$string_tok}{$FEATURE_tok}{$POS_tok}=1;
					$hash2{$string}{$FEATURE}{$POS}=1;
					#~ print "HELLO_WORLD_II:$string_tok\t$FEATURE_tok\t$POS_tok\n";
					#~ print "HELLO_WORLD_II:$counter2\t$string\t$FEATURE\t$POS\n";
				}
				elsif($string ne $string_tok)
				{
					foreach my $FEATURE_tok(sort keys %{$hash2{$string_tok}})
					{
						# Now we condense them using the condensing scheme developped for script 8_TRANSCRIPT_TABLE_condenser

							my @positions=();
							#~ print OUTPUT "la key es:$fields_tok\n";
							my @tmp_POS=sort{$a<=>$b}keys%{$hash2{$string_tok}{$FEATURE_tok}};
							#~ print OUTPUT "el array es:@tmp_POS\n";
							my $counter=0;
							my $max=scalar(@tmp_POS);
							#~ print OUTPUT "Máximo:$max\n";
							if(scalar(@tmp_POS) == 1)
							{
								print OUTPUT "\[$tmp_POS[0]\-$tmp_POS[0]\]\t$string_tok\t$FEATURE_tok\n";
							}
							elsif(scalar(@tmp_POS) > 1)
							{
								while(scalar(@tmp_POS) !=0)
								{
									my $contender=shift(@tmp_POS);
									push(@positions,$contender);
									$counter++;
									########
									#~ print OUTPUT "Counter:$counter\tEl array es:@positions\n";
									if($counter == 1)
									{
										print OUTPUT "\[$positions[$counter-1]\-";
									}
									elsif($counter != 1 && $counter != $max)
									{
										if($positions[$counter-1] == 1 + $positions[$counter-2])
										{
											# Do nothing
										}
										elsif($positions[$counter-1] != 1 + $positions[$counter-2])
										{
											print OUTPUT "$positions[$counter-2]\]\t$string_tok\t$FEATURE_tok\n";
											print OUTPUT "\[$positions[$counter-1]\-";
										}
									}
									elsif($counter == $max)
									{
										if($positions[$counter-1] == 1 + $positions[$counter-2])
										{
											print OUTPUT "$positions[$counter-1]\]\t$string_tok\t$FEATURE_tok\n";
										}
										elsif($positions[$counter-1] != 1 + $positions[$counter-2])
										{
											print OUTPUT "$positions[$counter-2]\]\t$string_tok\t$FEATURE_tok\n";
											print OUTPUT "\[$positions[$counter-1]\-$positions[$counter-1]\]\t$string_tok\t$FEATURE_tok\n";
										}
									}
								}
							}
					}
					#~ print "HELLO_WORLD_III:CHANGE:$counter2\t$string\t$FEATURE\t$POS\n";
					%hash1=();
					%hash2=();
					$hash1{$string}{$FEATURE}{$POS}=1;
					#~ print "HELLO_WORLD_I:$counter2\t$string\t$FEATURE\t$POS\n";
					$hash2{$string}{$FEATURE}{$POS}=1;
					#~ print "HELLO_WORLD_II:$counter2\t$string\t$FEATURE\t$POS\n";
				}
			}
		}
	}
}else {print "impossible to open INPUT1\n";die;}

# Here we open the previous file with genomic intervals and condense it further.

$time='['. timestamp(). ']'."\n";
print "Start charging & printing hash2:$time\n";

if (open(OUTPUT, $output1) && open(OUTPUT2, '>'.$output2))
{
		#OUTPUT=ALL_ISOFORMS_PROTEIN_table_midC.txt
	
	#INPUT1=ALL_ISOFORMS_PROTEIN_table.txt
	
#~ [12780885-12780948]     1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR013094**1
#~ [12785189-12785456]     1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR013094**1
#~ [12785667-12785879]     1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR013094**2
#~ [12780885-12780948]     1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR029058**1
#~ [12785189-12785618]     1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR029058**1
#~ [12785691-12785957]     1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR029058**2
#~ [12779651-12779693]     1       +       AADACL3 ENSG00000188984 ENST00000359318 IPR013094**1
#~ [12780885-12780948]     1       +       AADACL3 ENSG00000188984 ENST00000359318 IPR013094**1

	while(my $line=<OUTPUT>)
	{
		chomp($line);
		#~ print "AAAAA:$line:HHHH\n";
		if($line=~/([^\t]+)\t(.+)/)
		{
			#~ print "HELLO_WORLD_I:$line\n";
			my $interval=$1;
			my $string=$2;
			$hash3{$string}{$interval}=1;
		}
	}
	
	foreach my $string_tok(sort keys %hash3)
	{
		print OUTPUT2 "$string_tok\t";
	foreach my $interval_tok(sort keys %{$hash3{$string_tok}})
	{
		print OUTPUT2 "$interval_tok\t";	
	}
		print OUTPUT2 "\n";
	}
}
	
sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
