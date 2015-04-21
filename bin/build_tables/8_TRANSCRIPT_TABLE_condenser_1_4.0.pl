##	Script to condensate in genomic intervals the file with the trasncripts mapping to every genomic coordinate. 2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)


#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;

my %CHROM_accepted=();
my %hash_1=();
my %hash_2=();

my $input1=$ARGV[0];
my $output1=$ARGV[1];
my $output2=$ARGV[2];

my $time='['. timestamp(). ']'."\n";
print "Start charging hash_1:$time\n";

# We open the file with every genomic coordinate in a separate line with all the transcripts mapping to it

if (open (INPUT, $input1))
{
	#~ Input= TRANSCRIPTS_table.txt
	
	#~ 33772828        1       A3GALT2 -       ENST00000330379 ENST00000442999 
	#~ 33772829        1       A3GALT2 -       ENST00000330379 ENST00000442999 
	#~ 33772830        1       A3GALT2 -       ENST00000330379 ENST00000442999 
	#~ 33772831        1       A3GALT2 -       ENST00000330379 ENST00000442999 
	#~ 33772832        1       A3GALT2 -       ENST00000330379 ENST00000442999 
	#~ 33772833        1       A3GALT2 -       ENST00000330379 ENST00000442999 
	#~ 33772834        1       A3GALT2 -       ENST00000330379 ENST00000442999 
	#~ 33772835        1       A3GALT2 -       ENST00000330379 ENST00000442999 
	
	while(my $line=<INPUT>)
	{
		chomp $line;
		if ($line=~/(^[0-9]+)\s+(.+)/)
		{
			#print "hello_world\n";
			my $POS= $1;
			my $fields=$2;
			
			# We invert the order; we keep the common field of transcripts as a master key and as a dependent key the positions to which that combination of transcripts map
			# It should be noted that we can do it this way because genomic coordinates are unique; they are not repeated, Should there be rpetitions we should use a hash 
			# of array positions instead of a hash of hashes.
			
			$hash_1{$fields}{$POS}=1
		}
	}
}else{print "Unable to open INPUT\n";}

$time='['. timestamp(). ']'."\n";
print "Start printing:$time\n";

if(open(OUTPUT, '>'.$output1))
{
	foreach my $fields_tok(sort keys %hash_1)
	{
		# For every field key we obtain all the coding positions in an array, we initialize a counter variable an define a max parameter equal to the total number of positions 
		# in the array ($max)
		
		my @positions=();
		#~ print OUTPUT "la key es:$fields_tok\n";
		my @tmp_POS=sort{$a<=>$b}keys%{$hash_1{$fields_tok}};
		#~ print OUTPUT "el array es:@tmp_POS\n";
		my $counter=0;
		my $max=scalar(@tmp_POS);
		#~ print OUTPUT "Máximo:$max\n";
		
		# If there is only one position for the key then the interval begins and ends with that position.
		
		if(scalar(@tmp_POS) == 1)
		{
			print OUTPUT "\[$tmp_POS[0]\-$tmp_POS[0]\]\t$fields_tok\n";
		}
		
		# If there is more than one positionn in thea array
		
		elsif(scalar(@tmp_POS) > 1)
		{
			# WE extract the first element of the array while the array has elements
			
			while(scalar(@tmp_POS) !=0)
			{
				my $contender=shift(@tmp_POS);
				push(@positions,$contender);
				$counter++;
				#~ print OUTPUT "Counter:$counter\tEl array es:@positions\n";
				
				# If it is the very first element of the array ($counter == 1) then we open the first interval of positions with the same key
				
				if($counter == 1)
				{
					print OUTPUT "\[$positions[$counter-1]\-";
				}
				
				# If every other element unless the last
				
				elsif($counter != 1 && $counter != $max)
				{
					# If that elemente ($counter -1) is consecutive to the last (1 + positions[$counter-2]) do nothing
					
					if($positions[$counter-1] == 1 + $positions[$counter-2])
					{
						# Do nothing
					}
					
					# If that element is not consecutive; then close previous interval an open a new one
					
					elsif($positions[$counter-1] != 1 + $positions[$counter-2])
					{
						print OUTPUT "$positions[$counter-2]\]\t$fields_tok\n";
						print OUTPUT "\[$positions[$counter-1]\-";
					}
				}
				
				# If we are in the last position of the array; if it is consecutive to the previous, then close the interval. If the last position is not consecutive
				# to the previous position, then close the previous interval and create a last interval with the last position in both limits of the interval.
				elsif($counter == $max)
				{
					if($positions[$counter-1] == 1 + $positions[$counter-2])
					{
						print OUTPUT "$positions[$counter-1]\]\t$fields_tok\n";
					}
					elsif($positions[$counter-1] != 1 + $positions[$counter-2])
					{
						print OUTPUT "$positions[$counter-2]\]\t$fields_tok\n";
						print OUTPUT "\[$positions[$counter-1]\-$positions[$counter-1]\]\t$fields_tok\n";
					}
				}
			}
		}
	}
}else{print "Impossible to open OUTPUT\n";}

# Now we create a more condensed form of the intervals adding all the intervals of the same key in a per line basis.

if (open(OUTPUT, $output1) && open(OUTPUT2, '>'.$output2))
{
		#~ [180257501-180257652]   1       ACBD6   -       ENST00000367595 
		#~ [180283827-180283857]   1       ACBD6   -       ENST00000367595 
		#~ [180366651-180366740]   1       ACBD6   -       ENST00000367595 
		#~ [180382501-180382606]   1       ACBD6   -       ENST00000367595 
		#~ [180399315-180399397]   1       ACBD6   -       ENST00000367595 
		#~ [180461404-180461500]   1       ACBD6   -       ENST00000367595 
		#~ [180464596-180464660]   1       ACBD6   -       ENST00000367595 
		#~ [180471180-180471401]   1       ACBD6   -       ENST00000367595 
		
	while(my $line=<OUTPUT>)
	{
		chomp $line;
		if($line=~/^([^\t]+)\t(.+)\t$/)
		{

			my $interval=$1;
			my $key=$2;
			$hash_2{$key}{$interval}=1;
		}
		
	}
	foreach my $key_tok(sort keys %hash_2)
	{
		print OUTPUT2 "$key_tok\t";
	foreach my $interval_tok(sort keys %{$hash_2{$key_tok}})
	{
	
		print OUTPUT2 "$interval_tok\t";
	}
		print OUTPUT2 "\n";
	}
}else{print "Impossible to open INPUT2\n";}
$time='['. timestamp(). ']'."\n";
print "Fin del script:$time\n";
sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
