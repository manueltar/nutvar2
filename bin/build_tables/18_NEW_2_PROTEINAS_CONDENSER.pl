##	Script to condennse the protein features mapped to genomic positions in genomic intervals.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

#!/usr/bin/perl

# Da error al final del archivo input; ¿poner una línea más?

use strict;
use warnings;
use Time::localtime;

my %hash_1=();
my %hash_2=();

my $input1=$ARGV[0];
my $output1=$ARGV[1];
my $output2=$ARGV[2];

my $time='['. timestamp(). ']'."\n";
print "Start charging hash_1:$time\n";

# Here we parse the protein features on a genomic position per line basis

if (open (INPUT, $input1))
{
	# Input= ~/Escritorio/Proyecto_clasificador/Raw_Data/PROTEIN.txt
		
		#~ 12256540	19	-	ENST00000355738	ENSG00000257591	IPR007087**1	IPR013087**1	IPR015880**4	
		#~ 12256541	19	-	ENST00000355738	ENSG00000257591	IPR007087**1	IPR013087**1	IPR015880**4	
		#~ 12256542	19	-	ENST00000355738	ENSG00000257591	IPR007087**1	IPR013087**1	IPR015880**4	
		#~ 12256543	19	-	ENST00000355738	ENSG00000257591	IPR007087**1	IPR013087**1	IPR015880**4	
		#~ 12256544	19	-	ENST00000355738	ENSG00000257591	IPR007087**1	IPR013087**1	IPR015880**4	
		#~ 12256545	19	-	ENST00000355738	ENSG00000257591	IPR007087**1	IPR013087**1	IPR015880**4	
		#~ 12256546	19	-	ENST00000355738	ENSG00000257591	IPR007087**1	IPR013087**1	IPR015880**4	
		#~ 12256547	19	-	ENST00000355738	ENSG00000257591	IPR007087**1	IPR013087**1	IPR015880**4	
		#~ 12256548	19	-	ENST00000355738	ENSG00000257591	IPR007087**1	IPR013087**1	IPR015880**4	
		#~ 12256549	19	-	ENST00000355738	ENSG00000257591	IPR007087**1	IPR013087**1	IPR015880**4	
	
	while(my $line=<INPUT>)
	{
		chomp $line;
		if ($line=~/(^[0-9]+)\t(.+)\t$/)
		{
			#print "hello_world\n";
			my $POS= $1;
			my $fields=$2;
			$hash_1{$fields}{$POS}=1;
			#~ print "**$fields**\t$POS**\n";
		}
	}
}else{print "Unable to open INPUT\n";}

$time='['. timestamp(). ']'."\n";
print "Start printing:$time\n";

if(open(OUTPUT, '>'.$output1))
{
	foreach my $fields_tok(sort keys %hash_1)
	{
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
		
		# If there is more than one positionn in the array
		
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
					# If that element ($counter -1) is consecutive to the last (1 + positions[$counter-2]) do nothing
					
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


$time='['. timestamp(). ']'."\n";
print "Fin del script:$time\n";

sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
