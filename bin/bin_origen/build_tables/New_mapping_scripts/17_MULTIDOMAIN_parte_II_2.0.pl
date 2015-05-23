##	Script to condense in protein intervals the positions covered by each numbered domain prediction.2014.
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

my $time='['. timestamp(). ']'."\n";
print "Start charging hash_1:$time\n";

# We open the file with all the domain and their domain family predictors on a per line basis.

if (open (INPUT, $input1))
{
	# Input= multidomain.txt
		
		#~ P24821  IPR000742       186     PS50026__1      
		#~ P24821  IPR000742       187     PS50026__1      
		#~ P24821  IPR000742       188     PS50026__1      
		#~ P24821  IPR000742       189     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       190     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       191     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       192     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       193     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       194     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       195     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       196     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       197     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       198     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       199     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       200     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       201     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       202     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       203     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       204     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       205     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       206     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       207     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       208     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       209     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       210     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       211     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       212     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       213     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       214     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       215     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       216     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       217     PS50026__1      SM00181__1      
		#~ P24821  IPR000742       220     SM00181__2      

	
	while(my $line=<INPUT>)
	{
		chomp $line;
		#~ print "$line\n";
		if ($line=~/(^[^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)\t$/)
		{
			#~ print "hello_world:$line**\n";
			my $AC= $1;
			my $IPR=$2;
			my $POS=$3;
			my $fields=$4;
			
			$hash_1{$AC}{$IPR}{$fields}{$POS}=1;
			#~ print "**$AC\t$IPR\t$fields\t$POS**\n";
		}
	}
}else{print "Unable to open INPUT\n";}


$time='['. timestamp(). ']'."\n";
print "Start printing:$time\n";

if(open(OUTPUT, '>'.$output1))
{
	foreach my $AC_tok (sort keys %hash_1)
	{
	foreach my $IPR_tok (sort keys %{$hash_1{$AC_tok}})
	{
	foreach my $fields_tok(sort keys %{$hash_1{$AC_tok}{$IPR_tok}})
	{
		#~ print "$AC_tok\t$IPR_tok\t$fields_tok\n";
		my @positions=();
		#~ print "la key es:$fields_tok\n";
		my @tmp_POS=sort{$a<=>$b}keys%{$hash_1{$AC_tok}{$IPR_tok}{$fields_tok}};
		#~ print "el array es:@tmp_POS\n";
		my $counter=0;
		my $max=scalar(@tmp_POS);
		#~ print "Máximo:$max\n";
		
		# If there is only one position for the key then the interval begins and ends with that position.
		
		if(scalar(@tmp_POS) == 1)
		{
			print OUTPUT "$tmp_POS[0]\t$tmp_POS[0]\t$AC_tok\t$IPR_tok\t$fields_tok\n";
		}
		
		# If there is more than one positionnin the array
		
		elsif(scalar(@tmp_POS) > 1)
		{
			# WE extract the first element of the array while the arrya has elements
			
			while(scalar(@tmp_POS) !=0)
			{
				my $contender=shift(@tmp_POS);
				push(@positions,$contender);
				$counter++;
				#~ print "Counter:$counter\tEl array es:@positions\n";
				
				# If it is the very first element of the array ($counter == 1) then we open the first interval of positions with the same key
				
				if($counter == 1)
				{
					print OUTPUT "$positions[$counter-1]\t";
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
						print OUTPUT "$positions[$counter-2]\t$AC_tok\t$IPR_tok\t$fields_tok\n";
						print OUTPUT "$positions[$counter-1]\t";
					}
				}
				
				# If we are in the last position of the array; if it is consecutive to the previous, then close the interval. If the last position is not consecutive
				# to the previous position, then close the previous interval and create a last interval with the last position in both limits of the interval.
				elsif($counter == $max)
				{
					if($positions[$counter-1] == 1 + $positions[$counter-2])
					{
						print OUTPUT "$positions[$counter-1]\t$AC_tok\t$IPR_tok\t$fields_tok\n";
					}
					elsif($positions[$counter-1] != 1 + $positions[$counter-2])
					{
						print OUTPUT "$positions[$counter-2]\t$AC_tok\t$IPR_tok\t$fields_tok\n";
						print OUTPUT "$positions[$counter-1]\t$positions[$counter-1]\t$AC_tok\t$IPR_tok\t$fields_tok\n";
					}
				}
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
