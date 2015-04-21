##	Script to create the super-domains collapsed from overlapping boundaries of different domain predictors.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

#!/usr/bin/perl


use strict;
use warnings;
use Time::localtime;

my %hash_1=();
my %hash_2=();

my $input1=$ARGV[0];
my $output1=$ARGV[1];

my $line_counter=0;

my $time='['. timestamp(). ']'."\n";
print "Start charging hash_1:$time\n";

# We open the file that has the protein intervals with the domain predictions numbered and mapped

if (open (INPUT, $input1))
{
	# Input= multidomain_midC_sorted.txt
		
		#~ 1       987     A0AUZ9  IPR026180       PTHR22443__1
		#~ 795     915     A0AUZ9  IPR029332       PF15275__1
		#~ 44      411     A0AV02  IPR004841       PF00324__1
		#~ 71      71      A0AV96  IPR000504       PS50102__1
		#~ 72      72      A0AV96  IPR000504       PS50102__1      SM00360__1
		#~ 73      137     A0AV96  IPR000504       PF00076__1      PS50102__1      SM00360__1
	
	while(my $line=<INPUT>)
	{
		chomp $line;
		#~ print "$line\n";
		if ($line=~/(^[^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)/)
		{
			#~ print "hello_world:$line**\n";
			my $begin=$1;
			my $end=$2;
			my $AC= $3;
			my $IPR=$4;
			my $fields=$5;
			$hash_1{$AC}{$IPR}{$line_counter}{'BEGIN'}{$begin}=1;
			$hash_1{$AC}{$IPR}{$line_counter}{'END'}{$end}=1;
			$hash_1{$AC}{$IPR}{$line_counter}{'FAMILY'}{$fields}=1;
			#~ print "**$AC\t$IPR\t$line_counter\t$fields\t$begin\t$end**\n";
			$line_counter++;
		}
	}
}else{print "Unable to open INPUT\n";}


$time='['. timestamp(). ']'."\n";
print "Start printing:$time\n";

# Here we process the previous hash; we merge domains with overlapping boundaries for which a member of prediction family in the previuous is foun
# in the posterior

if(open(OUTPUT, '>'.$output1))
{
	foreach my $AC_tok (sort keys %hash_1)
	{
		my %hash_merge_definitive=();
	
	foreach my $IPR_tok (sort keys %{$hash_1{$AC_tok}})
	{
		
		my %hash_merge=();
		
		my @tmp_line_counter=sort{$a<=>$b}keys%{$hash_1{$AC_tok}{$IPR_tok}};
		for(my $i=0;$i<scalar(@tmp_line_counter);$i++)
		{
			my @begin=sort{$a<=>$b}keys%{$hash_1{$AC_tok}{$IPR_tok}{$tmp_line_counter[$i]}{'BEGIN'}};
			
			# This is just to treat begin as an scalar
			
			my $string_begin=join("",@begin);
			#~ print "contender:$IPR_tok\t$tmp_line_counter[$i]\t$string_begin\n";
						
			my @end=sort{$a<=>$b}keys%{$hash_1{$AC_tok}{$IPR_tok}{$tmp_line_counter[$i]}{'END'}};
			
			# This is just to treat end as an scalar
			
			my $string_end=join("",@end);
			#~ print "contender:$tmp_line_counter[$i]\t$string_end\n";
			my @family=sort keys%{$hash_1{$AC_tok}{$IPR_tok}{$tmp_line_counter[$i]}{'FAMILY'}};
			
			# Family is the predictor services (numbered) overlaping in that line of coordinates
			
			my $string_family=join("",@family);
			#~ print "contender:$tmp_line_counter[$i]\t$string_family\n";
			my $Flag1=0;
			my $Flag2=0;
			
			# With the start line we start charging the array
			
			if($i == 0)
			{
				#~ push(@{$hash_merge{'BEGIN'}},$string_begin);
				#~ push(@{$hash_merge{'END'}},$string_end);
				#~ push(@{$hash_merge{'FAMILY'}},$string_family);
				
				$hash_merge{'BEGIN'}{$tmp_line_counter[$i]}{$string_begin}=1;
				$hash_merge{'END'}{$tmp_line_counter[$i]}{$string_end}=1;
				$hash_merge{'FAMILY'}{$tmp_line_counter[$i]}{$string_family}=1;
				
				$hash_merge_definitive{$IPR_tok}{$string_begin}{$string_end}=1;
			}
			
			# As we advance in the lines of the file
			
			elsif($i > 0)
			{
				my @previous_begin=keys%{$hash_merge{'BEGIN'}{$tmp_line_counter[$i-1]}};
				my $string_previous_begin=join("",@previous_begin);
				my @previous_end=keys%{$hash_merge{'END'}{$tmp_line_counter[$i-1]}};
				my $string_previous_end=join("",@previous_end);
				#~ print "$string_previous_end\n";
				my @previous_family=keys%{$hash_merge{'FAMILY'}{$tmp_line_counter[$i-1]}};
				my $string_previous_family=join("",@previous_family);
				#~ print "$string_previous_family\n";
				
				# First condition, positions are consecutive
				
				if($string_begin - $string_previous_end <= 1)
				{
					$Flag1=1;
				}
				
				# Second condition, the same family is found across overlapping intervals
				my @tmp_family_splitted=split(/\t/,$string_family);
				my @tmp_family_previous_splitted=split(/\t/,$string_previous_family);
				foreach my $family_element_tok(@tmp_family_splitted)
				{
				foreach my $family_previous_element_tok(@tmp_family_previous_splitted)	
				{
					if($family_element_tok eq $family_previous_element_tok)
					{
						$Flag2=1;
					}
				}
				}
				
				# If both conditions are fulfilled then the new boundary for the IPR will be from previous begin to string end and the
				# family of predictors will be the last one
				
				if($Flag1 == 1 && $Flag2 == 1)
				{
					#~ print "FLAG1:$Flag1\tFLAG2:$Flag2\n";
					
					
					$hash_merge{'BEGIN'}{$tmp_line_counter[$i]}{$string_previous_begin}=1;
					$hash_merge{'END'}{$tmp_line_counter[$i]}{$string_end}=1;
					$hash_merge{'FAMILY'}{$tmp_line_counter[$i]}{$string_family}=1;
					
					#~ print"'BEGIN'\t$tmp_line_counter[$i]\t$string_previous_begin\n";
					#~ print"'END'\t$tmp_line_counter[$i]\t$string_end\n";
					#~ print"'FAMILY'\t$tmp_line_counter[$i]\t$string_family\n";
					#~ 
					$hash_merge_definitive{$IPR_tok}{$string_previous_begin}{$string_end}=1;
				}
				# If both conditions are not fulfilled then the hash maintains the begin, end and family of the new line
				
				else
				{
					#~ print "FLAG1:$Flag1\tFLAG2:$Flag2\n";
					
					
					$hash_merge{'BEGIN'}{$tmp_line_counter[$i]}{$string_begin}=1;
					$hash_merge{'END'}{$tmp_line_counter[$i]}{$string_end}=1;
					$hash_merge{'FAMILY'}{$tmp_line_counter[$i]}{$string_family}=1;
					
					#~ print"'BEGIN'\t$tmp_line_counter[$i]\t$string_begin\n";
					#~ print"'END'\t$tmp_line_counter[$i]\t$string_end\n";
					#~ print"'FAMILY'\t$tmp_line_counter[$i]\t$string_family\n";
					
					$hash_merge_definitive{$IPR_tok}{$string_begin}{$string_end}=1;
				}
			}
		}
	}#
	
	# Now we undo the hash; numbering domains depending on the number of begins remaining after the merging we ha just done
	
		foreach my $IPR_tok(sort keys %hash_merge_definitive)
		{
			my @begin_IPR_definitive=sort{$a<=>$b} keys%{$hash_merge_definitive{$IPR_tok}};
			
			for(my $j=0;$j<scalar(@begin_IPR_definitive);$j++)
			{
				my $domain_number=$j+1;
				
				my @end_IPR_definitive=sort{$a<=>$b}keys%{$hash_merge_definitive{$IPR_tok}{$begin_IPR_definitive[$j]}};
				
				my $end_max=$end_IPR_definitive[scalar(@end_IPR_definitive) -1];
				
				print OUTPUT "$AC_tok\t$IPR_tok"."__"."$domain_number\t$begin_IPR_definitive[$j]\t$end_max\n";
				
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
