##	Script to convert the diferent IPR positions to DOMAIN label and the diferent site positions to SITE label.
##  The final aim is to create a Pre-C table of intervals of DOMAIN or SITE positions .2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

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

# Here we parse the file with the genomic intervals of all the domains (IPR) and functional sites

if(open(INPUT1, $input1))
{
	
	while(my $line=<INPUT1>)
	{
		#~ Input= ALL_ISOFORMS_PROTEIN_table_full.txt
		
		#~ 1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR013094**1    [12780885-12780948]     [12785189-12785456]     
		#~ 1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR013094**2    [12785667-12785879]     
		#~ 1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR029058**1    [12780885-12780948]     [12785189-12785618]     
		#~ 1       +       AADACL3 ENSG00000188984 ENST00000332530 IPR029058**2    [12785691-12785957]     
		#~ 1       +       AADACL3 ENSG00000188984 ENST00000359318 IPR013094**1    [12779651-12779693]     [12780885-12780948]     [12785189-12785456]
			 

		chomp($line);
		#~ print "HelloI:$line**\n";
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)\t$/)
		{
			#~ print "HelloII:$line**\n";
			my $CHROM=$1;
			my $strand=$2;
			my $SYMBOL=$3;
			my $ENSG=$4;
			my $ENST=$5;
			my $FEATURE=$6;
			my $fields=$7;
			my $DOMAIN="NaN";
			
			# Here we transform IPR positions to DOMAIN and functional site positions to SITE
			
			if($FEATURE=~/^IPR\d+/){$DOMAIN="DOMAIN";} #print "$FEATURE\t$DOMAIN\n";}
			
			elsif($FEATURE !~/^IPR\d+/){$DOMAIN="SITE";} #print "$FEATURE\t$DOMAIN\n";}
			
			my $string=join("\t",$CHROM,$strand,$SYMBOL,$ENSG,$ENST,$DOMAIN);
			
			my @coordinates_tmp=split(/\t/,$fields);
			foreach my $coordinates_tmp_tok(@coordinates_tmp)
			{
				if($coordinates_tmp_tok=~/\[(.+)\-(.+)\]/)
				{
					my $begin=$1;
					my $end=$2;
					for(my $i=$begin;$i<=$end;$i++)
					{
						$hash_1{$string}{$i}=1;
					}
				}
			}
		}
	}
}

# Here we condense the previous hash without needing to build an intermediate file and print the final condensed file.

$time='['. timestamp(). ']'."\n";
print "Start printing_1:$time\n";

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
		
		# Here we use the same condensing scheme that we used in script 8_TRANSCRIPT_TABLE_condenser
		
		if(scalar(@tmp_POS) == 1)
		{
			print OUTPUT "\[$tmp_POS[0]\-$tmp_POS[0]\]\t$fields_tok\n";
		}
		elsif(scalar(@tmp_POS) > 1)
		{
			while(scalar(@tmp_POS) !=0)
			{
				my $contender=shift(@tmp_POS);
				push(@positions,$contender);
				$counter++;
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
						print OUTPUT "$positions[$counter-2]\]\t$fields_tok\n";
						print OUTPUT "\[$positions[$counter-1]\-";
					}
				}
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
print "Start printing_2:$time\n";

# Here we condense further the genomic intervals using the same approach that in script 23_ISOFORMAS_NO_DISPLAYED_CONDENSER

%hash_1=();

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
		#~ print "$line**\n";
		if($line=~/^([^\t]+)\t(.+)/)
		{
			#~ print "Hello_world:$line**\n";
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



sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
