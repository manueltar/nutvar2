##	Script to build a compressed Pre-C with DNA sequence and genomic positions of all the CDS positions in the transcript in ENSEMBL.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

use strict;
use warnings;
use Time::localtime;

my $input1=$ARGV[0];
my $output=$ARGV[1];

my $CHROM="NaN";
my $ENSG="NaN";
my $ENST="NaN";
my $SYMBOL="NaN";
my $strand="NaN";
my $seq="NaN";
my %position_hash=();

my $time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

# Here we parse the CDS nucleotide sequence and the associated genomic position

if(open (INPUT1, $input1) && open (OUTPUT,'>'.$output))
{
	
	#~ Input:CDS_genomic_coordinates.txt
	#~ 
	#~ >X      ENSG00000000003 TSPAN6  ENST00000373020 -
	#~ 99883779__A     99883780__A     99883781__T     99883782__G     99883783__T     99883784__G     99883785__A     99883786__T     99883787__A     99883788__G     99883789__A     99883790__G     99883791__T
	#~ //
	#~ >X      ENSG00000000005 TNMD    ENST00000373031 +
	#~ 99840016__A     99840017__T     99840018__G     99840019__G     99840020__C     99840021__A     99840022__A     99840023__A     99840024__G     99840025__A     99840026__A     99840027__T     99840028__C
	#~ //
	#~ >X      ENSG00000001497 LAS1L   ENST00000312391 -
	#~ 64732504__A     64732505__G     64732506__T     64732507__A     64732508__G     64732509__G     64732510__T     64732511__C     64732512__T     64732513__T     64732514__T     64732515__T     64732516__C
	#~ //

	while(my $line=<INPUT1>)
	{
		chomp $line;
		#~ print "$line\n";
		
		# If we are in the last line unload hashes and arrays and re-initialize them and variables
		
		if ($line=~/^\/\//)
		{	print OUTPUT ">$CHROM\t$ENSG\t$SYMBOL\t$ENST\t$strand\n";
			#~ print  ">$CHROM\t$ENSG\t$SYMBOL\t$ENST\t$strand\n";
			my @seq_tmp=@{$position_hash{'BASE'}};
			#~ print "***********@seq_tmp***************\n";
			#~ print "REVERSE**********************".reverse(@seq_tmp)."*************\n";
			my $seq="NaN";
			# The sequence is strand dependent. It is the join of the array if the strand is +
			if($strand eq '+')
			{
				$seq=join("",@seq_tmp);
			}
			# However if the strand is - then the sequence is the reverse array
			elsif($strand eq '-')
			{
				#~ print "****************HELLLLO\n";
				my @seq_tmp_reverse=reverse(@seq_tmp);
				#~ print "El array es:@seq_tmp_reverse\n";
				$seq=join("",@seq_tmp_reverse);
			}

			print OUTPUT "seq>$seq\n";
			#~ print  "seq>$seq\n";
			my @positions=();
			my @tmp_POS=sort {$a<=>$b}keys %{$position_hash{'COORDINATE'}};
					my $counter=0;
		my $max=scalar(@tmp_POS);
		#~ print OUTPUT "Máximo:$max\n";
		
		# Once the sequence is done we compress the genomic coordinates in intervals following the scheme developped in the script 8
		print OUTPUT "POS>";
		#~ print  "POS>";
		if(scalar(@tmp_POS) == 1)
		{
			print OUTPUT "\[$tmp_POS[0]\-$tmp_POS[0]\]\t";
			#~ print  "\[$tmp_POS[0]\-$tmp_POS[0]\]\t";
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
					#~ print  "\[$positions[$counter-1]\-";
				}
				elsif($counter != 1 && $counter != $max)
				{
					if($positions[$counter-1] == 1 + $positions[$counter-2])
					{
						# Do nothing
					}
					elsif($positions[$counter-1] != 1 + $positions[$counter-2])
					{
						print OUTPUT "$positions[$counter-2]\]\t";
						#~ print  "$positions[$counter-2]\]\t";
						print OUTPUT "\[$positions[$counter-1]\-";
						#~ print  "\[$positions[$counter-1]\-";
					}
				}
				elsif($counter == $max)
				{
					if($positions[$counter-1] == 1 + $positions[$counter-2])
					{
						print OUTPUT "$positions[$counter-1]\]";
						#~ print  "$positions[$counter-1]\]";
					}
					elsif($positions[$counter-1] != 1 + $positions[$counter-2])
					{
						print OUTPUT "$positions[$counter-2]\]\t";
						#~ print  "$positions[$counter-2]\]\t";
						print OUTPUT "\[$positions[$counter-1]\-$positions[$counter-1]\]";
						#~ print  "\[$positions[$counter-1]\-$positions[$counter-1]\]";
					}
				}
			}
		}
		print OUTPUT "\n";
		#~ print  "**\n";
		print OUTPUT "//\n";
		#~ print  "//\n";
			
			
			%position_hash=();	
		}
		# Here we parse the IDs
		elsif ($line=~/^>([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			$CHROM=$1;
			$ENSG=$2;
			$SYMBOL=$3;
			$ENST=$4;
			$strand=$5;
			#~ print "$CHROM\t$ENSG\t$SYMBOL\t$ENST\t$strand\n";
		}
		#Here we create the positions hash and the sequence array of the CDS
		if($line !~ /^>/ && $line !~ /^\/\//)
		{
			if ($line=~/^(.+)\t$/)
			{
				my $POS=$1;
				#~ print "****************************$POS\n";
				my @POS_tmp=split(/\t/,$POS);
				#~ print "****************************".join("**",@POS_tmp)."\n";
				foreach my $POS_tmp_tok(@POS_tmp)
				{
					my @tmp=split("__",$POS_tmp_tok);
					#~ print "****************************".join("**",@tmp)."\n";
					my $coordinate=$tmp[0];
					my $base=$tmp[1];
					$position_hash{'COORDINATE'}{$coordinate}=1;
					push(@{$position_hash{'BASE'}},$base);
					#~ print "$coordinate\t$base\n";
				}
			}
		}
	}
}else{print "Unable to open INPUT2\n";}

$time='['. timestamp(). ']'."\n";
print "Print FIN:$time\n";
	
sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
