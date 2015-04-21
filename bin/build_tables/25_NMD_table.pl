##	Script to build a Pre-C table of regions susceptible to Nonsense-Mediated-Decay.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

use strict;
use warnings;
use Time::localtime;

my %hash1=();
my %hash2=();
my %hash3=();
my %hash4=();

my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $input3=$ARGV[2];
my $output=$ARGV[3];

my $time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

# Here we obtain the equivalence SYMBOL-ENSG-strand

if (open (INPUT1, $input1))
{
	#INPUT1=gtf_output_ENSG.txt
	
	# ENSG00000000003 TSPAN6  -       99883667        99894988
	# ENSG00000000005 TNMD    +       99839799        99854882
	# ENSG00000000419 DPM1    -       49551404        49575092
	# ENSG00000000457 SCYL3   -       169818772       169863408
	# ENSG00000000460 C1orf112        +       169631245       169823221
	# ENSG00000000938 FGR     -       27938575        27961788
	# ENSG00000000971 CFH     +       196621008       196716634
	# ENSG00000001036 FUCA2   -       143815948       143832827
	# ENSG00000001084 GCLC    -       53362139        53481768
while(my $line=<INPUT1>)
	{
		chomp $line;
		#~ print "I:$line\n";
		if($line=~/^(ENSG[^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			#~ print "II:$line**\n";
			my $ENSG=$1;
			my $SYMBOL=$2;
			my $strand=$3;
			my $START=$4;
			my $END=$5;
			$hash1{$SYMBOL}{$strand}{$START}{$END}{$ENSG}=1;
			#~ print "$SYMBOL\t$strand\t$START\t$END\t$ENSG\n";
		}
	}
}else{print "Unable to open INPUT1\n";}

# Here we obtain all the protein coding trasncripts for a given ENSG 

$time='['. timestamp(). ']'."\n";
print "Start charging hash2:$time\n";

if (open(INPUT2, $input2))
{
	## Input file= gtf_output_ENST.txt

	#~ 1       ENSG00000186092 ENST00000335137 protein_coding  69091   70008
	#~ 1       ENSG00000237683 ENST00000423372 protein_coding  134901  139379
	#~ 1       ENSG00000235249 ENST00000426406 protein_coding  367640  368634
	#~ 1       ENSG00000185097 ENST00000332831 protein_coding  621059  622053


	while(my $line=<INPUT2>)
	{
		chomp($line);
		#~ print "I:$line\n";
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			#~ print "II:$line**\n";
			my $CHROM=$1;
			my $ENSG=$2;
			my $ENST=$3;
			my $BIOTYPE=$4;
			if($BIOTYPE eq 'protein_coding')
			{
				$hash2{$ENSG}{$ENST}{$CHROM}=1;
				#~ print "$ENSG\t$ENST\t$CHROM\n";
			}
		}
	}
}else {print "impossible to open INPUT2\n";die;}

# Here we obtain all the exon exon boundaries

$time='['. timestamp(). ']'."\n";
print "Start charging hash3:$time\n";

if (open(INPUT3, $input3))
{
	while(my $line=<INPUT3>)
	{
		chomp($line);
		# Input=gtf_ENSG_ENST_EXON.txt 
		# TSPAN6  -       99883667        99894988        ENST00000373020 99883667        99891803        8       [99891605]-[99891803]   [99890
		# TSPAN6  -       99883667        99894988        ENST00000494424 99888439        99894988        6       [99894942]-[99894988]   [99891
		# TSPAN6  -       99883667        99894988        ENST00000496771 99887538        99891686        6       [99891188]-[99891686]   [99890
		
		#~ print "I:$line\n";
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t[^\t]+\t[^\t]+\t([^\t]+)\t(.+)\t$/)
		{
			#~ print "II:$line**\n";
			my $SYMBOL=$1;
			my $strand=$2;
			my $START=$3;
			my $END=$4;
			my $ENST=$5;
			my $exon_number=$6;
			my $fields=$7;
			#print "$SYMBOL\t$strand\t$START\t$END\t$ENST\t$exon_number\t$fields\n";
			foreach my $ENSG_tok(sort keys%{$hash1{$SYMBOL}{$strand}{$START}{$END}})
			{
				foreach my $CHROM_tok(sort keys %{$hash2{$ENSG_tok}{$ENST}})
				{
					#my $string=join("\t",$CHROM_tok,$strand,$SYMBOL,$ENSG_tok);
					$hash3{$CHROM_tok}{$strand}{$SYMBOL}{$ENSG_tok}{$ENST}{$exon_number}{$fields}=1;
					#~ print "*********$CHROM_tok\t$strand\t$SYMBOL\t$ENSG_tok\t$ENST\t$exon_number\t$fields*****\n";
				}
			}
		}
	}
}else {print "impossible to open INPUT3\n";die;}

# Now we print all the exon-exon boundaries

$time='['. timestamp(). ']'."\n";
print "Start printing:$time\n";

if (open(OUTPUT, '>'.$output))
{
	foreach my $CHROM_tok(sort {$a<=>$b}keys%hash3)
	{
	foreach my $strand_tok(sort keys %{$hash3{$CHROM_tok}})
	{
	foreach my $SYMBOL_tok(sort keys %{$hash3{$CHROM_tok}{$strand_tok}})
	{
	foreach my $ENSG_tok(sort keys %{$hash3{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}})
	{
	foreach my $ENST_tok(sort keys %{$hash3{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENSG_tok}})
	{
	foreach my $exon_number_tok(sort{$a<=>$b} keys %{$hash3{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENSG_tok}{$ENST_tok}})
	{
	foreach my $fields_tok(sort keys %{$hash3{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENSG_tok}{$ENST_tok}{$exon_number_tok}})
	{
		print OUTPUT "$CHROM_tok\t$strand_tok\t$SYMBOL_tok\t$ENSG_tok\t$ENST_tok\t$exon_number_tok\t$fields_tok\n";
	}			
	}		
	}
	}		
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
