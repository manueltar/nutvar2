##	Script to create files with the SYMBOL, strand, ENST and exon coordinates on a per transcript basis of ENSEMBL. 2014.
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

# Here we open the file needed to obtain ENSG-SYMBOL-strand equivalence and parse the information

if (open (INPUT1, $input1))
{
	#~ INPUT1= gtf_output_ENSG.txt
	#~ 
#~ ENSG00000186092 OR4F5   +       69091   70008
#~ ENSG00000237683 AL627309.1      -       134901  139379

while(my $line=<INPUT1>)
	{
		chomp $line;
		#print "$line\n";
		if($line=~/^(ENSG[^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			my $ENSG=$1;
			my $SYMBOL=$2;
			my $strand=$3;
			my $begin=$4;
			my $end=$5;
			#print "hello_world:$SYMBOL\t$strand\t$ENSG\n";
			$hash1{$ENSG}{$SYMBOL}{$strand}{$begin}{$end}=1;
			#~ print "$ENSG\t$SYMBOL\t$strand\t$begin\t$end\n";
		}
	}
}else{print "Unable to open INPUT1\n";}

$time='['. timestamp(). ']'."\n";
print "Start charging hash2:$time\n";

# Here we obtain the protrin coding transcripts and their coordinates foreach ENSG

if (open(INPUT2, $input2))
{
#~ INPUT2=gtf_output_ENST.txt
#~ 
#~ 1       ENSG00000186092 ENST00000335137 protein_coding  69091   70008
#~ 1       ENSG00000237683 ENST00000423372 protein_coding  134901  139379
#~ 1       ENSG00000235249 ENST00000426406 protein_coding  367640  368634
#~ 1       ENSG00000185097 ENST00000332831 protein_coding  621059  622053
#~ 1       ENSG00000269831 ENST00000599533 protein_coding  738532  739137
#~ 1       ENSG00000269308 ENST00000594233 protein_coding  818043  819983
	
	while(my $line=<INPUT2>)
	{
		chomp($line);
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			my $CHROM=$1;
			my $ENSG=$2;
			my $ENST=$3;
			my $BIOTYPE=$4;
			my $begin=$5;
			my $end=$6;
			$hash2{$ENSG}{$ENST}{$begin}{$end}=1;
			#~ print "$ENSG\t$ENST\t$begin\t$end\n";
		
		}
	}
}else {print "impossible to open INPUT2\n";die;}

#~ exit;

$time='['. timestamp(). ']'."\n";
print "Start charging hash3:$time\n";

# Here we obtain the number of exons and exon coordinates for each transcripts

if (open(INPUT3, $input3))
{
			#~ INPUT3=gtf_output_EXON.txt
		#~ ENST00000335137 1       ENSE00002319515 69091   70008
		#~ ENST00000423372 1       ENSE00002221580 137621  139379
		#~ ENST00000423372 2       ENSE00002314092 134901  135802
	while(my $line=<INPUT3>)
	{
		chomp($line);
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			my $ENST=$1;
			my $EXON_number=$2;
			my $EXON=$3;
			my $begin=$4;
			my $end=$5;
			$hash3{$ENST}{$EXON_number}{$begin}{$end}=1;
			#~ print "$ENST\t$EXON_number\t$begin\t$end\n";
			$hash4{$ENST}{$begin}{$end}=1;
			#~ print "$ENST\t$begin\t$end\n";
		}
	}
}else {print "impossible to open INPUT3\n";die;}

$time='['. timestamp(). ']'."\n";
print "Start PRINTING:$time\n";

# Here we print all the information parsed before on a transcript per line basis

if (open(OUTPUT, '>'.$output))
{
	foreach my $ENSG_tok(sort keys%hash1)
	{
	foreach my $SYMBOL_tok(sort keys%{$hash1{$ENSG_tok}})
	{
	foreach my $strand_tok(sort keys%{$hash1{$ENSG_tok}{$SYMBOL_tok}})
	{
	foreach my $begin_tok(sort keys%{$hash1{$ENSG_tok}{$SYMBOL_tok}{$strand_tok}})
	{
	foreach my $end_tok(sort keys%{$hash1{$ENSG_tok}{$SYMBOL_tok}{$strand_tok}{$begin_tok}})
	{
		if($ENSG_tok eq 'ENSG00000175279' || $ENSG_tok eq 'ENSG00000188157' || $ENSG_tok eq 'ENSG00000196581'|| $ENSG_tok eq 'ENSG00000235098')
		{
			#~ print "$ENSG_tok\t$SYMBOL_tok\t$strand_tok\t$begin_tok\t$end_tok\n";
		}
		
		foreach my $ENST_tok(sort keys%{$hash2{$ENSG_tok}})
		{
			#~ print "**$ENST_tok**\n";
		foreach my $ENST_begin_tok(sort keys%{$hash2{$ENSG_tok}{$ENST_tok}})
		{
			#~ print "**$ENST_begin_tok**\n";
		foreach my $ENST_end_tok(sort keys%{$hash2{$ENSG_tok}{$ENST_tok}{$ENST_begin_tok}})
		{	
			
			#~ print "$SYMBOL_tok\t$strand_tok\t$begin_tok\t$end_tok\t$ENST_tok\t$ENST_begin_tok\t$ENST_end_tok\n";
			
			my @EXON_number_tmp=sort keys%{$hash3{$ENST_tok}};
			my $total_exon=scalar(@EXON_number_tmp);
			
			print OUTPUT "$SYMBOL_tok\t$strand_tok\t$begin_tok\t$end_tok\t$ENST_tok\t$ENST_begin_tok\t$ENST_end_tok\t$total_exon\t";
			
				
					my @EXON_begin_tmp=sort{$a<=>$b}keys%{$hash4{$ENST_tok}};
					#~ print "El array es:@EXON_begin_tmp\n";
					for(my $i=0;$i<scalar(@EXON_begin_tmp);$i++)
					{
						foreach my $EXON_end_tok(sort keys %{$hash4{$ENST_tok}{$EXON_begin_tmp[$i]}})
						{
							print OUTPUT "[$EXON_begin_tmp[$i]]"."-"."[$EXON_end_tok]\t";
						}
					}
			print OUTPUT "\n";
		}	
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
