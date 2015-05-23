use strict;
use warnings;
use Time::localtime;

my %hash1=();
my %hash2=();


my $input1=$ARGV[0];
my $output1="tabla_AFs_Beta_distribution.txt";
my $input2="OUT_plus_Intervals_of_confidence.txt";
my $output2=$ARGV[1];
my $output3=$ARGV[2];

my $time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

if (open (INPUT1, $input1) && open (OUTPUT1, '>'.$output1))
{
	#INPUT1=JOIN_vcf.vcf

#~ #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	
	
#~ 9	99017189	.	G	A	PASS	100	ANTONIO_ISTVAN;WT=2431;HET=1;HOM=0;AF_COUNTS=0.000205592105263158
#~ 9	99017189	.	G	A	PASS	100	ClinVar;5;NaN
#~ 6	132206079	.	C	T	PASS	100	ANTONIO_ISTVAN;WT=10562;HET=648;HOM=7;AF_COUNTS=0.0295087813140768
#~ 6	132206079	.	C	T	PASS	100	ClinVar;0;NaN
#~ 1	69261	.	C	A	PASS	100	ANTONIO_ISTVAN;WT=250;HET=0;HOM=1;AF_COUNTS=0.00398406374501992
#~ 1	69270	.	A	G	PASS	100	ANTONIO_ISTVAN;WT=31;HET=19;HOM=261;AF_COUNTS=0.869774919614148
#~ 1	69428	.	T	G	PASS	100	ANTONIO_ISTVAN;WT=5609;HET=158;HOM=99;AF_COUNTS=0.0303443573133311
#~ 1	69476	.	T	C	PASS	100	ANTONIO_ISTVAN;WT=5464;HET=0;HOM=1;AF_COUNTS=0.000182982616651418


#~ ## INFO=<ID=CLNSIG,Number=.,Type=String,Description="Variant Clinical Significance, 
#~ 0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 
#~ 4- Likely pathogenic, 5 - Pathogenic, 
 #~ 6 - drug response, 7 - histocompatibility, 255 - other">

while(my $line=<INPUT1>)
	{
		chomp $line;
		#print "$line\n";
		unless($line=~/^#CHROM/)
		{
			if($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
			{
				my $CHROM=$1;
				my $POS=$2;
				my $ID=$3;
				my $REF=$4;
				my $ALT=$5;
				my $FILTER=$6;
				my $QUAL=$7;
				my $INFO=$8;
				my $threshold=0;
				my @counter=();
				if($INFO=~/^ClinVar/)
				{
					my @INFO_tmp=split(";",$INFO);
					my $CLNSIG=$INFO_tmp[1];
					my @CLNSIG_tmp=split(/\|/,$CLNSIG);
					foreach my $CLNSIG_tmp_tok(@CLNSIG_tmp)
					{
						if($CLNSIG_tmp_tok == 2 || $CLNSIG_tmp_tok == 3){$threshold=1;}
						
						if($CLNSIG_tmp_tok == 5)
						{
							push(@counter,$CLNSIG_tmp_tok);
						}
						else
						{
							# Do nothing
						}
					}
					unless($threshold == 1)
					{
						if($CLNSIG=~/^5/||$CLNSIG=~/^4/)	
						{
							my $counter=scalar(@counter);
							#~ print "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$FILTER\t$QUAL\t$INFO\t$counter\n";
							$hash1{$CHROM}{$POS}{$ID}{$REF}{$ALT}{$FILTER}{$QUAL}{'INFO'}{$INFO}=1;
							$hash1{$CHROM}{$POS}{$ID}{$REF}{$ALT}{$FILTER}{$QUAL}{'COUNTER'}{$counter}=1;
						}
					}
				}	
				elsif($INFO=~/^ANTONIO_ISTVAN/)
				{
					my ($wt,$HET,$HOM,$AF)="NaN"x4;
					my @INFO_tmp=split(";",$INFO);
					foreach my $INFO_tmp_tok(@INFO_tmp)
					{
						if($INFO_tmp_tok=~/^WT=(.+)/)
						{
							$wt=$1;
						}
						elsif($INFO_tmp_tok=~/^HET=(.+)/)
						{
							$HET=$1;
						}
						elsif($INFO_tmp_tok=~/^HOM=(.+)/)
						{
							$HOM=$1;
						}
						elsif($INFO_tmp_tok=~/^AF_COUNTS=(.+)/)
						{
							$AF=$1;
						}
					}
					print OUTPUT1 "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$FILTER\t$QUAL\t$wt\t$HET\t$HOM\t$AF\n";
				}
			}
		}
	}
}else{print "Unable to open INPUT1\n";}

#~ $time='['. timestamp(). ']'."\n";
#~ print "Start R-script:$time\n";
#~ 
#~ system ("R --vanilla < Beta_distribution.R");
#~ 
#~ $time='['. timestamp(). ']'."\n";
#~ print "Finish R-script:$time\n";
#~ 
#~ $time='['. timestamp(). ']'."\n";
#~ print "creating Hash1 second part:$time\n";

if(open(INPUT2, $input2) && open(OUTPUT2, '>'.$output2))
{
	#~ Input=OUT_plus_Intervals_of_confidence.txt

#~ "V1"    "V2"    "V3"    "V4"    "V5"    "V6"    "V7"    "V8"    "V9"    "V10"   "V11"   "lower" "upper"
#~ "1"     1       69261   "."     "C"     "A"     "PASS"  100     250     0       1       0.00398406374501992     0.00122677281158299
     #~ 0.0142324971785862
#~ "2"     1       69270   "."     "A"     "G"     "PASS"  100     31      19      261     0.869774919614148       0.436668118369433
       #~ 0.493917043435119
while(my $line=<INPUT2>)
	{
		chomp $line;
		#print "$line\n";
		unless($line=~/^"V1"/)
		{
			#~ print "$line\n";
			if($line=~/^\"[^\t]+\"\t([^\t]+)\t([^\t]+)\t\"([^\t]+)\"\t\"([^\t]+)\"\t\"([^\t]+)\"\t\"([^\t]+)\"\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
			{
				#~ print "**$line\n";
				my $CHROM=$1;
				my $POS=$2;
				my $ID=$3;
				my $REF=$4;
				my $ALT=$5;
				my $FILTER=$6;
				my $QUAL=$7;
				my ($wt,$HET,$HOM,$AF,$lower_confidence,$upper_confidence)="NaN"x6;
				$wt=$8;
				$HET=$9;
				$HOM=$10;
				$AF=$11;
				$lower_confidence=$12;
				$upper_confidence=$13;
				#~ print "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$FILTER\t$QUAL\t$wt\t$HET\t$HOM\t$AF\t$lower_confidence\t$upper_confidence\n";
				my $tag=0;
				#print "*****$tag\n";
				if($upper_confidence=~/\d/)
				{
					# rare; tag -1; highCredible <= 0.01
					
					if(0.01 >= $upper_confidence)
					{
						$tag=-1;
					}
					elsif (0.01 <= $lower_confidence)
					{
						$tag=1;
					}
					else
					{
						# Do nothing
					}
				}
				print OUTPUT2 "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$FILTER\t$QUAL\tANTONIO_ISTVAN;WT=$wt;HET=$HET;HOM=$HOM;AF_COUNTS=$AF;lowCI=$lower_confidence;highCI=$upper_confidence\t$tag\n";
				#print "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$FILTER\t$QUAL\tANTONIO_ISTVAN;WT=$wt;HET=$HET;HOM=$HOM;AF_COUNTS=$AF;lowCI=$lower_confidence;highCI=$upper_confidence\t$tag\n";
			}
			
		}
	}
}

$time='['. timestamp(). ']'."\n";
print "creating Hash1 third part:$time\n";

if(open(OUTPUT2, $output2))
{
		#~ Input=OUT_plus_Intervals_of_confidence_plus_TAGS.txt

#~ 1       69261   .       C       A       PASS    100     ANTONIO_ISTVAN;WT=250;HET=0;HOM=1;AF_COUNTS=0.00398406374501992;lowCI=0.00122677281158299;highCI=0.0142324971785862
     #~ 0
#~ 1       69270   .       A       G       PASS    100     ANTONIO_ISTVAN;WT=31;HET=19;HOM=261;AF_COUNTS=0.869774919614148;lowCI=0.436668118369433;highCI=0.493917043435119
        #~ 0
while(my $line=<OUTPUT2>)
	{
		chomp $line;
		# print"$line\n";
		if($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
				# print"**$line\n";
				my $CHROM=$1;
				my $POS=$2;
				my $ID=$3;
				my $REF=$4;
				my $ALT=$5;
				my $FILTER=$6;
				my $QUAL=$7;
				my $INFO=$8;
				my $Tag=$9;
				## print"$CHROM\t$POS\t$ID\t$REF\t$ALT\t$FILTER\t$QUAL\t$INFO\t$Tag\n";
				$hash1{$CHROM}{$POS}{$ID}{$REF}{$ALT}{$FILTER}{$QUAL}{'INFO'}{$INFO}=1;
				# print"***$CHROM\t$POS\t$ID\t$REF\t$ALT\t$FILTER\t$QUAL\t'INFO'\t$INFO\n";
				$hash1{$CHROM}{$POS}{$ID}{$REF}{$ALT}{$FILTER}{$QUAL}{'TAG'}{$Tag}=1;
				# print"****$CHROM\t$POS\t$ID\t$REF\t$ALT\t$FILTER\t$QUAL\t'TAG'\t$Tag\n";
		}
	}
}

$time='['. timestamp(). ']'."\n";
print "Start PRINTING:$time\n";

if (open(OUTPUT3, '>'.$output3))
{
	print OUTPUT3 "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	Pathogenicity_Tag	AF_Tag	Localization_Tag\n";
	
	foreach my $CHROM_tok(sort{$a<=>$b} keys %hash1)
	{
	foreach my $POS_tok(sort{$a<=>$b} keys%{$hash1{$CHROM_tok}})
	{
	foreach my $ID_tok(sort keys%{$hash1{$CHROM_tok}{$POS_tok}})
	{
	foreach my $REF_tok(sort keys%{$hash1{$CHROM_tok}{$POS_tok}{$ID_tok}})
	{
	foreach my $ALT_tok(sort keys%{$hash1{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}})
	{
	foreach my $FILTER_tok(sort keys%{$hash1{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}})
	{
	foreach my $QUAL_tok(sort keys%{$hash1{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$FILTER_tok}})
	{
		my @INFO_tmp=sort keys%{$hash1{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$FILTER_tok}{$QUAL_tok}{'INFO'}};
		
		my $string=join(";",@INFO_tmp);
		
		my $Flag_localization="NaN";
		
		if($string =~ /ANTONIO_ISTVAN/ && $string =~ /ClinVar/)
		{
			$Flag_localization=3;
		}
		elsif($string =~ /ANTONIO_ISTVAN/ && $string !~ /ClinVar/)
		{
			$Flag_localization=1;
		}
		elsif($string !~ /ANTONIO_ISTVAN/ && $string =~ /ClinVar/)
		{
			$Flag_localization=2;	
		}
		
		my $Pathogenicity_Tag="NaN";
		my $AF_Tag="NaN";
		
		my @counter_tmp=sort keys%{$hash1{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$FILTER_tok}{$QUAL_tok}{'COUNTER'}};
		my @Tag_tmp=sort keys%{$hash1{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$FILTER_tok}{$QUAL_tok}{'TAG'}};
		
		if(scalar(@counter_tmp) >=1 )
		{
			my $string2=join(";",@counter_tmp);
			$Pathogenicity_Tag=$string2;
		}
		if(scalar(@Tag_tmp) >=1 )
		{
			my $string3=join(";",@Tag_tmp);
			$AF_Tag=$string3;
		}
		
		#print "$CHROM_tok\t$POS_tok\t$ID_tok\t$REF_tok\t$ALT_tok\t$FILTER_tok\t$QUAL_tok\t$string\t$Pathogenicity_Tag\t$AF_Tag\t$Flag_localization\n";
		print OUTPUT3 "$CHROM_tok\t$POS_tok\t$ID_tok\t$REF_tok\t$ALT_tok\t$QUAL_tok\t$FILTER_tok\t$string\t$Pathogenicity_Tag\t$AF_Tag\t$Flag_localization\n";
		
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
