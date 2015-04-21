use strict;
use warnings;
use Time::localtime;

my %hash1=();
my %hash2=();


my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $output=$ARGV[2];

my $time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

if (open (INPUT1, $input1))
{
	#INPUT1=OUTPUTVCF_CollapsedPerGene_Variants_PlusCountsFromIstvan.txt
	
#~ 12:6483977      A       G       7494    101     0
#~ 16:20335452     C       A       6502    1       0
#~ 14:74205652     T       C       10844   91      1
#~ 18:44424727     T       C       7588    2       0
#~ 2:191375186     TC      T       3567    1       0
#~ 10:3199163      T       C       1091    1       0


while(my $line=<INPUT1>)
	{
		chomp $line;
		#print "$line\n";
		unless($line=~/^base/)
		{
			if($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
			{
				my @INFO=();
				my $CHROM_POS=$1;
				my $REF=$2;
				my $ALT=$3;
				my $wild_type=$4;
				my $Heterozygous=$5;
				my $Homozygous=$6;
				#print "hello_world:$CHROM_POS\t$REF_ALT\t$SYMBOL\n";
				my @CHROM_POS_tmp=split(":",$CHROM_POS);
				my $CHROM=$CHROM_POS_tmp[0];
				my $POS=$CHROM_POS_tmp[1];
				my $ID=".";
				my $FILTER="PASS";
				my $QUAL=100;
				my $TOTAL=$wild_type+$Heterozygous+$Homozygous;
				my $AF=(2*$Homozygous + $Heterozygous)/(2*($TOTAL));
				push(@INFO,$wild_type,$Heterozygous,$Homozygous,$AF);
				my $INFO=join(";",@INFO);
				$hash1{$CHROM}{$POS}{$ID}{$REF}{$ALT}{$FILTER}{$QUAL}{$INFO}=1;
				#print "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$FILTER\t$QUAL\t$INFO\n";
			}
		}
	}
}else{print "Unable to open INPUT1\n";}

$time='['. timestamp(). ']'."\n";
print "Start charging hash2:$time\n";

if (open (INPUT2, $input2))
{
	#INPUT2=clinvar_Antonio_25_08_2014.vcf
	
#~ 1       100185102       rs149824246     G       A       100     PASS    FRRS1;No_modification_to_get_minimal_representation;NaN;255;RCV000060275.2
#~ 1       100316614       rs113994126     C       T       100     PASS    AGL;No_modification_to_get_minimal_representation;NaN;5|5;RCV000001153.2|RCV000020373.1
#~ 1       100316614       rs113994127     CAG     C       100     PASS    AGL;No_modification_to_get_minimal_representation;NaN;5|5;RCV000001155.2|RCV000020374.1

while(my $line=<INPUT2>)
	{
		chomp $line;
		#print "$line\n";
		unless($line=~/^base/)
		{
			if($line=~/^([^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
			{
				my @INFO=();
				my $CHROM=$1;
				my $POS=$2;
				my $ID="\.";
				my $REF=$3;
				my $ALT=$4;
				my $QUAL=$5;
				my $FILTER=$6;
				my $fields=$7;
				my @fields_tmp=split(";",$fields);
				my $CLNSIG=$fields_tmp[3];
				my $Mode_of_inheritance=$fields_tmp[2];
				$hash2{$CHROM}{$POS}{$ID}{$REF}{$ALT}{$FILTER}{$QUAL}{'CLNSIG'}{$CLNSIG}=1;
				$hash2{$CHROM}{$POS}{$ID}{$REF}{$ALT}{$FILTER}{$QUAL}{'MOI'}{$Mode_of_inheritance}=1;
				#~ print "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$FILTER\t$QUAL\t'CLNSIG'\t$CLNSIG\n";
				#~ print "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$FILTER\t$QUAL\t'MOI'\t$Mode_of_inheritance\n";
			}
		}
	}
}else{print "Unable to open INPUT2\n";}


$time='['. timestamp(). ']'."\n";
print "Start PRINTING:$time\n";

if (open(OUTPUT, '>'.$output))
{
	print OUTPUT "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	\n";
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
	foreach my $INFO_tok(sort keys%{$hash1{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$FILTER_tok}{$QUAL_tok}})
	{
		my @INFO_tmp=split(";",$INFO_tok);
		#~ $wild_type,$Heterozygous,$Homozygous,$AF
		my $wild_type=$INFO_tmp[0];
		my $Heterozygous=$INFO_tmp[1];
		my $Homozygous=$INFO_tmp[2];
		my $AF=$INFO_tmp[3];
		
		print OUTPUT "$CHROM_tok\t$POS_tok\t$ID_tok\t$REF_tok\t$ALT_tok\t$FILTER_tok\t$QUAL_tok\tANTONIO_ISTVAN;WT=$wild_type;HET=$Heterozygous;HOM=$Homozygous;AF_COUNTS=$AF\n";
		
	
	}
	}
	}
	}
	}
	}	
	}	
	}
	foreach my $CHROM_tok(sort{$a<=>$b} keys %hash2)
	{
	foreach my $POS_tok(sort{$a<=>$b} keys%{$hash2{$CHROM_tok}})
	{
	foreach my $ID_tok(sort keys%{$hash2{$CHROM_tok}{$POS_tok}})
	{
	foreach my $REF_tok(sort keys%{$hash2{$CHROM_tok}{$POS_tok}{$ID_tok}})
	{
	foreach my $ALT_tok(sort keys%{$hash2{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}})
	{
	foreach my $FILTER_tok(sort keys%{$hash2{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}})
	{
	foreach my $QUAL_tok(sort keys%{$hash2{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$FILTER_tok}})
	{
		my @CLNSIG_tmp=sort keys%{$hash2{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$FILTER_tok}{$QUAL_tok}{'CLNSIG'}};
		my $CLNSIG=join("|",@CLNSIG_tmp);
			
		my @ModeOfInheritance_tmp=sort keys%{$hash2{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$FILTER_tok}{$QUAL_tok}{'MOF'}};
		if(scalar(@ModeOfInheritance_tmp)==0){@ModeOfInheritance_tmp=("NaN");}
		my $inheritance=join("|",@ModeOfInheritance_tmp);
		
		print OUTPUT "$CHROM_tok\t$POS_tok\t$ID_tok\t$REF_tok\t$ALT_tok\t$FILTER_tok\t$QUAL_tok\tClinVar;$CLNSIG;$inheritance\n";
		
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
