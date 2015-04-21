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
	#~ baseNCBI37      Alleles_REF>ALT GeneNameHGNC    NumIsoformsInQueryGene  ratioIsoformsBearingTheVariant  ratioAffectedIsoforms_stop-gained       ratioAffectedIsoforms_frameshift        ratioAffectedIsoforms_splice    ratioAffectedIsoforms_coding-synonymous ratioAffectedIsoforms_missense  ratioAffectedIsoformsTargetedbyNMD      ratioAffectedIsoformsTargetedby_derived_NMD     IsPrincipalIsoformAffected      IsWithinLongestCCDS     IsWithinPervasiveIsoform        LongestCCDSLength       PercentagePrincipalOrLongestCCDSAffected        DomainINFOAvailable     PercentageOfDomainPositionsAffected     maxPercDomainAffected   NumberOfDomains100Damage        DomainMatched
   #~ SiteINFOAvailable       PercentageOfSitePositionsAffected       maxPercSiteAffected     NumberOfSites100Damage  SiteMatched

#~ #X:74334588     C>T     ABCB7   6       4       0.5     0       0       0       0.5     0.5     0       1       0       0       2324
     #~ 100     1       100     100     1       1       1       100     100     1       0

#~ 1:12779563      C>CT    AADACL3 2       0.5     1       1       0       0       0       1       0       1       1       NaN     1050
    #~ 92.0952380952381        1       94.8132780082988        100     3       1       0       NaN     NaN     NaN     NaN     

while(my $line=<INPUT1>)
	{
		chomp $line;
		#print "$line\n";
		unless($line=~/^base/ || $line=~/^#/)
		{
			if($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)\t$/)
			{
				my $CHROM_POS=$1;
				my $REF_ALT=$2;
				my $SYMBOL=$3;
				my $fields=$4;
				#print "hello_world:$CHROM_POS\t$REF_ALT\t$SYMBOL\t**$fields\n";
				#~ exit;
				my @CHROM_POS_tmp=split(":",$CHROM_POS);
				my @REF_ALT_tmp=split(">",$REF_ALT);
				my $CHROM=$CHROM_POS_tmp[0];
				my $POS=$CHROM_POS_tmp[1];
				my $REF=$REF_ALT_tmp[0];
				my $ALT=$REF_ALT_tmp[1];

				$hash1{$CHROM}{$POS}{$REF}{$ALT}{$SYMBOL}{$fields}=1;
				#~ print "$CHROM\t$POS\t$REF\t$ALT\t$SYMBOL\t$fields\n";
			}
		}
	}
}else{print "Unable to open INPUT1\n";}

$time='['. timestamp(). ']'."\n";
print "Start charging hash2:$time\n";

if (open(INPUT2, $input2))
{
	## Input file= Train\&Test.txt

#~ 1       21889712        .       G       A       100     PASS    ANTONIO_ISTVAN;WT=6909;HET=2;HOM=0;AF_COUNTS=0.000144696860078136;lowCI=4.47524797802744e-05;highCI=0.000522482649862658;ClinVar;5;NaN  1       1       3
#~ 1       21890587        .       G       A       100     PASS    ANTONIO_ISTVAN;WT=3521;HET=3;HOM=0;AF_COUNTS=0.000425652667423383;lowCI=0.000154567925541428;highCI=0.00124272386275024;ClinVar;5|5;NaN 2       1       3
#~ 1       21890632        .       G       A       100     PASS    ANTONIO_ISTVAN;WT=8445;HET=17;HOM=0;AF_COUNTS=0.00100449066414559;lowCI=0.000629791946004726;highCI=0.00160609361431182;ClinVar;5|5|5;NaN       3       1       3

	while(my $line=<INPUT2>)
	{
		chomp($line);
		unless($line=~/#CHROM/)
		{
			if($line=~/([^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
			{
				my $CHROM=$1;
				my $POS=$2;
				my $REF=$3;
				my $ALT=$4;
				my $Pathogenicity_Tag=$5;
				my $Confidence_Tag=$6;
				my $Location_Tag=$7;
				#~ print "$CHROM\t$POS\t$REF\t$ALT\t$Pathogenicity_Tag\t$Confidence_Tag\t$Location_Tag\n";
				if(exists($hash1{$CHROM}{$POS}{$REF}{$ALT}))
				{
					$hash2{$CHROM}{$POS}{$REF}{$ALT}{$Pathogenicity_Tag}{$Confidence_Tag}{$Location_Tag}=1;
				}
			}
		}
	}
}else {print "impossible to open INPUT2\n";die;}



$time='['. timestamp(). ']'."\n";
print "Start PRINTING:$time\n";

if (open(OUTPUT, '>'.$output))
{
	print OUTPUT "baseNCBI37\tAlleles_REF>ALT\tGeneNameHGNC\tNumIsoformsInQueryGene\tratioIsoformsBearingTheVariant\t".
				"ratioAffectedIsoforms_stop-gained\tratioAffectedIsoforms_frameshift\tratioAffectedIsoforms_splice\tratioAffectedIsoforms_coding-synonymous\tratioAffectedIsoforms_missense\t".
				"ratioAffectedIsoformsTargetedbyNMD\tratioAffectedIsoformsTargetedby_derived_NMD\t".
				"IsPrincipalIsoformAffected\tIsWithinLongestCCDS\tIsWithinPervasiveIsoform\t".
				"LongestCCDSLength\tPercentagePrincipalOrLongestCCDSAffected\t".
				"DomainINFOAvailable\tPercentageOfDomainPositionsAffected\tmaxPercDomainAffected\tNumberOfDomains100Damage\tDomainMatched\t".
				"SiteINFOAvailable\tPercentageOfSitePositionsAffected\tmaxPercSiteAffected\tNumberOfSites100Damage\tSiteMatched\t".
				#"pRDG_score\tRVIS_score\tIsInnateImmunity?\tIsAntiviral?\tIsISG?\tIsOMIMrecessive\n".
				"Pathogenicity_Tag\tCredible_Tag(-1,rare, 0,not_credible,1,common)\tLocation_Tag\n";
	
	print OUTPUT "#X:74334588\tC>T\tABCB7\t6\t4\t".
				"0.5\t0\t0\t0\t0.5\t".
				"0.5\t0\t".
				"1\t0\t0\t".
				"2324\t100\t".
				"1\t100\t100\t1\t1\t".
				"1\t100\t100\t1\t0\t".
				#"0.83\t0.5\t1\t0\t1\t0\n".
				"0\t1\t3\n";
				
	foreach my $CHROM_tok(sort{$a<=>$b}keys %hash1)
	{
	foreach my $POS_tok(sort{$a<=>$b}keys %{$hash1{$CHROM_tok}})
	{
	foreach my $REF_tok(sort keys %{$hash1{$CHROM_tok}{$POS_tok}})
	{
	foreach my $ALT_tok(sort keys %{$hash1{$CHROM_tok}{$POS_tok}{$REF_tok}})
	{
	foreach my $SYMBOL_tok(sort keys %{$hash1{$CHROM_tok}{$POS_tok}{$REF_tok}{$ALT_tok}})
	{		
	foreach my $fields_tok(sort keys %{$hash1{$CHROM_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$SYMBOL_tok}})
	{
		foreach my $Pathogenicity_Tag(sort keys %{$hash2{$CHROM_tok}{$POS_tok}{$REF_tok}{$ALT_tok}})
		{
		foreach my $Confidence_Tag(sort keys %{$hash2{$CHROM_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Pathogenicity_Tag}})
		{
		foreach my $Location_Tag(sort keys %{$hash2{$CHROM_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Pathogenicity_Tag}{$Confidence_Tag}})
		{
			print OUTPUT "$CHROM_tok:$POS_tok\t$REF_tok>$ALT_tok\t$SYMBOL_tok\t$fields_tok\tPEJMAN5\tPEJMAN6\t$Pathogenicity_Tag\t$Confidence_Tag\t$Location_Tag\n";
			print "$CHROM_tok:$POS_tok\t$REF_tok>$ALT_tok\t$SYMBOL_tok\t$fields_tok\t****$Pathogenicity_Tag\t$Confidence_Tag\t$Location_Tag\n";
		}	
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
