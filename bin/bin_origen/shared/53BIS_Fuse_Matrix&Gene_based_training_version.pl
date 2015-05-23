use strict;
use warnings;
use Time::localtime;

my %hash1=();
my %hash2=();
my %pRDG_hash=();
my %II_hash=();
my %Antiviral_hash=();
my %ISG_hash=();
my %OMIMrecessive_hash=();
my %RVIS_hash=();

my $input1=$ARGV[0];
#my $input2=$ARGV[1];
my $input3=$ARGV[1];
my $input4=$ARGV[2];
my $input5=$ARGV[3];
my $input6=$ARGV[4];
my $input7=$ARGV[5];
my $input8=$ARGV[6];
my $output=$ARGV[7];


my $time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

if (open (INPUT1, $input1))
{
	#~ Input=Matrix_snpeff.txt

#~ baseNCBI37	Alleles_REF>ALT	GeneNameHGNC	NumIsoformsInQueryGene	ratioIsoformsBearingTheVariant	ratioAffectedIsoforms_stop-gained	ratioAffectedIsoforms_frameshift	ratioAffectedIsoforms_splice	ratioAffectedIsoforms_coding-synonymous	ratioAffectedIsoforms_missense	ratioAffectedIsoformsTargetedbyNMD	ratioAffectedIsoformsTargetedby_derived_NMD	IsPrincipalIsoformAffected	IsWithinLongestCCDS	IsWithinPervasiveIsoform	LongestCCDSLength	PercentagePrincipalOrLongestCCDSAffected	DomainINFOAvailable	PercentageOfDomainPositionsAffected	maxPercDomainAffected	NumberOfDomains100Damage	DomainMatchedSiteINFOAvailable	PercentageOfSitePositionsAffected	maxPercSiteAffected	NumberOfSites100Damage	SiteMatched	Pathogenicity_Tag	Confidence_Tag	Location_Tag
#~ #X:74334588	C>T	ABCB7	6	4	0.5	0	0	0	0.5	0.5	0	1	0	0	2324	100	1	100	100	1	1	1	100	100	1	0	0	1	3
#~ 1:69620	TA>T	OR4F5	1	1	0	1	0	0	0	1	1	1	1	NaN	915	42.1857923497268	1	40.8371040723982	100	3	1	0	5.75539568345324	NaN	NaN	NaN		NaN	1	1

while(my $line=<INPUT1>)
	{
		chomp $line;
		#print "$line\n";
		unless($line=~/^base/ || $line=~/^#/)
		{
			if($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)/)
			{
				my $CHROM_POS=$1;
				my $REF_ALT=$2;
				my $SYMBOL=$3;
				my $fields=$4;

				$hash1{$CHROM_POS}{$REF_ALT}{$SYMBOL}{$fields}=1;
				if($SYMBOL eq 'UGT2A1')
				{print "$CHROM_POS\t$REF_ALT\t$SYMBOL\t************$fields**\n";}
			}
		}
	}
}else{print "Unable to open INPUT1\n";}


$time='['. timestamp(). ']'."\n";
print "Start charging hash3:$time\n";

if (open(INPUT3, $input3))
{
	## Input file= Genes_pRDG.txt


#~ Gene_Name       pRDG
#~ A1BG    0.31615396
#~ A1CF    0.359992897
#~ A2LD1   0.104827515
#~ A2M     0.796183485
#~ A2ML1   0.083306957


	while(my $line=<INPUT3>)
	{
		chomp($line);
		unless($line=~/^Gene_Name/)
		{
			if($line=~/([^\t]+)\t([^\t]+)/)
			{
				my $SYMBOL=$1;
				my $pRDG=$2;
				#~ print "$CHROM\t$POS\t$REF\t$ALT\t$Pathogenicity_Tag\t$Confidence_Tag\t$Location_Tag\n";
				$pRDG_hash{$SYMBOL}{$pRDG}=1;
				#~ print "$SYMBOL\t$pRDG\n";
			}
		}
	}
}else {print "impossible to open INPUT3\n";die;}


$time='['. timestamp(). ']'."\n";
print "Start charging hash4:$time\n";

if (open(INPUT4, $input4))
{
	## Input file= Genes_AllInnateImmunity.txt


#~ isInnateImmunity
#~ FGR
#~ CFH
#~ CFTR
#~ LAP3
#~ CASP10



	while(my $line=<INPUT4>)
	{
		chomp($line);
		unless($line=~/^isInnateImmunity/)
		{
			if($line=~/([^\t]+)/)
			{
				my $SYMBOL=$1;
				#~ print "$CHROM\t$POS\t$REF\t$ALT\t$Pathogenicity_Tag\t$Confidence_Tag\t$Location_Tag\n";
				$II_hash{$SYMBOL}=1;
				#~ print "$SYMBOL\n";
			}
		}
	}
}else {print "impossible to open INPUT4\n";die;}

$time='['. timestamp(). ']'."\n";
print "Start charging hash5:$time\n";

if (open(INPUT5, $input5))
{
	## Input file= Genes_Antiviral.txt


#~ isAntiViral
#~ AGPAT9
#~ ANKRD22
#~ APOL1
#~ APOL6
#~ B4GALT5
	while(my $line=<INPUT5>)
	{
		chomp($line);
		unless($line=~/^isAntiViral/)
		{
			if($line=~/([^\t]+)/)
			{
				my $SYMBOL=$1;
				$Antiviral_hash{$SYMBOL}=1;
				#~ print "$SYMBOL\n";
			}
		}
	}
}else {print "impossible to open INPUT5\n";die;}

$time='['. timestamp(). ']'."\n";
print "Start charging hash6:$time\n";

if (open(INPUT6, $input6))
{
	## Input file= Genes_ISGs.txt

#~ isISG
#~ ABCA9
#~ ABLIM3
#~ ABTB2
#~ ACSL1
#~ ADAMDEC1

	while(my $line=<INPUT6>)
	{
		chomp($line);
		unless($line=~/^isISG/)
		{
			if($line=~/([^\t]+)/)
			{
				my $SYMBOL=$1;
				$ISG_hash{$SYMBOL}=1;
				#~ print "$SYMBOL\n";
			}
		}
	}
}else {print "impossible to open INPUT6\n";die;}

$time='['. timestamp(). ']'."\n";
print "Start charging hash7:$time\n";

if (open(INPUT7, $input7))
{
	## Input file= Genes_OMIMrecessive.txt

#~ OMIM_Recessive
#~ AASS
#~ ABCA12
#~ ABCA4
#~ ABCB6
#~ ABCB7


	while(my $line=<INPUT7>)
	{
		chomp($line);
		unless($line=~/^OMIM_Recessive/)
		{
			if($line=~/([^\t]+)/)
			{
				my $SYMBOL=$1;
				$OMIMrecessive_hash{$SYMBOL}=1;
				#~ print "$SYMBOL\n";
			}
		}
	}
}else {print "impossible to open INPUT7\n";die;}

$time='['. timestamp(). ']'."\n";
print "Start charging hash8:$time\n";

if (open(INPUT8, $input8))
{
	#~ ## Input file= Genes_RVIS.txt 

#~ Gene_Name       RVIS
#~ A1BG    -0.466531444
#~ A1CF    -0.378346116
#~ A2LD1   0.189398298
#~ A2M     0.099178678


	while(my $line=<INPUT8>)
	{
		chomp($line);
		unless($line=~/^Gene_Name/)
		{
			if($line=~/([^\t]+)\t([^\t]+)/)
			{
				my $SYMBOL=$1;
				my $RVIS=$2;
				$RVIS_hash{$SYMBOL}{$RVIS}=1;
				#~ print "$SYMBOL\t$RVIS\n";
			}
		}
	}
}else {print "impossible to open INPUT8\n";die;}


$time='['. timestamp(). ']'."\n";
print "Start PRINTING:$time\n";

if (open(OUTPUT, '>'.$output))
{
		print OUTPUT "baseNCBI37\tAlleles_REF>ALT\tGeneNameHGNC\tNumIsoformsInQueryGene\tratioIsoformsBearingTheVariant\t".
				"ratioAffectedIsoforms_stop-gained\tratioAffectedIsoforms_frameshift\tratioAffectedIsoforms_splice\tratioAffectedIsoforms_coding-synonymous\tratioAffectedIsoforms_missense\t".
				"PEJMAN1\tPEJMAN2\tratioAffectedIsoformsTargetedbyNMD\tratioAffectedIsoformsTargetedby_derived_NMD\t".
				"IsPrincipalIsoformAffected\tIsWithinLongestCCDS\tIsWithinPervasiveIsoform\t".
				"PEJMAN3\tPEJMAN4\tLongestCCDSLength\tPercentagePrincipalOrLongestCCDSAffected\t".
				"DomainINFOAvailable\tPercentageOfDomainPositionsAffected\tmaxPercDomainAffected\tNumberOfDomains100Damage\tDomainMatched\t".
				"SiteINFOAvailable\tPercentageOfSitePositionsAffected\tmaxPercSiteAffected\tNumberOfSites100Damage\tSiteMatched\t".
				"PEJMAN5\tPEJMAN6\tPathogenicity_Tag\tCredible_Tag(-1,rare, 0,not_credible,1,common)\tLocation_Tag\t".
				"PEJMAN7\tPEJMAN8\tIsInnateImmunity\tIsAntiviral\tIsISG\tIsOMIMrecessive\tpRDG_score\tRVIS_score\n";
	
	print OUTPUT "#X:74334588\tC>T\tABCB7\t6\t4\t".
				"0.5\t0\t0\t0\t0.5\t".
				"PEJMAN1\tPEJMAN2\t0.5\t0\t".
				"1\t0\t0\t".
				"PEJMAN3\tPEJMAN4\t2324\t100\t".
				"1\t100\t100\t1\t1\t".
				"1\t100\t100\t1\t0\t".
				"PEJMAN5\tPEJMAN6\t1\t-1\t2\t".
				"PEJMAN7\tPEJMAN8\t1\t0\t1\t0\t0.83\t0.5\n";
				
	foreach my $CHROM_POS_tok(sort keys %hash1)
	{
	foreach my $REF_ALT_tok(sort keys %{$hash1{$CHROM_POS_tok}})
	{
	foreach my $SYMBOL_tok(sort keys %{$hash1{$CHROM_POS_tok}{$REF_ALT_tok}})
	{		
	foreach my $fields_tok(sort keys %{$hash1{$CHROM_POS_tok}{$REF_ALT_tok}{$SYMBOL_tok}})
	{
		print OUTPUT"$CHROM_POS_tok\t$REF_ALT_tok\t$SYMBOL_tok\t$fields_tok\tPEJMAN7\tPEJMAN8\t";
		my ($pRDG,$RVIS,$InnateImmunity,$Antiviral,$ISG,$OMIMrecessive)="NaN"x6;
		if(exists($II_hash{$SYMBOL_tok})){$InnateImmunity=1;}else{$InnateImmunity=0;}
		print OUTPUT"$InnateImmunity\t";
		if(exists($Antiviral_hash{$SYMBOL_tok})){$Antiviral=1;}else{$Antiviral=0;}
		print OUTPUT"$Antiviral\t";
		if(exists($ISG_hash{$SYMBOL_tok})){$ISG=1;}else{$ISG=0;}
		print OUTPUT"$ISG\t";
		if(exists($OMIMrecessive_hash{$SYMBOL_tok})){$OMIMrecessive=1;}else{$OMIMrecessive=0;}
		print OUTPUT"$OMIMrecessive\t";
		
					if(exists($pRDG_hash{$SYMBOL_tok}))
					{
						foreach my $pRDG_tok(sort keys %{$pRDG_hash{$SYMBOL_tok}})
						{
							$pRDG=$pRDG_tok;	
						}
					}
		else{$pRDG="NaN";}
		if(exists($RVIS_hash{$SYMBOL_tok}))
		{
			foreach my $RVIS_tok(sort keys %{$RVIS_hash{$SYMBOL_tok}})
			{
				$RVIS=$RVIS_tok;
				#~ print OUTPUT"$RVIS\t";
			}
		}
		else{$RVIS="NaN";}

		print OUTPUT"$pRDG\t";
		print OUTPUT"$RVIS\n";
		
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
