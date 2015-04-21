###############
###############

#!/usr/bin/perl

use strict;
# use warnings;
use Time::localtime;
use Bio::EnsEMBL::Registry;



my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $output1=$ARGV[2];
my $output2=$ARGV[3];

## Open input1 carrying protein information: equivalence ENST-AC_isoform, Domain_IDs, Domain_Coordinates, Postranslational sites and 'Active'sites and
## their coordinates.


# Print Time log

my $time='['. timestamp(). ']'."\n";
print "Start charging hash 1:$time\n";
my %hash1=();
my %hash2=();
my %hash3=();
my %hash4=();

if(open (INPUT1, $input1))
{
	while(my $line=<INPUT1>)
	{
		chomp($line);
		## Input1= UNIPROT_GLOBAL_POST_UNIQUE_PLUS_COORDINATES.txt

		#	A0A183  ENSG00000235942 ENST00000431011 ENSP00000411070 FEATURE:Late cornified envelope protein 6A.__0__240__CHAIN      CCDS:CCDS44227
		#	A0AUZ9  ENSG00000144445 ENST00000281772 ENSP00000281772 FEATURE:KAT8 regulatory NSL complex subunit 1-__0__2961__CHAIN  FEATURE:N6-ace
		#	A0AV02  ENSG00000221955 ENST00000393469 ENSP00000377112 FEATURE:Solute carrier family 12 member 8.__0__2142__CHAIN      IPR:IPR004841_
		
		my @tmp=split("\t",$line);
		my $ENST=$tmp[2];
		my $ENSG=$tmp[1];
		foreach my $tmp_tok(@tmp)
		{
			if ($tmp_tok=~/^FEATURE:(.+)/)
			{
				my $FEATURE_tok=$1;
				my @FEATURE_tmp=split("__",$FEATURE_tok);
				my $distance_begin=$FEATURE_tmp[1];
				my $distance_feature=$FEATURE_tmp[2];
				my $FEATURE_ID=$FEATURE_tmp[0];
				$hash1{$ENSG}{$ENST}=1;
				#print "$ENSG\t$ENST\n";
			}
		}
	}
}else {print "Unable to open INPUT1\n";}

$time='['. timestamp(). ']'."\n";
print "Tiempo de carga hash_2:$time\n";

#exit;

if (open(INPUT2, $input2))
{
	# input2= /home/bioinfo/Escritorio/Proyecto_clasificador/Documentos/UNIPROT/uniprot_sprot.dat
	
#	ID   CFAH_HUMAN              Reviewed;        1231 AA.
#	AC   P08603; A5PL14; P78435; Q14570; Q2TAZ5; Q38G77; Q5TFM3; Q8N708;
#	AC   Q9NU86;
#	CC   -!- ALTERNATIVE PRODUCTS:
#	CC       Event=Alternative splicing; Named isoforms=2;
#	CC       Name=1;
#	CC         IsoId=P08603-1; Sequence=Displayed;
#	CC       Name=2; Synonyms=FHL-1;
#	CC         IsoId=P08603-2; Sequence=VSP_001190, VSP_001191;
#	DR   Ensembl; ENST00000367429; ENSP00000356399; ENSG00000000971.
#	SQ   SEQUENCE   1231 AA;  139096 MW;  3C26D62A2BF9BFEE CRC64;
#     MRLLAKIICL MLWAICVAED CNELPPRRNT EILTGSWSDQ TYPEGTQAIY KCRPGYRSLG
#     NVIMVCRKGE WVALNPLRKC QKRPCGHPGD TPFGTFTLTG GNVFEYGVKA VYTCNEGYQL
#     LGEINYRECD TDGWTNDIPI CEVVKCLPVT APENGKIVSS AMEPDREYHF GQAVRFVCNS
#     GYKIEGDEEM HCSDDGFWSK EKPKCVEISC KSPDVINGSP ISQKIIYKEN ERFQYKCNMG
#     YEYSERGDAV CTESGWRPLP SCEEKSCDNP YIPNGDYSPL RIKHRTGDEI TYQCRNGFYP
#     ATRGNTAKCT STGWIPAPRC TLKPCDYPDI KHGGLYHENM RRPYFPVAVG KYYSYYCDEH
#     FETPSGSYWD HIHCTQDGWS PAVPCLRKCY FPYLENGYNQ NYGRKFVQGK SIDVACHPGY
#     ALPKAQTTVT CMENGWSPTP RCIRVKTCSK SSIDIENGFI SESQYTYALK EKAKYQCKLG
#     YVTADGETSG SITCGKDGWS AQPTCIKSCD IPVFMNARTK NDFTWFKLND TLDYECHDGY
#     ESNTGSTTGS IVCGYNGWSD LPICYERECE LPKIDVHLVP DRKKDQYKVG EVLK
	#my @tmp_START_positions=();
	
my @ENSG_ENST_tmp=();
my $FLAG=0;
my @sequence_tmp=();

while(my $line=<INPUT2>)
	{
		chomp($line);
		if ($line=~/^\/\//)
		{
			#print "hello_world_1\n";
			foreach my $ENSG_ENST_tok(@ENSG_ENST_tmp)
			{
				my @tmp=split("__",$ENSG_ENST_tok);
				my $ENSG=$tmp[0];
				my $ENST=$tmp[1];
				if(exists($hash1{$ENSG}{$ENST}))
				{
					#print "hello_world_1:$ENSG\t$ENST\n";
					my $seq=join("",@sequence_tmp);
					$hash2{$ENSG}{$ENST}{$seq}=1;
					#print "$ENSG\t$ENST\t$seq\n";
					
				}
			}
			@ENSG_ENST_tmp=();
			$FLAG=0;
			@sequence_tmp=();
		}
		elsif($line=~/^DR\s+Ensembl; (ENST[^\;]+); ENSP[^\;]+; (ENSG[^\.]+)\..*/)
		{
			my $ENSG=$2;
			my $ENST=$1;
			my $string=join("__",$ENSG,$ENST);
			#print "$string\n";
			push(@ENSG_ENST_tmp,$string);
			
		}
		
		if($line=~/^SQ\s+SEQUENCE/)
		{
			$FLAG=1;
		}
		if ($FLAG==1 && $line=~/^\s+(.+)\s*$/)
		{
			my $sequence_line=$1;
			#print "$sequence_line\n";
			$sequence_line=~s/\s+//g;
			push (@sequence_tmp, $sequence_line);
		}
	}
close (INPUT2);
}else {print "impossible to open INPUT2\n";die;}

#print "La posición es:$tmp_START_positions[0] \n";

$time='['. timestamp(). ']'."\n";
print "Tiempo de creación hash3:$time\n";

my %hash3=();
my $registry = "Bio::EnsEMBL::Registry";
$registry->load_registry_from_url('mysql://anonymous@ensembldb.ensembl.org');
my $gene_adaptor=$registry->get_adaptor( 'Human', 'Core', 'Gene' );

my @tmp_genes=sort keys %hash1;
my $ENSG_gene="NaN";
my $ENST_transcript="NaN";
my $seq="NaN";

foreach my $gene_tok(@tmp_genes)
{
	#my $gene=$gene_adaptor->fetch_by_stable_id($gene_tok);
	if(defined(my $gene=$gene_adaptor->fetch_by_stable_id($gene_tok)))
	{
		$ENSG_gene=$gene->stable_id();

		#print "AAAAA:$ENSG_gene:A\n";
		my @transcripts = @{ $gene->get_all_Transcripts };

		foreach my $transcript(@transcripts)
		{
			my $ENST_transcript=$transcript->stable_id();
			#print "$ENST_transcript\n";
			if(exists($hash1{$ENSG_gene}{$ENST_transcript}))
			{
				if(defined($transcript->translation()))
				{
					$seq=$transcript->translate()->seq();
					$seq=~s/\.|\s+//g;
					#print "CANDIDATE:$seq".":::::::::::\n";
					$hash3{$ENSG_gene}{$ENST_transcript}{$seq}=1;
					#print "$ENSG_gene\t$ENST_transcript\n";
				}
			}
		}
	}
}

my %hash4=();

if (open(OUTPUT,'>'.$output1))
{
	foreach my $ENSG_tok(sort keys %hash1)
	{
	foreach my $ENST_tok(sort keys %{$hash1{$ENSG_tok}})
	{
		foreach my $UNIPROT_seq_tok(sort keys %{$hash2{$ENSG_tok}{$ENST_tok}})
		{	
			foreach my $ENSEMBL_seq_tok(sort keys %{$hash3{$ENSG_tok}{$ENST_tok}})
			{
				if ($UNIPROT_seq_tok eq $ENSEMBL_seq_tok)
				{
					print OUTPUT ">$ENSG_tok\t$ENST_tok\n";
					print OUTPUT "$UNIPROT_seq_tok\n";
					print OUTPUT "$ENSEMBL_seq_tok\n";
				}
				else
				{
					$hash4{$ENSG_tok}{$ENST_tok}{$UNIPROT_seq_tok}{$ENSEMBL_seq_tok}=1;
				}
			}
		}
	}
	}
}

if (open(OUTPUT2, '>'.$output2))
{
	foreach my $ENSG_tok(sort keys %hash4)
	{
	foreach my $ENST_tok(sort keys %{$hash4{$ENSG_tok}})
	{
	foreach my $UNIPROT_seq_tok(sort keys %{$hash4{$ENSG_tok}{$ENST_tok}})
	{
	foreach my $ENSEMBL_seq_tok(sort keys %{$hash4{$ENSG_tok}{$ENST_tok}{$UNIPROT_seq_tok}})
	{
		print OUTPUT2 	">$ENSG_tok\t$ENST_tok\n";
		print OUTPUT2 "$UNIPROT_seq_tok\n";
		print OUTPUT2 "$ENSEMBL_seq_tok\n";
	}		
	}	
	}
	}
}

$time='['. timestamp(). ']'."\n";
print "Fin del script:$time\n";


sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
