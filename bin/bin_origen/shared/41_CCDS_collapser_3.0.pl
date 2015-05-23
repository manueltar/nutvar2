##	Script to convert all the files needed to create the matrix to CCDS isoforms from ENSEMBL transcripts. 2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

use strict;
use warnings;
use Time::localtime;

my $time='['. timestamp(). ']'."\n";
print "Start charging hash chromosomes:$time\n";

my %AcceptedChromosomes=();
my @AcceptChro=('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y');
foreach my $AcceptedChro_tok(@AcceptChro)
{
	$AcceptedChromosomes{$AcceptedChro_tok}=1;
}

my %hash1=();
my %hash1_reverse=();
my %hash5=();


my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $output1=$ARGV[2];
my $input3=$ARGV[3];
my $output2=$ARGV[4];
my $input4=$ARGV[5];
my $output3=$ARGV[6];
my $input5=$ARGV[7];
my $output4=$ARGV[8];
my $input6=$ARGV[9];
my $output5=$ARGV[10];
my $input7=$ARGV[11];
my $output6=$ARGV[12];
my $input8=$ARGV[13];
my $output7=$ARGV[14];
my $input9=$ARGV[15];
my $output8=$ARGV[16];

$time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

if (open (INPUT1, $input1))
{
	#Input= gtf_tabladef_sorted_by_SYMBOL.txt
	
#~ ENST00000263100 19      A1BG    -       protein_coding  CCDS12976       ENSP00000263100 UTR_5_prime:[58864804]-[58864865]       START:[58864801
#~ ]-[58864803]     CDS:[58864770]-[58864803]       CDS:[58864658]-[58864693]       CDS:[58864294]-[58864563]       CDS:[58863649]-[58863921]
       #~ CDS:[58862757]-[58863053]       CDS:[58861736]-[58862017]       CDS:[58858719]-[58859006]       CDS:[58858391]-[58858395]       STOP:[58
#~ 858388]-[58858390]      UTR_3_prime:[58858216]-[58858387]       
#~ ENST00000600966 19      A1BG    -       protein_coding  NaNein  ENSP00000470909 UTR_5_prime:[NaN]-[NaN] START:[NaN]-[NaN]       CDS:[58864294]-[58864495]       CDS:[58863649]-[58863921]       CDS:[58863464]-[58863550]       CDS:[58862757]-[58863053]       CDS:[58861960]-[58862017]
       #~ STOP:[NaN]-[NaN]        UTR_3_prime:[NaN]-[NaN] 
#~ ENST00000282641 10      A1CF    -       protein_coding  CCDS7242        ENSP00000282641 UTR_5_prime:[52645341]-[52645435]       UTR_5_prime:[52623793]-[52623840]       UTR_5_prime:[52619701]-[52619745]       START:[52619698]-[52619700]     CDS:[52619602]-[52619700]       CDS:[52603748]-[52603882]       CDS:[52601622]-[52601752]       CDS:[52595834]-[52596072]       CDS:[52587891]-[52588055]       CDS:[52580312]-[52580409]
       #~ CDS:[52575766]-[52576039]       CDS:[52573617]-[52573822]       CDS:[52570800]-[52570936]       CDS:[52569654]-[52569802]       CDS:[52566492]-[52566640]       STOP:[52566489]-[52566491]      UTR_3_prime:[52566327]-[52566488]       


while(my $line=<INPUT1>)
	{
		chomp $line;
		#print "AAAAAAAAA:$line\n";
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)\t.+/)
		{
			#print "BBBBBBBBBBBB:$line\n";
			my $CHROM=$2;
			my $SYMBOL=$3;
			my $BIOTYPE=$4;
			my $ENST=$1;
			if(exists($AcceptedChromosomes{$CHROM}))
			{
				if($BIOTYPE eq 'protein_coding')
				{
					my $CCDS=$5;
					$hash1{$CHROM}{$SYMBOL}{$CCDS}{$ENST}=1;
					$hash1_reverse{$CHROM}{$SYMBOL}{$ENST}{$CCDS}=1;
					#print "hello_world:$CHROM\t$SYMBOL\t$CCDS\t$ENST\n";
				}
			}
		}
	}
}else{print "Unable to open INPUT1\n";die;}


$time='['. timestamp(). ']'."\n";
print "Start charging hash2:$time\n";

my %hash_NMD=();

if (open(INPUT2, $input2))
{
	#Input= new_NMD_1.txt
		
		# X       ARAF    47428936        G       A       splice_acceptor_variant ENST00000377045 NMD_positive
		# X       CA5B    15802307        G       A       splice_acceptor_variant ENST00000454127 NaN
		# X       HEPH    65486284        G       T       splice_acceptor_variant ENST00000441993 NMD_negative

while(my $line=<INPUT2>)
	{
		chomp $line;
		#print "INICIO:$line\n";
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			#print "Hello_world_I:$line\n";
			my $CHROM=$1;
			my $SYMBOL=$2;
			my $POS=$3;
			my $REF=$4;
			my $ALT=$5;
			my $Effect=$6;
			my $ENST=$7;
			my $NMD=$8;
			
			# print "$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t$Effect\t$ENST\t$NMD\n";
			
			if(exists($AcceptedChromosomes{$CHROM}))
			{
				foreach my $CCDS(sort keys%{$hash1_reverse{$CHROM}{$SYMBOL}{$ENST}})
				{
					unless($CCDS eq 'NaNein')
					{
						$hash_NMD{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{$CCDS}{$NMD}=1;
					}
				}
			}
		}
	}
}

if(open(OUTPUT1, '>'.$output1))
{
foreach my $CHROM_tok(sort{$a<=>$b}keys%hash_NMD)
{
foreach my $SYMBOL_tok(sort keys %{$hash_NMD{$CHROM_tok}})	
{
foreach my $POS_tok(sort keys %{$hash_NMD{$CHROM_tok}{$SYMBOL_tok}})	
{
foreach my $REF_tok(sort keys %{$hash_NMD{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}})	
{
foreach my $ALT_tok(sort keys %{$hash_NMD{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}})	
{
foreach my $Effect_tok(sort keys %{$hash_NMD{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}})	
{
foreach my $CCDS_tok(sort keys %{$hash_NMD{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}})	
{
foreach my $NMD_tok(sort keys %{$hash_NMD{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$CCDS_tok}})	
{

	print OUTPUT1 "$CHROM_tok\t$SYMBOL_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$Effect_tok\t$CCDS_tok\t$NMD_tok\n";

}
}	
}	
}
}	
}	
}
}
}
%hash_NMD=();

$time='['. timestamp(). ']'."\n";
print "Start charging hash3:$time\n";

my %hash_NMD_derived=();


if (open(INPUT3, $input3))
{
while(my $line=<INPUT3>)
	{
		chomp($line);
		
		## Input file= 1GK_derived_NMD_new_version.txt

		# X       GRIA3   122336600       T       TG      frameshift_variant&feature_elongation   ENST00000371264 ->DERIVED:122338313     NMD_negative
		# X       GRIA3   122336600       T       TG      frameshift_variant&feature_elongation   ENST00000371266 ->DERIVED:122338313     NMD_negative
		# X       IL9R    155231142       GT      G       frameshift_variant&feature_truncation   ENST00000369423 ->DERIVED:155232595     NMD_positive
		# X       IL9R    155231142       GT      G       frameshift_variant&feature_truncation   ENST00000540897 ->DERIVED:155232595     NMD_positive
		# X       LCA10   153149707       C       CG      frameshift_variant&feature_elongation   ENST00000357566 ->DERIVED:153149717     NMD_positive



		if($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)/)
		{
			my $CHROM=$1;
			#my $strand=$2;
			my $SYMBOL=$2;
			my $POS=$3;
			my $REF=$4;
			my $ALT=$5;
			my $Effect=$6;
			my $ENST=$7;
			my $fields=$8;
			
			if(exists($AcceptedChromosomes{$CHROM}))
			{
				foreach my $CCDS(sort keys%{$hash1_reverse{$CHROM}{$SYMBOL}{$ENST}})
				{
					unless($CCDS eq 'NaNein')
					{
						$hash_NMD_derived{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{$CCDS}{$fields}=1;
						#print "hash_relleno:$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t$Effect\t$CCDS\t$fields\n";
					}
				}
			}
		}
	}
}

if(open(OUTPUT2, '>'.$output2))
{
foreach my $CHROM_tok(sort{$a<=>$b}keys%hash_NMD_derived)
{
foreach my $SYMBOL_tok(sort keys %{$hash_NMD_derived{$CHROM_tok}})	
{
foreach my $POS_tok(sort keys %{$hash_NMD_derived{$CHROM_tok}{$SYMBOL_tok}})	
{
foreach my $REF_tok(sort keys %{$hash_NMD_derived{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}})	
{
foreach my $ALT_tok(sort keys %{$hash_NMD_derived{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}})	
{
foreach my $Effect_tok(sort keys %{$hash_NMD_derived{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}})	
{
foreach my $CCDS_tok(sort keys %{$hash_NMD_derived{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}})	
{
foreach my $fields_tok(sort keys %{$hash_NMD_derived{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$CCDS_tok}})	
{

	print OUTPUT2 "$CHROM_tok\t$SYMBOL_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$Effect_tok\t$CCDS_tok\t$fields_tok\n";
	#print "$CHROM_tok\t$SYMBOL_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$Effect_tok\t$CCDS_tok\t$fields_tok\n";

}
}	
}	
}
}	
}	
}
}
}

%hash_NMD_derived=();
my %hash_ENSG=();

$time='['. timestamp(). ']'."\n";
print "Start charging hash4:$time\n";

if (open(INPUT4, $input4))
{

# Input=gtf_output_ENSG.txt

# ENSG00000000003 TSPAN6  -       99883667        99894988
# ENSG00000000005 TNMD    +       99839799        99854882
# ENSG00000000419 DPM1    -       49551404        49575092
# ENSG00000000457 SCYL3   -       169818772       169863408

while(my $line=<INPUT4>)
	{
		chomp($line);
		if($line=~/([^\t]+)\t([^\t]+)\t(.+)/)
		{
			my $ENSG=$1;
			my $SYMBOL=$2;
			my $fields=$3;
			
			$hash_ENSG{$ENSG}{$SYMBOL}{$fields}=1;
			#print "****$ENSG\t$SYMBOL\n";
		}
	}
}

if(open(OUTPUT3, '>'.$output3))
{
foreach my $ENSG_tok(sort keys%hash_ENSG)
{
foreach my $SYMBOL_tok(sort keys %{$hash_ENSG{$ENSG_tok}})	
{
foreach my $fields_tok(sort keys %{$hash_ENSG{$ENSG_tok}{$SYMBOL_tok}})	
{
	print OUTPUT3 "$ENSG_tok\t$SYMBOL_tok\t$fields_tok\n";
	#print "$ENSG_tok\t$SYMBOL_tok\t$fields_tok\n";
}
}	
}	
}

my %hash_ENST=();

$time='['. timestamp(). ']'."\n";
print "Start charging hash5:$time\n";


if (open(INPUT5, $input5))
{

# Input= gtf_output_ENST.txt

# X       ENSG00000000003 ENST00000373020 protein_coding  99883667        99891803
# X       ENSG00000000003 ENST00000494424 processed_transcript    99888439        99894988
# X       ENSG00000000003 ENST00000496771 processed_transcript    99887538        99891686
# X       ENSG00000000005 ENST00000373031 protein_coding  99839799        99854882
# X       ENSG00000000005 ENST00000485971 processed_transcript    99848621        99852528
while(my $line=<INPUT5>)
	{
		chomp($line);
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)/)
		{
			my $CHROM=$1;
			my $ENSG=$2;
			my $ENST=$3;
			my $BIOTYPE=$4;
			my $fields=$5;
			#print "$CHROM\t$ENSG\t$ENST\t$BIOTYPE\n";
			if(exists($AcceptedChromosomes{$CHROM}))
			{
				foreach my $SYMBOL_tok(sort keys %{$hash_ENSG{$ENSG}})
				{
					if($BIOTYPE eq 'protein_coding')
					{
						foreach my $CCDS(sort keys%{$hash1_reverse{$CHROM}{$SYMBOL_tok}{$ENST}})
						{
							unless($CCDS eq 'NaNein')
							{
								$hash_ENST{$CHROM}{$ENSG}{$CCDS}{$BIOTYPE}{$fields}=1;
								#print "hash_relleno:$CHROM\t$ENSG\t$CCDS\t$BIOTYPE\t$fields\n";
							}
						}
					}
				}
			}
		}
	}
}

if(open(OUTPUT4, '>'.$output4))
{
foreach my $CHROM_tok(sort{$a<=>$b}keys%hash_ENST)
{
foreach my $ENSG_tok(sort keys %{$hash_ENST{$CHROM_tok}})
{
foreach my $CCDS_tok(sort keys %{$hash_ENST{$CHROM_tok}{$ENSG_tok}})
{
foreach my $BIOTYPE_tok(sort keys %{$hash_ENST{$CHROM_tok}{$ENSG_tok}{$CCDS_tok}})
{
foreach my $fields_tok(sort keys %{$hash_ENST{$CHROM_tok}{$ENSG_tok}{$CCDS_tok}{$BIOTYPE_tok}})
{

	print OUTPUT4 "$CHROM_tok\t$ENSG_tok\t$CCDS_tok\t$BIOTYPE_tok\t$fields_tok\n";
	#print "$CHROM_tok\t$ENSG_tok\t$CCDS_tok\t$BIOTYPE_tok\t$fields_tok\n";
	
}
}	
}
}
}
}

$time='['. timestamp(). ']'."\n";
print "Start charging hash6:$time\n";

my %hash_ENST_full_condensed=();

if (open(INPUT6, $input6))
{
	# Input=ENST_table_full_condensed.txt
		
		#~ 1       A3GALT2 -       ENST00000330379 [33777791-33777822]     
		#~ 1       A3GALT2 -       ENST00000330379 ENST00000442999 [33772370-33773054]     [33777653-33777790]     
		#~ 1       A3GALT2 -       ENST00000442999 [33778102-33778191]     [33778408-33778491]     [33786677-33786699]     
		#~ 1       AADACL3 +       ENST00000332530 [12776344-12776347]     
		#~ 1       AADACL3 +       ENST00000332530 ENST00000359318 [12780885-12780948]     [12785189-12785960]     
		#~ 1       AADACL3 +       ENST00000359318 [12779480-12779693]
	while(my $line=<INPUT6>)
	{
		chomp($line);
		#~ print "Hello worldI:$line**\n";
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)\t$/)
		{
			my %order_hash=();
			#~ print "Hello worldII:$line**\n";
			my $CHROM=$1;
			my $SYMBOL=$2;
			my $strand=$3;
			my $fields=$4;
			#~ print "$CHROM\t$strand\t$SYMBOL\t**$fields**\n";
			
			my @fields_tmp=split(/\t/,$fields);
			foreach my $fields_tmp_tok(@fields_tmp)
			{
				if($fields_tmp_tok=~/^ENST\d+/)
				{
					#~ print "AA:$CHROM\t$strand\t$SYMBOL\t$fields_tmp_tok:AA\n";
					$order_hash{'ENST'}{$fields_tmp_tok}=1;
				}
				else
				{
					#~ print "CC:$CHROM\t$strand\t$SYMBOL\t$fields_tmp_tok:CC\n";
					$order_hash{'INTERVAL'}{$fields_tmp_tok}=1;
				}
			}
			foreach my $ENST_tok(sort keys%{$order_hash{'ENST'}})
			{
			foreach my $interval_tok(sort keys%{$order_hash{'INTERVAL'}})
			{
				
				$hash5{$CHROM}{$strand}{$SYMBOL}{$ENST_tok}{$interval_tok}=1;
				#~ print "$CHROM\t$strand\t$SYMBOL\t$ENST_tok\t$interval_tok**\n";
			}	
			}
		}
	}
}else {print "impossible to open INPUT6\n";die;}

$time='['. timestamp(). ']'."\n";
print "Start processing hash5:$time\n";

foreach my $CHROM_tok(sort keys%hash5)
{
foreach my $strand_tok(sort keys %{$hash5{$CHROM_tok}})
{
foreach my $SYMBOL_tok(sort keys %{$hash5{$CHROM_tok}{$strand_tok}})
{
foreach my $ENST_tok(sort keys %{$hash5{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}})
{
foreach my $interval_tok(sort keys %{$hash5{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENST_tok}})
{
	if(exists($AcceptedChromosomes{$CHROM_tok}))
	{
		# print "**$CHROM\t$strand\t$SYMBOL\t$ENST\t$FEATURE\t$fields\n";
		foreach my $CCDS_tok(sort keys%{$hash1_reverse{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}})
		{
			unless($CCDS_tok eq 'NaNein')
			{
				$hash_ENST_full_condensed{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$CCDS_tok}{$interval_tok}=1;
				# print "RELLENO:$CHROM\t$strand\t$SYMBOL\t$CCDS_tok\t$FEATURE\t$fields\n";
			}
		}
	}
}
}	
}	
}
}

%hash5=();

if(open(OUTPUT5, '>'.$output5))
{
foreach my $CHROM_tok(sort{$a<=>$b}keys%hash_ENST_full_condensed)
{
foreach my $strand_tok(sort keys %{$hash_ENST_full_condensed{$CHROM_tok}})
{
foreach my $SYMBOL_tok(sort keys %{$hash_ENST_full_condensed{$CHROM_tok}{$strand_tok}})
{
foreach my $CCDS_tok(sort keys %{$hash_ENST_full_condensed{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}})
{	
	print OUTPUT5 "$CHROM_tok\t$strand_tok\t$SYMBOL_tok\t$CCDS_tok\t";
	
foreach my $fields_tok(sort keys %{$hash_ENST_full_condensed{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$CCDS_tok}})
{

	print OUTPUT5 "$fields_tok\t";
	#print "$CHROM_tok\t$strand_tok\t$SYMBOL_tok\t$string\t$fields_tok\n";
}	
	print OUTPUT5 "\n";
}
}
}
}
}

%hash_ENST_full_condensed=();

$time='['. timestamp(). ']'."\n";
print "Start charging hash7:$time\n";

my %hash_APPRIS=();

if (open (INPUT7, $input7))
{
		#INPUT1=/Documentos/APPRIS/appris_principal_isoform_gencode_19_15_10_2014.txt

#~ METTL25 ENSG00000127720 ENST00000248306 CCDS9024.1      appris_principal
#~ OR13A1  ENSG00000256574 ENST00000374401 CCDS31188.1     appris_principal
#~ OR13A1  ENSG00000256574 ENST00000553795 CCDS31188.1     appris_principal
#~ OR13A1  ENSG00000256574 ENST00000536058 CCDS31188.1     appris_principal


while(my $line=<INPUT7>)
	{
		chomp $line;
		#print "$line\n";
		if($line=~/^([^\s]+)\s([^\s]+)\s([^\s]+)\s([^\s]+)\s([^\s]+)/)
		{
			my $SYMBOL=$1;
			my $ENSG=$2;
			my $ENST=$3;
			my $CCDS=$4;
			my $Appris=$5;
			#print "hello_world:$SYMBOL\t$ENST\t$CCDS\t$Appris\n";
			if($Appris eq 'appris_principal')
			{
				foreach my $CHROM_tok(sort keys %hash1_reverse)
				{
					if(exists($AcceptedChromosomes{$CHROM_tok}))
					{
						# print "**$CHROM\t$strand\t$SYMBOL\t$ENST\t$FEATURE\t$fields\n";
						foreach my $CCDS_tok(sort keys%{$hash1_reverse{$CHROM_tok}{$SYMBOL}{$ENST}})
						{
							unless($CCDS_tok eq 'NaNein')
							{
								$hash_APPRIS{$SYMBOL}{$ENSG}{$ENST}{$CCDS_tok}{$Appris}=1;
								#"hello_world2:$SYMBOL\t$ENSG\t$ENST\t$Appris\t$CCDS\n";
							}
						}
					}
				}
			}
		}
	}
}else{print "Unable to open INPUT7\n";}

if(open(OUTPUT6, '>'.$output6))
{
foreach my $SYMBOL_tok(sort keys %hash_APPRIS)
{
foreach my $ENSG_tok(sort keys %{$hash_APPRIS{$SYMBOL_tok}})
{
foreach my $ENST_tok(sort keys %{$hash_APPRIS{$SYMBOL_tok}{$ENSG_tok}})
{
foreach my $CCDS_tok(sort keys %{$hash_APPRIS{$SYMBOL_tok}{$ENSG_tok}{$ENST_tok}})
{
foreach my $Apppris_tok(sort keys %{$hash_APPRIS{$SYMBOL_tok}{$ENSG_tok}{$ENST_tok}{$CCDS_tok}})
{
	print OUTPUT6 "$SYMBOL_tok $ENSG_tok $ENST_tok $CCDS_tok $Apppris_tok\n";
}
}
}
}
}
}

$time='['. timestamp(). ']'."\n";
print "Start charging hash7:$time\n";

my %hash_Pervasive=();

if (open (INPUT8, $input8))
{
	#INPUT1=/Documentos/Pervasive.txt
	
# transcript_id gene_id transcript_biotype      number_of_samples_in_which_the_gene_is_expressed
# ENST00000000233 ENSG00000004059 protein_coding 15
# ENST00000002165 ENSG00000001036 protein_coding 16
# ENST00000003912 ENSG00000001461 protein_coding 16
# ENST00000004921 ENSG00000006074 protein_coding 12
# ENST00000004982 ENSG00000004776 protein_coding 15
# ENST00000005257 ENSG00000006451 protein_coding 16
# ENST00000005340 ENSG00000004975 protein_coding 16
# ENST00000006101 ENSG00000005243 protein_coding 16
# ENST00000006275 ENSG00000007255 protein_coding 16
# ENST00000006724 ENSG00000007306 protein_coding 1
# ENST00000007390 ENSG00000007520 protein_coding 16
# ENST00000007699 ENSG00000006047 protein_coding 12
# ENST00000008938 ENSG00000008438 protein_coding 5



while(my $line=<INPUT8>)
	{
		chomp $line;
		#print "$line\n";
		if($line=~/^(ENST[^\s]+)\s([^\s]+)\s([^\s]+)\s(.+)/)
		{
			my $ENSG=$2;
			my $ENST=$1;
			my $BIOTYPE=$3;
			my $samples=$4;
			
			foreach my $CHROM_tok(sort keys %hash1_reverse)
			{
				foreach my $SYMBOL_tok(sort keys %{$hash_ENSG{$ENSG}})
					{
						if($BIOTYPE eq 'protein_coding')
						{
							foreach my $CCDS(sort keys%{$hash1_reverse{$CHROM_tok}{$SYMBOL_tok}{$ENST}})
							{
								unless($CCDS eq 'NaNein')
								{
									$hash_Pervasive{$CCDS}{$ENSG}{$BIOTYPE}{$samples}=1;
									#print "hash_relleno:$CCDS\t$ENSG\t$BIOTYPE\t$samples\n";
								}
							}
						}
					}
			}
		}
	}
}else{print "Unable to open INPUT8\n";}

if(open(OUTPUT7, '>'.$output7))
{
	foreach my $CCDS_tok(sort keys %hash_Pervasive)
	{
		foreach my $ENSG_tok(sort keys %{$hash_Pervasive{$CCDS_tok}})
		{
		foreach my $BIOTYPE_tok(sort keys %{$hash_Pervasive{$CCDS_tok}{$ENSG_tok}})
		{
		foreach my $samples_tok(sort keys %{$hash_Pervasive{$CCDS_tok}{$ENSG_tok}{$BIOTYPE_tok}})
		{
			print OUTPUT7 "$CCDS_tok $ENSG_tok $BIOTYPE_tok $samples_tok\n";
		}
		}
		}
	}
}

$time='['. timestamp(). ']'."\n";
print "Start charging hash8:$time\n";

my %hash_first_table=();

if (open(INPUT9, $input9))
{
	#Input= 1GK_def_first_table.txt
		
		#1       A3GALT2 33772799        G       T       stop_gained     Constitutive|100;Principal_isoform;No_Pervasive_information     ENST00000330379 50.2923976608187        NMD_negative;NaN
		#1       A3GALT2 33772799        G       T       stop_gained     Constitutive|100;Principal_isoform;No_Pervasive_information     ENST00000442999 42.156862745098 NMD_negative;NaN
while(my $line=<INPUT9>)
	{
		chomp $line;
		#print "INICIO:$line\n";
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			my $CHROM=$1;
			my $SYMBOL=$2;
			my $POS=$3;
			my $REF=$4;
			my $ALT=$5;
			my $Effect=$6;
			my $Dispensable_group=$7;
			my $ENST=$8;
			my $Percentage_group=$9;
			#print "**$Percentage_group\n";
			my $NMD_group=$10;
			
			if(exists($AcceptedChromosomes{$CHROM}))
					{
						# print "**$CHROM\t$strand\t$SYMBOL\t$ENST\t$FEATURE\t$fields\n";
						foreach my $CCDS_tok(sort keys%{$hash1_reverse{$CHROM}{$SYMBOL}{$ENST}})
						{
							unless($CCDS_tok eq 'NaNein')
							{
								$hash_first_table{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{$Dispensable_group}{$CCDS_tok}{$Percentage_group}{$NMD_group}=1;
								#print "HASH_RELLENO:$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t$Effect\t$Dispensable_group\t$CCDS_tok\t$Percentage_group\t$NMD_group";
							}
						}
					}
		}
	}
}

if(open(OUTPUT8, '>'.$output8))
{
foreach my $CHROM_tok(sort{$a<=>$b} keys %hash_first_table)
{
foreach my $SYMBOL_tok(sort keys %{$hash_first_table{$CHROM_tok}})
{	
foreach my $POS_tok(sort keys %{$hash_first_table{$CHROM_tok}{$SYMBOL_tok}})
{	
foreach my $REF_tok(sort keys %{$hash_first_table{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}})
{	
foreach my $ALT_tok(sort keys %{$hash_first_table{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}})
{	
foreach my $Effect_tok(sort keys %{$hash_first_table{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}})
{	
foreach my $Dispensable_group_tok(sort keys %{$hash_first_table{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}})
{	
foreach my $CCDS_tok(sort keys %{$hash_first_table{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$Dispensable_group_tok}})
{	
foreach my $Percentage_tok(sort keys %{$hash_first_table{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$Dispensable_group_tok}{$CCDS_tok}})
{	
foreach my $NMD_tok(sort keys %{$hash_first_table{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$Dispensable_group_tok}{$CCDS_tok}{$Percentage_tok}})
{	
	print OUTPUT8 "$CHROM_tok\t$SYMBOL_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$Effect_tok\t".
					"$Dispensable_group_tok\t$CCDS_tok\t$Percentage_tok\t$NMD_tok\n";
	
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
}

	
sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
