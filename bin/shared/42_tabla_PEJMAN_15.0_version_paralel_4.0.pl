##	Script to create the matrix of numerical data using CCDS. 2014.
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
my %hash2=();
my %hash5=();
my %primer_hash=();
my %hash_CCDS=();
my %hash_ENSG=();
my %hash_ENSG2=();
my %hash_APPRIS=();
my %Pervasive_hash=();
my %Absolute_LONGEST_hash=();
my %hash_total_isoforms=();

my %RESULTS_hash=();


my %hash_distances=();
my %hash3=();
my %hash4=();



my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $input3=$ARGV[2];
my $input4=$ARGV[3];
my $input5=$ARGV[4];
my $input6=$ARGV[5];
my $input7=$ARGV[6];
my $input8=$ARGV[7];
my $output=$ARGV[8];

$time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

if (open(INPUT1, $input1))
{
	#Input= 1GK_def_first_table_CCDS.txt
		
		#~ X       AFF2    147743470       G       A       synonymous_variant      Alternative|80;No_Principal_isoform;No_Pervasive_information
    #~ CCDS14684       94.3808797355708;DOMAIN|94.4501018329939|1;SITE|100|0   NMD_positive;IPR007797**1|94.4501018329939|1;MOD_RES**1|10
#~ 0|0;MOD_RES**2|100|0;MOD_RES**3|100|0;NaN
#~ X       AFF2    147743470       G       A       synonymous_variant      Alternative|80;No_Principal_isoform;No_Pervasive_information
    #~ CCDS55520       94.5402298850575;DOMAIN|94.611561600837|1;SITE|100|0    NMD_positive;IPR007797**1|94.611561600837|1;MOD_RES**1|100
#~ |0;MOD_RES**2|100|0;MOD_RES**3|100|0;NaN
#~ X       AFF2    147743513       T       A       missense_variant        Alternative|80;No_Principal_isoform;No_Pervasive_information
    #~ CCDS14684       93.2875667429443;DOMAIN|93.3553971486762|1;SITE|100|0   NMD_positive;IPR007797**1|93.3553971486762|1;MOD_RES**1|10
#~ 0|0;MOD_RES**2|100|0;MOD_RES**3|100|0;NaN
#~ X       AFF2    147743513       T       A       missense_variant        Alternative|80;No_Principal_isoform;No_Pervasive_information
    #~ CCDS55520       93.4169278996865;DOMAIN|93.4867904786817|1;SITE|100|0   NMD_positive;IPR007797**1|93.4867904786817|1;MOD_RES**1|10
#~ 0|0;MOD_RES**2|100|0;MOD_RES**3|100|0;NaN

while(my $line=<INPUT1>)
	{
		chomp $line;
		#print "INICIO:$line\n";
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			#print "PARSEO:$line\n";
			my %domain_NMD_hash=();
			my %domain_NMD_hash2=();
			
			my $CHROM=$1;
			my $SYMBOL=$2;
			my $POS=$3;
			my $REF=$4;
			my $ALT=$5;
			my $Effect=$6;
			my $CCDS=$7;
			my $Percentage_group=$8;
			## # # print "**$Percentage_group\n";
			my $NMD_group=$9;
			
			
			if ($Effect=~/stop_gained/||$Effect=~/frameshift_variant/||$Effect=~/splice_donor_variant/||$Effect=~/splice_acceptor_variant/ || $Effect=~/inframe_deletion/ ||$Effect=~/disruptive_inframe_deletion/|| $Effect=~/missense_variant/ || $Effect=~/synonymous_variant/)
			{
				$hash1{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{'RATIO'}{$Effect}{$CCDS}=1;
				$hash1{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{'ISOFORM'}{$CCDS}{$Effect}=1;
				
				if ($Effect=~/stop_gained/||$Effect=~/frameshift_variant/||$Effect=~/splice_donor_variant/||$Effect=~/splice_acceptor_variant/ || $Effect=~/inframe_deletion/ ||$Effect=~/disruptive_inframe_deletion/)
				{#
					#print "*******************$NMD_group\n";
					$primer_hash{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}=1;
					
					my @Percentage_tmp=split(";",$Percentage_group);
					
					my $seq_percentage="NaN";
					my $domain_matched="NaN";
					my $domain_percentage="NaN";
					my $site_matched="NaN";
					my $site_percentage="NaN";
					
					foreach my $Percentage_tmp_tok(@Percentage_tmp)
					{
						if($Percentage_tmp_tok !~ /^DOMAIN/ && $Percentage_tmp_tok !~ /^SITE/)
						{
							$seq_percentage=$Percentage_tmp_tok;
						}
						elsif($Percentage_tmp_tok=~ /^DOMAIN/)
						{
							my @tmp_DOMAIN=split(/\|/,$Percentage_tmp_tok);
							
							$domain_percentage=$tmp_DOMAIN[1];
							$domain_matched=$tmp_DOMAIN[2];	
						}
						elsif($Percentage_tmp_tok=~ /^SITE/)
						{
							my @tmp_SITE=split(/\|/,$Percentage_tmp_tok);
							
							$site_percentage=$tmp_SITE[1];
							$site_matched=$tmp_SITE[2];	
						}
					}
					
					## # # print "La línea es:$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t$Effect\t$CCDS\t$Percentage_group\t$seq_percentage\t$domain_matched\t$domain_percentage\t$site_matched\t$site_percentage\n";
					$hash1{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{'EFFECT'}{$Effect}{$CCDS}{$seq_percentage}{$domain_percentage}{$domain_matched}{$site_percentage}{$site_matched}=1;
					if($SYMBOL eq 'PANK4')
					{
						#~ print "Coarse_domain_site:$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t'EFFECT'\t$Effect\t$CCDS\t$seq_percentage\t$domain_percentage\t$domain_matched\t$site_percentage\t$site_matched\n";
					}
					
					
					my @NMD_tmp=split(";",$NMD_group);
					
					my $Detail_domain_info=0;
					my $Detail_site_info=0;
					
					foreach my $NMD_tmp_tok(@NMD_tmp)
					{
						if($NMD_tmp_tok !~ /^NMD_/ && $NMD_tmp_tok !~ /^NaN/)
						{
							if($NMD_tmp_tok =~ /^IPR/)
							{
								$Detail_domain_info=1;
								my @IPR_tmp=split(/\|/,$NMD_tmp_tok);
								foreach my $IPR_tmp_tok(@IPR_tmp)
								{
									my $IPR=$IPR_tmp[0];
									my $Percentage_IPR=$IPR_tmp[1];
									my $matched_IPR=$IPR_tmp[2];
									$domain_NMD_hash{'IPR'}{$Detail_domain_info}{$IPR}{$Percentage_IPR}{$matched_IPR}=1;
									$domain_NMD_hash2{'IPR'}{$Detail_domain_info}{$Percentage_IPR}=1;
									# print "Hash_llenado_IPR:'IPR'\t'IPR'\t$Detail_domain_info\t$IPR\t$Percentage_IPR\t$matched_IPR\n";
								}
							
							}
							#################################################################################################3
							elsif($NMD_tmp_tok !~ /^IPR/) 
							{
								if($NMD_tmp_tok !~/REGION/ && $NMD_tmp_tok !~ /DOMAIN/ && $NMD_tmp_tok !~ /CHAIN/ && $NMD_tmp_tok !~ /CROSSLNK/ && $NMD_tmp_tok !~ /DISULFID/)
								{
									$Detail_site_info=1;
									my @SITE_tmp=split(/\|/,$NMD_tmp_tok);
										foreach my $SITE_tmp_tok(@SITE_tmp)
										{
											my $SITE=$SITE_tmp[0];
											my $Percentage_SITE=$SITE_tmp[1];
											my $matched_SITE=$SITE_tmp[2];
											$domain_NMD_hash{'SITE'}{$Detail_site_info}{$SITE}{$Percentage_SITE}{$matched_SITE}=1;
											$domain_NMD_hash2{'SITE'}{$Detail_site_info}{$Percentage_SITE}=1;
											# print "Hash_llenado_SITE:'SITE'\t$Detail_site_info\t$SITE\t$Percentage_SITE\t$matched_SITE\n";
										}
								}
							}
						}
					}

		#########################################################################################################################################################
					
					my @N_IPR_max_damage=();
					my %matched_IPR=();
					my $max_percentage_IPR="NaN";
					my $IPR_presence=0;
					
					foreach my $Detail_domain_info_tok(sort keys%{$domain_NMD_hash{'IPR'}})
					{
						# # # # # print "INICIO:$Detail_domain_info_tok\n";
						if($Detail_domain_info_tok ==1)
						{
							$IPR_presence=1;
							my @Percentage_IPR_tmp=reverse sort{$a<=>$b}keys%{$domain_NMD_hash2{'IPR'}{$Detail_domain_info_tok}};
							$max_percentage_IPR=$Percentage_IPR_tmp[0];
							foreach my $IPR_tok(sort keys %{$domain_NMD_hash{'IPR'}{$Detail_domain_info_tok}})
							{
								foreach my $Percentage_IPR_tmp_tok(sort keys %{$domain_NMD_hash{'IPR'}{$Detail_domain_info_tok}{$IPR_tok}})
								{
									# # # # # print "***** $Percentage_IPR_tmp_tok\n";
									if($Percentage_IPR_tmp_tok == 100)
									{
										push(@N_IPR_max_damage,$IPR_tok);
									}
									foreach my $matched_IPR_tmp_tok( sort keys %{$domain_NMD_hash{'IPR'}{$Detail_domain_info_tok}{$IPR_tok}{$Percentage_IPR_tmp_tok}})
									{
										# # # # # print "AAA:$matched_IPR_tmp_tok\n";
										
										$matched_IPR{$matched_IPR_tmp_tok}{$IPR_tok}=1;
									}
								}
							}
						}
						elsif($Detail_domain_info_tok ==0)
						{
							# Do nothing
						}
					}
					
					my $matched_IPR=0;
					my $number_of_100_damage_IPR=0;
					
					if($IPR_presence==1)
					{
						unless (scalar(@N_IPR_max_damage)==0){$number_of_100_damage_IPR=scalar(@N_IPR_max_damage);}
						foreach my $IPR_matched_tok(sort keys %matched_IPR)
						{
							if($IPR_matched_tok == 1)
							{
								$matched_IPR=1;
							}
						}
						$hash2{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{$CCDS}{'IPR'}{$IPR_presence}{$max_percentage_IPR}{$number_of_100_damage_IPR}{$matched_IPR}=1;
						if($SYMBOL eq 'PANK4')
						{
							#~ print "Hash de memoria:$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t'IPR'\t$Effect\t$CCDS\t$IPR_presence\t$max_percentage_IPR\t$number_of_100_damage_IPR\t$matched_IPR\n";
						}
					}
					elsif($IPR_presence==0)
					{
						$matched_IPR="NaN";
						$number_of_100_damage_IPR="NaN";
						$hash2{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{$CCDS}{'IPR'}{$IPR_presence}{$max_percentage_IPR}{$number_of_100_damage_IPR}{$matched_IPR}=1;
						if($SYMBOL eq 'PANK4')
						{
							#~ print "Hash de memoria:$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t'IPR'\t$Effect\t$CCDS\t$IPR_presence\t$max_percentage_IPR\t$number_of_100_damage_IPR\t$matched_IPR\n";
						}
					}


		########################################################################################################################################################
					my @N_SITE_max_damage=();
					my %matched_SITE=();
					my $max_percentage_SITE="NaN";
					my $SITE_presence=0;
					
					foreach my $Detail_site_info_tok(sort keys%{$domain_NMD_hash{'SITE'}})
					{
						# # # # # print "INICIO:$Detail_site_info_tok\n";
						if($Detail_site_info_tok ==1)
						{
							$SITE_presence=1;
							my @Percentage_SITE_tmp=reverse sort{$a<=>$b}keys%{$domain_NMD_hash2{'SITE'}{$Detail_site_info_tok}};
							$max_percentage_SITE=$Percentage_SITE_tmp[0];
							foreach my $SITE_tok(sort keys %{$domain_NMD_hash{'SITE'}{$Detail_site_info_tok}})
							{
								foreach my $Percentage_SITE_tmp_tok(sort keys %{$domain_NMD_hash{'SITE'}{$Detail_site_info_tok}{$SITE_tok}})
								{
									# # # # # print "***** $Percentage_SITE_tmp_tok\n";
									if($Percentage_SITE_tmp_tok == 100)
									{
										push(@N_SITE_max_damage,$SITE_tok);
									}
									foreach my $matched_SITE_tmp_tok( sort keys %{$domain_NMD_hash{'SITE'}{$Detail_site_info_tok}{$SITE_tok}{$Percentage_SITE_tmp_tok}})
									{
										# # # # # print "BBB:$matched_SITE_tmp_tok\n";
										
										$matched_SITE{$matched_SITE_tmp_tok}{$SITE_tok}=1;
									}
								}
							}
						}
						elsif($Detail_site_info_tok ==0)
						{
							# Do nothing
						}
					}
					
					my $matched_SITE=0;
					my $number_of_100_damage_SITE=0;
					
					if($SITE_presence==1)
					{
						unless (scalar(@N_SITE_max_damage)==0){$number_of_100_damage_SITE=scalar(@N_SITE_max_damage);}
						foreach my $SITE_matched_tok(sort keys %matched_SITE)
						{
							if($SITE_matched_tok == 1)
							{
								$matched_SITE=1;
							}
						}
						$hash2{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{$CCDS}{'SITE'}{$SITE_presence}{$max_percentage_SITE}{$number_of_100_damage_SITE}{$matched_SITE}=1;
						if($SYMBOL eq 'PANK4')
						{
							#~ print "Hash de memoria:$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t'SITE'\t$Effect\t$CCDS\t$SITE_presence\t$max_percentage_SITE\t$number_of_100_damage_SITE\t$matched_SITE\n";
						}
					}
					elsif($SITE_presence==0)
					{
						$matched_SITE="NaN";
						$number_of_100_damage_SITE="NaN";
						$hash2{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{$CCDS}{'SITE'}{$SITE_presence}{$max_percentage_SITE}{$number_of_100_damage_SITE}{$matched_SITE}=1;
						if($SYMBOL eq 'PANK4')
						{
							#~ print "Hash de memoria:$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t'SITE'\t$Effect\t$CCDS\t$SITE_presence\t$max_percentage_SITE\t$number_of_100_damage_SITE\t$matched_SITE\n";
						}
					}
				}
			}#
		}##
	}
}else {print "impossible to open INPUT1\n";die;}


$time='['. timestamp(). ']'."\n";
print "Start charging hash2:$time\n";


if (open(INPUT2, $input2))
{
	#Input= new_NMD_1.txt
		
		# X       ARAF    47428936        G       A       splice_acceptor_variant CCDS00000377045 NMD_positive
		# X       CA5B    15802307        G       A       splice_acceptor_variant CCDS00000454127 NaN
		# X       HEPH    65486284        G       T       splice_acceptor_variant CCDS00000441993 NMD_negative

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
			my $CCDS=$7;
			my $NMD=$8;
			
			# print "$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t$Effect\t$CCDS\t$NMD\n";
			if(exists($primer_hash{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}))
			{
				if($NMD eq 'NMD_positive')
				{
							my $NMD_binary=1;
							$hash2{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{$CCDS}{'NMD'}{$NMD_binary}=1;
							if($SYMBOL eq 'PANK4')
							{
								#~ print "Hash de memoria:$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t'NMD'\t$Effect\t$CCDS\t'NMD'\t$NMD\n";
							}
				}
				elsif($NMD eq 'NMD_negative')
				{
							my $NMD_binary=0;
							$hash2{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{$CCDS}{'NMD'}{$NMD_binary}=1;
							if($SYMBOL eq 'PANK4')
							{
								#~ print "Hash de memoria:$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t'NMD'\t$Effect\t$CCDS\t'NMD'\t$NMD\n";
							}
				}
				elsif($NMD eq 'NaN')
				{
							my $NMD_binary="NaN";
							$hash2{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{$CCDS}{'NMD'}{$NMD_binary}=1;
							if($SYMBOL eq 'PANK4')
							{
								#~ print "Hash de memoria:$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t'NMD'\t$Effect\t$CCDS\t'NMD'\t$NMD\n";
							}
				}
			}
		}
	}
}


$time='['. timestamp(). ']'."\n";
print "Start charging hash3:$time\n";


if (open(INPUT3, $input3))
{
while(my $line=<INPUT3>)
	{
		chomp($line);
		
		## Input file= 1GK_derived_NMD_new_version.txt

		# X       GRIA3   122336600       T       TG      frameshift_variant&feature_elongation   CCDS00000371264 ->DERIVED:122338313     NMD_negative
		# X       GRIA3   122336600       T       TG      frameshift_variant&feature_elongation   CCDS00000371266 ->DERIVED:122338313     NMD_negative
		# X       IL9R    155231142       GT      G       frameshift_variant&feature_truncation   CCDS00000369423 ->DERIVED:155232595     NMD_positive
		# X       IL9R    155231142       GT      G       frameshift_variant&feature_truncation   CCDS00000540897 ->DERIVED:155232595     NMD_positive
		# X       LCA10   153149707       C       CG      frameshift_variant&feature_elongation   CCDS00000357566 ->DERIVED:153149717     NMD_positive



		if($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t->DERIVED:([^\t]+)\t([^\t]+)/)
		{
			my $CHROM=$1;
			#my $strand=$2;
			my $SYMBOL=$2;
			my $POS=$3;
			my $REF=$4;
			my $ALT=$5;
			my $Effect=$6;
			my $CCDS=$7;
			my $POS_derived=$8;
			my $NMD_derived=$9;
			#print "hello_world:$CHROM\t$strand\t$SYMBOL\t$POS\t$REF\t$ALT\t$Pervasive\n";
			
			#print"$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t$Effect\t'TRANSCRIPT'\t$CCDS\t'DERIVED_PTC'\t$NMD_derived\n";
			if(exists($primer_hash{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}))
			{
					my $NMD_derived_binary="NaN";
					if($NMD_derived eq 'NMD_positive')
					{
						$NMD_derived_binary=1;
						$hash2{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{$CCDS}{'NMD_derived'}{$NMD_derived_binary}=1;
						if($SYMBOL eq 'PANK4')
						{
							#~ print "Hash de memoria:$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t'NMD_derived'\t$Effect\t$CCDS\t'NMD_derived'\t$NMD_derived_binary\n";
						}
					}
					elsif($NMD_derived eq 'NMD_negative')
					{
						$NMD_derived_binary=0;
						$hash2{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{$CCDS}{'NMD_derived'}{$NMD_derived_binary}=1;
						if($SYMBOL eq 'PANK4')
						{
							#~ print "Hash de memoria:$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t'NMD_derived'\t$Effect\t$CCDS\t'NMD_derived'\t$NMD_derived_binary\n";
						}

					}
					elsif($NMD_derived eq 'NaN')
					{
						$NMD_derived_binary="NaN";
						$hash2{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{$CCDS}{'NMD_derived'}{$NMD_derived_binary}=1;
						if($SYMBOL eq 'PANK4')
						{
							#~ print "Hash de memoria:$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t'NMD_derived'\t$Effect\t$CCDS\t'NMD_derived'\t$NMD_derived_binary\n";
						}
					}
			}
		}
	}
}

#~ exit;

$time='['. timestamp(). ']'."\n";
print "Start charging hash4:$time\n";


if (open(INPUT4, $input4))
{

# Input= gtf_output_CCDS.txt

# X       ENSG00000000003 CCDS00000373020 protein_coding  99883667        99891803
# X       ENSG00000000003 CCDS00000494424 processed_transcript    99888439        99894988
# X       ENSG00000000003 CCDS00000496771 processed_transcript    99887538        99891686
# X       ENSG00000000005 CCDS00000373031 protein_coding  99839799        99854882
# X       ENSG00000000005 CCDS00000485971 processed_transcript    99848621        99852528

#~ Y       ENSG00000012817 CCDS14794       protein_coding  21867301        21906809
#~ Y       ENSG00000012817 CCDS55554       protein_coding  21867306        21906647
#~ Y       ENSG00000012817 CCDS55555       protein_coding  21867303        21906825
#~ Y       ENSG00000067048 CCDS14782       protein_coding  15016019        15030451
#~ Y       ENSG00000067048 CCDS14782       protein_coding  15016742        15032390

while(my $line=<INPUT4>)
	{
		chomp($line);
		#~ print "$line\n";
		#~ exit;
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t.+/)
		{
			my $CHROM=$1;
			my $ENSG=$2;
			my $CCDS=$3;
			my $BIOTYPE=$4;
			#~ print "HELLO_WORLD:$CHROM\t$ENSG\t$CCDS\t$BIOTYPE\n";
			if(exists($AcceptedChromosomes{$CHROM}))
			{
				if($BIOTYPE eq 'protein_coding')
				{
					$hash_CCDS{$ENSG}{$CCDS}=1;
					#~ print "**$ENSG\t$CCDS\n";
				}
			}
		}
	}
}

$time='['. timestamp(). ']'."\n";
print "Start charging hash5:$time\n";

if (open(INPUT5, $input5))
{

# Input=gtf_output_ENSG.txt

# ENSG00000000003 TSPAN6  -       99883667        99894988
# ENSG00000000005 TNMD    +       99839799        99854882
# ENSG00000000419 DPM1    -       49551404        49575092
# ENSG00000000457 SCYL3   -       169818772       169863408

while(my $line=<INPUT5>)
	{
		chomp($line);
		#~ print "*********LINE:$line\n";
		if($line=~/([^\t]+)\t([^\t]+)\t.+/)
		{
			#~ print "*********LINE:$line\n";
			my $ENSG=$1;
			my $SYMBOL=$2;
			if(exists($hash_CCDS{$ENSG}))
			{
				$hash_ENSG{$SYMBOL}{$ENSG}=1;
				$hash_ENSG2{$ENSG}{$SYMBOL}=1;
				if($SYMBOL eq 'PANK4')
				{
					#~ print "****************$SYMBOL\t$ENSG\n";
				}
			}
		}
	}
}

$time='['. timestamp(). ']'."\n";
print "Start charging hash6:$time\n";

if (open(INPUT6, $input6))
{
	# Input=ENST_table_full_condensed.txt
		
#~ Y	+	BPY2	CCDS14800	[25138491-25138523]	[25138524-25138568]	[25140628-25140745]	[25143599-25143704]	[25144397-25144412]
		
	while(my $line=<INPUT6>)
	{
		chomp($line);
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)/)
		{
			my $CHROM=$1;
			my $strand=$2;
			my $SYMBOL=$3;
			my $CCDS=$4;
			my $fields=$5;

			my @fields_tmp=split("\t",$fields);
			foreach my $fields_tmp_tok(@fields_tmp)
			{
				if($fields_tmp_tok=~/\[(.+)\-(.+)\]/)
				{
					my $begin=$1;
					my $end=$2;
					my $distance=$end-$begin+1;
					if($SYMBOL eq 'PANK4')
					{
						#~ print "$CHROM\t$strand\t$SYMBOL\t$CCDS\t$end\t$begin\t$distance\n";
					}
					push(@{$hash_distances{$SYMBOL}{$CCDS}},$distance);
				}
			}
		}
	}
}else {print "impossible to open INPUT6\n";die;}


$time='['. timestamp(). ']'."\n";
print "Start charging hash7:$time\n";

if (open (INPUT7, $input7))
{
		#INPUT1=/Documentos/APPRIS/appris_principal_isoform_gencode_19_15_10_2014.txt

#~ A2M ENSG00000175899 ENST00000318602 CCDS44827 appris_principal
#~ A2ML1 ENSG00000166535 ENST00000299698 CCDS8596 appris_principal
#~ A3GALT2 ENSG00000184389 ENST00000442999 CCDS60080 appris_principal
#~ PANK4 ENSG00000128274 ENST00000249005 CCDS14041 appris_principal
#~ PANK4 ENSG00000128274 ENST00000381278 CCDS14041 appris_principal
#~ PANK4 ENSG00000128274 ENST00000401850 CCDS14041 appris_principal



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
			#print "hello_world:$SYMBOL\t$CCDS\t$CCDS\t$Appris\n";
			if($Appris eq 'appris_principal')
			{
				$hash_APPRIS{$SYMBOL}{$CCDS}=1;
				if($SYMBOL eq 'PANK4')
				{
					#~ print "Appris:$SYMBOL\t$ENSG\t$ENST\t$Appris\t$CCDS\n";
				}
			}
		}
	}
}else{print "Unable to open INPUT5\n";}

$time='['. timestamp(). ']'."\n";
print "Start charging hash8:$time\n";

if (open (INPUT8, $input8))
{
	#INPUT1=/Documentos/Pervasive.txt
	
# transcript_id gene_id transcript_biotype      number_of_samples_in_which_the_gene_is_expressed
# CCDS00000000233 ENSG00000004059 protein_coding 15
# CCDS00000002165 ENSG00000001036 protein_coding 16
# CCDS00000003912 ENSG00000001461 protein_coding 16
# CCDS00000004921 ENSG00000006074 protein_coding 12
# CCDS00000004982 ENSG00000004776 protein_coding 15
# CCDS00000005257 ENSG00000006451 protein_coding 16
# CCDS00000005340 ENSG00000004975 protein_coding 16
# CCDS00000006101 ENSG00000005243 protein_coding 16
# CCDS00000006275 ENSG00000007255 protein_coding 16
# CCDS00000006724 ENSG00000007306 protein_coding 1
# CCDS00000007390 ENSG00000007520 protein_coding 16
# CCDS00000007699 ENSG00000006047 protein_coding 12
# CCDS00000008938 ENSG00000008438 protein_coding 5



while(my $line=<INPUT8>)
	{
		chomp $line;
		#print "$line\n";
		if($line=~/^([^\s]+)\s([^\s]+)\s([^\s]+)\s.+/)
		{
			my $ENSG=$2;
			#my $ENSG=$2;
			my $CCDS=$1;
			#my $strand=$4;
			my $BIOTYPE=$3;
			#print "hello_world:$SYMBOL\t$CCDS\t$Appris\n";
			if($BIOTYPE eq 'protein_coding')
			{
				foreach my $SYMBOL_tok(sort keys %{$hash_ENSG2{$ENSG}})
				{
					$Pervasive_hash{$SYMBOL_tok}{$CCDS}=1;
					if($SYMBOL_tok eq 'PANK4')
					{
						#~ print "Pervasive:$SYMBOL_tok\t$CCDS\n";
					}
				}
			}
		}
	}
}else{print "Unable to open INPUT6\n";}

#~ exit;

$time='['. timestamp(). ']'."\n";
print "Start PROCESSING HASH_DISTANCES:$time\n";

foreach my $SYMBOL_tok(sort keys %hash_distances)
{
foreach my $CCDS_tok(sort keys %{$hash_distances{$SYMBOL_tok}})
{
	my @distance_tmp=@{$hash_distances{$SYMBOL_tok}{$CCDS_tok}};
	my $distance_CCDS = 0;
	for ( @distance_tmp )
	{
		$distance_CCDS += $_;
	}
	$hash3{$SYMBOL_tok}{$CCDS_tok}{$distance_CCDS}=1;
	$hash4{$SYMBOL_tok}{$distance_CCDS}{$CCDS_tok}=1;
	if($SYMBOL_tok eq 'PANK4')
	{
		#~ print "Hash_distancia:$SYMBOL_tok\t$CCDS_tok\t$distance_CCDS\n";
	}
}
}

$time='['. timestamp(). ']'."\n";
print "Start PROCESSING hash3:$time\n";

foreach my $SYMBOL_tok(sort keys %hash4)
{
	my @distances_tmp=reverse sort {$a<=>$b} keys%{$hash4{$SYMBOL_tok}};
	my $absolute_LONGEST_CCDS_distance=$distances_tmp[0];
	my @absolute_LONGEST_CCDS_tmp=sort keys %{$hash4{$SYMBOL_tok}{$distances_tmp[0]}};
	foreach my $absolute_LONGEST_CCDS_tmp_tok(@absolute_LONGEST_CCDS_tmp)
	{
		$Absolute_LONGEST_hash{$SYMBOL_tok}{$absolute_LONGEST_CCDS_tmp_tok}=1;
		if($SYMBOL_tok eq 'PANK4')
		{
			#~ print "$SYMBOL_tok\t'ABSOLUTE_LONGEST'\t$absolute_LONGEST_CCDS_tmp_tok\n";
		}
	}
}


$time='['. timestamp(). ']'."\n";
print "Start PROCESSING hash1:$time\n";

foreach my $SYMBOL_tok(sort keys %hash_ENSG)
{
	#~ print "hash_ENSG_symbol:$SYMBOL_tok\n";
foreach my $ENSG_tok(sort keys %{$hash_ENSG{$SYMBOL_tok}})
{
	#~ print "hash_ENSG_symbol:$ENSG_tok\n";
	my $total_CCDS=scalar(keys%{$hash_CCDS{$ENSG_tok}});
	push(@{$hash_total_isoforms{$SYMBOL_tok}},$total_CCDS);
	if($SYMBOL_tok eq 'PANK4')
	{
		my @tmp=@{$hash_total_isoforms{$SYMBOL_tok}};
		#~ print "El array es:@tmp\n";
	}
}
}

#~ exit;

foreach my $CHROM_tok(sort keys %hash1)
{
foreach my $SYMBOL_tok(sort keys %{$hash1{$CHROM_tok}})
{
	#
	my @total_CCDS=@{$hash_total_isoforms{$SYMBOL_tok}};
	my $total_CCDS = 0;
	for ( @total_CCDS )
	{
		$total_CCDS += $_;
	}
	
		if($SYMBOL_tok eq 'PANK4')
		{
			#~ print "**$CHROM_tok\t$SYMBOL_tok\t$total_CCDS**\n";
		}
	foreach my $POS_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}})
	{
	foreach my $REF_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}})
	{	
	foreach my $ALT_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}})
	{
		my $counter_missense=0;
		my $counter_synonimous=0;
		my $counter_stop=0;
		my $counter_frameshift=0;
		my $counter_splice=0;
		my $counter_inframe_del=0;
		my %hash_effects=();
		my %longest_hash=();
		my %NMD_ratio_hash=();
		my %Appris_isoform_hash=();
		if($SYMBOL_tok eq 'PANK4')
		{
			#~ print "$CHROM_tok\t$SYMBOL_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t'ISOFORM'\n";
		}
		my @CCDS_EFFECT_tmp=sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'ISOFORM'}};
		my $total_AFFECTED=scalar(@CCDS_EFFECT_tmp);
		if($SYMBOL_tok eq 'PANK4')
		{
			#~ print "El array es:@CCDS_EFFECT_tmp\tSCALAR:".scalar(@CCDS_EFFECT_tmp)."\n";
		}
		my $ratio_AFFECTED=$total_AFFECTED/$total_CCDS;
		if($SYMBOL_tok eq 'PANK4')
				{
					#~ print "$CHROM_tok\t$SYMBOL_tok\t$POS_tok\t$REF_tok\t$ALT_tok\tCCDS:$total_CCDS\tTOTAL:$total_AFFECTED\tratio:$ratio_AFFECTED\n";
				}
#####################################################	ratios	########################################################################################		
		foreach my $Effect_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'RATIO'}})
		{
				if($Effect_tok=~/missense_variant/)
				{
					my @tmp=sort keys%{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'RATIO'}{$Effect_tok}};
					$counter_missense=scalar(@tmp);
				}
				
				elsif($Effect_tok=~/synonymous_variant/)
				{
					my @tmp=sort keys%{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'RATIO'}{$Effect_tok}};
					$counter_synonimous=scalar(@tmp);			
				}
				
				elsif($Effect_tok=~/stop_gained/)
				{
					my @tmp=sort keys%{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'RATIO'}{$Effect_tok}};
					$counter_stop=scalar(@tmp);		
				}
				
				elsif($Effect_tok=~/frameshift_variant/)
				{
					my @tmp=sort keys%{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'RATIO'}{$Effect_tok}};
					$counter_frameshift=scalar(@tmp);
				}
				
				elsif($Effect_tok !~/frameshift_variant/ && $Effect_tok !~/stop_gained/)
				{
					if($Effect_tok=~/splice_donor_variant/||$Effect_tok=~/splice_acceptor_variant/)
					{
						my @tmp=sort keys%{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'RATIO'}{$Effect_tok}};
						$counter_splice=scalar(@tmp);	
					}
				}
		}
		my $ratio_stop_gained=$counter_stop/$total_AFFECTED;
		my $ratio_frameshift=$counter_frameshift/$total_AFFECTED;
		my $ratio_splice=$counter_splice/$total_AFFECTED;
		my $ratio_inframe_deletion=$counter_inframe_del/$total_AFFECTED;
		my $ratio_synonimous=$counter_synonimous/$total_AFFECTED;
		my $ratio_non_synonimous=$counter_missense/$total_AFFECTED;
##############################################################################################################################################################
		foreach my $Effect_tok_definitive(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'EFFECT'}})
		{
			my $Appris=0;
			my $Pervasive=0;
			my $absolute_LONGEST=0;
			$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{'ISOFORMS_SUMMARY'}{$total_CCDS}{$ratio_AFFECTED}{$ratio_stop_gained}{$ratio_frameshift}{$ratio_splice}{$ratio_synonimous}{$ratio_non_synonimous}=1;
			if($SYMBOL_tok eq 'PANK4')
			{
				#~ print "INICIO_ISOFORM:SUMMARY:\t$CHROM_tok\t$SYMBOL_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$total_CCDS\t$Effect_tok_definitive\t$ratio_AFFECTED\t$ratio_stop_gained\t$ratio_frameshift\t$ratio_splice\t$ratio_synonimous\t$ratio_non_synonimous\n";
			}
			
			foreach my $CCDS_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'EFFECT'}{$Effect_tok_definitive}})
			{
				if($SYMBOL_tok eq 'PANK4')
				{
					#~ print "Tránscritos:$CCDS_tok\n";
				}
				foreach my $distance_tok(sort keys %{$hash3{$SYMBOL_tok}{$CCDS_tok}})
				{
					if($SYMBOL_tok eq 'PANK4')
					{
						#~ print "Distance:$distance_tok\n";
					}
					foreach my $Percentage_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'EFFECT'}{$Effect_tok_definitive}{$CCDS_tok}})
					{
						if($SYMBOL_tok eq 'PANK4')
						{
							#~ print "Sequence_Percentage:$Percentage_tok\n";
						}
						#$longest_hash{$distance_tok}{$CCDS_tok}{$Percentage_tok}=1;
						foreach my $domain_percentage_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'EFFECT'}{$Effect_tok_definitive}{$CCDS_tok}{$Percentage_tok}})
						{
							if($SYMBOL_tok eq 'PANK4')
							{
								#~ print "Domain1:$domain_percentage_tok\n";
							}
						foreach my $domain_matched_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'EFFECT'}{$Effect_tok_definitive}{$CCDS_tok}{$Percentage_tok}{$domain_percentage_tok}})
						{
							if($SYMBOL_tok eq 'PANK4')
							{
								#~ print "Domain2:$domain_matched_tok\n";
							}
						foreach my $site_percentage_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'EFFECT'}{$Effect_tok_definitive}{$CCDS_tok}{$Percentage_tok}{$domain_percentage_tok}{$domain_matched_tok}})
						{
							if($SYMBOL_tok eq 'PANK4')
							{
								#~ print "SITE1:$site_percentage_tok\n";
							}
						foreach my $site_matched_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'EFFECT'}{$Effect_tok_definitive}{$CCDS_tok}{$Percentage_tok}{$domain_percentage_tok}{$domain_matched_tok}{$site_percentage_tok}})
						{
							if($SYMBOL_tok eq 'PANK4')
							{
								#~ print "SITE2:$site_matched_tok\n";
							}
							foreach my $NMD_tok(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'NMD'}})
							{
								if($SYMBOL_tok eq 'PANK4')
								{
									#~ print "NMD:$NMD_tok\n";
								}
								foreach my $NMD_derived_tok(sort keys%{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'NMD_derived'}})
								{
									if($SYMBOL_tok eq 'PANK4')
									{
										#~ print "NMD_derived:$NMD_derived_tok\n";
									}
									foreach my $IPR_presence_tok(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'IPR'}})
									{
										if($SYMBOL_tok eq 'PANK4')
										{
											#~ print "IPR1:$IPR_presence_tok\n";
										}
									foreach my $max_percentage_IPR_tok(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'IPR'}{$IPR_presence_tok}})
									{
										if($SYMBOL_tok eq 'PANK4')
										{
											#~ print "IPR2:$max_percentage_IPR_tok\n";
										}
									foreach my $number_of_100_damage_IPR(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'IPR'}{$IPR_presence_tok}{$max_percentage_IPR_tok}})
									{
										if($SYMBOL_tok eq 'PANK4')
										{
											#~ print "IPR3:$number_of_100_damage_IPR\n";
										}
									foreach my $matched_IPR(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'IPR'}{$IPR_presence_tok}{$max_percentage_IPR_tok}{$number_of_100_damage_IPR}})
									{
										if($SYMBOL_tok eq 'PANK4')
										{
											#~ print "IPR4:$matched_IPR\n";
										}
										foreach my $SITE_presence_tok(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'SITE'}})
										{
											if($SYMBOL_tok eq 'PANK4')
											{
												#~ print "SITE1:$SITE_presence_tok\n";
											}
										foreach my $max_percentage_SITE_tok(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'SITE'}{$SITE_presence_tok}})
										{
											if($SYMBOL_tok eq 'PANK4')
											{
												#~ print "SITE2:$max_percentage_SITE_tok\n";
											}
										foreach my $number_of_100_damage_SITE(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'SITE'}{$SITE_presence_tok}{$max_percentage_SITE_tok}})
										{
											if($SYMBOL_tok eq 'PANK4')
											{
												#~ print "SITE3:$number_of_100_damage_SITE\n";
											}
										foreach my $matched_SITE(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'SITE'}{$SITE_presence_tok}{$max_percentage_SITE_tok}{$number_of_100_damage_SITE}})
										{
											if($SYMBOL_tok eq 'PANK4')
											{
												#~ print "SITE4:$matched_SITE\n";
											}
##################################################################################################################################################################################################################################################################################################################
$longest_hash{$distance_tok}{$CCDS_tok}{$Percentage_tok}{$domain_percentage_tok}{$site_percentage_tok}{$IPR_presence_tok}{$max_percentage_IPR_tok}{$number_of_100_damage_IPR}{$matched_IPR}{$SITE_presence_tok}{$max_percentage_SITE_tok}{$number_of_100_damage_SITE}{$matched_SITE}=1;
$NMD_ratio_hash{$CCDS_tok}{$NMD_tok}{$NMD_derived_tok}=1;

if($SYMBOL_tok eq 'PANK4')
{
	#~ print "SEGUNDO_CCDS_a_CCDS:$distance_tok\t$CCDS_tok\t$Percentage_tok\t$domain_percentage_tok\t$site_percentage_tok\t$IPR_presence_tok\t$max_percentage_IPR_tok\t$number_of_100_damage_IPR\t$matched_IPR\t$SITE_presence_tok\t$max_percentage_SITE_tok\t$number_of_100_damage_SITE\t$matched_SITE\n";
}
#print "SEGUNDO_NMD_CCDS_a_CCDS:$CCDS_tok\t$NMD_tok\t$NMD_derived_tok\n";
##################################################################################################################################################################################################################################################################################################################
										}
										}	
										}
										}
									}	
									}
									}											
									}
								}
							}#
						}	
						}	
						}	
						}
					}
				}
#########################################################################################################################################################			
				if(exists($Pervasive_hash{$SYMBOL_tok}{$CCDS_tok})){$Pervasive=1;}
				
				$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{'PERVASIVE'}{$Pervasive}=1;
				
				if(exists($Absolute_LONGEST_hash{$SYMBOL_tok}{$CCDS_tok})){$absolute_LONGEST=1;}
				
				$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{'ABSOLUTE_LONGEST'}{$absolute_LONGEST}=1;
				
				if(exists($hash_APPRIS{$SYMBOL_tok}{$CCDS_tok}))
				{
					$Appris=1;

					foreach my $Percentage_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'EFFECT'}{$Effect_tok_definitive}{$CCDS_tok}})
					{
						foreach my $domain_percentage_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'EFFECT'}{$Effect_tok_definitive}{$CCDS_tok}{$Percentage_tok}})
						{
						foreach my $domain_matched_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'EFFECT'}{$Effect_tok_definitive}{$CCDS_tok}{$Percentage_tok}{$domain_percentage_tok}})
						{
						foreach my $site_percentage_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'EFFECT'}{$Effect_tok_definitive}{$CCDS_tok}{$Percentage_tok}{$domain_percentage_tok}{$domain_matched_tok}})
						{
						foreach my $site_matched_tok(sort keys %{$hash1{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{'EFFECT'}{$Effect_tok_definitive}{$CCDS_tok}{$Percentage_tok}{$domain_percentage_tok}{$domain_matched_tok}{$site_percentage_tok}})
						{
							# print "INICIO APPRIS:$CHROM_tok\t$SYMBOL_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t'EFFECT'\t$Effect_tok_definitive\t$CCDS_tok\t$Percentage_tok\t$domain_percentage_tok\t$domain_matched_tok\t$site_percentage_tok\n";
							foreach my $NMD_tok(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'NMD'}})
							{
								# print "NMD:$NMD_tok\n";
								foreach my $NMD_derived_tok(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'NMD_derived'}})
								{
									# print "NMD_derived:$NMD_derived_tok\n";
									foreach my $IPR_presence_tok(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'IPR'}})
									{
										# print "IPR1:$IPR_presence_tok\n";
									foreach my $max_percentage_IPR_tok(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'IPR'}{$IPR_presence_tok}})
									{
										# print "IPR2:$max_percentage_IPR_tok\n";
									foreach my $number_of_100_damage_IPR(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'IPR'}{$IPR_presence_tok}{$max_percentage_IPR_tok}})
									{
										# print "IPR3:$number_of_100_damage_IPR\n";
									foreach my $matched_IPR(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'IPR'}{$IPR_presence_tok}{$max_percentage_IPR_tok}{$number_of_100_damage_IPR}})
									{
										# print "IPR4:$matched_IPR\n";
										foreach my $SITE_presence_tok(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'SITE'}})
										{
											# print "SITE1:$SITE_presence_tok\n";
										foreach my $max_percentage_SITE_tok(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'SITE'}{$SITE_presence_tok}})
										{
											# print "SITE2:$max_percentage_SITE_tok\n";
										foreach my $number_of_100_damage_SITE(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'SITE'}{$SITE_presence_tok}{$max_percentage_SITE_tok}})
										{
											# print "SITE3:$number_of_100_damage_SITE\n";
										foreach my $matched_SITE(sort keys %{$hash2{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{$CCDS_tok}{'SITE'}{$SITE_presence_tok}{$max_percentage_SITE_tok}{$number_of_100_damage_SITE}})
										{
											# print "SITE4:$matched_SITE\n";
##################################################################################################################################################################################################################################################################################################################
$Appris_isoform_hash{$CCDS_tok}{$Percentage_tok}{$domain_percentage_tok}{$site_percentage_tok}{$IPR_presence_tok}{$max_percentage_IPR_tok}{$number_of_100_damage_IPR}{$matched_IPR}{$SITE_presence_tok}{$max_percentage_SITE_tok}{$number_of_100_damage_SITE}{$matched_SITE}=1;
# print "5_APPRIS:$CCDS_tok\t$Percentage_tok\t$domain_percentage_tok\t$site_percentage_tok\t$IPR_presence_tok\t$max_percentage_IPR_tok\t$number_of_100_damage_IPR\t$matched_IPR\t$SITE_presence_tok\t$max_percentage_SITE_tok\t$number_of_100_damage_SITE\t$matched_SITE\n";
##################################################################################################################################################################################################################################################################################################################
										}
										}	
										}
										}
									}	
									}
									}											
									}
								}
							}#							
						}	
						}	
						}	
						}
					}
					
				
				}
				
				else
				{
					# Do nothing
				}
				
				$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{'APPRIS_PRESENCE'}{$Appris}=1;	
			}		
##########################################################	Ratio NMD		###############################################################################
			#my @total_CCDS_affected_indep_NMD=();
			my @NMD_array=();
			my @NMD_derived_array=();
			 
			foreach my $CCDS_NMD_tok(sort keys %NMD_ratio_hash)
			{	
				 #push(@total_CCDS_affected_indep_NMD,$CCDS_NMD_tok);
				 foreach my $NMD_binary_tok(sort keys%{$NMD_ratio_hash{$CCDS_NMD_tok}})
				 {
					 if($NMD_binary_tok == 1){push(@NMD_array,$NMD_binary_tok);}
					 foreach my $NMD_derived_binary_tok(sort keys%{$NMD_ratio_hash{$CCDS_NMD_tok}{$NMD_binary_tok}})
					 {
						 if($NMD_derived_binary_tok == 1){push(@NMD_derived_array,$NMD_binary_tok);}	
					 }
				 }
			}
			 
			 
			 my $total_affected_CCDS_NMD=0;if(scalar(@NMD_array)>0){$total_affected_CCDS_NMD=scalar(@NMD_array);}
			 my $total_affected_CCDS_NMD_derived=0;if(scalar(@NMD_derived_array)>0){$total_affected_CCDS_NMD_derived=scalar(@NMD_derived_array);}
			 
			 # print "NMD:Ratio1:\t$CHROM_tok\t$SYMBOL_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$Effect_tok_definitive\t$total_affected_CCDS_NMD\t$total_AFFECTED\n";
			 # print "NMD_derived:Ratio2:\t$CHROM_tok\t$SYMBOL_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$Effect_tok_definitive\t$total_affected_CCDS_NMD_derived\t$total_AFFECTED\n";

			 
			 
			 my $ratio_NMD=$total_affected_CCDS_NMD/$total_AFFECTED;
			 my $ratio_NMD_derived=$total_affected_CCDS_NMD_derived/$total_AFFECTED;
			 
			 $RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{'NMD'}{$ratio_NMD}=1;
			 $RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{'NMD_derived'}{$ratio_NMD_derived}=1;
			
##########################################################  LONGEST->		  ##############################################################################
			
			my @longest_tmp=reverse sort{$a<=>$b}keys%longest_hash;
			# print "El array longest es:@longest_tmp\n";
			my @longest_CCDS=keys%{$longest_hash{$longest_tmp[0]}};
			# print "La longest es:$longest_tmp[0]\t**@longest_CCDS**\n";
			my @longest_CCDS_Percentage_def=();
			my @longest_CCDS_DOMAIN_percentage=();
			my @longest_CCDS_SITE_percentage=();
			
			my @longest_CCDS_IPR_presence=();my @longest_CCDS_IPR_max_percentage=();
			my @longest_CCDS_IPR_number_max=();my @longest_CCDS_IPR_matched=();

			my @longest_CCDS_SITE_presence=();my @longest_CCDS_SITE_max_percentage=();
			my @longest_CCDS_SITE_number_max=();my @longest_CCDS_SITE_matched=();
			
			foreach my $longest_CCDS_tok(@longest_CCDS)
			{
				foreach my $longest_CCDS_percentage_tok(sort keys%{$longest_hash{$longest_tmp[0]}{$longest_CCDS_tok}})
				{
					push(@longest_CCDS_Percentage_def,$longest_CCDS_percentage_tok);
					foreach my $DOMAIN_percentage_tok(sort keys%{$longest_hash{$longest_tmp[0]}{$longest_CCDS_tok}{$longest_CCDS_percentage_tok}})
					{
						push(@longest_CCDS_DOMAIN_percentage,$DOMAIN_percentage_tok);
						foreach my $SITE_percentage_tok(sort keys%{$longest_hash{$longest_tmp[0]}{$longest_CCDS_tok}{$longest_CCDS_percentage_tok}{$DOMAIN_percentage_tok}})
						{
							push(@longest_CCDS_SITE_percentage,$SITE_percentage_tok);
							foreach my $IPR_presence_tok(sort keys%{$longest_hash{$longest_tmp[0]}{$longest_CCDS_tok}{$longest_CCDS_percentage_tok}{$DOMAIN_percentage_tok}{$SITE_percentage_tok}})
							{#######################################################
								push(@longest_CCDS_IPR_presence,$IPR_presence_tok);
							foreach my $max_percentage_IPR_tok(sort keys%{$longest_hash{$longest_tmp[0]}{$longest_CCDS_tok}{$longest_CCDS_percentage_tok}{$DOMAIN_percentage_tok}{$SITE_percentage_tok}{$IPR_presence_tok}})
							{
								push(@longest_CCDS_IPR_max_percentage,$max_percentage_IPR_tok);
							foreach my $number_of_100_damage_IPR_tok(sort keys%{$longest_hash{$longest_tmp[0]}{$longest_CCDS_tok}{$longest_CCDS_percentage_tok}{$DOMAIN_percentage_tok}{$SITE_percentage_tok}{$IPR_presence_tok}{$max_percentage_IPR_tok}})
							{
								push(@longest_CCDS_IPR_number_max,$number_of_100_damage_IPR_tok);
							foreach my $matched_IPR_tok(sort keys%{$longest_hash{$longest_tmp[0]}{$longest_CCDS_tok}{$longest_CCDS_percentage_tok}{$DOMAIN_percentage_tok}{$SITE_percentage_tok}{$IPR_presence_tok}{$max_percentage_IPR_tok}{$number_of_100_damage_IPR_tok}})
							{
								push(@longest_CCDS_IPR_matched,$matched_IPR_tok);
							foreach my $SITE_presence_tok(sort keys%{$longest_hash{$longest_tmp[0]}{$longest_CCDS_tok}{$longest_CCDS_percentage_tok}{$DOMAIN_percentage_tok}{$SITE_percentage_tok}{$IPR_presence_tok}{$max_percentage_IPR_tok}{$number_of_100_damage_IPR_tok}{$matched_IPR_tok}})
							{	#################################################
								push(@longest_CCDS_SITE_presence,$SITE_presence_tok);
							foreach my $max_percentage_SITE_tok(sort keys%{$longest_hash{$longest_tmp[0]}{$longest_CCDS_tok}{$longest_CCDS_percentage_tok}{$DOMAIN_percentage_tok}{$SITE_percentage_tok}{$IPR_presence_tok}{$max_percentage_IPR_tok}{$number_of_100_damage_IPR_tok}{$matched_IPR_tok}{$SITE_presence_tok}})
							{
								push(@longest_CCDS_SITE_max_percentage,$max_percentage_SITE_tok);
							foreach my $number_of_100_damage_SITE_tok(sort keys%{$longest_hash{$longest_tmp[0]}{$longest_CCDS_tok}{$longest_CCDS_percentage_tok}{$DOMAIN_percentage_tok}{$SITE_percentage_tok}{$IPR_presence_tok}{$max_percentage_IPR_tok}{$number_of_100_damage_IPR_tok}{$matched_IPR_tok}{$SITE_presence_tok}{$max_percentage_SITE_tok}})
							{
								push(@longest_CCDS_SITE_number_max,$number_of_100_damage_SITE_tok);
							foreach my $matched_SITE_tok(sort keys%{$longest_hash{$longest_tmp[0]}{$longest_CCDS_tok}{$longest_CCDS_percentage_tok}{$DOMAIN_percentage_tok}{$SITE_percentage_tok}{$IPR_presence_tok}{$max_percentage_IPR_tok}{$number_of_100_damage_IPR_tok}{$matched_IPR_tok}{$SITE_presence_tok}{$max_percentage_SITE_tok}{$number_of_100_damage_SITE_tok}})
							{
								push(@longest_CCDS_SITE_matched,$matched_SITE_tok);
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
			}
			# my $string=join(";",@longest_CCDS);
			# my $string2=join(";",@longest_CCDS_Percentage_def);
			# my $string3=join(";",@longest_CCDS_DOMAIN_percentage);
			# my $string4=join(";",@longest_CCDS_SITE_percentage);
			# my $string5=join(";",@longest_CCDS_SITE_matched);
			
			my $string=$longest_CCDS[0];
			my $string2=$longest_CCDS_Percentage_def[0];
			my $string3=$longest_CCDS_DOMAIN_percentage[0];
			my $string4=$longest_CCDS_SITE_percentage[0];
			
			my $string6=$longest_CCDS_IPR_presence[0]; my $string7=$longest_CCDS_IPR_max_percentage[0];
			my $string8=$longest_CCDS_IPR_number_max[0]; my $string9=$longest_CCDS_IPR_matched[0];
			
			my $string10=$longest_CCDS_SITE_presence[0]; my $string11=$longest_CCDS_SITE_max_percentage[0];
			my $string12=$longest_CCDS_SITE_number_max[0]; my $string13=$longest_CCDS_SITE_matched[0];
			
			$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{'LONGEST'}{$string}{$longest_tmp[0]}{$string2}{$string3}{$string4}{$string6}{$string7}{$string8}{$string9}{$string10}{$string11}{$string12}{$string13}=1;
			#print "7_LONGEST:\t$CHROM_tok\t$SYMBOL_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$Effect_tok_definitive\t'LONGEST'\t$string\t$longest_tmp[0]\t$string2\t$string3\t$string4\t$string6\t$string7\t$string8\t$string9\t$string10\t$string11\t$string12\t$string13\n";

##########################################################  LONGEST->APPRIS  ##############################################################################

			my @Appris_CCDS=();
			my @Appris_CCDS_Percentage_def=();
			my @Appris_CCDS_DOMAIN_percentage=();
			my @Appris_CCDS_SITE_percentage=();
			
			my @Appris_CCDS_IPR_presence=();my @Appris_CCDS_IPR_max_percentage=();
			my @Appris_CCDS_IPR_number_max=();my @Appris_CCDS_IPR_matched=();

			my @Appris_CCDS_SITE_presence=();my @Appris_CCDS_SITE_max_percentage=();
			my @Appris_CCDS_SITE_number_max=();my @Appris_CCDS_SITE_matched=();
			
			foreach my $Appris_CCDS_tok(sort keys %Appris_isoform_hash)
			{
				# print "IMPRESIONES_DESGLOSADAS:$Appris_CCDS_tok\n";
				push(@Appris_CCDS,$Appris_CCDS_tok);
				foreach my $Appris_CCDS_percentage_tok(sort keys%{$Appris_isoform_hash{$Appris_CCDS_tok}})
				{
					# print "1_$Appris_CCDS_percentage_tok\n";
					push(@Appris_CCDS_Percentage_def,$Appris_CCDS_percentage_tok);
					foreach my $DOMAIN_percentage_tok(sort keys%{$Appris_isoform_hash{$Appris_CCDS_tok}{$Appris_CCDS_percentage_tok}})
					{
						# print "2_$DOMAIN_percentage_tok\n";
						push(@Appris_CCDS_DOMAIN_percentage,$DOMAIN_percentage_tok);
						foreach my $SITE_percentage_tok(sort keys%{$Appris_isoform_hash{$Appris_CCDS_tok}{$Appris_CCDS_percentage_tok}{$DOMAIN_percentage_tok}})
						{
							# print "3_$SITE_percentage_tok\n";
							push(@Appris_CCDS_SITE_percentage,$SITE_percentage_tok);	
							foreach my $IPR_presence_tok(sort keys%{$Appris_isoform_hash{$Appris_CCDS_tok}{$Appris_CCDS_percentage_tok}{$DOMAIN_percentage_tok}{$SITE_percentage_tok}})
							{
								# print "4_$IPR_presence_tok\n";
								push(@Appris_CCDS_IPR_presence,$IPR_presence_tok);
							foreach my $max_percentage_IPR_tok(sort keys%{$Appris_isoform_hash{$Appris_CCDS_tok}{$Appris_CCDS_percentage_tok}{$DOMAIN_percentage_tok}{$SITE_percentage_tok}{$IPR_presence_tok}})
							{
								# print "5_$max_percentage_IPR_tok\n";
								push(@Appris_CCDS_IPR_max_percentage,$max_percentage_IPR_tok);
							foreach my $number_of_100_damage_IPR_tok(sort keys%{$Appris_isoform_hash{$Appris_CCDS_tok}{$Appris_CCDS_percentage_tok}{$DOMAIN_percentage_tok}{$SITE_percentage_tok}{$IPR_presence_tok}{$max_percentage_IPR_tok}})
							{
								# print "6_$number_of_100_damage_IPR_tok\n";
								push(@Appris_CCDS_IPR_number_max,$number_of_100_damage_IPR_tok);
							foreach my $matched_IPR_tok(sort keys%{$Appris_isoform_hash{$Appris_CCDS_tok}{$Appris_CCDS_percentage_tok}{$DOMAIN_percentage_tok}{$SITE_percentage_tok}{$IPR_presence_tok}{$max_percentage_IPR_tok}{$number_of_100_damage_IPR_tok}})
							{
								# print "7_$matched_IPR_tok\n";
								push(@Appris_CCDS_IPR_matched,$matched_IPR_tok);
							foreach my $SITE_presence_tok(sort keys%{$Appris_isoform_hash{$Appris_CCDS_tok}{$Appris_CCDS_percentage_tok}{$DOMAIN_percentage_tok}{$SITE_percentage_tok}{$IPR_presence_tok}{$max_percentage_IPR_tok}{$number_of_100_damage_IPR_tok}{$matched_IPR_tok}})
							{	
								# print "8_$SITE_presence_tok\n";
								push(@Appris_CCDS_SITE_presence,$SITE_presence_tok);
							foreach my $max_percentage_SITE_tok(sort keys%{$Appris_isoform_hash{$Appris_CCDS_tok}{$Appris_CCDS_percentage_tok}{$DOMAIN_percentage_tok}{$SITE_percentage_tok}{$IPR_presence_tok}{$max_percentage_IPR_tok}{$number_of_100_damage_IPR_tok}{$matched_IPR_tok}{$SITE_presence_tok}})
							{
								# print "9_$max_percentage_SITE_tok\n";
								push(@Appris_CCDS_SITE_max_percentage,$max_percentage_SITE_tok);
							foreach my $number_of_100_damage_SITE_tok(sort keys%{$Appris_isoform_hash{$Appris_CCDS_tok}{$Appris_CCDS_percentage_tok}{$DOMAIN_percentage_tok}{$SITE_percentage_tok}{$IPR_presence_tok}{$max_percentage_IPR_tok}{$number_of_100_damage_IPR_tok}{$matched_IPR_tok}{$SITE_presence_tok}{$max_percentage_SITE_tok}})
							{
								# print "10_$number_of_100_damage_SITE_tok\n";
								push(@Appris_CCDS_SITE_number_max,$number_of_100_damage_SITE_tok);
							foreach my $matched_SITE_tok(sort keys%{$Appris_isoform_hash{$Appris_CCDS_tok}{$Appris_CCDS_percentage_tok}{$DOMAIN_percentage_tok}{$SITE_percentage_tok}{$IPR_presence_tok}{$max_percentage_IPR_tok}{$number_of_100_damage_IPR_tok}{$matched_IPR_tok}{$SITE_presence_tok}{$max_percentage_SITE_tok}{$number_of_100_damage_SITE_tok}})
							{
								# print "11_$matched_SITE_tok\n";
								push(@Appris_CCDS_SITE_matched,$matched_SITE_tok);
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
			}
			
			if($Appris==1)
			{
				my $string_Appris=$Appris_CCDS[0];
				#~ print "COTTON:@Appris_CCDS**\n";
				#~ print "COTTON:$Appris_CCDS[0]**\n";
				my $string_Appris2=$Appris_CCDS_Percentage_def[0];
				my $string_Appris3=$Appris_CCDS_DOMAIN_percentage[0];
				my $string_Appris4=$Appris_CCDS_SITE_percentage[0];
				my $string_Appris6=$Appris_CCDS_IPR_presence[0];
				my $string_Appris7=$Appris_CCDS_IPR_max_percentage[0];
				my $string_Appris8=$Appris_CCDS_IPR_number_max[0];
				my $string_Appris9=$Appris_CCDS_IPR_matched[0];
				my $string_Appris10=$Appris_CCDS_SITE_presence[0];  my $string_Appris11=$Appris_CCDS_SITE_max_percentage[0];
				my $string_Appris12=$Appris_CCDS_SITE_number_max[0];  my $string_Appris13=$Appris_CCDS_SITE_matched[0];
				
				
				
				$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok_definitive}{'APPRIS'}{$string_Appris}{$string_Appris2}{$string_Appris3}{$string_Appris4}{$string_Appris6}{$string_Appris7}{$string_Appris8}{$string_Appris9}{$string_Appris10}{$string_Appris11}{$string_Appris12}{$string_Appris13}=1;
				#~ print "8_APPRIS:\t$CHROM_tok\t$SYMBOL_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$Effect_tok_definitive\t'APPRIS'\t$string_Appris\t$string_Appris2\t$string_Appris3\t$string_Appris4\t$string_Appris6\t$string_Appris7\t$string_Appris8\t$string_Appris9\t$string_Appris10\t$string_Appris11\t$string_Appris12\t$string_Appris13\n";
			}
		}
	}#
	}
	}
}
}

#~ exit;

$time='['. timestamp(). ']'."\n";
print "Start PRINTING:$time\n";

if(open(OUTPUT, '>'.$output))
{
	print OUTPUT "baseNCBI37\tAlleles_REF>ALT\tGeneNameHGNC\tNumIsoformsInQueryGene\tratioIsoformsBearingTheVariant\t".
				"ratioAffectedIsoforms_stop-gained\tratioAffectedIsoforms_frameshift\tratioAffectedIsoforms_splice\tratioAffectedIsoforms_coding-synonymous\tratioAffectedIsoforms_missense\t".
				"PEJMAN1\tPEJMAN2\tratioAffectedIsoformsTargetedbyNMD\tratioAffectedIsoformsTargetedby_derived_NMD\t".
				"IsPrincipalIsoformAffected\tIsWithinLongestCCDS\tIsWithinPervasiveIsoform\t".
				"LongestCCDSLength\tPercentagePrincipalOrLongestCCDSAffected\t".
				"PEJMAN3\tPEJMAN4\tDomainINFOAvailable\tPercentageOfDomainPositionsAffected\tmaxPercDomainAffected\tNumberOfDomains100Damage\tDomainMatched\t".
				"SiteINFOAvailable\tPercentageOfSitePositionsAffected\tmaxPercSiteAffected\tNumberOfSites100Damage\tSiteMatched\n";
	
	print OUTPUT "#X:74334588\tC>T\tABCB7\t6\t4\t".
				"0.5\t0\t0\t0\t0.5\t".
				"PEJMAN1\tPEJMAN2\t0.5\t0\t".
				"1\t0\t0\t".
				"2324\t100\t".
				"PEJMAN3\tPEJMAN4\t1\t100\t100\t1\t1\t".
				"1\t100\t100\t1\t0\n";
	
	foreach my $CHROM_tok(sort {$a<=>$b}keys%RESULTS_hash)
	{
	foreach my $SYMBOL_tok(sort keys%{$RESULTS_hash{$CHROM_tok}})
	{
	foreach my $POS_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}})
	{
	foreach my $REF_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}})
	{
	foreach my $ALT_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}})
	{
	foreach my $Effect_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}})
	{
		
		if ($Effect_tok=~/stop_gained/||$Effect_tok=~/frameshift_variant/||$Effect_tok=~/splice_donor_variant/||$Effect_tok=~/splice_acceptor_variant/ || $Effect_tok=~/inframe_deletion/ ||$Effect_tok=~/disruptive_inframe_deletion/)
		{
			print OUTPUT "$CHROM_tok:$POS_tok\t$REF_tok>$ALT_tok\t$SYMBOL_tok\t";##############################################################
			#~ print "$CHROM_tok:$POS_tok\t$REF_tok>$ALT_tok\t$SYMBOL_tok\t";##############################################################

			foreach my $total_isoforms_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'ISOFORMS_SUMMARY'}})
			{
			foreach my $total_affected_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'ISOFORMS_SUMMARY'}{$total_isoforms_tok}})
			{
			foreach my $ratio_stop_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'ISOFORMS_SUMMARY'}{$total_isoforms_tok}{$total_affected_tok}})
			{
			foreach my $ratio_frameshift_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'ISOFORMS_SUMMARY'}{$total_isoforms_tok}{$total_affected_tok}{$ratio_stop_tok}})
			{
			foreach my $ratio_splice_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'ISOFORMS_SUMMARY'}{$total_isoforms_tok}{$total_affected_tok}{$ratio_stop_tok}{$ratio_frameshift_tok}})
			{
			#foreach my $ratio_inframe_del_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'ISOFORMS_SUMMARY'}{$total_isoforms_tok}{$total_affected_tok}{$ratio_stop_tok}{$ratio_frameshift_tok}{$ratio_splice_tok}})
			#{
			foreach my $ratio_syn_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'ISOFORMS_SUMMARY'}{$total_isoforms_tok}{$total_affected_tok}{$ratio_stop_tok}{$ratio_frameshift_tok}{$ratio_splice_tok}})
			{
			foreach my $ratio_non_syn_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'ISOFORMS_SUMMARY'}{$total_isoforms_tok}{$total_affected_tok}{$ratio_stop_tok}{$ratio_frameshift_tok}{$ratio_splice_tok}{$ratio_syn_tok}})
			{
				#print "$Effect_tok\n";
				print OUTPUT "$total_isoforms_tok\t$total_affected_tok\t$ratio_stop_tok\t$ratio_frameshift_tok\t$ratio_splice_tok\t$ratio_syn_tok\t$ratio_non_syn_tok\t";
				#~ print "Ratios$total_isoforms_tok\t$total_affected_tok\t$ratio_stop_tok\t$ratio_frameshift_tok\t$ratio_splice_tok\t$ratio_syn_tok\t$ratio_non_syn_tok\t";
			}	
			}	
			#}
			}	
			}	
			}	
			}	
			}
			
			foreach my $ratio_NMD_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'NMD'}})
			{
				print OUTPUT "PEJMAN1\tPEJMAN2\t$ratio_NMD_tok\t";
				#~ print "NMD:$ratio_NMD_tok\t";
				
			}

			foreach my $ratio_NMD_derived_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'NMD_derived'}})
			{
				print OUTPUT "$ratio_NMD_derived_tok\t";
				#~ print "NMD_derived:$ratio_NMD_derived_tok\t";
			}
			#############################################################################################################################################
			my $match_Appris="NaN";
			my $Flag=0;
			foreach my $CCDS_tok(sort keys %{$hash_APPRIS{$SYMBOL_tok}})
			{
				if(exists($hash_APPRIS{$SYMBOL_tok}{$CCDS_tok}))
				{
					$Flag++;
					#~ if($SYMBOL_tok eq 'UGT2A1')
					#~ {print "INICIO_ARRAY_Appris:$SYMBOL_tok\t$CCDS_tok\tFlag:$Flag\n";}
				}
			}
			if($Flag >= 1)
			{
					my @Appris_tmp=sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'APPRIS_PRESENCE'}};
					my $positive=1;
					if (grep /$positive/, @Appris_tmp) 
					{
						$match_Appris=1;
						print OUTPUT "$match_Appris\t";
						#~ print "Principal:$match_Appris\t";
					}
					else
					{
						$match_Appris=0;
						print OUTPUT "$match_Appris\t";
						#~ print "Principal:$match_Appris\t";
					}
			}
			elsif($Flag == 0)
			{
				$match_Appris="NaN";
				print OUTPUT "$match_Appris\t";
				#~ print "Principal:$match_Appris\t";
			}
			
			
			
			
			my @absolute_longest_tmp=sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'ABSOLUTE_LONGEST'}};
			my $match=1;
			if (grep /$match/, @absolute_longest_tmp) 
			{
				print OUTPUT "$match\t";
				#~ print "Longest:$match\t";
			}
			else
			{
				$match=0;
				print OUTPUT "$match\t";
				#~ print "Longest:$match\t";
			}
			
			my $match_Pervasive="NaN";
			my $Flag_Pervasive=0;
			foreach my $CCDS_tok(sort keys %{$Pervasive_hash{$SYMBOL_tok}})
			{
				if(exists($Pervasive_hash{$SYMBOL_tok}{$CCDS_tok}))
				{
					#print "INICIO_ARRAY_Pervasive:$SYMBOL_tok\n";
					$Flag_Pervasive++;
				}
			}
			if($Flag_Pervasive == 1)
			{
					my @Pervasive_tmp=sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'PERVASIVE'}};
					my $positive=1;
					if (grep /$positive/, @Pervasive_tmp) 
					{
						$match_Pervasive=1;
						print OUTPUT "$match_Pervasive\t";
						#~ print "Pervasive:$match_Pervasive\t";
					}
					else
					{
						$match_Pervasive=0;
						print OUTPUT "$match_Pervasive\t";
						#~ print "Pervasive:$match_Pervasive\t";
					}
			}
			elsif($Flag_Pervasive == 0)
			{
				$match_Pervasive="NaN";
				print OUTPUT "$match_Pervasive\t";
				#~ print "Pervasive:$match_Pervasive\t";
			}
			#~ print "EL MATCH APPRIS ES:$SYMBOL_tok\t$match_Appris\n";
			#################################################################################################################################################
			if($match_Appris != 1)
			{
				foreach my $CCDS_LONGEST(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'LONGEST'}})
				{
				foreach my $CCDS_LONGEST_length_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'LONGEST'}{$CCDS_LONGEST}})
				{
				foreach my $LONGEST_seq_percentage_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'LONGEST'}{$CCDS_LONGEST}{$CCDS_LONGEST_length_tok}})
				{
				foreach my $LONGEST_DOMAIN_percentage_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'LONGEST'}{$CCDS_LONGEST}{$CCDS_LONGEST_length_tok}{$LONGEST_seq_percentage_tok}})
				{
				foreach my $LONGEST_SITE_percentage_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'LONGEST'}{$CCDS_LONGEST}{$CCDS_LONGEST_length_tok}{$LONGEST_seq_percentage_tok}{$LONGEST_DOMAIN_percentage_tok}})
				{
				foreach my $IPR_presence_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'LONGEST'}{$CCDS_LONGEST}{$CCDS_LONGEST_length_tok}{$LONGEST_seq_percentage_tok}{$LONGEST_DOMAIN_percentage_tok}{$LONGEST_SITE_percentage_tok}})
				{
				foreach my $IPR_max_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'LONGEST'}{$CCDS_LONGEST}{$CCDS_LONGEST_length_tok}{$LONGEST_seq_percentage_tok}{$LONGEST_DOMAIN_percentage_tok}{$LONGEST_SITE_percentage_tok}{$IPR_presence_tok}})
				{
				foreach my $IPR_number_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'LONGEST'}{$CCDS_LONGEST}{$CCDS_LONGEST_length_tok}{$LONGEST_seq_percentage_tok}{$LONGEST_DOMAIN_percentage_tok}{$LONGEST_SITE_percentage_tok}{$IPR_presence_tok}{$IPR_max_tok}})
				{
				foreach my $IPR_matched_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'LONGEST'}{$CCDS_LONGEST}{$CCDS_LONGEST_length_tok}{$LONGEST_seq_percentage_tok}{$LONGEST_DOMAIN_percentage_tok}{$LONGEST_SITE_percentage_tok}{$IPR_presence_tok}{$IPR_max_tok}{$IPR_number_tok}})
				{
				foreach my $SITE_presence_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'LONGEST'}{$CCDS_LONGEST}{$CCDS_LONGEST_length_tok}{$LONGEST_seq_percentage_tok}{$LONGEST_DOMAIN_percentage_tok}{$LONGEST_SITE_percentage_tok}{$IPR_presence_tok}{$IPR_max_tok}{$IPR_number_tok}{$IPR_matched_tok}})
				{
				foreach my $SITE_max_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'LONGEST'}{$CCDS_LONGEST}{$CCDS_LONGEST_length_tok}{$LONGEST_seq_percentage_tok}{$LONGEST_DOMAIN_percentage_tok}{$LONGEST_SITE_percentage_tok}{$IPR_presence_tok}{$IPR_max_tok}{$IPR_number_tok}{$IPR_matched_tok}{$SITE_presence_tok}})
				{
				foreach my $SITE_number_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'LONGEST'}{$CCDS_LONGEST}{$CCDS_LONGEST_length_tok}{$LONGEST_seq_percentage_tok}{$LONGEST_DOMAIN_percentage_tok}{$LONGEST_SITE_percentage_tok}{$IPR_presence_tok}{$IPR_max_tok}{$IPR_number_tok}{$IPR_matched_tok}{$SITE_presence_tok}{$SITE_max_tok}})
				{
				foreach my $SITE_matched_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'LONGEST'}{$CCDS_LONGEST}{$CCDS_LONGEST_length_tok}{$LONGEST_seq_percentage_tok}{$LONGEST_DOMAIN_percentage_tok}{$LONGEST_SITE_percentage_tok}{$IPR_presence_tok}{$IPR_max_tok}{$IPR_number_tok}{$IPR_matched_tok}{$SITE_presence_tok}{$SITE_max_tok}{$SITE_number_tok}})
				{
					print OUTPUT "PEJMAN3\tPEJMAN4\t$CCDS_LONGEST_length_tok\t$LONGEST_seq_percentage_tok\t$IPR_presence_tok\t$LONGEST_DOMAIN_percentage_tok\t$IPR_max_tok\t$IPR_number_tok\t$IPR_matched_tok\t$SITE_presence_tok\t$LONGEST_SITE_percentage_tok\t$SITE_max_tok\t$SITE_number_tok\t$SITE_matched_tok\t";
					#~ print "Domain/site:$CCDS_LONGEST_length_tok\t$LONGEST_seq_percentage_tok\t$IPR_presence_tok\t$LONGEST_DOMAIN_percentage_tok\t$IPR_max_tok\t$IPR_number_tok\t$IPR_matched_tok\t$SITE_presence_tok\t$LONGEST_SITE_percentage_tok\t$SITE_max_tok\t$SITE_number_tok\t$SITE_matched_tok**\t";
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
				}
				}
			}
			elsif($match_Appris == 1)
			{
				foreach my $CCDS_LONGEST(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'LONGEST'}})
				{
				foreach my $CCDS_LONGEST_length_tok(sort keys%{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'LONGEST'}{$CCDS_LONGEST}})
				{
					foreach my $Appris_CCDS(sort keys %{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'APPRIS'}})
					{
					foreach my $Appris_CCDS_Percentage(sort keys %{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'APPRIS'}{$Appris_CCDS}})
					{
					foreach my $Appris_CCDS_DOMAIN_Percentage(sort keys %{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'APPRIS'}{$Appris_CCDS}{$Appris_CCDS_Percentage}})
					{
					foreach my $Appris_CCDS_SITE_Percentage(sort keys %{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'APPRIS'}{$Appris_CCDS}{$Appris_CCDS_Percentage}{$Appris_CCDS_DOMAIN_Percentage}})
					{
					foreach my $Appris_CCDS_IPR_presence(sort keys %{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'APPRIS'}{$Appris_CCDS}{$Appris_CCDS_Percentage}{$Appris_CCDS_DOMAIN_Percentage}{$Appris_CCDS_SITE_Percentage}})
					{
					foreach my $Appris_CCDS_IPR_max_tok(sort keys %{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'APPRIS'}{$Appris_CCDS}{$Appris_CCDS_Percentage}{$Appris_CCDS_DOMAIN_Percentage}{$Appris_CCDS_SITE_Percentage}{$Appris_CCDS_IPR_presence}})
					{
					foreach my $Appris_CCDS_IPR_number_tok(sort keys %{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'APPRIS'}{$Appris_CCDS}{$Appris_CCDS_Percentage}{$Appris_CCDS_DOMAIN_Percentage}{$Appris_CCDS_SITE_Percentage}{$Appris_CCDS_IPR_presence}{$Appris_CCDS_IPR_max_tok}})
					{
					foreach my $Appris_CCDS_IPR_matched_tok(sort keys %{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'APPRIS'}{$Appris_CCDS}{$Appris_CCDS_Percentage}{$Appris_CCDS_DOMAIN_Percentage}{$Appris_CCDS_SITE_Percentage}{$Appris_CCDS_IPR_presence}{$Appris_CCDS_IPR_max_tok}{$Appris_CCDS_IPR_number_tok}})
					{
					foreach my $Appris_CCDS_SITE_presence(sort keys %{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'APPRIS'}{$Appris_CCDS}{$Appris_CCDS_Percentage}{$Appris_CCDS_DOMAIN_Percentage}{$Appris_CCDS_SITE_Percentage}{$Appris_CCDS_IPR_presence}{$Appris_CCDS_IPR_max_tok}{$Appris_CCDS_IPR_number_tok}{$Appris_CCDS_IPR_matched_tok}})
					{
					foreach my $Appris_CCDS_SITE_max(sort keys %{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'APPRIS'}{$Appris_CCDS}{$Appris_CCDS_Percentage}{$Appris_CCDS_DOMAIN_Percentage}{$Appris_CCDS_SITE_Percentage}{$Appris_CCDS_IPR_presence}{$Appris_CCDS_IPR_max_tok}{$Appris_CCDS_IPR_number_tok}{$Appris_CCDS_IPR_matched_tok}{$Appris_CCDS_SITE_presence}})
					{
					foreach my $Appris_CCDS_SITE_number(sort keys %{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'APPRIS'}{$Appris_CCDS}{$Appris_CCDS_Percentage}{$Appris_CCDS_DOMAIN_Percentage}{$Appris_CCDS_SITE_Percentage}{$Appris_CCDS_IPR_presence}{$Appris_CCDS_IPR_max_tok}{$Appris_CCDS_IPR_number_tok}{$Appris_CCDS_IPR_matched_tok}{$Appris_CCDS_SITE_presence}{$Appris_CCDS_SITE_max}})
					{
					foreach my $Appris_CCDS_SITE_matched(sort keys %{$RESULTS_hash{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{'APPRIS'}{$Appris_CCDS}{$Appris_CCDS_Percentage}{$Appris_CCDS_DOMAIN_Percentage}{$Appris_CCDS_SITE_Percentage}{$Appris_CCDS_IPR_presence}{$Appris_CCDS_IPR_max_tok}{$Appris_CCDS_IPR_number_tok}{$Appris_CCDS_IPR_matched_tok}{$Appris_CCDS_SITE_presence}{$Appris_CCDS_SITE_max}{$Appris_CCDS_SITE_number}})
					{
						
						print OUTPUT "PEJMAN3\tPEJMAN4\t$CCDS_LONGEST_length_tok\t$Appris_CCDS_Percentage\t$Appris_CCDS_IPR_presence\t$Appris_CCDS_DOMAIN_Percentage\t$Appris_CCDS_IPR_max_tok\t$Appris_CCDS_IPR_number_tok\t$Appris_CCDS_IPR_matched_tok\t$Appris_CCDS_SITE_presence\t$Appris_CCDS_SITE_Percentage\t$Appris_CCDS_SITE_max\t$Appris_CCDS_SITE_number\t$Appris_CCDS_SITE_matched\t";
						#~ print "Domain/site:$CCDS_LONGEST_length_tok\t$Appris_CCDS_Percentage\t$Appris_CCDS_IPR_presence\t$Appris_CCDS_DOMAIN_Percentage\t$Appris_CCDS_IPR_max_tok\t$Appris_CCDS_IPR_number_tok\t$Appris_CCDS_IPR_matched_tok\t$Appris_CCDS_SITE_presence\t$Appris_CCDS_SITE_Percentage\t$Appris_CCDS_SITE_max\t$Appris_CCDS_SITE_number\t$Appris_CCDS_SITE_matched**\t";
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
					}
				}
				}
			}
			print OUTPUT "\n";
			#~ print "\n";
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
