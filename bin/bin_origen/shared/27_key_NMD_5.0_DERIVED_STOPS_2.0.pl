##	Script to calculate if the frameshift derived first stop-codon entails NMD. 2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne) 

use strict;
use warnings;
use Time::localtime;

my %hash2=();
my %hash3=();
my %hash_splice=();
my %hash_derived_PTC=();
my %hash_EXON=();
my %RESULT_hash=();
my %RESULT_hash_derived_PTC=();

my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $input3=$ARGV[2];
my $input4=$ARGV[3];
my $output=$ARGV[4];


my $time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

if (open (INPUT1, $input1))

{
	#INPUT1= 1GK_frameshift_derived_PTCs.txt
	
	#  X       GRIA3   122336600       T       TG      frameshift_variant&feature_elongation   ENST00000371264 ->DERIVED_PTC:122338313
#  X       GRIA3   122336600       T       TG      frameshift_variant&feature_elongation   ENST00000371266 ->DERIVED_PTC:122338313
#  X       IL9R    155231142       GT      G       frameshift_variant&feature_truncation   ENST00000369423 ->DERIVED_PTC:155232595
#  X       IL9R    155231142       GT      G       frameshift_variant&feature_truncation   ENST00000540897 ->DERIVED_PTC:155232595
#  X       LCA10   153149707       C       CG      frameshift_variant&feature_elongation   ENST00000357566 ->DERIVED_PTC:153149717
#  X       LCA10   153149708       G       GC      frameshift_variant&feature_elongation   ENST00000357566 ->DERIVED_PTC:153149717
#  X       LCA10   153151280       G       GCC     frameshift_variant&feature_elongation   ENST00000357566 ->DERIVED_PTC:153152479
#  X       LCA10   153151284       GT      G       frameshift_variant&feature_truncation   ENST00000357566 ->DERIVED_PTC:153152479

	
while(my $line=<INPUT1>)
	{
		chomp $line;
		#print "$line\n";
		#print "La línea es:$line\n";
		if($line=~/>DERIVED_PTC:/)
		{
			#print "La línea es:$line\n";
			my @tmp=split("\t->DERIVED_PTC:",$line);
			my $parental=$tmp[0];
			my $derived=$tmp[1];
			my ($CHROM,$SYMBOL,$POS,$REF,$ALT,$Effect,$ENST,$POS_derived)="NaN"x8;
			if($parental=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
			{
				$CHROM=$1;
				$SYMBOL=$2;
				$POS=$3;
				$REF=$4;
				$ALT=$5;
				$Effect=$6;
				$ENST=$7;
				
			}
			if($derived=~/(.+)/)
			{
				$POS_derived=$1;
			}
			
			$hash_derived_PTC{$CHROM}{$SYMBOL}{$ENST}{$POS}{$REF}{$ALT}{$Effect}{$POS_derived}=1;
			#~ print "$CHROM\t$SYMBOL\t$ENST\t$POS\t$REF\t$ALT\t$Effect\t$POS_derived\n";
		}
	}
}else{print "Unable to open INPUT1";}

$time='['. timestamp(). ']'."\n";
print "Start charging hash2:$time\n";

if (open (INPUT2, $input2))
{
	#INPUT1=XXX_out_vep_parsed.vcf
	
	#	13      32906565        .       CA      C       33.73   .       frameshift_variant&feature_truncation;BRCA2;ENST00000380152;protein_coding;;0.5;;;;0.5;;;;NaN;NaN;NaN;NaN;PIRSF_domain:PIRSF002397;;
	#	13      32906565        .       CA      C       33.73   .       frameshift_variant&feature_truncation;BRCA2;ENST00000530893;protein_coding;;0.5;;;;0.5;;;;NaN;NaN;NaN;NaN;NaN;;
	#	13      32906565        .       CA      C       33.73   .       frameshift_variant&feature_truncation;BRCA2;ENST00000544455;protein_coding;;0.5;;;;0.5;;;;NaN;NaN;NaN;NaN;PIRSF_domain:PIRSF002397;;
	#	13      32906576        .       CA      C       36.73   .       frameshift_variant&feature_truncation;BRCA2;ENST00000380152;protein_co
	
	#INPUT1b=XXX-eff_parsed.vcf
	
	#13      32906565        .       CA      C       33.73   .       frameshift_variant;BRCA2;ENST00000380152;protein_coding;T317;10;;CODING;aca/;3418;HIGH;1;0.5;NaN;NaN    LOF_positive_snpEff;BRCA2;ENSG00000139618;6;0.50        
	#13      32906565        .       CA      C       33.73   .       frameshift_variant;BRCA2;ENST00000530893;protein_coding;T194;10;;CODING;aca/;480;HIGH;1;0.5;WARNING_TRANSCRIPT_INCOMPLETE;NaN   LOF_positive_snpEff;BRCA2;ENSG00000139618;6;0.50        
	#13      32906565        .       CA      C       33.73   .       frameshift_variant;BRCA2;ENST00000544455;protein_coding;T317;10;;CODING;aca/;3418;HIGH;1;0.5;NaN;NaN    LOF_positive_snpEff;BRCA2;ENSG00000139618;6;0.50        
	#13      32906576        .       CA      C       36.73   .       frameshift_variant;BRCA2;ENST00000380152;protein_coding;Q321;10;;CODING;caa/;3418;HIGH;1;0.5;NaN;NaN    LOF_positive_snpEff;BRCA2;ENSG00000139618;6;0.50        

while(my $line=<INPUT2>)
	{
		chomp $line;
		#print "$line\n";
		if($line=~/^([^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)\t[^\t]+\t[^\t]+\t(.+)/)
		{
			my $CHROM=$1;
			my $POS=$2;
			my $REF=$3;
			my $ALT=$4;
			my $fields=$5;
			#print "hello_world:$CHROM\t$POS\t$fields\t$REF\t$ALT\n";
			my @fields_tmp=split(";",$fields);
			my $Effect=$fields_tmp[0];
			my $SYMBOL=$fields_tmp[1];
			my $ENST=$fields_tmp[2];
			unless (exists($hash_derived_PTC{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{$ENST}))
			{
				my $NMD_derived="NaN";
				$RESULT_hash_derived_PTC{$CHROM}{$SYMBOL}{$POS}{$REF}{$ALT}{$Effect}{$ENST}{$POS}{$NMD_derived}=1;
				#~ print "RESULT_hash_derived_PTC:$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t$Effect\t$ENST\t$POS\t$NMD_derived\n";
			}
			else
			{
				#~ print "$CHROM\t$SYMBOL\t$ENST\t$Effect\t$POS\t$REF\t$ALT\n";
			}
			
		}
	}
}else{print "Unable to open INPUT1\n";}

$time='['. timestamp(). ']'."\n";
print "Start charging hash3:$time\n";

if (open(INPUT3, $input3))
{
	## Input file= NMD_table.txt

# X       +       ABCD1   ENSG00000101986 ENST00000218104 10      [153008943]-[153010216] [152990323]-[152991621] [152994687]-[152994867] [153001566]-[153001708] [153001799]-[153001967] [153002611]-[153002705] [153005546]-[153005691] [153006028]-[153006173] [153008441]-[153008525] [153008675]-[153008800] 
# X       +       ABCD1   ENSG00000101986 ENST00000370129 2       [152991272]-[152991621] [152994687]-[152995352] 
# X       +       ABCD1   ENSG00000101986 ENST00000443684 6       [153000798]-[153000881] [153001566]-[153001708] [153001799]-[153001967] [153002611]-[153002705] [153005546]-[153005691] [153006028]-[153006058] 
# X       +       AC002365.1      ENSG00000234469 ENST00000445307 2       [9935392]-[9936005]     [9936107]-[9936134]     
# X       +       ACRC    ENSG00000147174 ENST00000373695 12      [70830531]-[70830669]   [70832205]-[70832409]   [70832712]-[70833432]
#   [70800135]-[70800730]   [70811972]-[70812018]   [70812387]-[70812420]   [70814181]-[70814227]   [70814596]-[70814629]   [70817800]-[70817888]   [70823438]-[70824526]   [70825513]-[70825579]   [70828823]-[70828967]   
# X       +       ACRC    ENSG00000147174 ENST00000373696 13      [70828823]-[70828967]   [70830531]-[70830669]   [70832205]-[70832409]
#   [70832712]-[70833433]   [70798261]-[70798372]   [70800670]-[70800730]   [70811972]-[70812018]   [70812387]-[70812420]   [70814181]-[70814227]   [70814596]-[70814629]   [70817800]-[70817888]   [70823438]-[70824526]   [70825513]-[70825579]   
# X       +       AFF2    ENSG00000155966 ENST00000286437 18      [148044245]-[148044467] [148048320]-[148048609] [148049159]-[148049222] [148055001]-[148055137] [148059463]-[148059534] [148059892]-[148059985] [148062268]-[148062320] [148068897]-[148069087] [148072741]-[148082193] [147800615]-[147800755] [147891400]-[147891444] [147924490]-[147924526] [147924906]-[147924957] [147967419]-[147967515] [147985751]-[147985788] [148035110]-[148035269] [148037133]-[148038143] [148039867]-[148039988] 
# X       +       AFF2    ENSG00000155966 ENST00000342251 20      [148037133]-[148038143] [148039867]-[148039988] [148044245]-[148044467] [148048320]-[148048609] [148049159]-[148049222] [148055001]-[148055137] [148059463]-[148059534] [148059892]-[148059985] [148062268]-[148062320] [148068897]-[148069087] [147582228]-[147582664] [148072741]-[148075954] [147733520]-[147733640] [147743429]-[147744289] [147891400]-[147891444] [147924490]-[147924526] [147924906]-[147924957] [147967419]-[147967515] [147985751]-[147985788] [148035110]-[148035269] 
# X       +       AFF2    ENSG00000155966 ENST00000370457 20      [148037133]-[148038143] [148039867]-[148039988] [148044245]-[148044461] [148048320]-[148048609] [148049159]-[148049222] [148055001]-[148055137] [148059463]-[148059534] [148059892]-[148059985] [148062268]-[148062320] [148068897]-[148069087] [147582139]-[147582664] [148072741]-[148082193] [147733520]-[147733640] [147743429]-[147744289] [147891400]-[147891444] [147924490]-[147924526] [147924906]-[147924957] [147967419]-[147967515] [147985751]-[147985788] [148035110]-[148035269] 

	while(my $line=<INPUT3>)
	{
		chomp($line);
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)/)
		{
			my $CHROM=$1;
			my $strand=$2;
			my $SYMBOL=$3;
			my $ENSG=$4;
			my $ENST=$5;
			my $exon_number=$6;
			my $fields=$7;
			
			my %hash_begin=();
			my %hash_end=();
				
			my @fields_tmp=split("\t",$fields);
			foreach my $fields_tmp_tok(@fields_tmp)
			{
				if($fields_tmp_tok=~/\[(.+)\]\-\[(.+)\]/)
				{
					my $begin=$1;
					my $end=$2;
					$hash_EXON{$CHROM}{$SYMBOL}{$ENST}{$strand}{'BEGIN'}{$begin}{$end}=1;
					$hash_EXON{$CHROM}{$SYMBOL}{$ENST}{$strand}{'END'}{$end}{$begin}=1;
				}
			}
			my @EXON_END_tmp=sort {$a<=>$b} keys %{$hash_EXON{$CHROM}{$SYMBOL}{$ENST}{$strand}{'END'}};
			my @EXON_BEGIN_tmp=sort {$a<=>$b} keys %{$hash_EXON{$CHROM}{$SYMBOL}{$ENST}{$strand}{'BEGIN'}};
			
			#~ print "INICIO:$CHROM\t$SYMBOL\t$ENST\t$strand\n";
			#~ print "El array de exones begin es:@EXON_BEGIN_tmp\n";
			#~ print "El array de exones end es:@EXON_END_tmp\n";
			if(scalar(@EXON_BEGIN_tmp > 1))
			{
				if($strand eq '+')
				{
					my %hash_index=();
					my %hash_position=();
							
					my @EXON_BEGIN_tmp_reversed=reverse(@EXON_BEGIN_tmp);
							
					for(my $i=1; $i<scalar(@EXON_BEGIN_tmp_reversed); $i++)
					{
						#~ print "$EXON_BEGIN_tmp_reversed[$i]\n";
						foreach my $end_tok(sort {$a<=>$b} keys %{$hash_EXON{$CHROM}{$SYMBOL}{$ENST}{$strand}{'BEGIN'}{$EXON_BEGIN_tmp_reversed[$i]}})
						{
							#~ print "$end_tok\n";
							my $distance=$end_tok-$EXON_BEGIN_tmp_reversed[$i]+1;
							#~ print "$EXON_BEGIN_tmp_reversed[$i]\t$end_tok\t$distance\n";
							$hash_index{$EXON_BEGIN_tmp_reversed[$i]}{$distance}=1;
							
						}
					}
							
					foreach my $begin_NMD_tok(sort keys %hash_index)
					{
						foreach my $distance_tok(sort keys %{$hash_index{$begin_NMD_tok}})
						{
							for(my $i=$begin_NMD_tok; $i<$begin_NMD_tok+$distance_tok;$i++)
							{
								$hash_position{$i}=1;
								#~ print "$i\n";
							}
						}
					}
					my @positions_tmp=reverse sort{$a<=>$b}keys%hash_position;
					#~ print "El array es:@positions_tmp\n";
					if(scalar(@positions_tmp <=50))
					{
						# Genes with less than 50 pb upstream last EJC are always NMD_negative
						
						my $threshold_NMD="NMD_negative";
						$hash_EXON{$CHROM}{$SYMBOL}{$ENST}{$strand}{'THRESHOLD'}{$threshold_NMD}=1;
						#~ print "WINDOW:$CHROM\t$SYMBOL\t$ENST\t$strand\t'THRESHOLD'\t$threshold_NMD\n";
					}
					elsif(scalar(@positions_tmp > 50))
					{
						my $threshold_NMD=$positions_tmp[49];
						$hash_EXON{$CHROM}{$SYMBOL}{$ENST}{$strand}{'THRESHOLD'}{$threshold_NMD}=1;
						#~ print "WINDOW:$CHROM\t$SYMBOL\t$ENST\t$strand\t'THRESHOLD'\t$threshold_NMD\n";
					}
				}
				if($strand eq '-')
				{
					my %hash_index=();
					my %hash_position=();
							
					for(my $i=1; $i<scalar(@EXON_BEGIN_tmp); $i++)
					{
						#~ print "$EXON_BEGIN_tmp[$i]\n";
						foreach my $end_tok(sort {$a<=>$b} keys %{$hash_EXON{$CHROM}{$SYMBOL}{$ENST}{$strand}{'BEGIN'}{$EXON_BEGIN_tmp[$i]}})
						{
							#~ print "$end_tok\n";
							my $distance=$end_tok-$EXON_BEGIN_tmp[$i]+1;
							#~ print "$EXON_BEGIN_tmp[$i]\t$end_tok\t$distance\n";
							$hash_index{$EXON_BEGIN_tmp[$i]}{$distance}=1;
							
						}
					}
							
					foreach my $begin_NMD_tok(sort keys %hash_index)
					{
						foreach my $distance_tok(sort keys %{$hash_index{$begin_NMD_tok}})
						{
							for(my $i=$begin_NMD_tok; $i<$begin_NMD_tok+$distance_tok;$i++)
							{
								$hash_position{$i}=1;
								#~ print "$i\n";
							}
						}
					}
					my @positions_tmp=sort{$a<=>$b}keys%hash_position;
					#~ print "El array es:@positions_tmp\n";
					if(scalar(@positions_tmp <=50))
					{
						# Genes with less than 50 pb downstream last EJC are always NMD_negative
						
						my $threshold_NMD="NMD_negative";
						$hash_EXON{$CHROM}{$SYMBOL}{$ENST}{$strand}{'THRESHOLD'}{$threshold_NMD}=1;
						#~ print "WINDOW:$CHROM\t$SYMBOL\t$ENST\t$strand\t'THRESHOLD'\t$threshold_NMD\n";
					}
					elsif(scalar(@positions_tmp > 50))
					{
						my $threshold_NMD=$positions_tmp[49];
						$hash_EXON{$CHROM}{$SYMBOL}{$ENST}{$strand}{'THRESHOLD'}{$threshold_NMD}=1;
						#~ print "WINDOW:$CHROM\t$SYMBOL\t$ENST\t$strand\t'THRESHOLD'\t$threshold_NMD\n";
					}
					
				}
			}
			else
			{
				# NMD does not apply to monoexonic genes
				
				my $threshold_NMD="NaN";
				$hash_EXON{$CHROM}{$SYMBOL}{$ENST}{$strand}{'THRESHOLD'}{$threshold_NMD}=1;
				#~ print "WINDOW:$CHROM\t$SYMBOL\t$ENST\t$strand\t'THRESHOLD'\t$threshold_NMD\n";
			}
		}
	}
}else{print "Impossible to open INPUT3\n";}

$time='['. timestamp(). ']'."\n";
print "Start charging hash4:$time\n";

if (open(INPUT4, $input4))
{    
	while(my $line=<INPUT4>)
	{
		chomp($line);
		# Input=ENST_table_full_condensed.txt
		
		#~ 1       A3GALT2 -       ENST00000330379 [33777791-33777822]     
		#~ 1       A3GALT2 -       ENST00000330379 ENST00000442999 [33772370-33773054]     [33777653-33777790]     
		#~ 1       A3GALT2 -       ENST00000442999 [33778102-33778191]     [33778408-33778491]     [33786677-33786699]     
		#~ 1       AADACL3 +       ENST00000332530 [12776344-12776347]     
		#~ 1       AADACL3 +       ENST00000332530 ENST00000359318 [12780885-12780948]     [12785189-12785960]     
		#~ 1       AADACL3 +       ENST00000359318 [12779480-12779693]
		
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
				
				$hash2{$CHROM}{$strand}{$SYMBOL}{$ENST_tok}{$interval_tok}=1;
				#~ print "$CHROM\t$strand\t$SYMBOL\t$ENST_tok\t$interval_tok**\n";
			}	
			}
		}
	}
}else {print "impossible to open INPUT4\n";die;}



$time='['. timestamp(). ']'."\n";
print "Start processing hash3:$time\n";

foreach my $CHROM_tok(sort keys%hash2)
{
foreach my $strand_tok(sort keys %{$hash2{$CHROM_tok}})
{
foreach my $SYMBOL_tok(sort keys %{$hash2{$CHROM_tok}{$strand_tok}})
{
foreach my $ENST_tok(sort keys %{$hash2{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}})
{
foreach my $interval_tok(sort keys %{$hash2{$CHROM_tok}{$strand_tok}{$SYMBOL_tok}{$ENST_tok}})
{
	if($interval_tok=~/\[(.+)\-(.+)\]/)
	{
		my $begin=$1;
		my $end=$2;
		
		# As there might be repeated distance values total distances should be kept in an array rahter than in a hash
		$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'BEGIN'}{$begin}{$end}=1;
		$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'END'}{$end}{$begin}=1;
	}
}
}	
}	
}
}

%hash2=();

$time='['. timestamp(). ']'."\n";
print "Start PRINTING:$time\n";

if(open (OUTPUT, '>'.$output))
{
foreach my $CHROM_tok(sort keys %hash3)
{
foreach my $SYMBOL_tok(sort keys %{$hash3{$CHROM_tok}})
{	
foreach my $ENST_tok(sort keys %{$hash3{$CHROM_tok}{$SYMBOL_tok}})
{	
foreach my $strand_tok(sort keys %{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}})
{
	#~ unless (exists($hash_derived_PTC{$CHROM}{$SYMBOL}{$ENST}{$Effect}{$POS}{$REF}{$ALT}))
					
	foreach my $POS_tok(sort keys %{$hash_derived_PTC{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}})
	{
	foreach my $REF_tok(sort keys %{$hash_derived_PTC{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}})
	{
	foreach my $ALT_tok(sort keys %{$hash_derived_PTC{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}{$REF_tok}})
	{
	foreach my $Effect_tok(sort keys %{$hash_derived_PTC{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}{$REF_tok}{$ALT_tok}})
	{
	foreach my $POS_derived_tok(sort keys %{$hash_derived_PTC{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}})
	{
		#~ print "INICIO:$CHROM_tok\t$SYMBOL_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$Effect_tok\t$ENST_tok\t$strand_tok\n";
		
		my @END_tmp=sort {$a<=>$b} keys %{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'END'}};
		my @BEGIN_tmp=sort {$a<=>$b} keys %{$hash3{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'BEGIN'}};
		
		#~ print "El array Begin codificante es:@BEGIN_tmp\n";
		#~ print "El array End codificante es:@END_tmp\n";
		
		my @EXON_END_tmp=sort {$a<=>$b} keys %{$hash_EXON{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'END'}};
		my @EXON_BEGIN_tmp=sort {$a<=>$b} keys %{$hash_EXON{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'BEGIN'}};
		
		#~ print "El array de exones begin es:@EXON_BEGIN_tmp\n";
		#~ print "El array de exones end es:@EXON_END_tmp\n";
		
		if($POS_derived_tok > $BEGIN_tmp[0] && $POS_derived_tok < $END_tmp[scalar(@END_tmp)-1])
		{
			if($strand_tok eq '+')
			{	
				foreach my $threshold_NMD_tok(sort keys %{$hash_EXON{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'THRESHOLD'}})
				{
					#~ print "Threshold_NMD:$threshold_NMD_tok\n";
							
					if($threshold_NMD_tok eq 'NMD_negative')
					{
						# Genes with less than 50 pb upstream last EJC
						
						my $NMD_derived="NMD_negative";
						$RESULT_hash_derived_PTC{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$ENST_tok}{$POS_tok}{$NMD_derived}=1;	
					}
					elsif($threshold_NMD_tok eq 'NaN')
					{
						# Monoexonic genes
						
						my $NMD_derived="NaN";
						$RESULT_hash_derived_PTC{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$ENST_tok}{$POS_tok}{$NMD_derived}=1;	
					}
					else
					{
						if($END_tmp[scalar(@END_tmp)-1] > $threshold_NMD_tok)
						{
							if($POS_derived_tok < $threshold_NMD_tok)
							{
								my $NMD_derived="NMD_positive";
								$RESULT_hash_derived_PTC{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$ENST_tok}{$POS_tok}{$NMD_derived}=1;	
							}
							elsif($POS_derived_tok >= $threshold_NMD_tok)
							{
								my $NMD_derived="NMD_negative";
								$RESULT_hash_derived_PTC{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$ENST_tok}{$POS_tok}{$NMD_derived}=1;	
							}
						}
						else
						{
							# NMD threshold in 5'UTR
							
							my $NMD_derived="NaN";
							$RESULT_hash_derived_PTC{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$ENST_tok}{$POS_tok}{$NMD_derived}=1;	
						}
					}
				}
			}
			
			elsif($strand_tok eq '-')
			{
				foreach my $threshold_NMD_tok(sort keys %{$hash_EXON{$CHROM_tok}{$SYMBOL_tok}{$ENST_tok}{$strand_tok}{'THRESHOLD'}})
				{
					#~ print "Threshold_NMD:$threshold_NMD_tok\n";
							
					if($threshold_NMD_tok eq 'NMD_negative')
					{
						# Genes with less than 50 pb upstream last EJC
						
						my $NMD_derived="NMD_negative";
						$RESULT_hash_derived_PTC{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$ENST_tok}{$POS_tok}{$NMD_derived}=1;	
					}
					elsif($threshold_NMD_tok eq 'NaN')
					{
						# Monoexonic genes
						
						my $NMD_derived="NaN";
						$RESULT_hash_derived_PTC{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$ENST_tok}{$POS_tok}{$NMD_derived}=1;	
					}
					else
					{
						if($BEGIN_tmp[0] < $threshold_NMD_tok)
						{
							if($POS_derived_tok > $threshold_NMD_tok)
							{
								my $NMD_derived="NMD_positive";
								$RESULT_hash_derived_PTC{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$ENST_tok}{$POS_tok}{$NMD_derived}=1;	
							}
							elsif($POS_derived_tok <= $threshold_NMD_tok)
							{
								my $NMD_derived="NMD_negative";
								$RESULT_hash_derived_PTC{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$ENST_tok}{$POS_tok}{$NMD_derived}=1;	
							}
						}
						else
						{
							# NMD threshold in 5'UTR
							
							my $NMD_derived="NaN";
							$RESULT_hash_derived_PTC{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$ENST_tok}{$POS_tok}{$NMD_derived}=1;	
						}
					}
				}
			}
		}
		else
		{
			# POS_tok outside coding region, NMD does not apply
			
			my $NMD_derived="NaN";
			$RESULT_hash_derived_PTC{$CHROM_tok}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$ENST_tok}{$POS_tok}{$NMD_derived}=1;	
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

foreach my $CHROM_tok2(sort{$a<=>$b}keys%RESULT_hash_derived_PTC)
	{
	foreach my $SYMBOL_tok(sort keys%{$RESULT_hash_derived_PTC{$CHROM_tok2}})
	{
	foreach my $POS_tok(sort keys%{$RESULT_hash_derived_PTC{$CHROM_tok2}{$SYMBOL_tok}})
	{
	foreach my $REF_tok(sort keys%{$RESULT_hash_derived_PTC{$CHROM_tok2}{$SYMBOL_tok}{$POS_tok}})
	{
	foreach my $ALT_tok(sort keys%{$RESULT_hash_derived_PTC{$CHROM_tok2}{$SYMBOL_tok}{$POS_tok}{$REF_tok}})
	{
	foreach my $Effect_tok(sort keys%{$RESULT_hash_derived_PTC{$CHROM_tok2}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}})
	{
	foreach my $ENST_tok(sort keys%{$RESULT_hash_derived_PTC{$CHROM_tok2}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}})
	{
	foreach my $POS_derived_tok(sort keys%{$RESULT_hash_derived_PTC{$CHROM_tok2}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$ENST_tok}})
	{	
	foreach my $NMD_tok(sort keys%{$RESULT_hash_derived_PTC{$CHROM_tok2}{$SYMBOL_tok}{$POS_tok}{$REF_tok}{$ALT_tok}{$Effect_tok}{$ENST_tok}{$POS_derived_tok}})
	{
		print OUTPUT"$CHROM_tok2\t$SYMBOL_tok\t$POS_tok\t$REF_tok\t$ALT_tok\t$Effect_tok\t$ENST_tok\t->DERIVED:$POS_derived_tok\t$NMD_tok\n";
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
