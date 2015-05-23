use strict;
use warnings;
use Time::localtime;

my %initial_hash=();
my %hash1=();
my %hash2=();
my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $output=$ARGV[2];

my $time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

if (open (INPUT1, $input1))
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

while(my $line=<INPUT1>)
	{
		chomp $line;
		#~ print "$line\n";
		if($line=~/^([^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)\t[^\t]+\t[^\t]+\t(.+)/)
		{
			my $CHROM=$1;
			my $POS=$2;
			my $REF=$3;
			my $ALT=$4;
			my $fields=$5;
			#~ print "hello_world:$CHROM\t$POS\t$fields\t$REF\t$ALT\n";
			my @fields_tmp=split(";",$fields);
			my $Effect=$fields_tmp[0];
			my $SYMBOL=$fields_tmp[1];
			my $ENST=$fields_tmp[2];
			
			if ($Effect =~/frameshift_variant/)
			{
				#~ print "hello_world:$CHROM\t$POS\t$fields\t$REF\t$ALT\n";
				my $length_REF=length($REF);
				my $length_ALT=length($ALT);
				if($length_REF > $length_ALT)
				{
					#~  ATG	->	AC ; 3-2; $length_ALT +1
					my $del_expand=$length_REF - $length_ALT;
					$hash1{$CHROM}{$SYMBOL}{$ENST}{$POS}{$REF}{$ALT}{$Effect}{'DEL'}{$length_ALT}{$del_expand}=1;
					#~ print "$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t$Effect\t$ENST\t'DEL'\t$length_ALT\t$del_expand\n";
					
				}
				elsif($length_REF < $length_ALT)
				{
					my $ins_expand=$length_ALT - $length_REF;
					$hash1{$CHROM}{$SYMBOL}{$ENST}{$POS}{$REF}{$ALT}{$Effect}{'INS'}{$length_REF}{$ins_expand}=1;
					#~ print "$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t$Effect\t$ENST\t'INS'\t$length_REF\t$ins_expand\n";
					
				}
				
			}
		}
	}
}else{print "Unable to open INPUT1\n";}

my $CHROM="NaN";
my $ENSG="NaN";
my $ENST="NaN";
my $SYMBOL="NaN";
my $strand="NaN";
my $seq="NaN";
my %position_hash=();

$time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

if(open (INPUT2, $input2) && open (OUTPUT,'>'.$output))
{
	
	#~ Input:CDS_genomic_coordinates_full_compresed.txt
	#~ 
	#~ >X      ENSG00000000003 TSPAN6  ENST00000373020 -
	#~ seq>NaN
	#~ POS>[99883779-99884516]
	#~ >X      ENSG00000000005 TNMD    ENST00000373031 +
	#~ seq>NaN
	#~ POS>[99840016-99840063] [99840228-99840359]     [99848892-99849032]     [99849258-99849359]     [99852501-99852654]     [99854013-99854179]     [99854505-99854714]

	while(my $line=<INPUT2>)
	{
		chomp $line;
		#~ print "$line\n";
		if ($line=~/^\/\//)
		{	
			my @seq_tmp=@{$position_hash{'BASE'}};
			if($strand eq '-')
			{
				@seq_tmp=reverse @seq_tmp;
			}
			
			my $seq=join("",@seq_tmp);
			my @pos_tmp=sort {$a<=>$b}keys %{$position_hash{'COORDINATE'}};
			
			
			%position_hash=();	
		}
		elsif ($line=~/^>([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			$CHROM=$1;
			$ENSG=$2;
			$SYMBOL=$3;
			$ENST=$4;
			$strand=$5;
			#~ print "$CHROM\t$ENSG\t$SYMBOL\t$ENST\t$strand\n";
		}
		if($line !~ /^>/ && $line !~ /^\/\//)
		{
			if ($line=~/^(.+)\t$/)
			{
				my $POS=$1;
				#~ print "****************************$POS\n";
				my @POS_tmp=split(/\t/,$POS);
				#~ print "****************************".join("**",@POS_tmp)."\n";
				foreach my $POS_tmp_tok(@POS_tmp)
				{
					my @tmp=split("__",$POS_tmp_tok);
					#~ print "****************************".join("**",@tmp)."\n";
					my $coordinate=$tmp[0];
					my $base=$tmp[1];
					$position_hash{'COORDINATE'}{$coordinate}=1;
					push(@{$position_hash{'BASE'}},$base);
					#~ print "$coordinate\t$base\n";
				}
			}
		}
	}
}else{print "Unable to open INPUT2\n";}

$time='['. timestamp(). ']'."\n";
print "Print FIN:$time\n";
	
sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
