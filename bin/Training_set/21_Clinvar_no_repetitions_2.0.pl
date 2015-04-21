##	Script to convert the clinVar.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

use strict;
use warnings;
use Time::localtime;

my %hash1=();
my %hash2=();

## Input file= prueba_clinvar2.txt

#	1	100185102	rs149824246	G	A	100	PASS	255;FRRS1;RCV000060275.2;autosomal dominant;No_modification_to_get_minimal_representation
#	1	100316614	rs113994126	C	T	100	PASS	5|5;AGL;RCV000001153.2|RCV000020373.1;autosomal recessive;No_modification_to_get_minimal_representation
#	1	100316614	rs113994127	CAG	C	100	PASS	5|5;AGL;RCV000001155.2|RCV000020374.1;autosomal recessive;No_modification_to_get_minimal_representation
#	1	10032156	rs387907294	G	A	100	PASS	5;NMNAT1;RCV000030771.1;autosomal dominant;No_modification_to_get_minimal_representation

my $input1=$ARGV[0];
my $output=$ARGV[1];

if (open(INPUT1, $input1))
{
	while(my $line=<INPUT1>)
	{
		chomp($line);
	#print "hello_world_1\n";
		if ($line=~/^[0-9]+/)
		{
			#print "$line";
			my @tmp=split("\t",$line);
			my $CHROM=$tmp[0];
			my $POS=$tmp[1];
			my $ID=$tmp[2];
			my $REF=$tmp[3];
			my $ALT=$tmp[4];
			my $FILTER=$tmp[5];
			my $PASS=$tmp[6];
			my $INFO=$tmp[7];
			my $string=join("\t",$CHROM,$POS,$ID,$REF,$ALT,$FILTER,$PASS);
			
			my @INFO_tmp=split(";",$INFO);
			
			my $CLNSIG=$INFO_tmp[0];
			my $SYMBOL=$INFO_tmp[1];
			my $CLNACC=$INFO_tmp[2];
			my $mode_of_inheritance=$INFO_tmp[3];
			my $TRIM=$INFO_tmp[4];
			
			$hash1{$string}{$SYMBOL}{$TRIM}{$mode_of_inheritance}{'CLNACC'}{$CLNACC}=1;
			$hash1{$string}{$SYMBOL}{$TRIM}{$mode_of_inheritance}{'CLNSIG'}{$CLNSIG}=1;
			#print "LA LÍNEA DE PARTIDA ES:$string\t$SYMBOL\t$TRIM\t$CLNACC\t$CLNSIG\n";
		}
	}
}

if (open(OUTPUT, '>'.$output))
{
	foreach my $string_tok(sort keys %hash1)
	{
		foreach my $SYMBOL_tok(sort keys %{$hash1{$string_tok}})
		{
			
			foreach my $TRIM_tok(sort keys %{$hash1{$string_tok}{$SYMBOL_tok}})
			{
				foreach my $mode_of_inheritance_tok(sort keys %{$hash1{$string_tok}{$SYMBOL_tok}{$TRIM_tok}})
				{
					
					print  OUTPUT "$string_tok\t$SYMBOL_tok;$TRIM_tok;$mode_of_inheritance_tok;";
					
					my @CLNSIG_tmp=sort keys %{$hash1{$string_tok}{$SYMBOL_tok}{$TRIM_tok}{$mode_of_inheritance_tok}{'CLNSIG'}};
					my @CLNSIG_tmp_joined=join("\|",@CLNSIG_tmp);
					print  OUTPUT "@CLNSIG_tmp_joined;";
					
					my @CLNACC_tmp=sort keys %{$hash1{$string_tok}{$SYMBOL_tok}{$TRIM_tok}{$mode_of_inheritance_tok}{'CLNACC'}};
					my @CLNACC_tmp_joined=join("\|",@CLNACC_tmp);
					print OUTPUT  "@CLNACC_tmp_joined\n";
					#print "El array1 es:@CLNACC_tmp\n";
					#print "El array2 es:@CLNSIG_tmp\n";
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
