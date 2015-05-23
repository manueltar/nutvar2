##	Script to create a file the transcripts mapping to each genomic coordinate. 2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;

my $input1=$ARGV[0];
my $output1=$ARGV[1];


my %features_hash=();

if (open(INPUT1, $input1))
{

#~ input1= gtf_tabladef_sorted_by_SYMBOL.txt
#~ 
#~ ENST00000263100 19      A1BG    -       protein_coding  CCDS12976       ENSP00000263100 UTR_5_prime:[58864804]-[58864865]       START:[58864801]-[58864803]     CDS:[58864770]-[58864803]       CDS:[58864658]-[58864693]       CDS:[58864294]-[58864563]       CDS:[58863649]-[58863921]
       #~ CDS:[58862757]-[58863053]       CDS:[58861736]-[58862017]       CDS:[58858719]-[58859006]       CDS:[58858391]-[58858395]       STOP:[58858388]-[58858390]      UTR_3_prime:[58858216]-[58858387]       
#~ ENST00000600966 19      A1BG    -       protein_coding  NaNein  ENSP00000470909 UTR_5_prime:[NaN]-[NaN] START:[NaN]-[NaN]       CDS:[58864294]-[58864495]       CDS:[58863649]-[58863921]       CDS:[58863464]-[58863550]       CDS:[58862757]-[58863053]       CDS:[58861960]-[58862017]
       #~ STOP:[NaN]-[NaN]        UTR_3_prime:[NaN]-[NaN] 
#~ ENST00000282641 10      A1CF    -       protein_coding  CCDS7242        ENSP00000282641 UTR_5_prime:[52645341]-[52645435]       UTR_5_prime:[52623793]-[52623840]       UTR_5_prime:[52619701]-[52619745]       START:[52619698]-[52619700]     CDS:[52619602]-[52619700]       CDS:[52603748]-[52603882]       CDS:[52601622]-[52601752]       CDS:[52595834]-[52596072]       CDS:[52587891]-[52588055]       CDS:[52580312]-[52580409]
       #~ CDS:[52575766]-[52576039]       CDS:[52573617]-[52573822]       CDS:[52570800]-[52570936]       CDS:[52569654]-[52569802]       CDS:[52566492]-[52566640] 

	while (my $line=<INPUT1>)
	{
		chomp $line;
		#~ print "$line\n";
		if ($line=~/^(ENST[^\t]*)\t([^\t]+)\t([^\t]+)\t([^\t]*)\t([^\t]*)\t(.+)\t$/)
		{
			#~ print "Hello_world:$line\n";
			my $CHROM=$2;
			my $SYMBOL=$3;
			# print "El symbol es:my $SYMBOL\n";
			my $ENST=$1;
			my $BIOTYPE=$5;
			my $strand=$4;
			my $tabladef_fields=$6;
			my @tabladef_fields_tmp=split(/\t/,$tabladef_fields);
			#~ my $array_string=join("**",@tabladef_fields_tmp);
			#~ print "El arrray es:$CHROM\t$SYMBOL\t$strand\t$array_string\n";
			
			# as we are interested only in the coding positions of transcripts (to map protein features). We only extract CDS postions from the table of transcript features.
			
			foreach my $tabladef_fields_tok(@tabladef_fields_tmp)
			{
				if($tabladef_fields_tok=~/CDS:\[(.+)\]\-\[(.+)\]/)
				{
					#~ print "hello_world_2:$tabladef_fields_tok\n";
					my $begin=$1;
					my $end=$2;
					
					# Now, we convert the begin/end coordinates of each CDDS segment into the begin coordinate plus a distance to the end to know how many bases
					# The CDS segment expands
					
					my $distance=$end-$begin+1;
					my $string=join("\t",$CHROM,$SYMBOL,$strand);
					
					# We store all in a hash
					
					$features_hash{$string}{$ENST}{$begin}{$distance}=1;
					#~ print "$SYMBOL**\t$string**\t$begin\t$distance**\n";
				}
				
			}
		}
	
	}
}else{print "Unable to open gtf_tabladef_sorted_by_SYMBOL.txt\n";}



if (open(OUTPUT, '>'.$output1))
{
	foreach my $string_tok(sort keys %features_hash)
	{
		# On a string basis (Chrom-symbol-strand, that is on a gene basis) we define a hash to be loaded with all the CDS positions of all the ENSTs in that gene
		
		my %position_hash=();
	
		foreach my $ENST_tok(sort keys %{$features_hash{$string_tok}})
		{
		foreach my $begin_tok(sort keys %{$features_hash{$string_tok}{$ENST_tok}})
		{
		foreach my $distance_tok(sort keys %{$features_hash{$string_tok}{$ENST_tok}{$begin_tok}})
		{
			# Now we transform the begin +distance into all the position the CDS segment expands;
			# So the positons start at the begin of the CDS and expand while $i< begin+distance
			# WE charge the array with all the trasncripts CDS positions and the trasncripts occupying them.
			
			for (my$i=$begin_tok;$i<($begin_tok+$distance_tok);$i++)
			{
				$position_hash{$i}{$string_tok}{$ENST_tok}=1;
			}
		}
		}
		}
		
		# Now we print the transcripts that map to genomic coordinates per gene on a per position basis
		
		foreach my $POS_tok(sort{$a<=>$b} keys %position_hash)
		{
			print OUTPUT "$POS_tok\t";
		foreach my $string_tok(sort keys %{$position_hash{$POS_tok}})
		{
			print OUTPUT "$string_tok\t";
		foreach my $ENST_tok(sort keys %{$position_hash{$POS_tok}{$string_tok}})
		{
			print OUTPUT "$ENST_tok\t";
		}
		}
			print OUTPUT "\n";	
		}
	}#
}else{print "Unable to open OUTPUT\n";}

sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
