##	Script to condense further the intervals of genomic positions mapping the protein features mapped.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;

my %hash1=();
my %hash2=();
my %hash3=();
my %hash4=();

my $input1=$ARGV[0];
my $output1=$ARGV[1];

my $time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

# Here we parse the features on interval per line basis

if (open (INPUT, $input1))
{
	# Input=  data/build_tables/PROTEIN_condensed.txt
	
#~ [209768395-209768412]   1       +       ENST00000009105 ENSG00000008118 IPR000719**1    IPR002290**1    IPR011009**1    IPR020636**1
#~ [209773404-209773456]   1       +       ENST00000009105 ENSG00000008118 IPR000719**1    IPR002290**1    IPR011009**1    IPR020636**1
#~ [209776559-209776633]   1       +       ENST00000009105 ENSG00000008118 IPR000719**1    IPR002290**1    IPR011009**1    IPR020636**1
#~ [209778881-209778998]   1       +       ENST00000009105 ENSG00000008118 IPR000719**1    IPR002290**1    IPR011009**1    IPR020636**1
#~ [209779683-209779788]   1       +       ENST00000009105 ENSG00000008118 IPR000719**1    IPR002290**1    IPR011009**1    IPR020636**1
#~ [209781203-209781278]   1       +       ENST00000009105 ENSG00000008118 IPR000719**1    IPR002290**1    IPR011009**1    IPR020636**1
#~ [209782325-209782437]   1       +       ENST00000009105 ENSG00000008118 IPR000719**1    IPR002290**1    IPR011009**1    IPR020636**1
#~ [209783196-209783274]   1       +       ENST00000009105 ENSG00000008118 IPR000719**1    IPR002290**1    IPR011009**1    IPR020636**1
#~ [209784810-209784813]   1       +       ENST00000009105 ENSG00000008118 IPR000719**1    IPR002290**1    IPR011009**1    IPR020636**1
#~ [209768368-209768394]   1       +       ENST00000009105 ENSG00000008118 IPR011009**1    IPR020636**1
#~ [209784814-209784897]   1       +       ENST00000009105 ENSG00000008118 IPR011009**1    IPR020636**1
#~ [209785137-209785157]   1       +       ENST00000009105 ENSG00000008118 IPR011009**1    IPR020636**1
#~ [209768329-209768367]   1       +       ENST00000009105 ENSG00000008118 IPR020636**1
     
	while(my $line=<INPUT>)
	{
		chomp $line;
		#~ print "hello_worldI:$line**\n";
		if ($line=~/(^\[.+\])\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)/)
		{
			#~ print "hello_worldII:$line**\n";
			my $interval= $1;
			my $CHROM=$2;
			my $strand=$3;
			my $ENST=$4;
			my $ENSG=$5;
			my $string=join("\t",$CHROM,$strand,$ENST,$ENSG);
			my $fields=$6;
			#print "$interval\t$CHROM\t$strand\t$ENST\t$ENSG*****************$fields\n";
			my @features_tmp=split(/\t/,$fields);
			foreach my $features_tmp_tok(@features_tmp)
			{
				$hash1{$string}{$features_tmp_tok}{$interval}=1;
				#~ print "$string\t$features_tmp_tok\t$interval\n";
			}
		}
	}
}else{print "Unable to open INPUT\n";}

# Here we print all the intervals belonging to the same key on a per line basis

$time='['. timestamp(). ']'."\n";
print "Start printing:$time\n";

if(open (OUTPUT, '>'.$output1))
{
	foreach my $string_tok(sort keys%hash1)
	{
	foreach my $feature_tok(sort keys%{$hash1{$string_tok}})
	{
		print OUTPUT "$string_tok\t$feature_tok\t";
	foreach my $interval_tok(sort keys %{$hash1{$string_tok}{$feature_tok}})
	{
		print OUTPUT "$interval_tok\t";	
	}
		print OUTPUT "\n";
	}
	}
}else {print "Impossible to open OUTPUT\n";}

$time='['. timestamp(). ']'."\n";
print "Fin del script:$time\n";
sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
