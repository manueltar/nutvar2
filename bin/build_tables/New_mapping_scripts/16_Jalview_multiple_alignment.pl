##	Script to extract the best UniPtor-ENSEMBL isoform pairs of aligment from the BLAST multiple alignment using as criteria the bit score. 
# The biggest the bit score, the best the alignment of the different queries for a given subject.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;

my $input=$ARGV[0];
my $output1=$ARGV[1];

my %hash1=();

my $time='['. timestamp(). ']'."\n";
print "Tiempo de incio de carga_hash1:$time\n";

# He open the file carrying the information of the alignment but without the sequences

if(open(INPUT,$input))
{
	#~ Input=resultsfmt7.out
#~ # BLASTP 2.2.28+
#~ # Query: sp|O43739|ENSG00000008256_ENST00000350796
#~ # Database: db.fasta
#~ # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
#~ # 2 hits found
#~ sp|O43739|ENSG00000008256_ENST00000350796       sp|O43739|ENSG00000008256_UNIPROT       99.73   374     0       1       27      399     27
      #~ 400     0.0      769
#~ sp|O43739|ENSG00000008256_ENST00000350796       sp|Q9Y4C0|ENSG00000021645_UNIPROT       30.19   53      35      2       172     224     863
     #~ 913     1.7     21.6
#~ # BLASTP 2.2.28+
#~ # Query: sp|O43739|ENSG00000008256_ENST00000396741
#~ # Database: db.fasta
#~ # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
#~ # 3 hits found
#~ sp|O43739|ENSG00000008256_ENST00000396741       sp|O43739|ENSG00000008256_UNIPROT       98.82   255     1       2       62      314     146
     #~ 400     0.0      526
#~ sp|O43739|ENSG00000008256_ENST00000396741       sp|Q9Y4C0|ENSG00000021645_UNIPROT       30.19   53      35      2       87      139     863
     #~ 913     1.8     21.2
#~ sp|O43739|ENSG00000008256_ENST00000396741       sp|Q9NXG0|ENSG00000044459_UNIPROT       40.00   20      12      0       54      73      1362
    #~ 1381    3.5     20.0


	while(my $line=<INPUT>)
		{
			chomp $line;
			unless ($line=~/^#/)
			{
				#~ print "**$line**\n";
				if($line=~/^([^\t]+)\t([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)/)
				{
					# We define query as the differnet trasncripts and subject as the UniProt displayed isoform 
					
					
					#~ print "Hello_world:$line\n";
					my $query=$1;
					my $subject=$2;
					my $bit_score=$3;
					my @tmp_query=split(/\|/,$query);
					my $AC_query=$tmp_query[1];
					#~ print "AC_query:$AC_query\n";
					my @tmp_subject=split(/\|/,$subject);
					my $AC_subject=$tmp_subject[1];
					#~ print "AC_subject:$AC_subject\n";
					if($AC_query eq $AC_subject)
					{
						# We store in a hash all the subjects, bit socres and queries.
						$hash1{$subject}{$bit_score}{$query}=1;
						#~ print "$subject\t$bit_score\t$query\n";
					}
				}
				
			}
			
		}
}

$time='['. timestamp(). ']'."\n";
print "Tiempo de IMPRESIÓN:$time\n";

if(open(OUTPUT, '>'.$output1))
{
foreach my $subject_tok(sort keys %hash1)
{
	# For every subject (UniProt isoform) que recover the query with the highest bit score.
	# The biggest the bit score, the best the alignment of the different queries for a given subject.
	
	my @bit_score_tmp=reverse sort{$a<=>$b} keys%{$hash1{$subject_tok}};
	my @query_tmp=sort keys%{$hash1{$subject_tok}{$bit_score_tmp[0]}};
	
	print OUTPUT "$query_tmp[0]\t$subject_tok\t$bit_score_tmp[0]\n";
	
}
}

sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}



