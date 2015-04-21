use strict;
# use warnings;
use Time::localtime;

use Bio::EnsEMBL::Registry;
use Bio::SeqIO;



my %hash1=();
my %hash2=();
my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $output=$ARGV[2];

if (open(INPUT1, $input1))
{
	# Input file= UNIPROT_recovery_genes_global_2.txt

# ENSG00000000005 Q9H2S6  ENST00000373031 ENSP00000362122 NaN
# ENSG00000002746 Q76N89  ENST00000395891 ENSP00000379228 NaN
# ENSG00000002746 Q76N89  ENST00000453890 ENSP00000407774 NaN
# ENSG00000004948 P30988  ENST00000415529 ENSP00000413179 P30988-6
# ENSG00000004948 P30988  ENST00000423724 ENSP00000391369 P30988-5
# ENSG00000006377 P56179  ENST00000007660 ENSP00000007660 P56179-2
# ENSG00000006377 P56179  ENST00000518156 ENSP00000428480 P56179-3
# ENSG00000008256 O43739  ENST00000350796 ENSP00000297044 O43739-2
# ENSG00000008735 Q13387  ENST00000008876 ENSP00000008876 Q13387-3
# ENSG00000008853 Q9BYZ6  ENST00000251822 ENSP00000251822 NaN
# ENSG00000008853 Q9BYZ6  ENST00000519685 ENSP00000427926 NaN
# ENSG00000008853 Q9BYZ6  ENST00000522948 ENSP00000429141 NaN
# ENSG00000010671 Q06187  ENST00000308731 ENSP00000308176 NaN
# ENSG00000010932 Q01740  ENST00000354841 ENSP00000346901 NaN
# ENSG00000010932 Q01740  ENST00000367750 ENSP00000356724 NaN
# ENSG00000010932 Q01740  ENST00000402921 ENSP00000385543 NaN

	while(my $line=<INPUT1>)
	{
		chomp($line);
		
		#	First, we create a hash with the genes with no equivalence ENST-displayed isoform
		
		my @tmp=split("\t",$line);
		#my $string=join("***",@tmp);
		#print "$string\n";
		my $ENSG=$tmp[0];
		my $AC=$tmp[1];
		my $ENST=$tmp[2];
		my $ENSP=$tmp[3];
		my $AC_isoform=$tmp[4];
		$hash1{$AC}{$ENSG}{$ENST}{$ENSP}{$AC_isoform}=1;
		#print "$ENSG\t$ENST\t$ENSP\t$AC\t$AC_isoform\n";
	}
close(INPUT1);
}else {print "impossible to open INPUT1\n";die;}

if (open(INPUT2, $input2))
{
	# Input2=  uniprot_sprot.dat
	#	ID   1A_CMVO                 Reviewed;         993 AA.
	#AC   P20122;
	#DT   01-FEB-1991, integrated into UniProtKB/Swiss-Prot.
	#FT   CHAIN         1    993       Replication protein 1a.
	#FT                                /FTId=PRO_0000083262.
	#FT   DOMAIN      687    838       (+)RNA virus helicase ATP-binding.
	#FT   DOMAIN      839    993       (+)RNA virus helicase C-terminal.
	#FT   NP_BIND     714    721       ATP (Potential).
	#FT   REGION       50    409       Methyltransferase.
	#FT   REGION      712    975       ATP-dependent helicase.
	#SQ   SEQUENCE   993 AA;  111266 MW;  8FDEC1F3C66EBB4C CRC64;
	#     MATSSFNINE LVASHGDKGL LATALVDKAA HEQLEEQLQH QRRGRKVYVR NVLSVKDSEV
	#     IRNRYGGKYD LHLTQQEFAP HGLAGALRLC ETLDCLDSFP SSGLRQDLVL DFGGSWVTHY
	#     LRGHNVHCCS PCLGIRDKMR HTERLMNMRK IILNDPQQFD GRQPDFCTHP AADCKVQAHF
	#//

my @AC_tmp=();
my $FLAG=0;
my @sequence_tmp=();

while(my $line=<INPUT2>)
	{
		chomp($line);
		if ($line=~/^\/\//)
		{
			#print "hello_world_1\n";
			foreach my $AC_tok(@AC_tmp)
			{
				if(exists($hash1{$AC_tok}))
				{
					#print "hello_world_1:$AC_tok\n";
					foreach my $ENSG_tok(sort keys %{$hash1{$AC_tok}})
					{
						my $seq=join("",@sequence_tmp);
						$hash2{$ENSG_tok}{$AC_tok}{$seq}=1;
						#print "$ENSG_tok\t$AC_tok\t$seq\n";
					}
				}
			}
			@AC_tmp=();
			$FLAG=0;
			@sequence_tmp=();
		}
		elsif($line=~/^AC\s+([^\;]+);/)
		{
			my $AC=$1;
			#print "$AC\n";
			my @tmp=split(";",$AC);
			foreach my $tmp_tok(@tmp)
			{
				$tmp_tok=~s/\s+//g;
				push(@AC_tmp,$tmp_tok);
			}
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

my %hash3=();
my $registry = "Bio::EnsEMBL::Registry";
$registry->load_registry_from_url('mysql://anonymous@ensembldb.ensembl.org');
my $gene_adaptor=$registry->get_adaptor( 'Human', 'Core', 'Gene' );




my @tmp_genes=sort keys(%hash2);
print "El scalar es:".scalar(@tmp_genes)."\n";
my $ENSG_gene="NaN";
my $ENST_transcript="NaN";
my $seq="NaN";

if(open(OUTPUT, '>'.$output))
{
foreach my $gene_tok(@tmp_genes)
{
	#my $gene=$gene_adaptor->fetch_by_stable_id($gene_tok);
	if(defined(my $gene=$gene_adaptor->fetch_by_stable_id($gene_tok)))
	
	{
		$ENSG_gene=$gene->stable_id();

		# print "AAAAA:$ENSG_gene:A\n";
		my @transcripts = @{ $gene->get_all_Transcripts };

		foreach my $transcript(@transcripts)
		{
			my $ENST_transcript=$transcript->stable_id();
			# print "$ENST_transcript\n";
			if(defined($transcript->translation()))
			{
				my $CDS=$transcript->translateable_seq();
				## print "$CDS\n";
				my $ENSP=$transcript->translation()->stable_id();
				## print "$ENSP\n";
				$seq=$transcript->translate()->seq();
				$seq=~s/\.|\s+//g;
				# print "CANDIDATE:$seq".":::::::::::\n";
				$hash3{$ENSG_gene}{$ENST_transcript}{$seq}=1;
				foreach my $AC_tok(sort keys %{$hash2{$gene_tok}})
				{
					# print "HELLO_WORLD_IV:$AC_tok\n";
					foreach my $seq_displayed_tok(sort keys %{$hash2{$gene_tok}{$AC_tok}})
					{
						# print "OPPONENTS:$seq_displayed_tok"."****************\n";
						if ($seq eq $seq_displayed_tok)
						{
							print OUTPUT ">$AC_tok\t$ENSG_gene\t$ENST_transcript\t$ENSP\n";
							print OUTPUT "$seq\n";
						}
						
						elsif($seq=~/$seq_displayed_tok/)
						{
							my @seq_tmp=split('',$seq);
							my @seq_displayed_tok_tmp=split('',$seq_displayed_tok);
							my $POS_trimmed=0;
							if($seq_tmp[scalar(@seq_tmp)-1] eq $seq_displayed_tok_tmp[scalar(@seq_displayed_tok_tmp)-1])
							{
								while ($seq_tmp[scalar(@seq_tmp)-1] eq $seq_displayed_tok_tmp[scalar(@seq_displayed_tok_tmp)-1] && scalar(@seq_displayed_tok_tmp)>=0)
								{
									## Extract last element of the array.
									#print OUTPUT "Hello_world_1\n";
									pop (@seq_tmp);
									pop (@seq_displayed_tok_tmp);
								}
								
							$POS_trimmed=scalar(@seq_tmp);
							print OUTPUT "ALTERED_COORDINATES"."("."upstream:$POS_trimmed".")>"."$AC_tok\t$ENSG_gene\t$ENST_transcript\t$ENSP\n";
							print OUTPUT "$seq\n";
								
							}
							elsif($seq_tmp[0] eq $seq_displayed_tok_tmp[0])
							{													
								while ($seq_tmp[0] eq $seq_displayed_tok_tmp[0] && scalar(@seq_displayed_tok_tmp)>=0)
								{
									#print OUTPUT "Hello_world_2\n";
									shift (@seq_tmp);
									shift (@seq_displayed_tok_tmp);
								}
							$POS_trimmed=scalar(@seq_tmp);
							print OUTPUT "ALTERED_COORDINATES"."("."downstream:$POS_trimmed".")>"."$AC_tok\t$ENSG_gene\t$ENST_transcript\t$ENSP\n";
							print OUTPUT "$seq\n";
							}
							elsif($seq_tmp[0] ne $seq_displayed_tok_tmp[0] && $seq_tmp[scalar(@seq_tmp)-1] ne $seq_displayed_tok_tmp[scalar(@seq_displayed_tok_tmp)-1]&& scalar(@seq_displayed_tok_tmp)>=0)
							{	
								my @number=();												
								while ($seq_tmp[0] ne $seq_displayed_tok_tmp[0])
								{
									#print OUTPUT "Hello_world_2\n";
									my $number=shift (@seq_tmp);
									push (@number, $number);
								}
							$POS_trimmed=scalar(@number);
							print OUTPUT "ALTERED_COORDINATES"."("."upstream:$POS_trimmed".")>"."$AC_tok\t$ENSG_gene\t$ENST_transcript\t$ENSP\n";
							print OUTPUT "$seq\n";
							}
						}
						elsif ($seq ne $seq_displayed_tok)
						{
							my @seq_tmp=split('',$seq);
							my @seq_displayed_tok_tmp=split('',$seq_displayed_tok);
							my $POS_trimmed=0;
							while ($seq_tmp[scalar(@seq_tmp)-1] eq $seq_displayed_tok_tmp[scalar(@seq_displayed_tok_tmp)-1])
								{
									## Extract last element of the array.
									#print OUTPUT "Hello_world_1\n";
									pop (@seq_tmp);
									pop (@seq_displayed_tok_tmp);
								}
								while ($seq_tmp[0] eq $seq_displayed_tok_tmp[0])
								{
									#print OUTPUT "Hello_world_2\n";
									shift (@seq_tmp);
									shift (@seq_displayed_tok_tmp);
								}
								
							if(scalar(@seq_tmp)==1 && scalar(@seq_displayed_tok_tmp)==1)
							{
								print OUTPUT "MISMATCH"."("."1".")>"."$AC_tok\t$ENSG_gene\t$ENST_transcript\t$ENSP\n";
								print OUTPUT "$seq\n";
							}
						}
					}
				}
			}
			else{}#print "$ENST_transcript\tnon_protein_coding\n";}
		}
	}
	else{print "$gene_tok\tno_stable_id_in_ENSEMBL_release_75\n";}
}
close(OUTPUT);
}else{print "impossible to open OUTPUT\n";die;}
sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
