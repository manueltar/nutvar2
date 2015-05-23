##	Script to calculate the different regions not to be displayed and offsets for the selected alignments regarding their bit scores.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

use strict;
use warnings;
use Time::localtime;

my %hash1=();
my %hash2=();
my %check_length_query=();
my %check_length_subject=();



my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $input3=$ARGV[2];
my $input4=$ARGV[3];
my $output=$ARGV[4];

my $time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

# We extract the selected alignments; query= ENSEMBL trasncript; suject= UniProt isoform displayed

if (open (INPUT1, $input1))
{
	#~ Input=selected_alignments.txt
	#~ 
#~ sp|O43739|ENSG00000008256_ENST00000350796       sp|O43739|ENSG00000008256_UNIPROT        769
#~ sp|P00450|ENSG00000047457_ENST00000264613       sp|P00450|ENSG00000047457_UNIPROT       2236
#~ sp|P07358|ENSG00000021852_ENST00000371237       sp|P07358|ENSG00000021852_UNIPROT       1237
#~ sp|P07686|ENSG00000049860_ENST00000261416       sp|P07686|ENSG00000049860_UNIPROT       1155
#~ sp|P08588|ENSG00000043591_ENST00000369295       sp|P08588|ENSG00000043591_UNIPROT        949
#~ sp|P08603|ENSG00000000971_ENST00000367429       sp|P08603|ENSG00000000971_UNIPROT       2545
#~ sp|P15502|ENSG00000049540_ENST00000358929       sp|P15502|ENSG00000049540_UNIPROT        316
#~ sp|P30988|ENSG00000004948_ENST00000359558       sp|P30988|ENSG00000004948_UNIPROT       1064
#~ sp|P41180|ENSG00000036828_ENST00000296154       sp|P41180|ENSG00000036828_UNIPROT       2227
#~ sp|Q13387|ENSG00000008735_ENST00000329492       sp|Q13387|ENSG00000008735_UNIPROT       1355
#~ sp|Q5BKX6|ENSG00000022567_ENST00000517878       sp|Q5BKX6|ENSG00000022567_UNIPROT       1584
#~ sp|Q99594|ENSG00000007866_ENST00000338863       sp|Q99594|ENSG00000007866_UNIPROT        898
#~ sp|Q9C0H9|ENSG00000017373_ENST00000264659       sp|Q9C0H9|ENSG00000017373_UNIPROT       2001
#~ sp|Q9NQ90|ENSG00000047617_ENST00000356134       sp|Q9NQ90|ENSG00000047617_UNIPROT       2080
#~ sp|Q9NUL3|ENSG00000040341_ENST00000524300       sp|Q9NUL3|ENSG00000040341_UNIPROT       1166
#~ sp|Q9NXG0|ENSG00000044459_ENST00000262360       sp|Q9NXG0|ENSG00000044459_UNIPROT       2840
#~ sp|Q9P299|ENSG00000005243_ENST00000006101       sp|Q9P299|ENSG00000005243_UNIPROT        382

while(my $line=<INPUT1>)
	{
		chomp $line;
		#~ print "$line\n";
		if($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			#~ print "Hello_world:$line\n";
			my $query=$1;
			my $subject=$2;
			my $bit_score=$3;
			$hash1{$query}{$subject}{$bit_score}=1;
			#~ print "$query\t$subject\t$bit_score\n";
		}
	}
}else{print "Unable to open INPUT1\n";}

my $query="NaN";

$time='['. timestamp(). ']'."\n";
print "Start charging hash2:$time\n";

my $query_check="NaN";

# To check the original length of the sequences we extract them from the files introduced in the aligment. Here the queries.

if (open (INPUT2, $input2))
{
	#~ Input=no_aligned2.fasta
	#~ >sp|A0FGR8|ENSG00000117868_ENST00000251527
	#~ MTPPSRAEAGVRRSRVPSEGRWRGAEPPGISASTQPASAGRAARHCGAMSGARGEGPEAGAGGAGGRAAP
	#~ ENPGGVLSVELPGLLAQLARSFALLLPVYALGYLGLSFSWVLLALALLAWCRRSRGLKALRLCRALALLE
	#~ DEERVVRLGVRACDLPAWVHFPDTERAEWLNKTVKHMWPFICQFIEKLFRETIEPAVRGANTHLSTFSFT
	#~ KVDVGQQPLRINGVKVYTENVDKRQIILDLQISFVGNCEIDLEIKRYFCRAGVKSIQIHGTMRVILEPLI
	#~ GDMPLVGALSIFFLRKPLLEINWTGLTNLLDVPGLNGLSDTIILDIISNYLVLPNRITVPLVSEVQIAQL
	#~ RFPVPKGVLRIHFIEAQDLQGKDTYLKGLVKGKSDPYGIIRVGNQIFQSRVIKENLSPKWNEVYEALVYE
	#~ HPGQELEIELFDEDPDKDDFLGSLMIDLIEVEKERLLDEWFTLDEVPKGKLHLRLEWLTLMPNASNLDKV
	#~ LTDIKADKDQANDGLSSALLILYLDSARNLPSGKKISSNPNPVVQMSVGHKAQESKIRYKTNEPVWEENF
	#~ TFFIHNPKRQDLEVEVRDEQHQCSLGNLKVPLSQLLTSEDMTVSQRFQLSNSGPNSTIKMKIALRVLHLE
	#~ KRERPPDHQHSAQVKRPSVSKEGRKTSIKSHMSGSPGPGGSNTAPSTPVIGGSDKPGMEEKAQPPEAGPQ
	#~ GLHDLGRSSSSLLASPGHISVKEPTPSIASDISLPIATQELRQRLRQLENGTTLGQSPLGQIQLTIRHSS
	#~ QRNKLIVVVHACRNLIAFSEDGSDPYVRMYLLPDKRRSGRRKTHVSKKTLNPVFDQSFDFSVSLPEVQRR
	#~ TLDVAVKNSGGFLSKDKGLLGKVLVALASEELAKGWTQWYDLTEDGTRPQAMT
	#~ >sp|A0FGR8|ENSG00000117868_ENST00000275418
	#~ SASTQPASAGRAARHCGAMSGARGEGPEAGAGGAGGRAAPENPGGVLSVELPGLLAQLARSFALLLPVYA
	#~ LGYLGLSFSWVLLALALLAWCRRSRGLKALRLCRALALLEDEERVVRLGVRACDLPAWVHFPDTERAEWL
	#~ NKTVKHMWPFICQFIEKLFRETIEPAVRGANTHLSTFSFTKVDVGQQPLRINGVKVYTENVDKRQIILDL
	#~ QISFVGNCEIDLEIKRYFCRAGVKSIQIHGTMRVILEPLIGDMPLVGALSIFFLRKPLLEINWTGLTNLL
	#~ DVPGLNGLSDTIILDIISNYLVLPNRITVPLVSEVQIAQLRFPVPKGVLRIHFIEAQDLQGKDTYLKGLV
	#~ KGKSDPYGIIRVGNQIFQSRVIKENLSPKWNEVYEALVYEHPGQELEIELFDEDPDKDDFLGSLMIDLIE
	#~ VEKERLLDEWFTLDEVPKGKLHLRLEWLTLMPNASNLDKVLTDIKADKDQANDGLSSALLILYLDSARNL
	#~ PSNPLEFNPDVLKKTAVQRALKSGKKISSNPNPVVQMSVGHKAQESKIRYKTNEPVWEENFTFFIHNPKR
	#~ QDLEVEVRDEQHQCSLGNLKVPLSQLLTSEDMTVSQRFQLSNSGPNSTIKMKIALRVLHLEKRERPPDHQ
	#~ HSAQVKRPSVSKEGRKTSIKSHMSGSPGPGGSNTAPSTPVIGGSDKPGMEEKAQPPEAGPQGLHDLGRSS
	#~ SSLLASPGHISVKEPTPSIASDISLPIATQELRQRLRQLENGTTLGQSPLGQIQLTIRHSSQRNKLIVVV
	#~ HACRNLIAFSEDGSDPYVRMYLLPDKRRSGRRKTHVSKKTLNPVFDQSFDFSVSLPEVQRRTLDVAVKNS
	#~ GGFLSKDKGLLGKVLVALASEELAKGWTQWYDLTEDGTRPQAMT
	while(my $line=<INPUT2>)
	{
		chomp($line);
		if($line=~/^>(.+)/)
		{
			$query_check=$1;
		}
		if(exists($hash1{$query_check}))
		{
			unless($line=~/^>/)
			{
				#~ print "query_check $query_check:$line**\n";
				push(@{$check_length_query{$query_check}},$line);
			}
			
		}
	}
}

$time='['. timestamp(). ']'."\n";
print "Start charging hash3:$time\n";

# To check the original length of the sequences we extract them from the files introduced in the aligment. Here the subjects.

my $subject_check=$1;
my $AC_check="NaN";

if (open (INPUT3, $input3))
{
	#~ Input=db.fasta
	#~ 
	#~ >sp|A0FGR8|ENSG00000117868_UNIPROT
	#~ MTANRDAALSSHRHPGCAQRPRTPTFASSSQRRSAFGFDDGNFPGLGERSHAPGSRLGARRRAKTARGLR
	#~ GHRQRGAGAGLSRPGSARAPSPPRPGGPENPGGVLSVELPGLLAQLARSFALLLPVYALGYLGLSFSWVL
	#~ LALALLAWCRRSRGLKALRLCRALALLEDEERVVRLGVRACDLPAWVHFPDTERAEWLNKTVKHMWPFIC
	#~ QFIEKLFRETIEPAVRGANTHLSTFSFTKVDVGQQPLRINGVKVYTENVDKRQIILDLQISFVGNCEIDL
	#~ EIKRYFCRAGVKSIQIHGTMRVILEPLIGDMPLVGALSIFFLRKPLLEINWTGLTNLLDVPGLNGLSDTI
	#~ ILDIISNYLVLPNRITVPLVSEVQIAQLRFPVPKGVLRIHFIEAQDLQGKDTYLKGLVKGKSDPYGIIRV
	#~ GNQIFQSRVIKENLSPKWNEVYEALVYEHPGQELEIELFDEDPDKDDFLGSLMIDLIEVEKERLLDEWFT
	#~ LDEVPKGKLHLRLEWLTLMPNASNLDKVLTDIKADKDQANDGLSSALLILYLDSARNLPSGKKISSNPNP
	#~ VVQMSVGHKAQESKIRYKTNEPVWEENFTFFIHNPKRQDLEVEVRDEQHQCSLGNLKVPLSQLLTSEDMT
	#~ VSQRFQLSNSGPNSTIKMKIALRVLHLEKRERPPDHQHSAQVKRPSVSKEGRKTSIKSHMSGSPGPGGSN
	#~ TAPSTPVIGGSDKPGMEEKAQPPEAGPQGLHDLGRSSSSLLASPGHISVKEPTPSIASDISLPIATQELR
	#~ QRLRQLENGTTLGQSPLGQIQLTIRHSSQRNKLIVVVHACRNLIAFSEDGSDPYVRMYLLPDKRRSGRRK
	#~ THVSKKTLNPVFDQSFDFSVSLPEVQRRTLDVAVKNSGGFLSKDKGLLGKVLVALASEELAKGWTQWYDL
	#~ TEDGTRPQAMT
	
	while(my $line=<INPUT3>)
	{
		chomp($line);
		#~ print "db.fasta:$line\n";
		if($line=~/^>(.+)/)
		{
			$subject_check=$1;
		}
		
		my @tmp=split(/\|/,$subject_check);
		$AC_check=$tmp[1];
		#~ print "$AC_check\n";
		unless($line=~/^>/)
		{
			#~ print "subject_check:$subject_check\t$line**\n";
			push(@{$check_length_subject{$AC_check}},$line);
		}
	}
}

my %hash3=();

$time='['. timestamp(). ']'."\n";
print "Start charging hash4:$time\n";

my $FLAG_line="NaN";

# We open the file with alignments and the sequence mismatches/ loops

if (open(INPUT4, $input4) && open(OUTPUT, '>'.$output))
{
	## Input file= ~/SOFTWARE/BLAST/ncbi-blast-2.2.30+/db/resultsfmt3.out

#~ Query= sp|O43739|ENSG00000008256_ENST00000350796
#~ 
#~ Length=399
                                                                      #~ Score     E
#~ Sequences producing significant alignments:                          (Bits)  Value
#~ 
#~ sp|O43739|ENSG00000008256_UNIPROT  unnamed protein product             769    0.0  
#~ sp|Q9Y4C0|ENSG00000021645_UNIPROT  unnamed protein product            21.6    1.7  
#~ 
#~ 
#~ 
#~ Query_1  27   IRRRKKELIDDIERLKYEIAEVMTEIDNLTSVEESKTTQRNKQIAMGRKKFNMDPKKGIQ  86
#~ O43739   27   ............................................................  86
#~ 
#~ Query_1  87   FLIENDLLQSSPEDVAQFLYKGEGLNKTVIGDYLGERDEFNIKVLQAFVELHEFADLNLV  146
#~ O43739   87   ............................................................  146
#~ 
#~ Query_1  147  QALRQFLWSFRLPGEAQKIDRMMEAFASRYCLCNPGVFQSTDTCYVLSFAIIMLNTSLHN  206
#~ O43739   147  ............................................................  206
#~ Q9Y4C0   863                           .GL.NIIAD.VT.K-.KSS.-..L.TLQAY..M.L  895
#~ 
#~ Query_1  207  HNVRDKPTAERFIAMNRGINEGGDLPEELLRNLYESIKNEPFKIPEDDGNDLTHTFFNPD  266
#~ O43739   207  ............................................................  266
#~ Q9Y4C0   896  FFQFKTTSPDG..LF.S.                                            913
#~ 
#~ Query_1  267  REGWLLKL-GGRVKTWKRRWFILTDNCLYYFEYTTDKEPRGIIPLENLSIREVEDPRKPN  325
#~ O43739   267  ........G...................................................  326
#~ 
#~ Query_1  326  CFELYNPSHKGQVIKACKTEADGRVVEGNHVVYRISAPSPEEKEEWMKSIKASISRDPFY  385
#~ O43739   327  ............................................................  386


	while(my $line=<INPUT4>)
	{
		chomp($line);
		unless($line !~ /\w/)
		{
			#~ print "LINE:$line**\n";
			if($line=~/^Query= (.+)/)
			{
				$query=$1;
			}
			
			# If the alignment corresponds to the selected pair query-subject
			
			if(exists($hash1{$query}))
			{
				#~ print OUTPUT "Query:$query\n";
				$hash3{$query}=1;
				my @tmp=split(/\|/,$query);
				my $AC=$tmp[1];
				
				# Recover the values of length of query and subject from the check hashes
				
				my @check_length_query_tmp=@{$check_length_query{$query}};
				my $check_length_query_tmp_string=join("",@check_length_query_tmp);
				my @check_length_query_tmp_string_residues_tmp=split("",$check_length_query_tmp_string);
				my $length_query=scalar(@check_length_query_tmp_string_residues_tmp);
				
				my @check_length_subject_tmp=@{$check_length_subject{$AC}};
				my $check_length_subject_tmp_string=join("",@check_length_subject_tmp);
				my @check_length_subject_tmp_string_residues_tmp=split("",$check_length_subject_tmp_string);
				my $length_subject=scalar(@check_length_subject_tmp_string_residues_tmp);
				
				
				# If we are in the last line (Effective space search) or in lines displaying the sequence of the alignment (AC(->subject) and query)
				if($line=~/^Effective search space used/ ||$line=~/^$AC/ || $line=~/^Query_/)
				{
					#~ print "SELECTED_ALIGNMENT\n";
					#~ print OUTPUT "LINE:$line**\n";
					#~ print "LINE:$line**\n";
					
					# I we are in the last line of each alignment block of lines
					
					if($line=~/^Effective search space used/)
					{
						my @positions_query=();
						my @positions_subject=();
						
						foreach my $query_tok(sort keys %hash3)
						{
							print OUTPUT ">$query_tok\n";
						}
						foreach my $AC_tok(sort keys%hash2)
						{
							my @index_tmp=sort{$a<=>$b}keys%{$hash2{$AC_tok}};
							#~ print "El array es:@index_tmp\n";
							# We select 1 and 2 as the line indexes for the alignment query-sequence we are going to use (see the problem with multiple repetitions in alignment query-subject in the subject line)
							my $index_query=$index_tmp[0];
							my $index_subject=$index_tmp[1];
							#~ print "Los índices son:$index_query\t$index_subject\n";
							
							# For each of these indexes we extract the begin and end positions of each aligment FASTA line
							
							my @query_begin_tmp=sort{$a<=>$b}keys%{$hash2{$AC_tok}{$index_query}{'begin'}};
							my @subject_begin_tmp=sort{$a<=>$b}keys%{$hash2{$AC_tok}{$index_subject}{'begin'}};
							
							
							# We start navigating the aligment using all the start and end positions for each FASTA line
							
							my $start_subject=$subject_begin_tmp[0];
							my $start_query=$query_begin_tmp[0];
							
							my @query_end_tmp=sort{$a<=>$b}keys%{$hash2{$AC_tok}{$index_query}{'begin'}{$query_begin_tmp[scalar(@query_begin_tmp)-1]}};
							my @subject_end_tmp=sort{$a<=>$b}keys%{$hash2{$AC_tok}{$index_subject}{'begin'}{$subject_begin_tmp[scalar(@subject_begin_tmp)-1]}};
							
							my $finish_subject=$subject_end_tmp[0];
							my $finish_query=$query_end_tmp[0];
							
							#~ if($subject_begin_tmp[0] !=1)
							#~ {
								#~ print "WARNING\t$AC_tok: alignment does not start at the begginning of UNIPROT\n";
							#~ }
							my $counter=$subject_begin_tmp[0];
							#~ print "INICIO_counter=$counter\n";
							
							my @query_tmp=@{$hash2{$AC}{$index_query}{'seq'}};
							my $query_seq=join("",@query_tmp);
							#~ print "query_SEQ:$AC_tok\t$query_seq**\n";
							my @subject_tmp=@{$hash2{$AC}{$index_subject}{'seq'}};
							my $subject_seq=join("",@subject_tmp);
							#~ print "subject_SEQ:$AC_tok\t$subject_seq**\n";
							
							# We store all the aminoacids of each line in a query array and a subject array. As we are going to shift them we define auxiliary hashes (_shifted)
							
							my @aminoacid_query_tmp=split("",$query_seq);
							my @aminoacid_query_tmp_shifted=@aminoacid_query_tmp;
							
							my @aminoacid_subject_tmp=split("",$subject_seq);
							my @aminoacid_subject_tmp_shifted=@aminoacid_subject_tmp;
							
							my $length_query_alignment=scalar(@aminoacid_query_tmp);
							my $length_subject_alignment=scalar(@aminoacid_subject_tmp);
							
							if($AC eq 'P58107'){print "Lonfitud de alineamientos:\tquery:$length_query_alignment\tsubject:$length_subject_alignment\n";}
							
								while(scalar(@aminoacid_subject_tmp_shifted) !=0)
								{
									my $query_contender=shift(@aminoacid_query_tmp_shifted);
									my $subject_contender=shift(@aminoacid_subject_tmp_shifted);
									
									# We extract a compare aminoacids from each hash
									
									if($subject_contender eq '.')
									{
										#~ print OUTPUT "subject eq . :$counter\n";
										# If aa subject eq . means that there is correspondnce query-subject. We move one position in the counter
										$counter++;
										push(@positions_query,$query_contender);
										push(@positions_subject,$subject_contender);
									}
									elsif($subject_contender eq '-' && $query_contender eq '-')
									{
										# Do nothing, these postions are the result of Multiple Alignment and do not exist in query and subject
										
										#~ print OUTPUT "MAL gap:$counter\n";
									}
									elsif($subject_contender ne '-' && $query_contender eq '-')
									{
										# GAP in ENST, record position 
										
										push(@{$hash3{$AC_tok}{'GAP_ENST'}{$counter}},$subject_contender);
										print OUTPUT "GAP_ENST:$counter\t$subject_contender\n";
										$counter++;
										push(@positions_query,$query_contender);
										push(@positions_subject,$subject_contender);
									}
									elsif($subject_contender eq '-' && $query_contender ne '-')
									{
										# GAP in UNIPROT, record position
										
										push(@{$hash3{$AC_tok}{'GAP_UNIPROT'}{$counter}},$query_contender);
										print OUTPUT "GAP_UNIPROT:$counter\t$query_contender\n";
										push(@positions_query,$query_contender);
										push(@positions_subject,$subject_contender);
									}
									
									elsif($subject_contender ne '.' && $subject_contender ne '-' && $query_contender ne '-')
									{
										# 'Missense' between ENSEMBL and UNIPROT order; subject->query
										my $missense=join("__",$subject_contender,$query_contender);
										push(@{$hash3{$AC_tok}{'MISSENSE'}{$counter}},$missense);
										print OUTPUT "Missense:$counter\t$missense\n";
										$counter++;
										push(@positions_query,$query_contender);
										push(@positions_subject,$subject_contender);
									}
									
								}#while2
							if($length_query_alignment == $length_subject_alignment)
							{
								# Do nothing
							}
							else
							{
								# Sometimes alignments are incomplete and we are not going to use those
								
								print OUTPUT "WARNING\t$AC_tok: incomplete_alignment_of_query_or_subject\n";
							}
							my $length_query_alignment_refined=scalar(@positions_query);
							my $length_subject_alignment_refined=scalar(@positions_subject);
							
							# we print the lines with all the fields informative of the alignment.
							
							print OUTPUT "Length:\tquery:$length_query_alignment_refined($length_query;$start_query|$finish_query)\tsubject:$length_subject_alignment_refined($length_subject;$start_subject|$finish_subject)\n";
						}
						
						%hash2=();
						%hash3=();
					}
					
					# If we are in the query line we store in a hash the beggining and end of each row of alignment and the sequence in an array 
					# (for repetitions of lines in fasta format may be destroyed in a hash of hashes)
					elsif($line=~/^Query_[^\s]+\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)/)
					{
						# We use a line code for every line of the alignment, the line starting with the query is always 1
						$FLAG_line=1;
						my $begin=$1;
						my $seq=$2;
						my $end=$3;
						#~ print "$begin\t$seq\t$end\n";
						$hash2{$AC}{$FLAG_line}{'begin'}{$begin}{$end}=1;
						#~ print "$AC\t$FLAG_line\t'begin'\t$begin\t$end\n";
						push(@{$hash2{$AC}{$FLAG_line}{'seq'}},$seq);
						#~ my @tmp2=@{$hash2{$AC}{'query'}{'seq'}};
						#~ print "El array es:@tmp2\n";
					}
					# If we are in the subject line we store in a hash the beggining and end of each row of alignment and the sequence in an array 
					# (for repetitions of lines in fasta format may be destroyed in a hash of hashes)
					elsif($line=~/^$AC+\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)/)
					{
						# We use a line code for every line of the alignment, the line starting with the subject is 1+the query. This is due to the possibility that the same subject
						# aligns in different parts of the same query (for example in porteins with multiple repetitions). The best alignment fot the same subect is the one closest to
						# the query line (in this case line 2)
						
						$FLAG_line++;
						my $begin=$1;
						my $seq=$2;
						my $end=$3;
						#~ print "$begin\t$seq\t$end\n";
						$hash2{$AC}{$FLAG_line}{'begin'}{$begin}{$end}=1;
						#~ print "$AC\t$FLAG_line\t'begin'\t$begin\t$end\n";
						push(@{$hash2{$AC}{$FLAG_line}{'seq'}},$seq);
						#~ my @tmp2=@{$hash2{$AC}{'subject'}{'seq'}};
						#~ print "El array es:@tmp2\n";
					}
					
					
				}
			}
		}
	}
}else {print "impossible to open INPUT4\n";die;}

	
sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
