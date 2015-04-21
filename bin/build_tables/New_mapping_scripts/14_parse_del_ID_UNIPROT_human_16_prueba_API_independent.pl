##	Script to ascertain the correspondence between the UNIPROT displayed sequence and the different ENSEMBL polypetides for the same gene. 2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)


#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;


my $input=$ARGV[0];
my $input2=$ARGV[1];
my $input3=$ARGV[2];
my $input4=$ARGV[3];
my $output1=$ARGV[4];

my %initial_hash=();
my %ID_hash=();

my $time='['. timestamp(). ']'."\n";
	print "Tiempo de incio de carga_hash1:$time\n";

my %check_hash=();

# We open this file to obtain the correspondence SYMBOL-ENSG-strand

if(open(INPUT,$input))
{
	# Input= gtf_output_ENSG.txt
	
	#~ ENSG00000186092 OR4F5   +       69091   70008
	#~ ENSG00000237683 AL627309.1      -       134901  139379
	
	while(my $line=<INPUT>)
	{
		chomp $line;
		my @tmp=split(/\t/,$line);
		my $ENSG=$tmp[0];
		my $SYMBOL=$tmp[1];
		$check_hash{$SYMBOL}{$ENSG}=1;
		#~ print "$SYMBOL\t$ENSG\n";
	}	
}else{print "Unable to open $input\n";}

# We open this file to extract all the UniProt identifiers (AC number) and the corresponding symbol for all the human isoforms displayed in UniProt.

$time='['. timestamp(). ']'."\n";
	print "Tiempo de incio de carga_hash2:$time\n";
if(open(INPUT2,$input2))#&& open (OUTPUT2, '>'.$output2))
{
	# Input= HUMAN.fa (fasta de UNIPROT para human del 28/07/2014) o paralelizado FASTA_chunk_1.fa
	
	#>tr|A0A024R3B9|A0A024R3B9_HUMAN Crystallin, alpha B, isoform CRA_c OS=Homo sapiens GN=CRYAB PE=4 SV=1
	#MRLEKDRFSVNLDVKHFSPEELKVKVLGDVIEVHGKHEERQDEHGFISREFHRKYRIPAD
	#VDPLTITSSLSSDGVLTVNGPRKQVSGPERTIPITREEKPAVTAAPKK
	#>tr|A0A024R7E8|A0A024R7E8_HUMAN Elongation factor 1 homolog (S. cerevisiae), isoform CRA_a OS=Homo sapiens GN=ELOF1 PE=4 SV=1
	#MVRSRLTAVSASWVQAHPPADMGRRKSKRKPPPKKKMTGTLETQFTCPFCNHEKSCDVKM
	
	
	while(my $line=<INPUT2>)
	{
		chomp $line;
		if ($line=~/^>/)
		{
			$line=~s/>//ig;
			## print "$line\n";
			my @tmp=split('\|',$line);
			my $AC=$tmp[1];
			my $description=$tmp[2];
			my $SYMBOL="NaN";
			my @description_tmp=split(/\s+/,$description);
			foreach my $description_tmp_tok(@description_tmp)
			{
				if($description_tmp_tok =~ /GN=(.+)/)
				{
					$SYMBOL=$1;
				}
			}
			foreach my $ENSG_tok(sort keys %{$check_hash{$SYMBOL}})
			{
				$initial_hash{$AC}{$SYMBOL}{$ENSG_tok}=1;
				#~ print "$AC\t$SYMBOL\t$ENSG_tok\n";
			}
		}
	
	}	
}else{print "Unable to open $input2\n";}


my $ID="NaN";
my @AC_tmp=();
my $status="Full";
my %transcript_hash=();
my $FLAG=0;
my @sequence_tmp=();
my %UNIPROT_hash=();


# We open the file carrying all the swiss prot features to extract the sequence of aminoacids corresponding to each displayed isoform of every gene

$time='['. timestamp(). ']'."\n";
	print "Tiempo de incio de carga_hash3:$time\n";
if(open(INPUT3,$input3))
{
	
	# Input2=  uniprot_sprot_human.dat
	
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
#     AISIHGGYDM GFRGLCEAMN AHGTTILKGT MMFDGAMMFD DQGIIPELNC QWRKIRNAFS
#     ETEDVTPLVG KLNSTVFSRV RKFKTLVAFD FINESTMSYV HDWENIKSFL TDQTYSYKGM
#     TYGIERCVIN AGIMTYKIIG VPGMCPPELI RHCIWFPSIK DYVGLKIPAS QDLVEWKTVR
#     ILTSTLRETE EIAMRCYNDK KAWMEQFKVI LGVLSAKSST IVINGMSMQS GERIDINDYH
#     YIGFAILLHT KMKYEQLGKM YDMWNASSIS KWFAALTRPV RVFFSSAVHA LFPTLRPREE
#     KEFLIKLSTF VTFNEECSFD GGEEWDVISS AAYVATQAVT DGKVLAAQKA EKLAEKLAQP
#     VDEVSDSPEV PSSTPDDTAD VCGKEQEVSE LDSLSAQTRS PITRVAERAT AMLEYAAYEK
#     QLHDTTVSNL KRIWNMAGGD DKRNSLEGNL KFVFDTYFTV DPMVNIHFST GRWMRPVPEG
#     IVYSVGYNER GLGPKSDGEL FIVNSECVIC NSESLSAVTR SLQAPTGTIS QVDGVAGCGK
#     TTAIKSIFEP STDMIVTANK KSAQDVRMAL FKSSDSKEAC AFVRTADSVL LNECPTVSRV
#     LVDEVVLLHF GQLCAVMSKL KAVRAICFGD SEQIAFSSRD ASFDMRFSKI IPDETSDADT
#     TFRSPQDVVP LVRLMATKAL PKGTHSKYTK WVSQSKVKRS VTSRSIASVT LVDLDSSRFY
#     ITMTQADKAS LISRAKEMNL PKTFWNERIK TVHESQGISE DHVTLVRLKS TKCDLFKQFS
#     YCLVALTRHK VTFRYEYCGV LNGDLIAECI ARA
#//
	while(my $line=<INPUT3>)
	{
		#print "Hello_world1\n";
		chomp $line;
		if ($line=~/^\/\//)
		{
			my $seq=join("",@sequence_tmp);
			#~ print "$seq\n";
			foreach my $AC_tmp_tok(@AC_tmp)
			{
				if(exists($initial_hash{$AC_tmp_tok}))
				{
					#~ print "candidate:$AC_tmp_tok\n";
					foreach my $SYMBOL_tok (sort keys%{$initial_hash{$AC_tmp_tok}})
					{
						foreach my $ENSG_tok(sort keys %{$initial_hash{$AC_tmp_tok}{$SYMBOL_tok}})
						{
							$UNIPROT_hash{$ENSG_tok}{$AC_tmp_tok}{$seq}=1;
							#~ print"$ENSG_tok\t$AC_tmp_tok\t$seq\n";
						}
					}
				}
			}
		
			$ID="NaN";
			@AC_tmp=();
			$status="Full";
			%transcript_hash=();
			@sequence_tmp=();
			$FLAG=0;
		}
			
		if($line=~/^ID\s+([^\s]+)/)
		{
			$ID=$1;
			#~ print "*********************$ID\n";
		}
		elsif($line=~/^AC\s+(.+)/)
		{
			my $AC=$1;
			my @tmp=split(";",$AC);
			foreach my $tmp_tok(@tmp)
			{
				$tmp_tok=~s/\s+//g;
				push(@AC_tmp,$tmp_tok);
			}
		}
		#~ print "El_array_es:@AC_tmp\n";
		foreach my $AC_tmp_tok(@AC_tmp)
		{
			if(exists($initial_hash{$AC_tmp_tok}))
			{
				if ($line=~/^DE\s+Flags: ([^\;]+)/)
				{
					$status=$1;
					#~ print "STATUS:$status\n";
				}
				unless($status eq 'Fragment' |$status eq 'Fragments')
				{
						# There are proteins with various isoforms and only one of them is displayed.
						if($line=~/^SQ\s+SEQUENCE/)
						{
							$FLAG=1;
							#~ print "I:$line\n";
						}
						if ($FLAG==1 && $line=~/^\s+(.+)\s*$/)
						{
							#~ print "II:$line\n";
							my $sequence_line=$1;
							#~ print "$sequence_line\n";
							$sequence_line=~s/\s+//g;
							push (@sequence_tmp, $sequence_line);
						}
				}
			}
		}
	}# while	
}else{print "Unable to open $input3\n";}

$time='['. timestamp(). ']'."\n";
print "Tiempo de carga del HASH ENSEMBL:$time\n";

my $CHROM="NaN";
my $ENSG="NaN";
my $ENST="NaN";
my $ENSG_biotype="NaN";
my $ENST_biotype="NaN";
my %sequence_hash=();

# We open this file to obtain all the ENSEMBL protein sequences correspondding to all the protein coding trasncripts of each gene

if(open (INPUT4, $input4))
{
	#~ Input=Homo_sapiens.GRCh37.75.pep.all.fa
	
	#~ >ENSP00000455919 pep:known chromosome:GRCh37:HG185_PATCH:39868578:39872221:1 gene:ENSG00000260334 transcript:ENST00000561776 gene_biotype:protein_coding transcript_biotype:protein_coding
	#~ MQRLCVYVLIFALALAAFSEASWKPRSQQPDAPLGTGANRDLELPWLEQQGPASHHRRQL
	#~ GPQGPPHLVADPSKKQGPWLEEEEEAYGWMDFGRRSAEDEN
	#~ >ENSP00000456403 pep:putative chromosome:GRCh37:HG185_PATCH:39873994:39881333:-1 gene:ENSG00000261125 transcript:ENST00000567371 gene_biotype:protein_coding transcript_biotype:protein_coding
	#~ ADIMRGEDFTPAEEFVPQEELGAAKKVPAEEGVMEEAELVSEETEGWEEVELELDEATRM
	#~ NVVTSALEASGLGPSHLDMNYVLQQLANWQDAHYRRQLRWKMLQKEPQPRT
	
	while(my $line=<INPUT4>)
	{
		chomp $line;
		if ($line=~/^>/)
		{
			$line=~s/>//ig;
			#~ print "$line\n";
			my @tmp=split(/\s+/,$line);
			foreach my $tmp_tok(@tmp)
			{
				if($tmp_tok =~ /chromosome:([^:]+)/)
				{
					$CHROM=$1;
					#~ print "El CHROM es:$CHROM\n";
				}
				
				elsif($tmp_tok =~ /gene:(.+)/)
				{
					$ENSG=$1;
				}
				elsif($tmp_tok =~ /transcript:(.+)/)
				{
					$ENST=$1;
				}
				elsif($tmp_tok =~ /gene_biotype:(.+)/)
				{
					$ENSG_biotype=$1;
				}
				elsif($tmp_tok =~ /transcript_biotype:(.+)/)
				{
					$ENST_biotype=$1;
				}
			}
		}
		else
		{
			if($CHROM eq 'GRCh37')
			{
				if(exists($UNIPROT_hash{$ENSG}))
				{
					if($ENSG_biotype eq 'protein_coding' && $ENST_biotype eq 'protein_coding')
					{
						push(@{$sequence_hash{$ENSG}{$ENST}},$line);
					}
				}
			}
		}
	}
}

# Now we match sequences coming from UniProt and sequences coming from ENSEMBL

if(open(OUTPUT, '>'.$output1))
{
	#~ $UNIPROT_hash{$ENSG_tok}{$AC_tmp_tok}{$seq}=1;
	
	foreach my $ENSG_tok(sort keys %UNIPROT_hash)
	{
	foreach my $AC_tok(sort keys %{$UNIPROT_hash{$ENSG_tok}})
	{
	foreach my $UNIPROT_seq_tok(sort keys %{$UNIPROT_hash{$ENSG_tok}{$AC_tok}})
	{
		# We first print the UniProt polypetide sequence
		
		print OUTPUT "UNIPROT>$AC_tok\t$ENSG_tok\n";
		#~ print "UNIPROT>$AC_tok\t$ENSG_tok\n";
		print OUTPUT "$UNIPROT_seq_tok\n";
		#~ print "$UNIPROT_seq_tok\n";
		
		foreach my $ENST_tok(sort keys %{$sequence_hash{$ENSG_tok}})
		{
			my @sequence_tmp=@{$sequence_hash{$ENSG_tok}{$ENST_tok}};
			my $ENST_seq=join("",@sequence_tmp);
			
			# If the any of the ENSEMBL polypetide sequences for the same gene contains the UniProt displayed isoform
			# (in many cases many polypetides of the same gene contain the isoform displayed)
			
			if($ENST_seq =~ /$UNIPROT_seq_tok/)
			{
				# If the sequence not only contains but is exactly equal to the UniProt
				
				if($ENST_seq eq $UNIPROT_seq_tok)
				{
					print OUTPUT "UPSTREAM=0;DOWNSTREAM=0>$AC_tok\t$ENSG_tok\t$ENST_tok\n";
					print OUTPUT "$ENST_seq\n";
							
				}
				
				# If ENSEMBL polypeptide is not equal
				
				elsif ($ENST_seq ne $UNIPROT_seq_tok)
				{
					# If it has a N-ter offset
					
					if($ENST_seq =~ /(.+)$UNIPROT_seq_tok/)
					{
						my $upstream=$1;
						my @upstream_tmp=split("",$upstream);
						print OUTPUT "UPSTREAM=".scalar(@upstream_tmp).";DOWNSTREAM=0>$AC_tok\t$ENSG_tok\t$ENST_tok\n";
						print OUTPUT "$ENST_seq\n";
					}
					
					# If it has N-ter and C-ter offsets
					
					elsif($ENST_seq =~ /(.+)$UNIPROT_seq_tok(.+)/)
					{
						my $upstream=$1;
						my @upstream_tmp=split("",$upstream);
						my $downstream=$2;
						my @downstream_tmp=split("",$downstream);
						print OUTPUT "UPSTREAM=".scalar(@upstream_tmp).";DOWNSTREAM=".scalar(@downstream_tmp).">$AC_tok\t$ENSG_tok\t$ENST_tok\n";
						print OUTPUT "$ENST_seq\n";
					}
					
					# If it has a C-ter offset
					
					elsif($ENST_seq =~ /$UNIPROT_seq_tok(.+)/)
					{
						my $downstream=$1;
						my @downstream_tmp=split("",$downstream);
						print OUTPUT "UPSTREAM=0;DOWNSTREAM=".scalar(@downstream_tmp).">$AC_tok\t$ENSG_tok\t$ENST_tok\n";
						print OUTPUT "$ENST_seq\n";
								
					}
				}
			}
			
			# If the ENSEMBL polypetide does not contain the UniProt isoform
			
			elsif($ENST_seq !~ /$UNIPROT_seq_tok/)
			{
				print OUTPUT "UPSTREAM=NO_ALIGN;DOWNSTREAM=NO_ALIGN>$AC_tok\t$ENSG_tok\t$ENST_tok\n";
				print OUTPUT "$ENST_seq\n";
						
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
