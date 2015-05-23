###############
###############

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;
use Bio::EnsEMBL::Registry;
use Bio::SeqIO;


my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $input3=$ARGV[2];
my $output1=$ARGV[3];
my $output2=$ARGV[4];

my %check_hash=();
my %initial_hash=();
my %ID_hash=();
my $time="NaN";

$time='['. timestamp(). ']'."\n";
	print "Tiempo de incio de carga_hash1:$time\n";
if(open(INPUT,$input1))
{
	# Input= checkII.txt
	
	#~ ENSG00000004139 SARM1
	#~ ENSG00000004142 POLDIP2
	#~ ENSG00000005189 AC004381.6
	#~ ENSG00000006534 ALDH3B1
	#~ ENSG00000026036 RTEL1-TNFRSF6B
	#~ ENSG00000064489 MEF2BNB-MEF2B
	#~ ENSG00000068781 STON1-GTF2A1L
	#~ ENSG00000073169 SELO
	#~ ENSG00000083842 ZNF8
	#~ ENSG00000085365 SCAMP1
	#~ ENSG00000088899 LZTS3
	#~ ENSG00000091436 MLTK
	
	while(my $line=<INPUT>)
	{
		chomp $line;
		my @tmp=split(/\t/,$line);
		my $ENSG=$tmp[0];
		my $SYMBOL=$tmp[1];
		$check_hash{$SYMBOL}{$ENSG}=1;
		print "$SYMBOL\t$ENSG\n";
	}	
}else{print "Unable to open $input1\n";}



$time='['. timestamp(). ']'."\n";
	print "Tiempo de incio de carga_hash2:$time\n";
if(open(INPUT2,$input2) && open (OUTPUT2, '>'.$output2))
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
			my $db=$tmp[0];
			## print "$AC\t$db\t$description\n";
			if(exists($check_hash{$SYMBOL}))
			{
			foreach my $ENSG_tok(sort keys %{$check_hash{$SYMBOL}})
			{
				$initial_hash{$AC}{$SYMBOL}{$ENSG_tok}=1;
				print OUTPUT2 "$AC\t$SYMBOL\t$ENSG_tok\n";
			}
				
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


$time='['. timestamp(). ']'."\n";
	print "Tiempo de incio de carga_hash3:$time\n";
if(open(INPUT3,$input3))
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
print "Tiempo de carga del Registry:$time\n";

my %hash3=();
my $registry = "Bio::EnsEMBL::Registry";
$registry->load_registry_from_url('mysql://anonymous@ensembldb.ensembl.org');
my $gene_adaptor=$registry->get_adaptor( 'Human', 'Core', 'Gene' );

$time='['. timestamp(). ']'."\n";
print "Tiempo de IMPRESIÓN:$time\n";

my $iteration=1;

if(open(OUTPUT, '>'.$output1))
{
	foreach my $ENSG_tok(sort keys %UNIPROT_hash)
	{
	foreach my $AC_tok(sort keys %{$UNIPROT_hash{$ENSG_tok}})
	{	
		print "UNIPROT>$AC_tok\t$ENSG_tok\t$iteration\n";
		$iteration++;
	foreach my $seq_tok(sort keys %{$UNIPROT_hash{$ENSG_tok}{$AC_tok}})
	{
		print OUTPUT "UNIPROT>$AC_tok\t$ENSG_tok\n";
		print OUTPUT "$seq_tok\n";
		if(defined(my $gene=$gene_adaptor->fetch_by_stable_id($ENSG_tok)))
		{
			my $ENSG_gene=$gene->stable_id();
			#~ print "AAAAA:$ENSG_gene:A\n";
			
			my @transcripts = @{ $gene->get_all_Transcripts };
			foreach my $transcript(@transcripts)
			{
				my $ENST_transcript=$transcript->stable_id();
				#~ print "$ENST_transcript\n";
				if(defined($transcript->translation()))
				{
					my $seq=$transcript->translate()->seq();
					$seq=~s/\.|\s+//g;
					#~ print "CANDIDATE:$seq".":::::::::::\n";
					if($seq =~ /$seq_tok/)
					{
						if($seq eq $seq_tok)
						{
							print OUTPUT "UPSTREAM=0;DOWNSTREAM=0>$AC_tok\t$ENSG_gene\t$ENST_transcript\n";
							print OUTPUT "$seq\n";
							
						}
						elsif ($seq ne $seq_tok)
						{
							if($seq =~ /(.+)$seq_tok/)
							{
								my $upstream=$1;
								my @upstream_tmp=split("",$upstream);
								print OUTPUT "UPSTREAM=".scalar(@upstream_tmp).";DOWNSTREAM=0>$AC_tok\t$ENSG_gene\t$ENST_transcript\n";
								print OUTPUT "$seq\n";
							}
							elsif($seq =~ /(.+)$seq_tok(.+)/)
							{
								my $upstream=$1;
								my @upstream_tmp=split("",$upstream);
								my $downstream=$2;
								my @downstream_tmp=split("",$downstream);
								print OUTPUT "UPSTREAM=".scalar(@upstream_tmp).";DOWNSTREAM=".scalar(@downstream_tmp).">$AC_tok\t$ENSG_gene\t$ENST_transcript\n";
								print OUTPUT "$seq\n";
							}
							elsif($seq =~ /$seq_tok(.+)/)
							{
								my $downstream=$1;
								my @downstream_tmp=split("",$downstream);
								print OUTPUT "UPSTREAM=0;DOWNSTREAM=".scalar(@downstream_tmp).">$AC_tok\t$ENSG_gene\t$ENST_transcript\n";
								print OUTPUT "$seq\n";
								
							}
						}
					}
					elsif($seq !~ /$seq_tok/)
					{
						print OUTPUT "UPSTREAM=NO_ALIGN;DOWNSTREAM=NO_ALIGN>$AC_tok\t$ENSG_gene\t$ENST_transcript\n";
						print OUTPUT "$seq\n";
						
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
