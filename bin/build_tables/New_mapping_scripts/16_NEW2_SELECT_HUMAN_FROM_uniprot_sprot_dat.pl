###############
###############

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;


my $input=$ARGV[0];
my $input2=$ARGV[1];
my $input3=$ARGV[2];
my $output1=$ARGV[3];

my %initial_hash=();
my %ID_hash=();
my $time="NaN";

$time='['. timestamp(). ']'."\n";
	print "Tiempo de incio de carga_hash1:$time\n";

my %check_hash=();

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



my @AC_tmp=();
my %hash_line=();
my $counter=0;


$time='['. timestamp(). ']'."\n";
	print "Tiempo de incio de carga_hash3:$time\n";
if(open(INPUT3,$input3) && open(OUTPUT, '>'.$output1))
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
			foreach my $AC_tmp_tok(@AC_tmp)
			{
				if(exists($initial_hash{$AC_tmp_tok}))
				{
					foreach my $counter_tok(sort {$a<=>$b} keys %hash_line)
					{
					foreach my $line_tok(sort keys %{$hash_line{$counter_tok}})
					{
						print OUTPUT "$line_tok\n";
						#~ print "$line_tok\n";
					}
					}
					print OUTPUT "//\n";
					#~ print "//\n";
				}
			}
			@AC_tmp=();
			%hash_line=();
			$counter=0;
		}
			
		else
		{
			$counter++;
			$hash_line{$counter}{$line}=1;
			if($line=~/^AC\s+(.+)/)
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
		}
	}# while	
}else{print "Unable to open $input3\n";}


sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
