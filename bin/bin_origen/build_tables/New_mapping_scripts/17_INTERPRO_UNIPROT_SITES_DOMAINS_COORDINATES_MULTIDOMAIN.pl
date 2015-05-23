##	Script to select protein features (domain IDs and coordinates and type of functional sites) from UniProt.2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;


my $input=$ARGV[0];
my $input2=$ARGV[1];
my $input3=$ARGV[2];
my $output1=$ARGV[3];
my $output2=$ARGV[4];

my %initial_hash=();
my %ID_hash=();
my $time="NaN";

# From this file we obtain all the UniProt identifiers of human proteins


$time='['. timestamp(). ']'."\n";
	print "Tiempo de incio de carga_hash1:$time\n";
if(open(INPUT,$input))#&& open (OUTPUT2, '>'.$output2))
{
	# Input= HUMAN.fa (fasta de UNIPROT para human del 28/07/2014) o paralelizado FASTA_chunk_1.fa
	
	#>tr|A0A024R3B9|A0A024R3B9_HUMAN Crystallin, alpha B, isoform CRA_c OS=Homo sapiens GN=CRYAB PE=4 SV=1
	#MRLEKDRFSVNLDVKHFSPEELKVKVLGDVIEVHGKHEERQDEHGFISREFHRKYRIPAD
	#VDPLTITSSLSSDGVLTVNGPRKQVSGPERTIPITREEKPAVTAAPKK
	#>tr|A0A024R7E8|A0A024R7E8_HUMAN Elongation factor 1 homolog (S. cerevisiae), isoform CRA_a OS=Homo sapiens GN=ELOF1 PE=4 SV=1
	#MVRSRLTAVSASWVQAHPPADMGRRKSKRKPPPKKKMTGTLETQFTCPFCNHEKSCDVKM
	
	
	while(my $line=<INPUT>)
	{
		chomp $line;
		if ($line=~/^>/)
		{
			$line=~s/>//ig;
			## print "$line\n";
			my @tmp=split('\|',$line);
			my $AC=$tmp[1];
			my $description=$tmp[2];
			my $db=$tmp[0];
			## print "$AC\t$db\t$description\n";
			$initial_hash{$AC}{$db}{$description}=1;
			#print OUTPUT2 "$AC\n";
		}
	
	}	
}else{print "Unable to open $input\n";}


my $ID="NaN";
my @AC_tmp=();
my $status="Full";
my %protein_hash=();
my $non_experimental_qualifier="NaN";
my %hash_displayed=();


# From this file we extract the protein features

$time='['. timestamp(). ']'."\n";
	print "Tiempo de incio de carga_hash2:$time\n";
if(open(INPUT2,$input2) && open (OUTPUT1,'>'.$output1))
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
	while(my $line=<INPUT2>)
	{
		#print "Hello_world1\n";
		chomp $line;
		
		# Print & reinitialize variables when we get to the last line of every record
		
		if ($line=~/^\/\//)
		{
			foreach my $AC_tmp_tok(@AC_tmp)
			{
				if(exists($initial_hash{$AC_tmp_tok}))
				{
					foreach my $IPR_tok(sort keys %{$protein_hash{$AC_tmp_tok}{'IPR'}})
					{
						$hash_displayed{$AC_tmp_tok}{'IPR'}{$IPR_tok}=1;
						#~ print "$AC_tmp_tok\t'IPR'\t$IPR_tok\n";
					}
					foreach my $FT_tok(sort keys %{$protein_hash{$AC_tmp_tok}{'FEATURES'}})
					{
					foreach my $begin_FT_tok(sort keys %{$protein_hash{$AC_tmp_tok}{'FEATURES'}{$FT_tok}})
					{
					foreach my $end_FT_tok(sort keys %{$protein_hash{$AC_tmp_tok}{'FEATURES'}{$FT_tok}{$begin_FT_tok}})
					{
					foreach my $name_FT_tok(sort keys %{$protein_hash{$AC_tmp_tok}{'FEATURES'}{$FT_tok}{$begin_FT_tok}{$end_FT_tok}})
					{
						print OUTPUT1 "$AC_tmp_tok\t$FT_tok\t$begin_FT_tok\t$end_FT_tok\t$name_FT_tok\n";
						#~ print "$AC_tmp_tok\t$FT_tok\t$begin_FT_tok\t$end_FT_tok\t$name_FT_tok\n";
					
					}
					}
					}	
					}
				}
			}
			$ID="NaN";
			@AC_tmp=();
			$status="Full";
			%protein_hash=();
			$non_experimental_qualifier="NaN";
		}
		
		# We obtain the ID
			
		if($line=~/^ID\s+([^\s]+)/)
		{
			$ID=$1;
			## print "$ID\n";
		}
		
		# We obtain UniProt AC number
		
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
		foreach my $AC_tmp_tok(@AC_tmp)
		{
			if(exists($initial_hash{$AC_tmp_tok}))
			{
				if ($line=~/^DE\s+Flags: ([^\;]+)/)
				{
					$status=$1;
					## print "$status\n";
				}
				# We discards records of non-complete proteins
				
				unless($status eq 'Fragment' |$status eq 'Fragments')
				{
					# There are proteins with various isoforms and only one of them is displayed.
							
					# We parse the InterPro domains present in the record
					
					if($line=~/^DR\s+(.+)/)
					{
						my $DR=$1;
						if ($DR=~/InterPro; ([^;]+); ([^.]+)/)
						{
							my $IPR=$1;
							my $description_IPR=$2;
							# # print "$IPR\t$description_IPR\n";
							$protein_hash{$AC_tmp_tok}{'IPR'}{$IPR}=1;
						}
					}
					
					# We parse the level of evidence of protein existence
					
					elsif($line=~/^PE\s+(.+)/)
					{
						my $protein_existence= $1;
						## print "$protein_existence\n";
						$protein_hash{$AC_tmp_tok}{'PROTEIN_EVIDENCE'}{$protein_existence}=1;
					}
					elsif($line=~/^FT\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+(.*)/)
					{
						# We obtain the site features and coordinates
								
						my $FT=$1;
						# # print "$FT\n";
						my $begin=$2;
						# # print "$begin\n";
						my $end=$3;
						# # print "$end\n";
						my $name=$4;if ($name !~/\w/){$name="NaN";}
						
						# We parse the following funtional sites
						
						if ($FT=~/REGION/|$FT=~/REPEAT/|$FT=~/NP_BIND/|$FT=~/BINDING/|$FT=~/NON_STD/|$FT=~/METAL/|$FT=~/DNA_BIND/|$FT=~/MOD_RES/|$FT=~/LIPID/|$FT=~/CARBOHYD/|$FT=~/DISULFID/|$FT=~/CROSSLNK/|$FT=~/CA_BIND/|$FT=~/ACT_SITE/)
						{
							# For some features coordinates swissprot curation includes uncertainty (? < or >) we discard those coordinates and their features that match non-digit character
							unless ($begin=~/\D/|$end=~/\D/)
							{
								# The confidence in the existence of the site is indicated by an indicator. Absence of the indicator equals maximum
								# confidence. When present the indicator varies between potential, probable, and by similarity.
								
								if($name=~/.+\((Pootential)\).+/|$name=~/(^Pootential)\./|$name=~/.+\((Potential)\).+/|$name=~/(^Potential)\./|$name=~/(^Probable)\./|$name=~/.+\((Probable)\).+/|$name=~/(^By similarity)\./|$name=~/.+\((By similarity)\).+/)
								{
									## print "Hello_world_1:$1\n";
									$non_experimental_qualifier=$1;
									## print "AA:$non_experimental_qualifier:AAA\n";
								}else{$non_experimental_qualifier="literature";}
								
								# We discard features whose level of confidence is labelled as potential
											
								if ($non_experimental_qualifier eq 'literature'|$non_experimental_qualifier eq 'Probable'|$non_experimental_qualifier eq 'By similarity')
								{
									if ($end !~ /\d/){print "The line is*****************************************:$line\n"; print"Theparsed fields are:$FT\t"."AA:$begin:AA\tAA$end"."AA\n";}
									$protein_hash{$AC_tmp_tok}{'FEATURES'}{$FT}{$begin}{$end}{$name}=1;
									#~ print "PARSEO:$AC_tmp_tok\t'FEATURES'\t$FT\t$begin\t$end\t$name\n";
								}
							}
						}
					}	
				}
			}
		}
	}# while	
}else{print "Unable to open $input2\n";}

$time='['. timestamp(). ']'."\n";
print "Carga del hash3:$time\n";

my %hash_IPR_coordinates=();

# We open the file with the IPR domain coordinates and the predictor service they come from

if(open(INPUT3,$input3))
{	
	#~ Input=protein2ipr_human.dat
#~ 
	#~ X6RM45	IPR006020	PF00640	19	118
	#~ X6RM45	IPR006020	PS01179	12	140
	#~ X6RM45	IPR011993	G3DSA:2.30.29.30	12	118
	#~ X6RM45	IPR029586	PTHR10337:SF2	19	383
	#~ X6RM59	IPR006434	PF05822	86	331
	#~ X6RM59	IPR006434	PTHR13045	4	331
	#~ X6RM59	IPR006434	TIGR01544	50	330
	#~ X6RM59	IPR023214	G3DSA:3.40.50.1000	80	95
	#~ X6RM59	IPR023214	G3DSA:3.40.50.1000	173	281

while(my $line=<INPUT3>)
	{
		#print "Hello_world1\n";
		chomp $line;
		#~ print "**$line**\n";
		if ($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			#~ print "Hello_world:$line**\n";
			my $AC_INTERPRO=$1;
			my $IPR_INTERPRO=$2;
			my $family=$3;
			my $begin=$4;
			my $end=$5;
			if(exists($hash_displayed{$AC_INTERPRO}{'IPR'}{$IPR_INTERPRO}))
			{
				# We extract all these data
				
				$hash_IPR_coordinates{$AC_INTERPRO}{$IPR_INTERPRO}{$family}{$begin}{$end}=1;
				#~ print "$AC_INTERPRO\t'IPR'\t$IPR_INTERPRO\t$family\t$begin\t$end\n";
			}
		}
	}
}
#~ print "\n";

$time='['. timestamp(). ']'."\n";
print "Merging overlaping domains:$time\n";

if(open(OUTPUT2, '>'.$output2))
{
foreach my $AC_tok(sort keys %hash_IPR_coordinates)
{
	#~ print "1:$AC_tok\n";
foreach my $IPR_INTERPRO_tok(sort keys %{$hash_IPR_coordinates{$AC_tok}})
{
	#~ print "2:$IPR_INTERPRO_tok\n";
	my %hash_merge=();
	foreach my $family_tok(sort keys %{$hash_IPR_coordinates{$AC_tok}{$IPR_INTERPRO_tok}})
	{
		#~ print "3:$family_tok\n";
		
		# For a given domain predicted by a given service of prediction, there is as many repetitions as begin coordinates.
		
		my @begin_tmp=sort {$a<=>$b} keys%{$hash_IPR_coordinates{$AC_tok}{$IPR_INTERPRO_tok}{$family_tok}};
		#~ print "El array es:@begin_tmp\n";
		for(my $i=0; $i<scalar(@begin_tmp);$i++)
		{
			#~ print "4:$begin_tmp[$i]\n";
			
			# Here we number the repetitions of a domain for a given service of prediction
			 
			my $index_family=$i+1;
			my $string=join("__",$family_tok,$index_family);
			
			#~ print "5:$string\n";
			
			foreach my $end_tok(sort keys %{$hash_IPR_coordinates{$AC_tok}{$IPR_INTERPRO_tok}{$family_tok}{$begin_tmp[$i]}})
			{
				#~ print "6:$end_tok\n";
				
				# We define a distance parameter to obtain the positions in the protein covered by each repetition of the domain
				my $distance=$end_tok-$begin_tmp[$i]+1;
				for(my $j=$begin_tmp[$i];$j < $begin_tmp[$i]+$distance;$j++)
				{
					# We create a hash with all the positions
					
					$hash_merge{$j}{$string}=1;
					#~ print "Hash_merge:$IPR_INTERPRO_tok\t$j\t$string\n";
					
				}
			}
		}
	}
	
	# We print a file where we list on a protein position basis the domain positions covered for each domain repetition of each predictor service
	
	foreach my $POS_tok(sort {$a<=>$b} keys %hash_merge)
	{
		print OUTPUT2 "$AC_tok\t$IPR_INTERPRO_tok\t$POS_tok\t";
	foreach my $string_tok(sort keys %{$hash_merge{$POS_tok}})
	{
		print OUTPUT2 "$string_tok\t";
	}	
		print OUTPUT2 "\n";
	}	
}	
}
}



$time='['. timestamp(). ']'."\n";
print "Fin del script:$time\n";

sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
