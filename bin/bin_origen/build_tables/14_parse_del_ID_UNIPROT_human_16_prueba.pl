###############
###############

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;


my $input=$ARGV[0];
my $input2=$ARGV[1];
my $output1=$ARGV[2];
my $output2=$ARGV[3];

my %initial_hash=();
my %ID_hash=();
my $time="NaN";

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
my %isoform_hash=();
my %transcript_hash=();
my %protein_hash=();
my $FLAG=0;
my $non_experimental_qualifier="NaN";
my %hash_displayed=();
my %hash_non_displayed=();


$time='['. timestamp(). ']'."\n";
	print "Tiempo de incio de carga_hash2:$time\n";
if(open(INPUT2,$input2))
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
		if ($line=~/^\/\//)
		{
			foreach my $AC_tmp_tok(@AC_tmp)
			{
				if(exists($initial_hash{$AC_tmp_tok}))
				{
					foreach my $ENSG_tok (sort keys%{$transcript_hash{$AC_tmp_tok}})
					{
						foreach my $ENST_tok(sort keys %{$transcript_hash{$AC_tmp_tok}{$ENSG_tok}})
						{
						foreach my $ENSP_tok(sort keys %{$transcript_hash{$AC_tmp_tok}{$ENSG_tok}{$ENST_tok}})
						{
							my @AC_isoform_tmp=sort keys %{$transcript_hash{$AC_tmp_tok}{$ENSG_tok}{$ENST_tok}{$ENSP_tok}};
							if ($FLAG==0)
							{
								foreach my $AC_tmp_tok_isoform_tok(@AC_isoform_tmp)
								{
							
									# print "INICIO_UNIQUE:$AC_tmp_tok\t$ENSG_tok\t$ENST_tok\t$ENSP_tok\t$AC_tmp_tok_isoform_tok\n";
									if ($AC_tmp_tok_isoform_tok eq 'NaN'){$AC_tmp_tok_isoform_tok="Unique";}
									my $Iso_displayed="Displayed";
									# print "HELLO WORLD_I:$AC_tmp_tok\t$ENSG_tok\t$ENST_tok\t$ENSP_tok\t$AC_tmp_tok_isoform_tok\t$Iso_displayed:I\n";
									
									foreach my $multiple_keys(sort keys%{$protein_hash{$AC_tmp_tok}})
									{
										if ($multiple_keys eq 'FEATURES')
										{
											foreach my $FT_tok(sort keys%{$protein_hash{$AC_tmp_tok}{'FEATURES'}})
											{
												my @tmp_repeated_feature=sort{ $a <=> $b } keys%{$protein_hash{$AC_tmp_tok}{'FEATURES'}{$FT_tok}};
												for (my $i=0;$i<scalar(@tmp_repeated_feature);$i++)
												{
													foreach my $distance_feature(sort keys%{$protein_hash{$AC_tmp_tok}{'FEATURES'}{$FT_tok}{$tmp_repeated_feature[$i]}})
													{
													foreach my $name(sort keys%{$protein_hash{$AC_tmp_tok}{'FEATURES'}{$FT_tok}{$tmp_repeated_feature[$i]}{$distance_feature}})
													{
														my $counter=$i+1;
														my $FT_counted=join("**",$FT_tok,$counter);
														my $feature=join("__",$name,$tmp_repeated_feature[$i],$distance_feature,$FT_counted);
														$hash_displayed{$ENSG_tok}{$AC_tmp_tok}{$ENST_tok}{$ENSP_tok}{$AC_tmp_tok_isoform_tok}{$Iso_displayed}{'FEATURES'}{$feature}=1;
													}
													}
												}												
											}
										}
										elsif ($multiple_keys eq 'IPR')
										{	
											foreach my $IPR_tok(sort keys%{$protein_hash{$AC_tmp_tok}{'IPR'}})
											{
												foreach my $IPR_description_tok(sort keys %{$protein_hash{$AC_tmp_tok}{'IPR'}{$IPR_tok}})
												{
													# print "IPR=$IPR_tok\t$IPR_description_tok\t";
													$hash_displayed{$ENSG_tok}{$AC_tmp_tok}{$ENST_tok}{$ENSP_tok}{$AC_tmp_tok_isoform_tok}{$Iso_displayed}{'IPR'}{$IPR_tok}{$IPR_description_tok}=1;
												}
											}
										}
										if ($multiple_keys eq 'PFAM')
										{
											foreach my $PFAM_tok(sort keys%{$protein_hash{$AC_tmp_tok}{'PFAM'}})
											{
												foreach my $PFAM_description_tok (sort keys%{$protein_hash{$AC_tmp_tok}{'PFAM'}{$PFAM_tok}})
												{
													foreach my $number_tok (sort keys%{$protein_hash{$AC_tmp_tok}{'PFAM'}{$PFAM_tok}{$PFAM_description_tok}})
													{
														# print "PFAM=$PFAM_tok\t$PFAM_description_tok\t$number_tok\t";
														$hash_displayed{$ENSG_tok}{$AC_tmp_tok}{$ENST_tok}{$ENSP_tok}{$AC_tmp_tok_isoform_tok}{$Iso_displayed}{'PFAM'}{$PFAM_tok}{$PFAM_description_tok}{$number_tok}=1;

													}
												}
											}
										}
										if ($multiple_keys eq 'CCDS')
										{
											foreach my $CCDS_tok(sort keys%{$protein_hash{$AC_tmp_tok}{'CCDS'}})
											{
												foreach my $equivalence_tok (sort keys%{$protein_hash{$AC_tmp_tok}{'CCDS'}{$CCDS_tok}})
												{
													# print "CCDS=$CCDS_tok\t$equivalence_tok\t";
													$hash_displayed{$ENSG_tok}{$AC_tmp_tok}{$ENST_tok}{$ENSP_tok}{$AC_tmp_tok_isoform_tok}{$Iso_displayed}{'CCDS'}{$CCDS_tok}=1;
												}
											}
										}
										if ($multiple_keys eq 'GO')
										{
											foreach my $GO_tok(sort keys%{$protein_hash{$AC_tmp_tok}{'GO'}})
											{
												# print "GO=$GO_tok\t";
												$hash_displayed{$ENSG_tok}{$AC_tmp_tok}{$ENST_tok}{$ENSP_tok}{$AC_tmp_tok_isoform_tok}{$Iso_displayed}{'GO'}{$GO_tok}=1;
											}
										}
										if ($multiple_keys eq 'PROTEIN_EVIDENCE')
										{
											foreach my $protein_existence_tok(sort keys%{$protein_hash{$AC_tmp_tok}{'PROTEIN_EVIDENCE'}})
											{
												# print "PROTEIN_EVIDENCE=$protein_existence_tok\t";
												$hash_displayed{$ENSG_tok}{$AC_tmp_tok}{$ENST_tok}{$ENSP_tok}{$AC_tmp_tok_isoform_tok}{$Iso_displayed}{'PROTEIN_EVIDENCE'}{$protein_existence_tok}=1;
											}
										}
									}
								# print "\n";
								}
							}
							elsif ($FLAG==1)
							{
								foreach my $AC_tmp_tok_isoform_tok(@AC_isoform_tmp)
								{
									# print "INICIO_ISOFORM:$AC_tmp_tok\t$ENSG_tok\t$ENST_tok\t$ENSP_tok\t$AC_tmp_tok_isoform_tok\n";
											if (exists($isoform_hash{$AC_tmp_tok}{$AC_tmp_tok_isoform_tok}))
											{
												foreach my $Iso_displayed_tok(sort keys %{$isoform_hash{$AC_tmp_tok}{$AC_tmp_tok_isoform_tok}})
												{
													if ($Iso_displayed_tok=~/^Displayed/)
													{
														# print "HELLO WORLD_III:$AC_tmp_tok\t$ENSG_tok\t$ENST_tok\t$ENSP_tok\t$AC_tmp_tok_isoform_tok\t$Iso_displayed_tok\n";
														foreach my $multiple_keys(sort keys%{$protein_hash{$AC_tmp_tok}})
														{
															if ($multiple_keys eq 'FEATURES')
															{
																foreach my $FT_tok(sort keys%{$protein_hash{$AC_tmp_tok}{'FEATURES'}})
																{
																	my @tmp_repeated_feature=sort{ $a <=> $b }keys%{$protein_hash{$AC_tmp_tok}{'FEATURES'}{$FT_tok}};
																	for (my $i=0;$i<scalar(@tmp_repeated_feature);$i++)
																	{
																		foreach my $distance_feature(sort keys%{$protein_hash{$AC_tmp_tok}{'FEATURES'}{$FT_tok}{$tmp_repeated_feature[$i]}})
																		{
																		foreach my $name(sort keys%{$protein_hash{$AC_tmp_tok}{'FEATURES'}{$FT_tok}{$tmp_repeated_feature[$i]}{$distance_feature}})
																		{
																			my $counter=$i+1;
																			my $FT_counted=join("**",$FT_tok,$counter);
																			my $feature=join("__",$name,$tmp_repeated_feature[$i],$distance_feature,$FT_counted);
																			$hash_displayed{$ENSG_tok}{$AC_tmp_tok}{$ENST_tok}{$ENSP_tok}{$AC_tmp_tok_isoform_tok}{$Iso_displayed_tok}{'FEATURES'}{$feature}=1;
																		}
																		}
																	}												
																}
															}
															if ($multiple_keys eq 'IPR')
															{	
																foreach my $IPR_tok(sort keys%{$protein_hash{$AC_tmp_tok}{'IPR'}})
																{
																	foreach my $IPR_description_tok(sort keys %{$protein_hash{$AC_tmp_tok}{'IPR'}{$IPR_tok}})
																	{
																		# print "IPR=$IPR_tok\t$IPR_description_tok\t";
																		$hash_displayed{$ENSG_tok}{$AC_tmp_tok}{$ENST_tok}{$ENSP_tok}{$AC_tmp_tok_isoform_tok}{$Iso_displayed_tok}{'IPR'}{$IPR_tok}{$IPR_description_tok}=1;
																		
																	}
																}
															}
															if ($multiple_keys eq 'PFAM')
															{
																foreach my $PFAM_tok(sort keys%{$protein_hash{$AC_tmp_tok}{'PFAM'}})
																{
																	foreach my $PFAM_description_tok (sort keys%{$protein_hash{$AC_tmp_tok}{'PFAM'}{$PFAM_tok}})
																	{
																		foreach my $number_tok (sort keys%{$protein_hash{$AC_tmp_tok}{'PFAM'}{$PFAM_tok}{$PFAM_description_tok}})
																		{
																			# print "PFAM=$PFAM_tok\t$PFAM_description_tok\t$number_tok\t";
																			$hash_displayed{$ENSG_tok}{$AC_tmp_tok}{$ENST_tok}{$ENSP_tok}{$AC_tmp_tok_isoform_tok}{$Iso_displayed_tok}{'PFAM'}{$PFAM_tok}{$PFAM_description_tok}{$number_tok}=1;

																		}
																	}
																}
															}
															if ($multiple_keys eq 'CCDS')
															{
																foreach my $CCDS_tok(sort keys%{$protein_hash{$AC_tmp_tok}{'CCDS'}})
																{
																	foreach my $equivalence_tok (sort keys%{$protein_hash{$AC_tmp_tok}{'CCDS'}{$CCDS_tok}})
																	{
																		# print "CCDS=$CCDS_tok\t$equivalence_tok\t";
																		$hash_displayed{$ENSG_tok}{$AC_tmp_tok}{$ENST_tok}{$ENSP_tok}{$AC_tmp_tok_isoform_tok}{$Iso_displayed_tok}{'CCDS'}{$CCDS_tok}=1;
																	}
																}
															}
															if ($multiple_keys eq 'GO')
															{
																foreach my $GO_tok(sort keys%{$protein_hash{$AC_tmp_tok}{'GO'}})
																{
																	# print "GO=$GO_tok\t";
																	$hash_displayed{$ENSG_tok}{$AC_tmp_tok}{$ENST_tok}{$ENSP_tok}{$AC_tmp_tok_isoform_tok}{$Iso_displayed_tok}{'GO'}{$GO_tok}=1;
																}
															}
															if ($multiple_keys eq 'PROTEIN_EVIDENCE')
															{
																foreach my $protein_existence_tok(sort keys%{$protein_hash{$AC_tmp_tok}{'PROTEIN_EVIDENCE'}})
																{
																	# print "PROTEIN_EVIDENCE=$protein_existence_tok\t";
																	$hash_displayed{$ENSG_tok}{$AC_tmp_tok}{$ENST_tok}{$ENSP_tok}{$AC_tmp_tok_isoform_tok}{$Iso_displayed_tok}{'PROTEIN_EVIDENCE'}{$protein_existence_tok}=1;

																}
															}
														}
													# print "\n";
													}
													elsif ($Iso_displayed_tok !~ /^Displayed/)
													{
														# print "HELLO WORLD_IV:$AC_tmp_tok\t$ENSG_tok\t$ENST_tok\t$ENSP_tok\t$AC_tmp_tok_isoform_tok\t$Iso_displayed_tok\n";
														$hash_non_displayed{$ENSG_tok}{$AC_tmp_tok}{$ENST_tok}{$ENSP_tok}{$AC_tmp_tok_isoform_tok}{$Iso_displayed_tok}=1;												
													}
												}
											}
											else
											{
												my $Iso_displayed_tok="NaN";
												# print "HELLO WORLD_V:$AC_tmp_tok\t$ENSG_tok\t$ENST_tok\t$ENSP_tok\t$AC_tmp_tok_isoform_tok\n";
												$hash_non_displayed{$ENSG_tok}{$AC_tmp_tok}{$ENST_tok}{$ENSP_tok}{$AC_tmp_tok_isoform_tok}{$Iso_displayed_tok}=1;
											}
										}
										
							}
						}
					}
					}
				}
			}
			$ID="NaN";
			@AC_tmp=();
			$status="Full";
			%isoform_hash=();
			%transcript_hash=();
			%protein_hash=();
			$FLAG=0;
			$non_experimental_qualifier="NaN";
		}
			
		if($line=~/^ID\s+([^\s]+)/)
		{
			$ID=$1;
			## print "$ID\n";
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
		foreach my $AC_tmp_tok(@AC_tmp)
		{
			if(exists($initial_hash{$AC_tmp_tok}))
			{
				if ($line=~/^DE\s+Flags: ([^\;]+)/)
				{
					$status=$1;
					## print "$status\n";
				}
				unless($status eq 'Fragment' |$status eq 'Fragments')
				{
						# There are proteins with various isoforms and only one of them is displayed.
							if ($line=~/^CC.+/)
							{
								if($line=~/^.+ALTERNATIVE PRODUCTS.+/)
								{
									# print "$line\n";
									$FLAG=1; 
									## print "La proteína:$AC_tmp_tok\tTIENE ISOFORMAS\n";
								}
								elsif ($line=~/^.+IsoId=([^;]+); Sequence=([^;]+)/ && $FLAG==1)
								{
									my $Iso_ID=$1;
									my $Iso_features=$2;
									$isoform_hash{$AC_tmp_tok}{$Iso_ID}{$Iso_features}=1;
									# print "Has isoforms:$AC_tmp_tok\t$Iso_ID\t$Iso_features\n";
								}
							}
							elsif($line=~/^DR\s+(.+)/)
									{
										my $DR=$1;
										## print "$DR\n";
										if ($DR=~/CCDS; (CCDS[^;]+); (.+)/)
										{
											my $CCDS=$1;
											# # print "$CCDS\n";
											my $equivalence_CCDS_AC=$2;
											$equivalence_CCDS_AC=~s/^\-|\s|\.|\[|\]//g;
											## print "$CCDS\t$equivalence_CCDS_AC\n";
											$protein_hash{$AC_tmp_tok}{'CCDS'}{$CCDS}{$equivalence_CCDS_AC}=1;
										}
										elsif ($DR=~/Ensembl; (ENST[^;]+); (ENSP[^;]+); (ENSG[^\.]+)(.*)/)
										{
												# print "$DR\n";
												my $ENST=$1;
												my $ENSP=$2;
												my $ENSG=$3;
												my $AC_tmp_tok_isoform=$4;
												if ($AC_tmp_tok_isoform !~/\w/){$AC_tmp_tok_isoform="NaN";}
												else{$AC_tmp_tok_isoform=~s/\[|\]|\s|\.//g;}
												
												$transcript_hash{$AC_tmp_tok}{$ENSG}{$ENST}{$ENSP}{$AC_tmp_tok_isoform}=1;
												# print "LINE:$AC_tmp_tok\t$ENSG\t$ENST\t$ENSP\t$AC_tmp_tok_isoform\n";
												
											#}
										}
										elsif ($DR=~/GO; (GO[^;]+)/)
										{
											my $GO_term=$1;
											# # print "$GO_term\n";
											$protein_hash{$AC_tmp_tok}{'GO'}{$GO_term}=1;
										}
										elsif ($DR=~/InterPro; ([^;]+); ([^.]+)/)
										{
											my $IPR=$1;
											my $description_IPR=$2;
											# # print "$IPR\t$description_IPR\n";
											$protein_hash{$AC_tmp_tok}{'IPR'}{$IPR}{$description_IPR}=1;

										}
										
										elsif ($DR=~/Pfam; ([^;]+); ([^;]+); ([^.]+)/)
										{
											my $PFAM=$1;
											my $description_PFAM=$2;
											my $number=$3;
											# # print "$PFAM\t$description_PFAM\t$number\n";
											$protein_hash{$AC_tmp_tok}{'PFAM'}{$PFAM}{$description_PFAM}{$number}=1;

										}
									}
									# We parse this field in order to be able to use it as a feature of selection in the future.
									elsif($line=~/^PE\s+(.+)/)
									{
										my $protein_existence= $1;
										## print "$protein_existence\n";
										$protein_hash{$AC_tmp_tok}{'PROTEIN_EVIDENCE'}{$protein_existence}=1;
									}
									elsif($line=~/^FT\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+(.*)/)
									{
									
										my $FT=$1;
										# # print "$FT\n";
										my $begin=$2;
										# # print "$begin\n";
										my $end=$3;
										# # print "$end\n";
										my $name=$4;if ($name !~/\w/){$name="NaN";}
										if ($FT=~/REGION/|$FT=~/REPEAT/|$FT=~/NP_BIND/|$FT=~/BINDING/|$FT=~/NON_STD/|$FT=~/METAL/|$FT=~/DNA_BIND/|$FT=~/MOD_RES/|$FT=~/LIPID/|$FT=~/CARBOHYD/|$FT=~/DISULFID/|$FT=~/CROSSLNK/|$FT=~/CA_BIND/|$FT=~/ACT_SITE/)
										{
											# For some features coordinates swissprot curation includes uncertainty (? < or >) we discard those coordinates and their features that match non-digit character
											unless ($begin=~/\D/|$end=~/\D/)
											{
												if($name=~/.+\((Pootential)\).+/|$name=~/(^Pootential)\./|$name=~/.+\((Potential)\).+/|$name=~/(^Potential)\./|$name=~/(^Probable)\./|$name=~/.+\((Probable)\).+/|$name=~/(^By similarity)\./|$name=~/.+\((By similarity)\).+/)
												{
													## print "Hello_world_1:$1\n";
													$non_experimental_qualifier=$1;
													## print "AA:$non_experimental_qualifier:AAA\n";
												}else{$non_experimental_qualifier="literature";}
												
												if ($non_experimental_qualifier eq 'literature'|$non_experimental_qualifier eq 'Probable'|$non_experimental_qualifier eq 'By similarity')
													{
														if ($end !~ /\d/){print "The line is*****************************************:$line\n"; print"Theparsed fields are:$FT\t"."AA:$begin:AA\tAA$end"."AA\n";}
														#else{# print "The line is: $FT\t"."HH:$begin:HH\tHH$end"."HH\n";}
														# transform it in distance to the first base of the start codon of the transcript
														my $distance_begin=($begin-1)*3;
														# transform it in distance to the last base of the transcript containing the feature
														my $distance_feature=($end-$begin +1)*3;
														#my $feature=join("__",$name,$distance_begin,$distance_feature,$FT);
														$protein_hash{$AC_tmp_tok}{'FEATURES'}{$FT}{$distance_begin}{$distance_feature}{$name}=1;
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
print "Tiempo de incio de impresión:$time\n";

if(open(OUTPUT, '>'.$output1))
{
foreach my $ENSG_tok(sort keys %hash_displayed)
{
foreach my $AC_tok(sort keys %{$hash_displayed{$ENSG_tok}})
{
foreach my $ENST_tok(sort keys %{$hash_displayed{$ENSG_tok}{$AC_tok}})
{
foreach my $ENSP_tok(sort keys %{$hash_displayed{$ENSG_tok}{$AC_tok}{$ENST_tok}})
{
foreach my $AC_isoform_tok(sort keys %{$hash_displayed{$ENSG_tok}{$AC_tok}{$ENST_tok}{$ENSP_tok}})
{
foreach my $Iso_displayed_tok(sort keys %{$hash_displayed{$ENSG_tok}{$AC_tok}{$ENST_tok}{$ENSP_tok}{$AC_isoform_tok}})
{
print OUTPUT "$ENSG_tok\t$AC_tok\t$ENST_tok\t$ENSP_tok\t$AC_isoform_tok\t$Iso_displayed_tok\t";
	foreach my $feature_tok(sort keys %{$hash_displayed{$ENSG_tok}{$AC_tok}{$ENST_tok}{$ENSP_tok}{$AC_isoform_tok}{$Iso_displayed_tok}{'FEATURES'}})
	{
		print OUTPUT "FEATURE:$feature_tok\t";
	}
	foreach my $IPR_tok(sort keys %{$hash_displayed{$ENSG_tok}{$AC_tok}{$ENST_tok}{$ENSP_tok}{$AC_isoform_tok}{$Iso_displayed_tok}{'IPR'}})
	{
		foreach my $IPR_description_tok(sort keys %{$hash_displayed{$ENSG_tok}{$AC_tok}{$ENST_tok}{$ENSP_tok}{$AC_isoform_tok}{$Iso_displayed_tok}{'IPR'}{$IPR_tok}})
		{
			print OUTPUT "IPR:$IPR_tok"."__"."$IPR_description_tok\t";
		}
	}
	foreach my $PFAM_tok(sort keys %{$hash_displayed{$ENSG_tok}{$AC_tok}{$ENST_tok}{$ENSP_tok}{$AC_isoform_tok}{$Iso_displayed_tok}{'PFAM'}})
	{
		foreach my $PFAM_description_tok(sort keys %{$hash_displayed{$ENSG_tok}{$AC_tok}{$ENST_tok}{$ENSP_tok}{$AC_isoform_tok}{$Iso_displayed_tok}{'PFAM'}{$PFAM_tok}})
		{
			foreach my $PFAM_number_tok(sort keys %{$hash_displayed{$ENSG_tok}{$AC_tok}{$ENST_tok}{$ENSP_tok}{$AC_isoform_tok}{$Iso_displayed_tok}{'PFAM'}{$PFAM_tok}{$PFAM_description_tok}})
			{
				print OUTPUT "PFAM:$PFAM_tok"."__"."$PFAM_description_tok"."__"."$PFAM_number_tok\t";
			}
		}
	}
	foreach my $protein_evidence_tok(sort keys %{$hash_displayed{$ENSG_tok}{$AC_tok}{$ENST_tok}{$ENSP_tok}{$AC_isoform_tok}{$Iso_displayed_tok}{'PROTEIN_EVIDENCE'}})
	{
		print OUTPUT "PROTEIN_EVIDENCE:$protein_evidence_tok\t";
	}
	
	foreach my $GO_tok(sort keys %{$hash_displayed{$ENSG_tok}{$AC_tok}{$ENST_tok}{$ENSP_tok}{$AC_isoform_tok}{$Iso_displayed_tok}{'GO'}})
	{
		print OUTPUT "GO:$GO_tok\t";
	}
	foreach my $CCDS_tok(sort keys %{$hash_displayed{$ENSG_tok}{$AC_tok}{$ENST_tok}{$ENSP_tok}{$AC_isoform_tok}{$Iso_displayed_tok}{'CCDS'}})
	{
		print OUTPUT "CCDS:$CCDS_tok\t";
	}
print OUTPUT "\n";
}
}	
}		
}	
}	
}
}

## Now we print ENSG that do not have a displayed isoform; either the transcripts map to isoforms different from the displayed isoform or transcripts map to no isform (all NaN) although there is one displayed isoform.

if (open(OUTPUT2, '>'.$output2))
{
	foreach my $ENSG_tok(sort keys %hash_non_displayed)
	{
		unless (exists($hash_displayed{$ENSG_tok}))
		{
			foreach my $AC_tok(sort keys %{$hash_non_displayed{$ENSG_tok}})
			{
			foreach my $ENST_tok(sort keys %{$hash_non_displayed{$ENSG_tok}{$AC_tok}})
			{
			foreach my $ENSP_tok(sort keys %{$hash_non_displayed{$ENSG_tok}{$AC_tok}{$ENST_tok}})
			{
			foreach my $AC_isoform_tok(sort keys %{$hash_non_displayed{$ENSG_tok}{$AC_tok}{$ENST_tok}{$ENSP_tok}})
			{
			foreach my $Iso_displayed_tok(sort keys %{$hash_non_displayed{$ENSG_tok}{$AC_tok}{$ENST_tok}{$ENSP_tok}{$AC_isoform_tok}})
			{
				print OUTPUT2 "$ENSG_tok\t$AC_tok\t$ENST_tok\t$ENSP_tok\t$AC_isoform_tok\n";
			}
			}	
			}	
			}	
			}	
		}
	}
}

#$hash_non_displayed{$ENSG_tok}{$AC}{$ENST_tok}{$ENSP_tok}{$AC_isoform_tok}{$Iso_displayed_tok}=1;


$time='['. timestamp(). ']'."\n";
print "Fin del script:$time\n";

sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
