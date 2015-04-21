##	
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;




my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $output_minimal_representation= $ARGV[2];
#my $output2=$ARGV[3];
my $time="NaN";
my %hash_XML=();
## Orden:  perl ~/Escritorio/Proyecto_clasificador/Scripts_Ad_Hoc/Perl/14_parse_del_ID_UNIPROT_human_3.0.pl HUMAN_CCDC180.fa id_mapping_CCD180.txt sprot_CCD180.txt prueba.txt
$time='['. timestamp(). ']'."\n";
	print "Tiempo de incio de carga_hash1:$time\n";

## Here we parse the XML file of clinvar, we obtain the mode of inheritance, and the identifiers of the variable

if (open(INPUT1, $input1))
{
	##	Input1: ClinVarFullRelease_2014-08.xml
	
	# <?xml version="1.0" encoding="UTF-8" standalone="yes"?>
# <ReleaseSet Dated="2014-08-07" Type="full">

# <ClinVarSet ID="1738117">
  # <RecordStatus>current# </RecordStatus>
  # <Title>NM_015120.4:c.11869+18G&gt;A AND Cardiomyopathy# </Title>
  # <ReferenceClinVarAssertion DateCreated="2012-08-13" DateLastUpdated="2014-04-18" ID="86824">
    # <ClinVarAccession Acc="RCV000029321" Version="1" Type="RCV" DateUpdated="2014-04-19"/>
    # <RecordStatus>current# </RecordStatus>
    # <ClinicalSignificance DateLastEvaluated="2011-08-18">
      # <ReviewStatus>classified by single submitter# </ReviewStatus>
      # <Description>Uncertain significance# </Description>
    # </ClinicalSignificance>
    # <Assertion Type="variation to disease"/>
    # <AttributeSet>
      # <Attribute Type="ModeOfInheritance" integerValue="483">Autosomal unknown </Attribute>
                #<Attribute Type="ModeOfInheritance" integerValue="262">Autosomal dominant inheritance</Attribute>

    # </AttributeSet>
    # <ObservedIn>
      # <Sample>
        # <Origin>germline# </Origin>

my $Title="NaN"; 
my @accesion_tmp=();
my $mode_of_inheritance="NaN";
my @inheritance_tmp=();
my $CHROM="NaN";
my $POS="NaN";
my $REF="NaN";
my $ALT="NaN";

while (my $line = <INPUT1>)
	{
		chomp ($line);
		# Last line of every XML registry. Introduce the elements in the hash
		if ($line=~/\<\/ClinVarSet\>/)
		{
			my $GENE_def="NaN";
			my $NM_transcript_def="NaN";
			foreach my $accesion_def_tok(@accesion_tmp)
			{
				print "HELLO_WORLD_I:@accesion_tmp\n";	
				if($mode_of_inheritance eq 'NaN')
				{
					
					print "hello_world_1:$mode_of_inheritance\n";
					my @title_tmp=split(/\:/,$Title);
					foreach my $title_tmp_tok(@title_tmp)
					{
						print "@title_tmp\n";
						if ($title_tmp_tok=~/(NM_[^\(]+)\((.+)\)/)
						{
							$GENE_def=$2;
							$NM_transcript_def=$1;	
							$hash_XML{$accesion_def_tok}{$GENE_def}{$NM_transcript_def}{$mode_of_inheritance}=1;
							print "$accesion_def_tok\t$GENE_def\t$NM_transcript_def\t$CHROM\t$POS\t$REF\t$ALT\t$mode_of_inheritance\n";
									
						}
						else{$hash_XML{$accesion_def_tok}{$Title}{$mode_of_inheritance}=1;print "$accesion_def_tok\t$Title\t$CHROM\t$POS\t$REF\t$ALT\t$mode_of_inheritance\n";}		
					}
				}
				else
				{
					print "hello_world_2:$mode_of_inheritance\n";
					my $mode_of_inheritance_def=join("|",@inheritance_tmp);
					my @title_tmp=split(/\:/,$Title);
					foreach my $title_tmp_tok(@title_tmp)
					{
						if ($title_tmp_tok=~/(NM_[^\(]+)\((.+)\)/)
							{
								$GENE_def=$2;
								$NM_transcript_def=$1;
								foreach my $accesion_def_tok(@accesion_tmp)
								{
									$hash_XML{$accesion_def_tok}{$GENE_def}{$NM_transcript_def}{$mode_of_inheritance_def}=1;
									print "$accesion_def_tok\t$GENE_def\t$NM_transcript_def\t$CHROM\t$POS\t$REF\t$ALT\t$mode_of_inheritance\n";
								}
							}
							else{$hash_XML{$accesion_def_tok}{$Title}{$mode_of_inheritance_def}=1;print "$accesion_def_tok\t$Title\t$CHROM\t$POS\t$REF\t$ALT\t$mode_of_inheritance_def\n";}
					}
				}
			}
			
			$Title="NaN";
			@accesion_tmp=();
			$mode_of_inheritance="NaN";
			@inheritance_tmp=();
			$CHROM="NaN";
			$POS="NaN";
			$REF="NaN";
			$ALT="NaN";
		}
		# First line of XML. Parse information
		elsif ($line=~/\<Title\>(.+)\<\/Title\>/)
		{
			 #print "La línea es:$line:HH\n";
			 $Title=$1;
			 #print "El título es:$Title:HH\n"; 
		}
		elsif($line=~/\<ClinVarAccession Acc=(.+)/)
		{
			#print "La línea es:$line:AA\n";
			my $Clinvar_accesion=$1;
			my @tmp=split("\"",$Clinvar_accesion);
			#print "###################################################$accesion_tmp[1]\t$accesion_tmp[3]\n";
			my $Clinvar_accesion_def=$tmp[1];
			my $version=$tmp[3];
			my $Accesion_def=join("\.",$Clinvar_accesion_def,$version);
			print "$Accesion_def\n";
			push(@accesion_tmp,$Accesion_def);
			print "El array es:@accesion_tmp\n";
		}
		# Get CHROM, POS, REF and ALT for the assembly 37
		elsif($line=~/\<SequenceLocation Assembly="GRCh37" Chr="(.+)" Accession=".+" start="(.+)" stop=".+" variantLength=".+" referenceAllele="(.+)" alternateAllele="(.+)"\/>/)
		{
			#print "*************************************************************$line\n";
			$CHROM=$1;
			$POS=$2;
			$REF=$3;
			$ALT=$4;
			#print "La variante está en:$CHROM\t$POS\t$REF\t$ALT\n";
		}
		elsif($line=~/\<Attribute Type="ModeOfInheritance".+\>(.+)\<\/Attribute\>/)
		{
			$mode_of_inheritance=$1;
			#print "Hello_world_IMPORTANT:$line\t$mode_of_inheritance\n";
			$mode_of_inheritance=lc($mode_of_inheritance);
			#print "######################################################33$mode_of_inheritance\n";
			if($mode_of_inheritance=~/(.+)inheritance/){$mode_of_inheritance=$1;$mode_of_inheritance=~s/\s+$//g;}
			#print "*****************************************************$mode_of_inheritance\n";
			push(@inheritance_tmp,$mode_of_inheritance);
		}
		
	}# while
}else {print "Impossible to open INPUT1\n";}

$time='['. timestamp(). ']'."\n";
	print "Tiempo de incio de carga_hash2:$time\n";

## Here we parse the vcf file of clinvar, we obtain the mode of the Clinical Significance and all the vcf parameters (Chrom, POS,ID, REF,ALT and INFO)
## WE combine all the information and apply the scheme to get minimal representation applied in the scrip 2_Script_minimal_representation


my %HI0=();
my %IH0=();
my %hash_definitive=();
my $REF_post_trimming="NaN";
my $ALT_post_trimming="NaN";

if(open(INPUT2, $input2) && open(OUTPUT, '>'.$output_minimal_representation))
{
while (my $line = <INPUT2>)
	{
		#input 2: clinvar_20140807.vcf
		# 1       8835# 16  rs267598747     G       A       .       .       RS=267598747;RSPOS=8835# 16;dbSNPBuildID=# 137;SSR=0;SAO=3;VP=0x050060000305000002# 100# 120;GENEINFO=NOC2L:26# 155;WGT=# 1;VC=SNV;PM;REF;SYN;ASP;OTHERKG;LSD;CLNALLE=# 1;CLNHGVS=NC_00000# 1.# 10:g.8835# 16G>A;CLNSRC=ClinVar;CLNORIGIN=2;CLNSRCID=NM_0# 15658.3:c.# 1654C>T;CLNSIG=255;CLNDSDB=MedGen:SNOMED_CT;CLNDSDBID=C0025202:2092003;CLNDBN=Malignant_melanoma;CLNACC=RCV000064926.2
		# 1       89# 1344  rs267598748     G       A       .       .       RS=267598748;RSPOS=89# 1344;dbSNPBuildID=# 137;SSR=0;SAO=3;VP=0x050060000305000002# 100# 120;GENEINFO=NOC2L:26# 155;WGT=# 1;VC=SNV;PM;REF;SYN;ASP;OTHERKG;LSD;CLNALLE=# 1;CLNHGVS=NC_00000# 1.# 10:g.89# 1344G>A;CLNSRC=ClinVar;CLNORIGIN=2;CLNSRCID=NM_0# 15658.3:c.657C>T;CLNSIG=255;CLNDSDB=MedGen:SNOMED_CT;CLNDSDBID=C0025202:2092003;CLNDBN=Malignant_melanoma;CLNACC=RCV000064927.2
		# 1       906# 168  rs267598759     G       A       .       .       RS=267598759;RSPOS=906# 168;dbSNPBuildID=# 137;SSR=0;SAO=3;VP=0x050060080a05000002# 100# 120;GENEINFO=PLEKHN# 1:84069;WGT=# 1;VC=SNV;PM;NSM;REF;INT;ASP;OTHERKG;LSD;CLNALLE=# 1;CLNHGVS=NC_00000# 1.# 10:g.906# 168G>A;CLNSRC=ClinVar;CLNORIGIN=2;CLNSRCID=NM_00# 1# 160# 184.# 1:c.484+30G>A;CLNSIG=255;CLNDSDB=MedGen:SNOMED_CT;CLNDSDBID=C0025202:2092003;CLNDBN=Malignant_melanoma;CLNACC=RCV000064940.2
		# 1       955597  rs115173026     G       T       .       .       RS=115173026;RSPOS=955597;dbSNPBuildID=132;SSR=0;SAO=0;VP=0x050000000305170016000100;GENEINFO=AGRN:375790;WGT=1;VC=SNV;REF;SYN;ASP;VLD;G5A;G5;KGPhase1;KGPROD;OTHERKG;CLNALLE=1;CLNHGVS=NC_000001.10:g.955597G>T;CLNSRC=ClinVar|.|University_of_Chicago;CLNORIGIN=1;CLNSRCID=NM_198576.3:c.45G>T|.|NM_198576.2(AGRN):c.45G>;CLNSIG=3;CLNDSDB=.;CLNDSDBID=.;CLNDBN=AllHighlyPenetrant;CLNACC=RCV000116272.2;CAF=[0.6938,0.3062];COMMON=1


##	From input2 file we load HI & IH hashes with the heather fields.

		chomp ($line);
		if (index($line,'#CHROM')==0) 
				
			{
				$line=~/^#(.+)/;
				my @tmp_HI0 = split ("\t",$1);
				for(my $i=0;$i<scalar(@tmp_HI0);$i++)
										
				{
					$HI0{$tmp_HI0[$i]}=$i;
					$IH0{$i}=$tmp_HI0[$i];
					#print "$HI0{$tmp_HI0[$i]} \n";
					#print "$IH0{$i} \n";
				}
			}

##	From lines carrying data we store every column in an array
						
		elsif ($line !~ /^#/)
			{	
				##	Here we define and re-initialize the arrays on a per-line basis.
				my %duplicate_genes=();
				my %hash_REF_ALT=();
				my %hash_transcripts=();
				my %hash_1=();
				my %hash_2=();
				my %hash_3=();
				my %hash_INFO=();
				my $trim="NaN";
				my $pos_trimmed="NaN";
				my $clinvar_accesion="NaN";
				my $version="NaN";

				#print "$line\n";
				my @tmp= split(" ",$line);
				
				my $CHROM="NaN";if (exists($HI0{"CHROM"})){$CHROM=$tmp[$HI0{"CHROM"}];}unless(defined($CHROM)){$CHROM="NaNein"; print "ERROR in CHROM\n";}
				my $POS="NaN";if (exists($HI0{"POS"})){$POS=$tmp[$HI0{"POS"}];}unless(defined($POS)){$POS="NaNein"; print "ERROR in POS\n";}
				my $ID="NaN";if (exists($HI0{"ID"})){$ID=$tmp[$HI0{"ID"}];}unless(defined($ID)){$ID="NaNein"; print "ERROR in ID\n";}
				my $REF="NaN";if (exists($HI0{"REF"})){$REF=$tmp[$HI0{"REF"}];}unless(defined($REF)){$REF="NaNein"; print "ERROR in REF\n";}
				my $ALT="NaN";if (exists($HI0{"ALT"})){$ALT=$tmp[$HI0{"ALT"}];}unless(defined($ALT)){$ALT="NaNein"; print "ERROR in ALT\n";}
				my $QUAL="NaN";if (exists($HI0{"QUAL"})){$QUAL=$tmp[$HI0{"QUAL"}];}if($QUAL=~/\./){$QUAL=100;}
				my $FILTER="NaN";if (exists($HI0{"FILTER"})){$FILTER=$tmp[$HI0{"FILTER"}];}if($FILTER=~/\./){$FILTER="PASS";}
				my $INFO="NaN";if (exists($HI0{"INFO"})){$INFO=$tmp[$HI0{"INFO"}];}unless(defined($INFO)){$INFO="NaNein"; print "ERROR in INFO\n";}
				my $FORMAT="NaN";if (exists($HI0{"FORMAT"})){$FORMAT=$tmp[$HI0{"FORMAT"}];}unless(defined($FORMAT)){$FORMAT="NaNein"; print "ERROR in FORMAT\n";}
				my $NA12878="NaN";if (exists($HI0{"NA12878"})){$NA12878=$tmp[$HI0{"NA12878"}];}unless(defined($NA12878)){$NA12878="NaNein"; print "ERROR in NA12878\n";}
				
				my @INFO_tmp=split(";",$INFO);
				my @clinalle_tmp=();
				my @clinacc_tmp=();
				my @clinacc_tmp_2=();
				my @genes_tmp=();
				my @transcripts_tmp=();
				my @transcripts_tmp2=();
				my @clinsig_tmp=();
				{
					foreach my $INFO_tmp_tok(@INFO_tmp)
					{
						#if ($INFO_tmp_tok=~/GENEINFO=([^:]+)\:.+/)
						if ($INFO_tmp_tok=~/GENEINFO=(.+)/)
						{
							my $SYMBOL=$1;
							# print"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIEL GEN ES:$SYMBOL\n";
							push(@genes_tmp,$SYMBOL);
							# print"EL ARRAY DE GENES ES:@genes_tmp\n";
						}
						elsif($INFO_tmp_tok=~/CLNSRCID=(.+)/)
						{
							my $field=$1;
							## print"El tránscrito es:$field\n";
							if ($field=~/\,/)
							{
								#my @fields_derived_tmp=split("\,",$field);
								@transcripts_tmp=split("\,",$field);
							}
							else
							{
								push(@transcripts_tmp,$field);
							}
						}
						elsif($INFO_tmp_tok=~/CLNACC=(.+)/)
						{
							$clinvar_accesion=$1;
							## print"#####################################$clinvar_accesion\n";
							if($clinvar_accesion=~/\,/)
							{
								@clinacc_tmp=split("\,",$clinvar_accesion);
							}
							else
							{
								push(@clinacc_tmp,$clinvar_accesion);
							}
						}
						elsif ($INFO_tmp_tok=~/CLNALLE=(.+)/)
						{
							my $clinvar_accesion=$1;
							## print"EL CLNALLE ES:$clinvar_accesion\n";
							if($clinvar_accesion=~/\,/)
							{
								@clinalle_tmp=split("\,",$clinvar_accesion);
							}
							else
							{
								push(@clinalle_tmp,$clinvar_accesion);
							}
							#if ($INFO_tmp_tok=~/\,/){# print"CLINALLE:$CHROM\t$POS\t$ID\t$REF\t$ALT\t$FILTER\t@clinalle_tmp\n";}
						}
						elsif ($INFO_tmp_tok=~/CLNSIG=(.+)/)
						{
							my $CLNSIG=$1;
							# print"EL CLNSIG ES:$CLNSIG\n";
							if($CLNSIG=~/\,/)
							{
								@clinsig_tmp=split("\,",$CLNSIG);
							}
							else
							{
								push(@clinsig_tmp,$CLNSIG);
							}
							#if ($INFO_tmp_tok=~/\,/){# print"CLINSIG:$CHROM\t$POS\t$ID\t$REF\t$ALT\t$FILTER\t@clinsig_tmp\n";}
						}
					}
				}
				my @ALT_tmp=split(',',$ALT);
				
				foreach my $GENE_tmp_tok(@genes_tmp)
				{
					## print"Hello_world_1:$GENE_tmp_tok\n";
					foreach my $clinalle_tmp_tok(@clinalle_tmp)
					{
						$hash_1{$GENE_tmp_tok}{$clinalle_tmp_tok}=1;
						# print"Hello_world_2\t$GENE_tmp_tok\t$clinalle_tmp_tok\n";
					}
				}
				for(my$i=0;$i<scalar(@clinalle_tmp);$i++)
				{
					$hash_2{$clinalle_tmp[$i]}{$clinsig_tmp[$i]}=1;
					# print"Hello_world3:$clinalle_tmp[$i]\t$clinsig_tmp[$i]\n";
				}
				for(my$i=0;$i<scalar(@clinalle_tmp);$i++)
				{
					$hash_3{$clinalle_tmp[$i]}{$clinacc_tmp[$i]}=1;
					# print"Hello_world17:$clinalle_tmp[$i]\t$clinacc_tmp[$i]\n";
				}
				for(my$i=0;$i<scalar(@clinalle_tmp);$i++)
				{
					$hash_transcripts{$clinalle_tmp[$i]}{$transcripts_tmp[$i]}=1;
					# print"Hello_world18:$clinalle_tmp[$i]\t$transcripts_tmp[$i]\n";
				}
				for(my$i=0;$i<scalar(@ALT_tmp);$i++)
				{
					my $index=$i+1;
					$hash_REF_ALT{$index}{$ALT_tmp[$i]}=1;
					# print"Hello_world4:$index\t$ALT_tmp[$i]\n";
				}
				
				foreach my $GENE_tok(sort keys %hash_1)
				{
					print "IIIIIIIIIIIIIIIIIII:$GENE_tok\t";
					foreach my $clinalle_tok(sort keys %{$hash_1{$GENE_tok}})
						{
							print "$clinalle_tok\t";
							my @clinacc_tmp=();
							foreach my $transcript_tok(sort keys %{$hash_transcripts{$clinalle_tok}})
							{
								print "$transcript_tok\t";
								if($transcript_tok=~/\|/)
									{
										@transcripts_tmp2=split(/\|/,$transcript_tok);
									}
								else{push(@transcripts_tmp2,$transcript_tok);}
								foreach my $clinacc_tok(sort keys%{$hash_3{$clinalle_tok}})
								{
									print "$clinacc_tok\t";
									foreach my $clinsig_tok(sort keys%{$hash_2{$clinalle_tok}})
									{
										print "$clinsig_tok\t";
										foreach my $ALT_tok(sort keys%{$hash_REF_ALT{$clinalle_tok}})
										{
											print "$ALT_tok\n";
											if ($GENE_tok=~/\|/)
											{
												print"BANDIDO:*******************************************************$clinacc_tok\t$GENE_tok\n";
												my @GENE_tok_tmp=split(/\|/,$GENE_tok);
												foreach my $GENE_tok_tmp_tok(@GENE_tok_tmp)
												{
													$GENE_tok_tmp_tok=~s/:.+//g;
													
													foreach my $transcript_tmp2_tok(@transcripts_tmp2)
													{
														if($transcript_tmp2_tok=~/^(NM_[^\:]+).+/)
														{
															$transcript_tmp2_tok=~s/:.+//g;
															if($clinacc_tok=~/\|/)
															{
																@clinacc_tmp_2=split(/\|/,$clinacc_tok);
																#print"El array es:@clinacc_tmp\n";
															}
															else{push(@clinacc_tmp_2,$clinacc_tok);}
															print "OOOOOOOOOOOOOOOOOOOOOOOOO$clinacc_tok\t$GENE_tok_tmp_tok\t$transcript_tmp2_tok\n";
															foreach my $clinacc_tmp_tok(@clinacc_tmp_2)
															{
																print"BAIT macheo:$clinacc_tmp_tok\t$GENE_tok_tmp_tok\t$transcript_tmp2_tok\n";
																if (exists($hash_XML{$clinacc_tmp_tok}{$GENE_tok_tmp_tok}{$transcript_tmp2_tok}))
																{
																		foreach my $mode_of_inheritance_tok(sort keys %{$hash_XML{$clinacc_tmp_tok}{$GENE_tok_tmp_tok}{$transcript_tmp2_tok}})
																		{
																		print "Hello_world_25:$clinacc_tmp_tok\t$GENE_tok_tmp_tok\t$transcript_tok\n";
																			if (length($REF)==1 && length($ALT_tok)==1)
																			{
																										#print  "El alelo partido por comas es y que cumple las NOMODIF es:$ALT_tok\n";
																										$trim= "No_modification_to_get_minimal_representation";
																										my $INFO_clinvar=join(";",$clinsig_tok,$GENE_tok_tmp_tok,$clinacc_tok,$mode_of_inheritance_tok,$trim);
																										$hash_definitive{$CHROM}{$POS}{$ID}{$REF}{$ALT_tok}{$QUAL}{$FILTER}{$INFO_clinvar}=1;
																										print  "$CHROM\t$POS\t$ID\t$REF\t$ALT_tok\t$QUAL\t$FILTER\t$INFO_clinvar\n";	
																			}
																			elsif (length($REF) == 1 && length($ALT_tok) != 1)
																			{
																										#print  "El alelo partido por comas es y que cumple las NOMODIF es:$ALT_tok\n";
																										$trim= "No_modification_to_get_minimal_representation";
																										my $INFO_clinvar=join(";",$clinsig_tok,$GENE_tok_tmp_tok,$clinacc_tok,$mode_of_inheritance_tok,$trim);
																										$hash_definitive{$CHROM}{$POS}{$ID}{$REF}{$ALT_tok}{$QUAL}{$FILTER}{$INFO_clinvar}=1;
																										print  "$CHROM\t$POS\t$ID\t$REF\t$ALT_tok\t$QUAL\t$FILTER\t$INFO_clinvar\n";	
																			}
																			elsif (length($REF) != 1 && length($ALT_tok) == 1)
																			{
																										#print  "El alelo partido por comas es y que cumple las NOMODIF es:$ALT_tok\n";
																										$trim= "No_modification_to_get_minimal_representation";
																										my $INFO_clinvar=join(";",$clinsig_tok,$GENE_tok_tmp_tok,$clinacc_tok,$mode_of_inheritance_tok,$trim);
																										$hash_definitive{$CHROM}{$POS}{$ID}{$REF}{$ALT_tok}{$QUAL}{$FILTER}{$INFO_clinvar}=1;
																										print  "$CHROM\t$POS\t$ID\t$REF\t$ALT_tok\t$QUAL\t$FILTER\t$INFO_clinvar\n";	
																			}
																			elsif(length($REF) > 1 && length($ALT_tok) > 1)
																			{
																										##	We split by letters the components of REF and ALT
																										
																										#print  "$REF\t$ALT_tok\n";
																										my @REF_tmp=split('',$REF);
																										#print  "@REF_tmp\n";
																										my @ALT_tmp=split('',$ALT_tok);
																										#print  "@ALT_tmp\n";
																										my $POS_trimmed="NaN";
																										
																										## Conditions of trimming downstream; same element in last position of REF and ALT and REF and ALT more than one in length. 
																										
																										while (($REF_tmp[scalar(@REF_tmp)-1] eq $ALT_tmp[scalar(@ALT_tmp)-1] && scalar(@REF_tmp)>1 && scalar(@ALT_tmp)>1))
																											{
																												## Extract last element of the array.
																												#print  "Hello_world_1\n";
																												pop (@REF_tmp);
																												pop (@ALT_tmp);
																											}
																										
																										## Conditions of trimming upstream; same element in last position of REF and ALT and REF and ALT more than one in length.
																										
																										while ($REF_tmp[0] eq $ALT_tmp[0] && scalar(@REF_tmp)>1 && scalar(@ALT_tmp)>1)
																											{
																												#print  "Hello_world_2\n";
																												shift (@REF_tmp);
																												shift (@ALT_tmp);
																												#print  "$POS\n";
																												$POS++;
																												#print  "$POS\n";
																											}
																										
																											$REF_post_trimming=join('',@REF_tmp);
																											$ALT_post_trimming=join('',@ALT_tmp);
																											$trim= "This_line_has_been_modified_to_get_minimal_representation";
																											my $INFO_clinvar=join(";",$clinsig_tok,$GENE_tok_tmp_tok,$clinacc_tok,$mode_of_inheritance_tok,$trim);
																											$hash_definitive{$CHROM}{$POS}{$ID}{$REF_post_trimming}{$ALT_post_trimming}{$QUAL}{$FILTER}{$INFO_clinvar}=1;
																											print  "$CHROM\t$POS\t$ID\t$REF_post_trimming\t$ALT_post_trimming\t$QUAL\t$FILTER\t$INFO_clinvar\n";
																					}
																					
																				}	
																			}
																		}
																	}
																}
																
															
														}
											}
											else
											{
												$GENE_tok=~s/:.+//g;
												foreach my $transcript_tmp2_tok(@transcripts_tmp2)
													{
														if($transcript_tmp2_tok=~/^(NM_[^\:]+).+/)
														{
															$transcript_tmp2_tok=~s/:.+//g;
															if($clinacc_tok=~/\|/)
															{
																@clinacc_tmp_2=split(/\|/,$clinacc_tok);
																#print"El array es:@clinacc_tmp\n";
															}
															else{push(@clinacc_tmp_2,$clinacc_tok);}
															print "OOOOOOOOOOOOOOOOOOOOOOOOO$clinacc_tok\t$GENE_tok\t$transcript_tmp2_tok\n";
															foreach my $clinacc_tmp_tok(@clinacc_tmp_2)
															{
																print"BAIT macheo:$clinacc_tmp_tok\t$GENE_tok\t$transcript_tmp2_tok\n";
																if (exists($hash_XML{$clinacc_tmp_tok}{$GENE_tok}{$transcript_tmp2_tok}))
																{
																		foreach my $mode_of_inheritance_tok(sort keys %{$hash_XML{$clinacc_tmp_tok}{$GENE_tok}{$transcript_tmp2_tok}})
																		{
																		print "Hello_world_26:$clinacc_tmp_tok\t$GENE_tok\t$transcript_tok\n";	
																		if (length($REF)==1 && length($ALT_tok)==1)
																		{
																			#print  "El alelo partido por comas es y que cumple las NOMODIF es:$ALT_tok\n";
																			$trim= "No_modification_to_get_minimal_representation";
																			my $INFO_clinvar=join(";",$clinsig_tok,$GENE_tok,$clinacc_tok,$mode_of_inheritance_tok,$trim);
																			$hash_definitive{$CHROM}{$POS}{$ID}{$REF}{$ALT_tok}{$QUAL}{$FILTER}{$INFO_clinvar}=1;
																			print  "$CHROM\t$POS\t$ID\t$REF\t$ALT_tok\t$QUAL\t$FILTER\t$INFO_clinvar\n";	
																		}
																		elsif (length($REF) == 1 && length($ALT_tok) != 1)
																		{
																			#print  "El alelo partido por comas es y que cumple las NOMODIF es:$ALT_tok\n";
																			$trim= "No_modification_to_get_minimal_representation";
																			my $INFO_clinvar=join(";",$clinsig_tok,$GENE_tok,$clinacc_tok,$mode_of_inheritance_tok,$trim);
																			$hash_definitive{$CHROM}{$POS}{$ID}{$REF}{$ALT_tok}{$QUAL}{$FILTER}{$INFO_clinvar}=1;
																			print  "$CHROM\t$POS\t$ID\t$REF\t$ALT_tok\t$QUAL\t$FILTER\t$INFO_clinvar\n";	
																		}
																		elsif (length($REF) != 1 && length($ALT_tok) == 1)
																		{
																			#print  "El alelo partido por comas es y que cumple las NOMODIF es:$ALT_tok\n";
																			$trim= "No_modification_to_get_minimal_representation";
																			my $INFO_clinvar=join(";",$clinsig_tok,$GENE_tok,$clinacc_tok,$mode_of_inheritance_tok,$trim);
																			$hash_definitive{$CHROM}{$POS}{$ID}{$REF}{$ALT_tok}{$QUAL}{$FILTER}{$INFO_clinvar}=1;
																			print  "$CHROM\t$POS\t$ID\t$REF\t$ALT_tok\t$QUAL\t$FILTER\t$INFO_clinvar\n";	
																		}
																		elsif(length($REF) > 1 && length($ALT_tok) > 1)
																		{
																			##	We split by letters the components of REF and ALT
																													
																			#print  "$REF\t$ALT_tok\n";
																			my @REF_tmp=split('',$REF);
																			#print  "@REF_tmp\n";
																			my @ALT_tmp=split('',$ALT_tok);
																			#print  "@ALT_tmp\n";
																			my $POS_trimmed="NaN";
																													
																			## Conditions of trimming downstream; same element in last position of REF and ALT and REF and ALT more than one in length. 
																													
																			while (($REF_tmp[scalar(@REF_tmp)-1] eq $ALT_tmp[scalar(@ALT_tmp)-1] && scalar(@REF_tmp)>1 && scalar(@ALT_tmp)>1))
																			{
																				## Extract last element of the array.
																				#print  "Hello_world_1\n";
																				pop (@REF_tmp);
																				pop (@ALT_tmp);
																			}										
																			## Conditions of trimming upstream; same element in last position of REF and ALT and REF and ALT more than one in length.
																													
																			while ($REF_tmp[0] eq $ALT_tmp[0] && scalar(@REF_tmp)>1 && scalar(@ALT_tmp)>1)
																			{
																				#print  "Hello_world_2\n";
																				shift (@REF_tmp);
																				shift (@ALT_tmp);
																				#print  "$POS\n";
																				$POS++;
																				#print  "$POS\n";
																			}
																													
																			$REF_post_trimming=join('',@REF_tmp);
																			$ALT_post_trimming=join('',@ALT_tmp);
																			$trim= "This_line_has_been_modified_to_get_minimal_representation";
																			my $INFO_clinvar=join(";",$clinsig_tok,$GENE_tok,$clinacc_tok,$mode_of_inheritance_tok,$trim);
																			$hash_definitive{$CHROM}{$POS}{$ID}{$REF_post_trimming}{$ALT_post_trimming}{$QUAL}{$FILTER}{$INFO_clinvar}=1;
																			print  "$CHROM\t$POS\t$ID\t$REF_post_trimming\t$ALT_post_trimming\t$QUAL\t$FILTER\t$INFO_clinvar\n";	
																		}
																	}
																}
															}
														}
													}
											}
										}
									}
								}
							}
						}
				}
			}
	}#WHILE
		
}else {print "unable to open $input2 or $output_minimal_representation\n";}

my %hash_no_repetitions=();
#my %hash_no_repetitions_2=();


foreach my $CHROM_tok(sort keys%hash_definitive)
{
foreach my $POS_tok(sort keys%{$hash_definitive{$CHROM_tok}})
{
foreach my $ID_tok(sort keys%{$hash_definitive{$CHROM_tok}{$POS_tok}})
{
foreach my $REF_tok(sort keys%{$hash_definitive{$CHROM_tok}{$POS_tok}{$ID_tok}})
{
foreach my $ALT_tok(sort keys%{$hash_definitive{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}})
{
foreach my $QUAL_tok(sort keys%{$hash_definitive{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}})
{
foreach my $FILTER_tok(sort keys%{$hash_definitive{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$QUAL_tok}})
{
foreach my $INFO_tok(sort keys%{$hash_definitive{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$QUAL_tok}{$FILTER_tok}})
{
	print OUTPUT "$CHROM_tok\t$POS_tok\t$ID_tok\t$REF_tok\t$ALT_tok\t$QUAL_tok\t$FILTER_tok\t$INFO_tok\n";
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
