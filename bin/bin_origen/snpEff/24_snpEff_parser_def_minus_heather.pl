###############
###############

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;

## Declare the array of accepted chromosomes. We won't accept variants in patches.

my @AcceptedChrArray=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");


## Create the hash for accepted chromosomes.

my %AcceptedChr=();
foreach my $AcceptedChrTmp(@AcceptedChrArray){$AcceptedChr{"$AcceptedChrTmp"}=1;}


## Create the hash with the effects we are going to accept. These are: missense_variant, stop_gained, frameshift_variant, stop_retained_variant, start_retained_variant,
## synonymous_variant,splice_acceptor_variant, splice_donor_variant, splice_region_variant, inframe_deletion, disruptive_inframe_deletion, exon_loss_variant,
## and start_lost.


my %AcceptedEffects=();
$AcceptedEffects{"splice_acceptor_variant"}=1;
$AcceptedEffects{"splice_donor_variant"}=1;
#$AcceptedEffects{"splice_region_variant"}=1;

#~ $AcceptedEffects{"exon_loss_variant"}=1;
#~ $AcceptedEffects{"transcript_ablation"}=1;

#~ $AcceptedEffects{"inframe_deletion"}=1;
#~ $AcceptedEffects{"disruptive_inframe_deletion"}=1;

$AcceptedEffects{"missense_variant"}=1; # missense_variant    (CURRENT_SVN) SO Accession:	SO:0001583 ANNOVAR:nonsynonymous SNV, VAAST:non_synonymous_codon, missense, missense codon, vep:NON_SYSNONYMOUS_CODING, SO:0001584, SO:0001783

$AcceptedEffects{"stop_gained"}=1; # término SO: stop_gained    (CURRENT_SVN) SO Accession:	SO:0001587       ANNOVAR:stopgain, nonsense, nonsense codon, vep:STOP_GAINED, stop gained, VAAST:stop_gained

$AcceptedEffects{"frameshift_variant"}=1; # término SO: frameshift_variant    (CURRENT_SVN) SO Accession:	SO:0001589  ANNOVAR:frameshift block substitution, frameshift variant, frameshift_, frameshift_coding, vep:FRAME_SHIFT, VAAST:frameshift_variant


#~ $AcceptedEffects{"stop_retained_variant"}=1; # término SO: stop_retained_variant  SO Accession:	SO:0001567  (CURRENT_SVN) vep:SYNONYMOUS_STOP, stop retained variant, VAAST:stop_retained
#~ $AcceptedEffects{"start_retained_variant"}=1; # término SO: start_retained_variant    (CURRENT_SVN) SO Accession:	SO:0002019 synonymous: vep??????????
$AcceptedEffects{"synonymous_variant"}=1; # término SO: synonymous_variant    SO Accesion: SO:0001819		ANNOVAR:synonymous SNV, silent mutation, silent substitution, silent_mutation, coding-synon, vep:SYNONYMOUS_CODING, synonymous codon, synonymous_coding, synonymous_codon, VAAST:synonymous_codon, SO:0001588, SO:0001588

#~ $AcceptedEffects{"start_lost"}=1;

## In the case o VEP this process is paralelized in 10 files so input and output files vary from 01 to 10, see the BASH SCRIPT to run it in a cluster.
## As the input file has been splited in 10 files and it is time consuming to copy the heather lines to each of them they have 
## been copied in a separate file (heathers_VEP.txt) which is used to load the hashes.
## Even though this is the parser of snpeff results I have decided to follow the same approach of reading the heathers from a separate file just in case
## snpeff needs to be paralelized sometime.

## With this command line we create the heather file using only lines starting with #

my $list_1= $ARGV[0];
my $OUT_snpeff_parsed=$ARGV[1];
my %hash1=();
my %hash2=();
my %hash3=();

my $time='['. timestamp(). ']'."\n";
print "Start reading INFILE:$time\n";
	
if(open (LIST_1, $list_1) && open (OUTPUT,'>'.$OUT_snpeff_parsed))
{

## Here we declare the hashes we are going to use to load the fields from the heather lines from the input file that contains only the lines starting 
## with # 

my %HI0_snpeff=();
my %IH0_snpeff=();
my %HI1_snpeff=();
my %IH1_snpeff=();
my %HI2_snpeff=();
my %IH2_snpeff=();
my %HI3_snpeff=();
my %IH3_snpeff=();

	while (my $line = <LIST_1>)
	{
	chomp ($line);
			
		if($line=~/^#/)
		{
			if($line=~/^##INFO=<ID=EFF.+Format: '(.+)' ">$/)	
			{	
				my $parse_field=$1;
				$parse_field=~s/(\(|\))/\|/g;
				$parse_field=~s/(\s|\]|\[)//g;
				#~ print "H:$parse_field:H\n";		
				my @tmp_HI1_snpeff = split (/\|/,$parse_field);		
				
				for(my $i=0;$i<scalar(@tmp_HI1_snpeff);$i++)			
				{
					$HI1_snpeff{$tmp_HI1_snpeff[$i]}=$i;
					$IH1_snpeff{$i}=$tmp_HI1_snpeff[$i];
					#~ print "$HI1_snpeff{$tmp_HI1_snpeff[$i]} \n";
					#~ print "$IH1_snpeff{$i} \n";			 
				}
				#~ my @dumkeys=keys (%HI1_snpeff);
				#~ foreach my $dumkeys_tok(@dumkeys){ print "$dumkeys_tok\n";}			
			}	
			elsif($line=~/^##INFO=<ID=LOF.+ Format: '(.+)' ">$/)	
			{
				my $parse_field=$1;
				$parse_field=~s/(\(|\)|\]|\[)/\|/g;
				$parse_field=~s/\s//g;
				#~ print "H:$parse_field:H\n";
				#~ print "$parse_field\n";
				my @tmp_HI2_snpeff = split (/\|/,$parse_field);
							
				for(my $i=0;$i<scalar(@tmp_HI2_snpeff);$i++)			
				{
					$HI2_snpeff{$tmp_HI2_snpeff[$i]}=$i;
					$IH2_snpeff{$i}=$tmp_HI2_snpeff[$i];
					#~ print "H:$HI2_snpeff{$tmp_HI2_snpeff[$i]}:H\n";
					#~ print "H:$IH2_snpeff{$i}:H\n";
				}
				#~ my @dumkeys=keys (%HI2_snpeff);
				#~ foreach my $dumkeys_tok(@dumkeys){print "$dumkeys_tok\n";}	
			}	
			elsif($line=~/^##INFO=<ID=NMD.+ Format: '(.+)' ">$/)					
			{	
				my $parse_field=$1;
				$parse_field=~s/(\(|\)|\]|\[)/\|/g;
				$parse_field=~s/\s//g;
				#~ print "H:$parse_field:H\n";
				#~ print "$parse_field\n";				
				my @tmp_HI3_snpeff = split (/\|/,$parse_field);
						
				for(my $i=0;$i<scalar(@tmp_HI3_snpeff);$i++)			
				{
					$HI3_snpeff{$tmp_HI3_snpeff[$i]}=$i;
					$IH3_snpeff{$i}=$tmp_HI3_snpeff[$i];
					#~ print "H:$HI3_snpeff{$tmp_HI3_snpeff[$i]}:H\n";
					#~ print "H:$IH3_snpeff{$i}:H\n";
				}
				#~ my @dumkeys=keys (%HI3_snpeff);
				#~ foreach my $dumkeys_tok(@dumkeys){ print "$dumkeys_tok\n";}		
			}	
			elsif ($line=~/^#CHROM/) 	
			{
				$line=~/^#(.+)/;
				my @tmp_HI0_snpeff = split (/\t/,$1);
				for(my $i=0;$i<scalar(@tmp_HI0_snpeff);$i++)						
				{
					$HI0_snpeff{$tmp_HI0_snpeff[$i]}=$i;
					$IH0_snpeff{$i}=$tmp_HI0_snpeff[$i];
					#~ print "$HI0_snpeff{$tmp_HI0_snpeff[$i]} \n";
					#~ print "$IH0_snpeff{$i} \n";
				}							
			}
		}			
		elsif ($line !~ /^#/)
		{
					#print "INICIO:$line_2\n";
						my @tmp = split (/\t/,$line);
					  		
					  		## Here we define the fields we are going to collect. Of note snpeff warns of bad transcripts (i.e. not known start_site) and
					  		## errors. We are collecting these two fields as separate features in case we want to filter out affected trannscripts.
					  					  
							my $CHROM="NaN";if (exists($HI0_snpeff{"CHROM"})){$CHROM=$tmp[$HI0_snpeff{"CHROM"}];}unless(defined($CHROM)){ print "ERROR in CHROM_SNPEFF";}
							my $FILTER="NaN";if (exists($HI0_snpeff{"FILTER"})){$FILTER=$tmp[$HI0_snpeff{"FILTER"}];}unless(defined($FILTER)){ print "ERROR in FILTER_SNPEFF";}
							my $INFO="NaN";if (exists($HI0_snpeff{"INFO"})){$INFO=$tmp[$HI0_snpeff{"INFO"}];}unless(defined($INFO)){ print "ERROR in INFO_SNPEFF";}
							
							#print "*************************$CHROM\t$FILTER\t$INFO\n";
						
						if((exists($AcceptedChr{"$CHROM"})) and ($INFO =~ /EFF\=/))
						
							{	
								#print "Hello_world_1:$CHROM\t$FILTER\t$INFO\n";
								
								my $POS="NaN";if (exists($HI0_snpeff{"POS"})){$POS=$tmp[$HI0_snpeff{"POS"}];}unless(defined($POS)){ print "ERROR in POS_SNPEFF";}
								my $ID="NaN";if (exists($HI0_snpeff{"ID"})){$ID=$tmp[$HI0_snpeff{"ID"}];}unless(defined($ID)){ print "ERROR in ID_SNPEFF";}
								my $REF="NaN";if (exists($HI0_snpeff{"REF"})){$REF=$tmp[$HI0_snpeff{"REF"}];}unless(defined($REF)){ print "ERROR in REF_SNPEFF";}
								my $ALT="NaN";if (exists($HI0_snpeff{"ALT"})){$ALT=$tmp[$HI0_snpeff{"ALT"}];}unless(defined($ALT)){ print "ERROR in ALT_SNPEFF";}
								my $QUAL="NaN";if (exists($HI0_snpeff{"QUAL"})){$QUAL=$tmp[$HI0_snpeff{"QUAL"}];}unless(defined($QUAL)){ print "ERROR in QUAL_SNPEFF";}
								my $AC="NaN";
								my $AN="NaN";
								my $AF="NaN";
								my $AFassessed="NaN";

								my $Effect="NaN"; 
								my $Effect_Impact="NaN";
								my $Exon_Rank="NaN"; 
								my $Amino_Acid_Change="NaN";
								my $Transcript_BioType="NaN";
								my $Functional_Class="NaN";
								my $Gene_Coding="NaN";
								my $Codon_Change="NaN";
								my $Amino_Acid_length="NaN";
								my $Transcript_ID="NaN";
								my $Genotype_Number="NaN";
								my $Gene_Name="NaN";
								my $ERRORS="NaN";
								my $WARNINGS="NaN";
								
								## Snpeff predicts LOF and LOF-based-on-NMD. Criteria according to P.Cingolani are:

								my $LOF="LOF_negative_snpEff";
								my $Gene_Name_LOF="NaN";
								my $Gene_ID_LOF="NaN";
								my $Number_of_transcripts_in_gene_LOF="NaN";
								my $Percent_of_transcripts_affected_LOF="NaN";

								my $NMD="NMD_negative_snpEff";
								my $Gene_Name_NMD="NaN";
								my $Gene_ID_NMD="NaN";
								my $Number_of_transcripts_in_gene_NMD="NaN";
								my $Percent_of_transcripts_affected_NMD="NaN";
								my @ INFO_derived_AFs= split ("\;",$INFO);
								
								
							
									foreach my $INFO_derived_AFs_tok (@ INFO_derived_AFs)
									
										{	
											#print "Estoy entrando en el búcle AF\n";
											if($INFO_derived_AFs_tok=~/^AC\=/){$INFO_derived_AFs_tok=~s/^AC\=//;$AC=$INFO_derived_AFs_tok;}
											if($INFO_derived_AFs_tok=~/^AN\=/){$INFO_derived_AFs_tok=~s/^AN\=//;$AN=$INFO_derived_AFs_tok;}
											if($INFO_derived_AFs_tok=~/^AF\=/){$INFO_derived_AFs_tok=~s/^AF\=//;$AF=$INFO_derived_AFs_tok;}
											
											if(($AC=~/[0-9]/)and($AN=~/[0-9]/)){$AFassessed=$AC/$AN;}
											
											if($INFO_derived_AFs_tok=~/LOF\=/)
															{
																#print "ORIGEN: $INFO_derived_AFs_tok\n";
																$LOF= "LOF_positive_snpEff";
																#print "$LOF\n";
																my @ INFO_derived_AFs_LOF_tmp= split ("LOF\=", $INFO_derived_AFs_tok);
																my $INFO_derived_AFs_LOF_tmp_tok=$INFO_derived_AFs_LOF_tmp[1];
																$INFO_derived_AFs_LOF_tmp_tok=~s/(\(|\)|\]|\[)|\s//g;
																#print "split_LOF: $INFO_derived_AFs_LOF_tmp_tok:H\n";
																my @parse_fields_LOF= split (/\|/, $INFO_derived_AFs_LOF_tmp_tok);
																#print "El array es:".join("***",@parse_fields_LOF)."\n";
																if (exists($HI2_snpeff{"Gene_Name"})){$Gene_Name_LOF=$parse_fields_LOF[$HI2_snpeff{"Gene_Name"}];if ($Gene_Name_LOF!~/\w/){$Gene_Name_LOF="NaN";}}else{ print "NaN in Gene_Name_LOF";}
																if (exists($HI2_snpeff{"Gene_ID"})){$Gene_ID_LOF=$parse_fields_LOF[$HI2_snpeff{"Gene_ID"}];if ($Gene_ID_LOF!~/\w/){$Gene_ID_LOF="NaN";}}else{ print "NaN in Gene_ID_LOF";}
																if (exists($HI2_snpeff{"Number_of_transcripts_in_gene"})){$Number_of_transcripts_in_gene_LOF=$parse_fields_LOF[$HI2_snpeff{"Number_of_transcripts_in_gene"}];if ($Number_of_transcripts_in_gene_LOF!~/\w/){$Number_of_transcripts_in_gene_LOF="NaN";}}else{ print "NaN  in Number_of_transcripts_in_gene_LOF";}
																if (exists($HI2_snpeff{"Percent_of_transcripts_affected"})){$Percent_of_transcripts_affected_LOF=$parse_fields_LOF[$HI2_snpeff{"Percent_of_transcripts_affected"}];if ($Percent_of_transcripts_affected_LOF!~/\w/){$Percent_of_transcripts_affected_LOF="NaN";}}else{ print "NaN in Percent_of_transcripts_affected_LOF";}
																#print "AAA:$Gene_Name_LOF\t$Gene_ID_LOF\t$Number_of_transcripts_in_gene_LOF\t$Percent_of_transcripts_affected_LOF:AAA\n";
																my $string=join(";",$LOF,$Gene_Name_LOF,$Gene_ID_LOF,$Number_of_transcripts_in_gene_LOF,$Percent_of_transcripts_affected_LOF);
																$hash1{$CHROM}{$POS}{$ID}{$REF}{$ALT}{$Gene_Name_LOF}{$string}=1;
															}
																#print "HHH:$Gene_Name_LOF\t$Gene_ID_LOF\t$Number_of_transcripts_in_gene_LOF\t$Percent_of_transcripts_affected_LOF:HHH\n";

											if($INFO_derived_AFs_tok=~/NMD\=/)
															{#
																$NMD= "NMD_positive_snpEff";
																#print "$NMD\n";
																my @ INFO_derived_AFs_NMD_tmp= split ("NMD\=", $INFO_derived_AFs_tok);	
																my $INFO_derived_AFs_NMD_tmp_tok=$INFO_derived_AFs_NMD_tmp[1];
																$INFO_derived_AFs_NMD_tmp_tok=~s/(\(|\)|\]|\[)|\s//g;
																#print "THE LINE TO PARSE NMD IS:$INFO_derived_AFs_NMD_tmp_tok\n";
																my @ parse_fields_NMD= split (/\|/, $INFO_derived_AFs_NMD_tmp_tok);
																		
																if (exists($HI3_snpeff{"Gene_Name"})){$Gene_Name_NMD=$parse_fields_NMD[$HI3_snpeff{"Gene_Name"}];if ($Gene_Name_NMD!~/\w/){$Gene_Name_NMD="NaN";}}else{ print "NaN in Gene_Name_NMD";}
																if (exists($HI3_snpeff{"Gene_ID"})){$Gene_ID_NMD=$parse_fields_NMD[$HI3_snpeff{"Gene_ID"}];if ($Gene_ID_NMD!~/\w/){$Gene_ID_NMD="NaN";}}else{ print "NaN in Gene_ID_NMD";}
																if (exists($HI3_snpeff{"Number_of_transcripts_in_gene"})){$Number_of_transcripts_in_gene_NMD=$parse_fields_NMD[$HI3_snpeff{"Number_of_transcripts_in_gene"}];if ($Number_of_transcripts_in_gene_NMD!~/\w/){$Number_of_transcripts_in_gene_NMD="NaN";}}else{ print "NaN  in Number_of_transcripts_in_gene_NMD";}
																if (exists($HI3_snpeff{"Percent_of_transcripts_affected"})){$Percent_of_transcripts_affected_NMD=$parse_fields_NMD[$HI3_snpeff{"Percent_of_transcripts_affected"}];if ($Percent_of_transcripts_affected_NMD!~/\w/){$Percent_of_transcripts_affected_NMD="NaN";}}else{ print "NaN in Percent_of_transcripts_affected_NMD";}
																my $string=join(";",$NMD,$Gene_Name_NMD,$Gene_ID_NMD,$Number_of_transcripts_in_gene_NMD,$Percent_of_transcripts_affected_NMD);
																$hash2{$CHROM}{$POS}{$ID}{$REF}{$ALT}{$Gene_Name_NMD}{$string}=1;
															}
											
											if($INFO_derived_AFs_tok=~/EFF\=/)
											{
												my @ INFO_derived_AFs_EFFs_tmp= split ("EFF\=", $INFO_derived_AFs_tok);
												foreach my $INFO_derived_AFs_EFFs_tmp_tok(@ INFO_derived_AFs_EFFs_tmp)
												{
													my @ INFO_derived_AFs_EFFs_tmp_tok_tmp= split ("\,", $INFO_derived_AFs_EFFs_tmp_tok);
													foreach my $INFO_derived_AFs_EFFs_tmp_tok_tmp_tok(@ INFO_derived_AFs_EFFs_tmp_tok_tmp)
													
													{
														$INFO_derived_AFs_EFFs_tmp_tok_tmp_tok=~s/(\(|\)|\]|\[)/\|/g;
														$INFO_derived_AFs_EFFs_tmp_tok_tmp_tok=~s/\s//g;
														my @parse_fields=split(/\|/,$INFO_derived_AFs_EFFs_tmp_tok_tmp_tok);														
														
														
														if(exists($HI1_snpeff{"Effect"})){$Effect=$parse_fields[$HI1_snpeff{"Effect"}];unless (defined($Effect)){$Effect="NaN";}} else {print "Error in Effect\n";}
														if(exists($HI1_snpeff{"Effect_Impact"})){$Effect_Impact=$parse_fields[$HI1_snpeff{"Effect_Impact"}];unless (defined($Effect_Impact)){$Effect_Impact="NaN";}} else {print "Error in Effect_Impact\n";}
														if(exists($HI1_snpeff{"Exon_Rank"})){$Exon_Rank=$parse_fields[$HI1_snpeff{"Exon_Rank"}];unless (defined($Exon_Rank)){$Exon_Rank="NaN";}} else {print "Error in Exon_Rank\n";}
														if(exists($HI1_snpeff{"Amino_Acid_Change"})){$Amino_Acid_Change=$parse_fields[$HI1_snpeff{"Amino_Acid_Change"}];unless (defined($Amino_Acid_Change)){$Amino_Acid_Change="NaN";}} else {print "Error in Amino_Acid_Change\n";}
														if(exists($HI1_snpeff{"Transcript_BioType"})){$Transcript_BioType=$parse_fields[$HI1_snpeff{"Transcript_BioType"}];unless (defined($Transcript_BioType)){$Transcript_BioType="NaN";}} else {print "Error in Transcript_BioType\n";}
														if(exists($HI1_snpeff{"Functional_Class"})){$Functional_Class=$parse_fields[$HI1_snpeff{"Functional_Class"}];unless (defined($Functional_Class)){$Functional_Class="NaN";}} else {print "Error in Functional_Class\n";}
														if(exists($HI1_snpeff{"Gene_Coding"})){$Gene_Coding=$parse_fields[$HI1_snpeff{"Gene_Coding"}];unless (defined($Gene_Coding)){$Gene_Coding="NaN";}} else {print "Error in Gene_Coding\n";}
														if(exists($HI1_snpeff{"Codon_Change"})){$Codon_Change=$parse_fields[$HI1_snpeff{"Codon_Change"}];unless (defined($Codon_Change)){$Codon_Change="NaN";}} else {print "Error in Codon_Change\n";}
														if(exists($HI1_snpeff{"Amino_Acid_length"})){$Amino_Acid_length=$parse_fields[$HI1_snpeff{"Amino_Acid_length"}];unless (defined($Amino_Acid_length)){$Amino_Acid_length="NaN";}} else {print "Error in Amino_Acid_length\n";}
														if(exists($HI1_snpeff{"Transcript_ID"})){$Transcript_ID=$parse_fields[$HI1_snpeff{"Transcript_ID"}];unless (defined($Transcript_ID)){$Transcript_ID="NaN";}} else {print "Error in Transcript_ID\n";}
														if(exists($HI1_snpeff{"Genotype_Number"})){$Genotype_Number=$parse_fields[$HI1_snpeff{"Genotype_Number"}];unless (defined($Genotype_Number)){$Genotype_Number="NaN";}} else {print "Error in Genotype_Number\n";}
														if(exists($HI1_snpeff{"Gene_Name"})){$Gene_Name=$parse_fields[$HI1_snpeff{"Gene_Name"}];unless (defined($Gene_Name)){$Gene_Name="NaN";}} else {print "Error in Gene_Name\n";}
														if(exists($HI1_snpeff{"ERRORS"})){$ERRORS=$parse_fields[$HI1_snpeff{"ERRORS"}];unless (defined($ERRORS)){$ERRORS="NaN";}} else {print "Error in ERRORS\n";}
														if(exists($HI1_snpeff{"WARNINGS"})){$WARNINGS=$parse_fields[$HI1_snpeff{"WARNINGS"}];unless (defined($WARNINGS)){$WARNINGS="NaN";}} else {print "Error in WARNINGS\n";}
														
														#print "Hello_world_2:$Effect\n";
																
														if(exists($AcceptedEffects{$Effect}))
															{
																	
																	#print "La línea es:$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$Effect\t$Gene_Name\t$Transcript_ID\t$Transcript_BioType\t$Amino_Acid_Change\t$Exon_Rank\t$Functional_Class\t$Gene_Coding\t$Codon_Change\t$Amino_Acid_length\t$Effect_Impact\t$Genotype_Number\t$AFassessed\t$ERRORS\t$WARNINGS\t$LOF\t$Gene_Name_LOF\t$Gene_ID_LOF\t$Number_of_transcripts_in_gene_LOF\t$Percent_of_transcripts_affected_LOF\t$NMD\t$Gene_Name_NMD\t$Gene_ID_NMD\t$Number_of_transcripts_in_gene_NMD\t$Percent_of_transcripts_affected_NMD\n";
																	#print "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$Effect\t$Gene_Name\t$Transcript_ID\t$Transcript_BioType\t$Amino_Acid_Change\t$Exon_Rank\t$Functional_Class\t$Gene_Coding\t$Codon_Change\t$Amino_Acid_length\t$Effect_Impact\t$Genotype_Number\t$AFassessed\t$ERRORS\t$WARNINGS\n";
																	my $string1=join(";",$Effect,$Gene_Name,$Transcript_ID,$Transcript_BioType,$Amino_Acid_Change,$Exon_Rank,$Functional_Class,$Gene_Coding,$Codon_Change,$Amino_Acid_length,$Effect_Impact,$Genotype_Number,$AFassessed,$ERRORS,$WARNINGS);
																	my $string2=join("\t",$CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER);
																	my $string3=join("\t",$string2,$string1);
																	#~ $hash3{$CHROM}{$POS}{$ID}{$REF}{$ALT}{$Gene_Name}{$string3}=1;
																	print OUTPUT "$string3\n";
															}
																
													
													}
												
												}
											
											}
									}# 	iNFOFIELDS AF					
					}# Filter Chrom, Pass Effect
					
				} #line !=#	
	}# while

}else {print "unable to open FILES: $list_1 or $OUT_snpeff_parsed\n";}

close (LIST_1);

#~ $time='['. timestamp(). ']'."\n";
#~ print "Start PRINTING:$time\n";
#~ 
#~ if()
#~ {
	#~ foreach my $CHROM_tok(sort keys%hash3)
	#~ {
	#~ foreach my $POS_tok(sort keys%{$hash3{$CHROM_tok}})
	#~ {
	#~ foreach my $ID_tok(sort keys%{$hash3{$CHROM_tok}{$POS_tok}})
	#~ {
	#~ foreach my $REF_tok(sort keys%{$hash3{$CHROM_tok}{$POS_tok}{$ID_tok}})
	#~ {
	#~ foreach my $ALT_tok(sort keys%{$hash3{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}})
	#~ {
	#~ foreach my $GENE_tok(sort keys%{$hash3{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}})
	#~ {
		#~ foreach my $string_tok(sort keys%{$hash3{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$GENE_tok}})
		#~ {
			#~ #print "$CHROM_tok\t$POS_tok\t$ID_tok\t$REF_tok\t$ALT_tok\t$GENE_tok\n";
			#~ print OUTPUT "$string_tok\t";
			#~ 
			#~ if (exists($hash1{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$GENE_tok}))
			#~ {
				#~ #print OUTPUT "HELLO_WORLD_III:\n";
				#~ foreach my $LOF_tok(sort keys%{$hash1{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$GENE_tok}})
				#~ {
					#~ print OUTPUT "$LOF_tok\t";
				#~ }
			#~ }
			#~ if (exists($hash2{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$GENE_tok}))
			#~ {
				#~ #print OUTPUT "HELLO_WORLD_IV:\n";
				#~ foreach my $NMD_tok(sort keys%{$hash1{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$GENE_tok}})
				#~ {
					#~ print OUTPUT "$NMD_tok\t";
				#~ }
			#~ }
			#~ print OUTPUT "\n";
		#~ }
#~ 
	#~ }	
		#~ 
	#~ }	
	#~ }	
	#~ }	
	#~ }	
	#~ }
#~ }

sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}


