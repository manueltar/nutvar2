##	Script to convert .vcf into .vcf with joint calls expressed as minimal representation and undo cryptic alleles. 2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne) 

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;

# Here we define the chromosomes we are going to accept; autosomes and sexual chromosomes.

my $time='['. timestamp(). ']'."\n";
print "Start charging Chromosome_hash:$time\n";

my %AcceptedChromosomes=();
my @AcceptChro=('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y');
foreach my $AcceptedChro_tok(@AcceptChro)
{
	$AcceptedChromosomes{$AcceptedChro_tok}=1;
}

# Here we declare the files we are going to use

my $input=$ARGV[0];
my $output=$ARGV[1];
my $output_2=$ARGV[2];
my $output_3=$ARGV[3];
my $output_4=$ARGV[4];
my $output_5=$ARGV[5];
my $output_6=$ARGV[6];
my $output_7=$ARGV[7];
my $output_8=$ARGV[8];

# Here we declare two variables needed indifferent parts of the script

my $unique_ID_gene="NaN";
my $unique_ID_transcript="NaN";

$time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

if(open(INPUT,$input) && open(OUTPUT1, '>'.$output) && open(OUTPUT2, '>'.$output_2)&& open(OUTPUT3, '>'.$output_3)&& open(OUTPUT4, '>'.$output_4)&& open(OUTPUT5, '>'.$output_5) && open(OUTPUT6, '>'.$output_6)&& open(OUTPUT7, '>'.$output_7)&& open(OUTPUT8, '>'.$output_8))
{
	while (my $line=<INPUT>)
	{
		# The first input is the gtf of ENSEMBL GRCh37.75
		
		#~ Input=  Homo_sapiens.GRCh37.75.gtf
		
		#~ #!genome-build GRCh37.p13
		#~ #!genome-version GRCh37
		#~ #!genome-date 2009-02
		#~ #!genome-build-accession NCBI:GCA_000001405.14
		#~ #!genebuild-last-updated 2013-09
		#~ 1       pseudogene      gene    11869   14412   .       +       .       gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene";
		#~ 1       processed_transcript    transcript      11869   14409   .       +       .       gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene"; transcript_name "DDX11L1-002"; transcript_source "havana";
		#~ 1       processed_transcript    exon    11869   12227   .       +       .       gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "1"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene"; transcript_name "DDX11L1-002"; transcript_source "havana"; exon_id "ENSE00002234944";

		chomp ($line);
		unless ($line=~/^#/)
		{
			# Here we declare the hashes of the different features we are going to extract from the file
			
			my % gene_hash=();
			my % gene_hash_rev=();
			my % transcript_hash=();
			my % transcript_hash_rev=();
			my % exon_hash=();
			my % exon_hash_rev=();
			my % CDS_hash=();
			my % CDS_hash_rev=();
			my % start_codon_hash=();
			my % start_codon_hash_rev=();
			my % stop_codon_hash=();
			my % stop_codon_hash_rev=();
			my % UTR_hash=();
			my % UTR_hash_rev=();
			my % Selenocysteine_hash=();
			my % Selenocysteine_hash_rev=();
			
			my($unique_ID_gene,$ENSG,$SYMBOL,$ENST,$NUMBER,$ENSE,$CCDS,$ENSP)="NaN"x8;
			
			my @tmp=split(/\t/,$line);
			
			# Here we accept only transcripts belonging to autosomes and sexual chromosmes that are 'protein coding'
			
			if($tmp[1] eq 'protein_coding' && exists($AcceptedChromosomes{$tmp[0]}))
			{
				# We obtain gene information
			
				if ($tmp[2]=~/gene/)
				{
					#print "H:$line:H\n";
					my @gene_tmp=split (/\t/,$line);
					my @gene_tmp_info_tmp=split (";",$gene_tmp[8]);
					
					foreach my $gene_tmp_info_tmp_tok(@gene_tmp_info_tmp)
						{
							if ($gene_tmp_info_tmp_tok=~/gene_id/){my @gene_tmp_info_tmp_unique_ID_gene_tmp=split ("\"", $gene_tmp_info_tmp_tok);$unique_ID_gene=$gene_tmp_info_tmp_unique_ID_gene_tmp[1];}
							elsif ($gene_tmp_info_tmp_tok=~/gene_name/){my @gene_tmp_info_tmp_ENST_tmp=split ("\"", $gene_tmp_info_tmp_tok);$SYMBOL=$gene_tmp_info_tmp_ENST_tmp[1];}
						}
					
					my $gene_coordinates=join("\t",$SYMBOL,$gene_tmp[6],$gene_tmp[3],$gene_tmp[4]);
					$gene_hash{$unique_ID_gene}{$gene_coordinates}{$tmp[0]}=1;
					
				}
				
				# We obtain transcript information
				
				elsif ($tmp[2]=~/transcript/)
				{
					#print "A:$line:A\n";
					my @transcript_tmp=split (/\t/,$line);
					my @transcript_tmp_info_tmp=split (";",$transcript_tmp[8]);
					foreach my $transcript_tmp_info_tmp_tok(@transcript_tmp_info_tmp)
						{
							if ($transcript_tmp_info_tmp_tok=~/gene_id/){my @transcript_tmp_info_tmp_ENSG_tmp=split ("\"", $transcript_tmp_info_tmp_tok);$ENSG=$transcript_tmp_info_tmp_ENSG_tmp[1];}
							elsif ($transcript_tmp_info_tmp_tok=~/transcript_id/){my @transcript_tmp_info_tmp_ENST_tmp=split ("\"", $transcript_tmp_info_tmp_tok);$ENST=$transcript_tmp_info_tmp_ENST_tmp[1];}
						}
					
					my $transcript_coordinates=join("\t",$transcript_tmp[3],$transcript_tmp[4]);
					$transcript_hash{$tmp[0]}{$ENSG}{$ENST}{$tmp[1]}{$transcript_coordinates}=1;
				}
				
				# We obtain exon information
				
				elsif ($tmp[2]=~/exon/)
				{
					#print "E:$line:E\n";
					my @exon_tmp=split (/\t/,$line);
					my @exon_tmp_info_tmp=split (";",$exon_tmp[8]);
					foreach my $exon_tmp_info_tmp_tok(@exon_tmp_info_tmp)
						{
							if ($exon_tmp_info_tmp_tok=~/transcript_id/){my @exon_tmp_info_tmp_ENST_tmp=split ("\"", $exon_tmp_info_tmp_tok);$ENST=$exon_tmp_info_tmp_ENST_tmp[1];}
							elsif ($exon_tmp_info_tmp_tok=~/exon_number/){my @exon_tmp_info_tmp_NUMBER_tmp=split ("\"", $exon_tmp_info_tmp_tok);$NUMBER=$exon_tmp_info_tmp_NUMBER_tmp[1];}
							elsif ($exon_tmp_info_tmp_tok=~/exon_id/){my @exon_tmp_info_tmp_ENSE_tmp=split ("\"", $exon_tmp_info_tmp_tok);$ENSE=$exon_tmp_info_tmp_ENSE_tmp[1];}
						}
					
					my $exon_coordinates=join("\t",$NUMBER,$ENSE,$exon_tmp[3],$exon_tmp[4]);
					$exon_hash{$ENST}{$exon_coordinates}=1;

				}
				
				# We obtain codifying sequence information
				
				elsif ($tmp[2]=~/CDS/)
				{
					#print "CDS:$line:CDS\n";
					my @CDS_tmp=split (/\t/,$line);
					my @CDS_tmp_info_tmp=split (";",$CDS_tmp[8]);
					foreach my $CDS_tmp_info_tmp_tok(@CDS_tmp_info_tmp)
						{
							if ($CDS_tmp_info_tmp_tok=~/transcript_id/){my @CDS_tmp_info_tmp_ENST_tmp=split ("\"", $CDS_tmp_info_tmp_tok);$ENST=$CDS_tmp_info_tmp_ENST_tmp[1];}
							elsif ($CDS_tmp_info_tmp_tok=~/ccds_id/){my @CDS_tmp_info_tmp_CCDS_tmp=split ("\"", $CDS_tmp_info_tmp_tok);$CCDS=$CDS_tmp_info_tmp_CCDS_tmp[1];}
							elsif ($CDS_tmp_info_tmp_tok=~/protein_id/){my @CDS_tmp_info_tmp_ENSP_tmp=split ("\"", $CDS_tmp_info_tmp_tok);$ENSP=$CDS_tmp_info_tmp_ENSP_tmp[1];}
						}
					
					unless(defined($CCDS)){$CCDS="NaNein";}
					#print "El valor es:$CCDS\n";
					my $CDS_coordinates=join("\t",$CCDS,$ENSP,$CDS_tmp[3],$CDS_tmp[4]);
					$CDS_hash{$ENST}{$CDS_coordinates}=1;
				}
			
				#	We obtain start codon information
			
				elsif ($tmp[2]=~/start_codon/)
				{
					#print "start_codon:$line:start_codon\n";
					my @start_codon_tmp=split (/\t/,$line);
					my @start_codon_tmp_info_tmp=split (";",$start_codon_tmp[8]);
					foreach my $start_codon_tmp_info_tmp_tok(@start_codon_tmp_info_tmp)
						{
							if ($start_codon_tmp_info_tmp_tok=~/transcript_id/){my @start_codon_tmp_info_tmp_ENST_tmp=split ("\"", $start_codon_tmp_info_tmp_tok);$ENST=$start_codon_tmp_info_tmp_ENST_tmp[1];}
							elsif ($start_codon_tmp_info_tmp_tok=~/ccds_id/){my @start_codon_tmp_info_tmp_CCDS_tmp=split ("\"", $start_codon_tmp_info_tmp_tok);$CCDS=$start_codon_tmp_info_tmp_CCDS_tmp[1];}
							elsif ($start_codon_tmp_info_tmp_tok=~/exon_number/){my @start_codon_tmp_info_tmp_NUMBER_tmp=split ("\"", $start_codon_tmp_info_tmp_tok);$NUMBER=$start_codon_tmp_info_tmp_NUMBER_tmp[1];}
						}
					
					unless(defined($CCDS)){$CCDS="NaNein";}
					my $start_codon_coordinates=join("\t",$NUMBER,$CCDS,$start_codon_tmp[3],$start_codon_tmp[4]);
					$start_codon_hash{$ENST}{$start_codon_coordinates}=1;
				}
				
				#	We obtain stop codon information
			
				elsif ($tmp[2]=~/stop_codon/)
				{
					#print "stop_codon:$line:stop_codon\n";
					my @stop_codon_tmp=split (/\t/,$line);
					my @stop_codon_tmp_info_tmp=split (";",$stop_codon_tmp[8]);
					foreach my $stop_codon_tmp_info_tmp_tok(@stop_codon_tmp_info_tmp)
						{
							if ($stop_codon_tmp_info_tmp_tok=~/transcript_id/){my @stop_codon_tmp_info_tmp_ENST_tmp=split ("\"", $stop_codon_tmp_info_tmp_tok);$ENST=$stop_codon_tmp_info_tmp_ENST_tmp[1];}
							elsif ($stop_codon_tmp_info_tmp_tok=~/ccds_id/){my @stop_codon_tmp_info_tmp_CCDS_tmp=split ("\"", $stop_codon_tmp_info_tmp_tok);$CCDS=$stop_codon_tmp_info_tmp_CCDS_tmp[1];}
							elsif ($stop_codon_tmp_info_tmp_tok=~/exon_number/){my @stop_codon_tmp_info_tmp_NUMBER_tmp=split ("\"", $stop_codon_tmp_info_tmp_tok);$NUMBER=$stop_codon_tmp_info_tmp_NUMBER_tmp[1];}
						}
					
					unless(defined($CCDS)){$CCDS="NaNein";}
					my $stop_codon_coordinates=join("\t",$NUMBER,$CCDS,$stop_codon_tmp[3],$stop_codon_tmp[4]);
					$stop_codon_hash{$ENST}{$stop_codon_coordinates}=1;
				}			
			
				#	We obtain UTR information
			
				elsif ($tmp[2]=~/UTR/)
				{
					#print "UTR:$line:UTR\n";
					my @UTR_tmp=split (/\t/,$line);
					my @UTR_tmp_info_tmp=split (";",$UTR_tmp[8]);
					foreach my $UTR_tmp_info_tmp_tok(@UTR_tmp_info_tmp)
						{
							if ($UTR_tmp_info_tmp_tok=~/transcript_id/){my @UTR_tmp_info_tmp_ENST_tmp=split ("\"", $UTR_tmp_info_tmp_tok);$ENST=$UTR_tmp_info_tmp_ENST_tmp[1];}
						}
					my $UTR_coordinates=join("\t",$UTR_tmp[3],$UTR_tmp[4]);
					$UTR_hash{$ENST}{$UTR_coordinates}=1;
				}
				
				#	We obtain Selenocysteine residues information
							
				elsif ($tmp[2]=~/Selenocysteine/)
				{
					#print "Selenocysteine:$line:Selenocysteine\n";
					my @Selenocysteine_tmp=split (/\t/,$line);
					my @Selenocysteine_tmp_info_tmp=split (";",$Selenocysteine_tmp[8]);
					foreach my $Selenocysteine_tmp_info_tmp_tok(@Selenocysteine_tmp_info_tmp)
						{
							if ($Selenocysteine_tmp_info_tmp_tok=~/transcript_id/){my @Selenocysteine_tmp_info_tmp_ENST_tmp=split ("\"", $Selenocysteine_tmp_info_tmp_tok);$ENST=$Selenocysteine_tmp_info_tmp_ENST_tmp[1];}
							elsif ($Selenocysteine_tmp_info_tmp_tok=~/ccds_id/){my @Selenocysteine_tmp_info_tmp_CCDS_tmp=split ("\"", $Selenocysteine_tmp_info_tmp_tok);$CCDS=$Selenocysteine_tmp_info_tmp_CCDS_tmp[1];}
						}
					unless(defined($CCDS)){$CCDS="NaNein";}
					my $Selenocysteine_coordinates=join("\t",$CCDS,$Selenocysteine_tmp[3],$Selenocysteine_tmp[4]);
					$Selenocysteine_hash{$ENST}{$Selenocysteine_coordinates}=1;
				}
		
			}
			
			# From now onwards we start printing different files with the information we have stored in the hashes on a per line basis
			
			foreach my $gene_hash_tok (sort keys %gene_hash)
			{
				#print "$gene_hash_tok\n";
				
					#print "$gene_hash_tok\t$symbol_tok\n";
					foreach my $coordinates_tok (sort keys %{$gene_hash{$gene_hash_tok}})
					{
						foreach my $CHROM_tok (sort keys %{$gene_hash{$gene_hash_tok}{$coordinates_tok}})
						{
							print OUTPUT1 "$CHROM_tok\t$gene_hash_tok\t$coordinates_tok\n";
						}
					}	
			}
			
			foreach my $CHROM_tok (sort{$a<=>$b} keys %transcript_hash)
			{
				foreach my $ENSG_tok (sort keys %{$transcript_hash{$CHROM_tok}})
				{
					#print "$transcript_hash_tok\n";
					foreach my $ENST_tok (sort keys %{$transcript_hash{$CHROM_tok}{$ENSG_tok}})
					{
						foreach my $BIOTYPE_tok (sort keys %{$transcript_hash{$CHROM_tok}{$ENSG_tok}{$ENST_tok}})
						{
							#print "$transcript_hash_tok\t$symbol_tok\t$strand_tok\n";
							foreach my $transcript_coordinates_tok (sort keys %{$transcript_hash{$CHROM_tok}{$ENSG_tok}{$ENST_tok}{$BIOTYPE_tok}})
							{
								
								print OUTPUT2 "$CHROM_tok\t$ENSG_tok\t$ENST_tok\t$BIOTYPE_tok\t$transcript_coordinates_tok\n";
								
							}
						}
						
					}
					
				}#DECONV Transcript
			}

		# $exon_hash{$ENST}{$exon_coordinates}=1;
		
		foreach my $ENST_tok (sort keys %exon_hash)
				{
					#print "$exon_hash_tok\n";
					foreach my $exon_coordinates (sort keys %{$exon_hash{$ENST_tok}})
					{		
						print OUTPUT3 "$ENST_tok"."\t"."$exon_coordinates\n";		
					}
					
				}#DECONV Exon
				
		# $CDS_hash{$ENST}{$CDS_coordinates}=1;
		
		foreach my $ENST_tok (sort keys %CDS_hash)
				{
					#print "$CDS_hash_tok\n";
					foreach my $CDS_coordinates(sort keys %{$CDS_hash{$ENST_tok}})
					{		
						print OUTPUT4 "$ENST_tok"."\t"."$CDS_coordinates\n";		
					}
					
				}#DECONV CDS
				
		foreach my $ENST_tok (sort keys %start_codon_hash)
				{
					#print "$start_codon_hash_tok\n";
					foreach my $start_codon_coordinates(sort keys %{$start_codon_hash{$ENST_tok}})
					{		
						print OUTPUT5 "$ENST_tok"."\t"."$start_codon_coordinates\n";		
					}
					
				}#DECONV $start_codon_hash{$ENST}{$start_codon_coordinates}=1;
				
		foreach my $ENST_tok (sort keys %stop_codon_hash)
				{
					#print "$stop_codon_hash_tok\n";
					foreach my $stop_codon_coordinates(sort keys %{$stop_codon_hash{$ENST_tok}})
					{		
						print OUTPUT6 "$ENST_tok"."\t"."$stop_codon_coordinates\n";		
					}
					
				}#DECONV $stop_codon_hash{$ENST}{$stop_codon_coordinates}=1;
				
		foreach my $ENST_tok (sort keys %UTR_hash)
				{
					#print "$UTR_hash_tok\n";
					foreach my $UTR_coordinates(sort keys %{$UTR_hash{$ENST_tok}})
					{		
						print OUTPUT7 "$ENST_tok"."\t"."$UTR_coordinates\n";		
					}
					
				}#DECONV $UTR_hash{$ENST}{$UTR_coordinates}=1;
		
		foreach my $ENST_tok (sort keys %Selenocysteine_hash)
				{
					#print "$Selenocysteine_hash_tok\n";
					foreach my $Selenocysteine_coordinates(sort keys %{$Selenocysteine_hash{$ENST_tok}})
					{		
						print OUTPUT8 "$ENST_tok"."\t"."$Selenocysteine_coordinates\n";		
					}
					
				}#DECONV $Selenocysteine_hash{$ENST}{$Selenocysteine_coordinates}=1;
		}
	}
	
}else {print "unable to open $output or $input\n";}


sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}


