##	Script to parse transcript features from the gtf file of ENSEMBL. 2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne) 

#!/usr/bin/perl

use strict;
use warnings;

my $input=$ARGV[0];
my $input_2=$ARGV[1];
my $input_3=$ARGV[2];
my $input_4=$ARGV[3];
my $input_5=$ARGV[4];
my $output1=$ARGV[5];
my $output2=$ARGV[6];
my $output3=$ARGV[7];
my $output4=$ARGV[8];
my $output5=$ARGV[9];

my % UTR_hash=();
my % START_hash=();
my % CDS_hash=();
my % Selenocysteine_hash=();
my % STOP_hash=();
	


if(open(INPUT,$input) && open(INPUT2, $input_2)&& open(INPUT3, $input_3)&& open(INPUT4, $input_4)&& open(INPUT5, $input_5) && open(OUTPUT1, '>'.$output1)&& open(OUTPUT2, '>'.$output2)&& open(OUTPUT3, '>'.$output3)&& open(OUTPUT4, '>'.$output4)&& open(OUTPUT5, '>'.$output5))
{
	while (my $line=<INPUT>)
	{
		chomp ($line);
		#~ print "$line\n";
		if($line=~/(ENST[^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			
			#~ Input= gtf_output_UTR.txt
			#~ 
			#~ ENST00000423372 139310  139379
			#~ ENST00000423372 137621  138529
			#~ ENST00000423372 134901  135802
			#~ ENST00000426406 367640  367658
			#~ ENST00000426406 368598  368634
			
			#~ print "****$line******\n";
			my $unique_ID_UTR=$1;
			my $UTR_begin=$2;
			my $UTR_end=$3;
			
			# Here we store all the UTR coordinates on a per transcript basis
			
			# # print "Lo que saca el array es:$unique_ID_UTR\t$UTR_fields\n";
			$UTR_hash{$unique_ID_UTR}{$UTR_begin}{$UTR_end}=1;
			#~ print "AAAAA:$unique_ID_UTR\t$UTR_begin\t$UTR_end:AAAA\n";	
		}
	}
close (INPUT);
foreach my $unique_ID_UTR_tok(sort keys %UTR_hash)
					{
						# Here we print all the coordinates of the UTR for each transcript collapsed in the same line
						
						#~ my @UTR_fields_tok_tmp=();
						#~ my @UTR_begin_tmp=();
						#~ my @UTR_end_tmp=();
						print OUTPUT1 "$unique_ID_UTR_tok\t";
						#~ print "UTR"."$unique_ID_UTR_tok\t";
						
						my @UTR_begin_tmp=sort{$a<=>$b}keys%{$UTR_hash{$unique_ID_UTR_tok}};
						for (my $i=0; $i<scalar(@UTR_begin_tmp);$i++)
						{
							foreach my $UTR_end_tok(sort{$a<=>$b}keys%{$UTR_hash{$unique_ID_UTR_tok}{$UTR_begin_tmp[$i]}})
							{
								print OUTPUT1 "UTR:["."$UTR_begin_tmp[$i]"."]"."-"."["."$UTR_end_tok"."]"."\t";
								#~ print "UTR:["."$UTR_begin_tmp[$i]"."]"."-"."["."$UTR_end_tok"."]"."\t";
							}
						}
						print OUTPUT1 "\n";
						#~ print "\n";
					}	

	while (my $line2=<INPUT2>)
		{
			#~ Input= gtf_output_start_codons.txt
			#~ 
			#~ ENST00000335137 1       CCDS30547       69091   69093
			#~ ENST00000426406 1       CCDS41220       367659  367661
			#~ ENST00000332831 1       CCDS41221       622032  622034
			#~ ENST00000342066 2       CCDS2   861322  861324
			#~ ENST00000598827 1       NaNein  866443  866445

			chomp ($line2);
				if($line2=~/(ENST[^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)/)
			{
				my $unique_ID_START=$1;
				my $number_exon_START=$2;
				my $begin_START=$3;
				my $end_START=$4;
				
				# Here we extract the coordinates and number of exon where it lies of the START CODON
				
				$START_hash{$unique_ID_START}{$number_exon_START}{$begin_START}{$end_START}=1;
				#~ print "Lo que saca el array es:$unique_ID_START\t$number_exon_START\t$begin_START\t$end_START\n";
			}
		}
close (INPUT2);

foreach my $unique_ID_START_tok(sort keys %START_hash)
					{
						
						# Here we print all the coordinates and exon number of the START for each transcript collapsed in the same line
						
						print OUTPUT2 "$unique_ID_START_tok\t";
						#~ print "START"."$unique_ID_START_tok\t";
						
						my @START_exon_tmp=sort{$a<=>$b}keys%{$START_hash{$unique_ID_START_tok}};
						for (my $i=0; $i<scalar(@START_exon_tmp);$i++)
						{
							my @START_begin_tmp=sort{$a<=>$b}keys%{$START_hash{$unique_ID_START_tok}{$START_exon_tmp[$i]}};
							for (my $j=0; $j<scalar(@START_begin_tmp);$j++)
							{
								foreach my $START_end_tok(sort{$a<=>$b}keys%{$START_hash{$unique_ID_START_tok}{$START_exon_tmp[$i]}{$START_begin_tmp[$j]}})
								{
									print OUTPUT2 "START:$START_exon_tmp[$i]:"."["."$START_begin_tmp[$j]"."]"."-"."["."$START_end_tok"."]"."\t";
									#~ print "START:$START_exon_tmp[$i]:"."["."$START_begin_tmp[$j]"."]"."-"."["."$START_end_tok"."]"."\t";
								}
							}
						}
						print OUTPUT2 "\n";
						#~ print "\n";
					}
	
	while (my $line3=<INPUT3>)
			{
				#~ Input=gtf_output_CDS.txt
				#~ 
				#~ ENST00000335137 CCDS30547       ENSP00000334393 69091   70005
				#~ ENST00000423372 NaNein  ENSP00000473460 138533  139309
				#~ ENST00000426406 CCDS41220       ENSP00000409316 367659  368594
				#~ ENST00000332831 CCDS41221       ENSP00000329982 621099  622034
				
				chomp ($line3);
				
				if($line3=~/(ENST[^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
				{
					# Here we extract the coordinates and CCDS of the different CDS intervals
					
					my $unique_ID_CDS=$1;
					my $CCDS=$2;
					my $ENSP=$3;
					my $begin_CDS=$4;
					my $end_CDS=$5;
					$CDS_hash{$unique_ID_CDS}{$CCDS}{$ENSP}{$begin_CDS}{$end_CDS}=1;
					#~ print "$unique_ID_CDS\t$CCDS\t$ENSP\t$begin_CDS\t$end_CDS\n";
				}
			}
close (INPUT3);		
foreach my $unique_ID_CDS_tok(sort keys %CDS_hash)
{
foreach my $CCDS_tok(sort keys %{$CDS_hash{$unique_ID_CDS_tok}})
{
foreach my $ENSP_tok(sort keys %{$CDS_hash{$unique_ID_CDS_tok}{$CCDS_tok}})
{
	
	# Here we print all the coordinates and CCDS of the CDS for each transcript collapsed in the same line

	print OUTPUT3 "$unique_ID_CDS_tok\tCDS:$CCDS_tok\tCDS:$ENSP_tok\t";
	#~ print "$unique_ID_CDS_tok\t$CCDS_tok\t$ENSP_tok\t";
	my @CDS_begin_tmp=sort{$a<=>$b}keys%{$CDS_hash{$unique_ID_CDS_tok}{$CCDS_tok}{$ENSP_tok}};
	for (my $i=0; $i<scalar(@CDS_begin_tmp);$i++)
						{
							foreach my $CDS_end_tok(sort{$a<=>$b}keys%{$CDS_hash{$unique_ID_CDS_tok}{$CCDS_tok}{$ENSP_tok}{$CDS_begin_tmp[$i]}})
							{
								print OUTPUT3 "CDS:["."$CDS_begin_tmp[$i]"."]"."-"."["."$CDS_end_tok"."]"."\t";
								#~ print "CDS:["."$CDS_begin_tmp[$i]"."]"."-"."["."$CDS_end_tok"."]"."\t";
							}
						}
						print OUTPUT3 "\n";
						#~ print "\n";
}	
}	
}

	while (my $line4=<INPUT4>)
			{
				#~ Input=gtf_output_Selenocysteine.txt
				#~ 
				#~ ENST00000361547 CCDS41282       26128584        26128586
				#~ ENST00000361547 CCDS41282       26139280        26139282
				#~ ENST00000374315 CCDS41283       26139280        26139282

				
				chomp ($line4);
				#~ print "$line\n";
				if($line4=~/(ENST[^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)/)
				{
					# Here we extract the coordinates of the different Selenocysteine residues
					
					#~ print "****$line******\n";
					my $unique_ID_Selenocysteine=$1;
					my $Selenocysteine_begin=$2;
					my $Selenocysteine_end=$3;
					# # print "Lo que saca el array es:$unique_ID_Selenocysteine\t$Selenocysteine_fields\n";
					$Selenocysteine_hash{$unique_ID_Selenocysteine}{$Selenocysteine_begin}{$Selenocysteine_end}=1;
					#~ print "AAAAA:$unique_ID_Selenocysteine\t$Selenocysteine_begin\t$Selenocysteine_end:AAAA\n";	
				}
			}
close (INPUT4);	
foreach my $unique_ID_Selenocysteine_tok(sort keys %Selenocysteine_hash)
					{
						# Here we print all the coordinates of the Selenocysteine residues for each transcript collapsed in the same line

						print OUTPUT4 "$unique_ID_Selenocysteine_tok\t";
						#~ print "Selenocysteine"."$unique_ID_Selenocysteine_tok\t";
						
						my @Selenocysteine_begin_tmp=sort{$a<=>$b}keys%{$Selenocysteine_hash{$unique_ID_Selenocysteine_tok}};
						for (my $i=0; $i<scalar(@Selenocysteine_begin_tmp);$i++)
						{
							foreach my $Selenocysteine_end_tok(sort{$a<=>$b}keys%{$Selenocysteine_hash{$unique_ID_Selenocysteine_tok}{$Selenocysteine_begin_tmp[$i]}})
							{
								print OUTPUT4 "Seleno:["."$Selenocysteine_begin_tmp[$i]"."]"."-"."["."$Selenocysteine_end_tok"."]"."\t";
								#~ print "Selenocysteine:["."$Selenocysteine_begin_tmp[$i]"."]"."-"."["."$Selenocysteine_end_tok"."]"."\t";
							}
						}
						print OUTPUT4 "\n";
						#~ print "\n";
					}
	while (my $line5=<INPUT5>)
		{
			#~ Input=gtf_output_stop_codons.txt
			#~ 
			#~ ENST00000335137 1       CCDS30547       70006   70008
			#~ ENST00000423372 1       NaNein  138530  138532
			#~ ENST00000426406 1       CCDS41220       368595  368597
			#~ ENST00000332831 1       CCDS41221       621096  621098
			
			chomp ($line5);
				if($line5=~/(ENST[^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)/)
			{
				# Here we extract the coordinates and number of exon where it lies of the STOP CODON
				
				my $unique_ID_STOP=$1;
				my $number_exon_STOP=$2;
				my $begin_STOP=$3;
				my $end_STOP=$4;
				$STOP_hash{$unique_ID_STOP}{$number_exon_STOP}{$begin_STOP}{$end_STOP}=1;
				#~ print "Lo que saca el array es:$unique_ID_STOP\t$number_exon_STOP\t$begin_STOP\t$end_STOP\n";
			}
		}
close (INPUT5);

foreach my $unique_ID_STOP_tok(sort keys %STOP_hash)
					{
						print OUTPUT5 "$unique_ID_STOP_tok\t";
						#~ print "STOP"."$unique_ID_STOP_tok\t";
						
						# Here we print all the coordinates of the stop codon and the exon number where it lies for each transcript collapsed in the same line
						
						my @STOP_exon_tmp=sort{$a<=>$b}keys%{$STOP_hash{$unique_ID_STOP_tok}};
						for (my $i=0; $i<scalar(@STOP_exon_tmp);$i++)
						{
							my @STOP_begin_tmp=sort{$a<=>$b}keys%{$STOP_hash{$unique_ID_STOP_tok}{$STOP_exon_tmp[$i]}};
							for (my $j=0; $j<scalar(@STOP_begin_tmp);$j++)
							{
								foreach my $STOP_end_tok(sort{$a<=>$b}keys%{$STOP_hash{$unique_ID_STOP_tok}{$STOP_exon_tmp[$i]}{$STOP_begin_tmp[$j]}})
								{
									print OUTPUT5 "STOP:$STOP_exon_tmp[$i]:"."["."$STOP_begin_tmp[$j]"."]"."-"."["."$STOP_end_tok"."]"."\t";
									#~ print "STOP:$STOP_exon_tmp[$i]:"."["."$STOP_begin_tmp[$j]"."]"."-"."["."$STOP_end_tok"."]"."\t";

								}
							}
						}
						print OUTPUT5 "\n";
						#~ print "\n";
					}
						
}else {print "unable to open $input\n";}

