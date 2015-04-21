##	Script to create files with the SYMBOL, strand, ENST and transcript features (UTR 5 and 3 prime, START, STOP and Selenocys) on a per transcript basis of ENSEMBL. 2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)


use strict;
use warnings;
use Time::localtime;

my %hash1=();
my %hash2=();
my %hash3=();
my %hash4=();
my %hash5=();
my %hash6=();
my %ENST_hash=();


my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $input3=$ARGV[2];
my $input4=$ARGV[3];
my $input5=$ARGV[4];
my $input6=$ARGV[5];
my $output=$ARGV[6];

my $time='['. timestamp(). ']'."\n";
print "Start charging hash1:$time\n";

# Here we obtain ENSG - SYMBOL (HGNC nomenclature) - strand equivalence

if (open (INPUT1, $input1))
{
	#~ INPUT1= gtf_output_ENSG.txt
	#~ 
#~ ENSG00000186092 OR4F5   +       69091   70008
#~ ENSG00000237683 AL627309.1      -       134901  139379

while(my $line=<INPUT1>)
	{
		chomp $line;
		#print "$line\n";
		if($line=~/^(ENSG[^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			my $ENSG=$1;
			my $SYMBOL=$2;
			my $strand=$3;
			my $begin=$4;
			my $end=$5;
			#print "hello_world:$SYMBOL\t$strand\t$ENSG\n";
			$hash1{$ENSG}{$SYMBOL}{$strand}=1;
			#~ print "$ENSG\t$SYMBOL\t$strand\n";
		}
	}
}else{print "Unable to open INPUT1\n";}

$time='['. timestamp(). ']'."\n";
print "Start charging hash2:$time\n";

# Here we obtain transcript - ENSG equivalence

if (open(INPUT2, $input2))
{
#~ INPUT2=gtf_output_ENST.txt
#~ 
#~ 1       ENSG00000186092 ENST00000335137 protein_coding  69091   70008
#~ 1       ENSG00000237683 ENST00000423372 protein_coding  134901  139379
#~ 1       ENSG00000235249 ENST00000426406 protein_coding  367640  368634
#~ 1       ENSG00000185097 ENST00000332831 protein_coding  621059  622053
#~ 1       ENSG00000269831 ENST00000599533 protein_coding  738532  739137
#~ 1       ENSG00000269308 ENST00000594233 protein_coding  818043  819983
	
	while(my $line=<INPUT2>)
	{
		chomp($line);
		if($line=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
			my $CHROM=$1;
			my $ENSG=$2;
			my $ENST=$3;
			my $BIOTYPE=$4;
			my $begin=$5;
			my $end=$6;
			$hash2{$ENSG}{$CHROM}{$ENST}{$BIOTYPE}=1;
			#~ print "$ENSG\t$ENST\n";
		
		}
	}
}else {print "impossible to open INPUT2\n";die;}

# From now onwards we parse each of the files carrying the features of the transcripts and store each in a hash

$time='['. timestamp(). ']'."\n";
print "Start charging hash3:$time\n";

if (open(INPUT3, $input3))
{
			#~ INPUT3=gtf_output_JOINED_CDS.txt
		
		#~ ENST00000000233 CDS:CCDS34745   CDS:ENSP00000000233     CDS:[127228553]-[127228619]     CDS:[127229137]-[127229217]     CDS:[127229539]-[127229648]     CDS:[127230120]-[127230191]     CDS:[127231017]-[127231142]     CDS:[127231267]-[127231350]     
		#~ ENST00000000412 CDS:CCDS8598    CDS:ENSP00000000412     CDS:[9094417]-[9094536] CDS:[9095012]-[9095138] CDS:[9096001]-[9096131] CDS:[9096397]-[9096506] CDS:[9098014]-[9098180] CDS:[9098825]-[9099000] 

	while(my $line=<INPUT3>)
	{
		chomp($line);
		if($line=~/([^\t]+)\tCDS:([^\t]+)\tCDS:([^\t]+)\t(.+)/)
		{
			my $ENST=$1;
			my $CCDS=$2;
			my $ENSP=$3;
			my $fields=$4;
			$hash3{$ENST}{$CCDS}{$ENSP}{$fields}=1;
			#~ print "$ENST\t$CCDS\t$ENSP\t**$fields**\n";
		}
	}
}else {print "impossible to open INPUT3\n";die;}

$time='['. timestamp(). ']'."\n";
print "Start charging hash4:$time\n";

if (open(INPUT4, $input4))
{
			#~ INPUT4=gtf_output_JOINED_UTR.txt 
		
	#~ ENST00000000233 UTR:[127228399]-[127228552]     UTR:[127231354]-[127231759]     
	#~ ENST00000000412 UTR:[9092961]-[9094413] UTR:[9099001]-[9099001] UTR:[9102084]-[9102551] 
	#~ ENST00000000442 UTR:[64073050]-[64073208]       UTR:[64074640]-[64074651]       UTR:[64083439]-[64084210]   
	
	while(my $line=<INPUT4>)
	{
		chomp($line);
		if($line=~/([^\t]+)\t(.+)/)
		{
			my $ENST=$1;
			my $fields=$2;
			$hash4{$ENST}{$fields}=1;
			#~ print "$ENST\t**$fields**\n";
		}
	}
}else {print "impossible to open INPUT3\n";die;}

$time='['. timestamp(). ']'."\n";
print "Start charging hash5:$time\n";

if (open(INPUT5, $input5))
{
			#~ INPUT5=gtf_output_JOINED_START.txt 
		
	#~ ENST00000000233 START:1:[127228553]-[127228555] 
	#~ ENST00000000412 START:2:[9098998]-[9099000]  
	
	while(my $line=<INPUT5>)
	{
		chomp($line);
		if($line=~/([^\t]+)\t(.+)/)
		{
			my $ENST=$1;
			my $fields=$2;
			$hash5{$ENST}{$fields}=1;
			#~ print "$ENST\t**$fields**\n";
		}
	}
}else {print "impossible to open INPUT6\n";die;}

if (open(INPUT6, $input6))
{
			#~ INPUT6=gtf_output_JOINED_STOP.txt 
		
		#~ ENST00000000233 STOP:6:[127231351]-[127231353]  
		#~ ENST00000000412 STOP:7:[9094414]-[9094416]      
		#~ ENST00000000442 STOP:7:[64083436]-[64083438]  
	
	while(my $line=<INPUT6>)
	{
		chomp($line);
		if($line=~/([^\t]+)\t(.+)/)
		{
			my $ENST=$1;
			my $fields=$2;
			$hash6{$ENST}{$fields}=1;
			#~ print "$ENST\t**$fields**\n";
		}
	}
}else {print "impossible to open INPUT6\n";die;}

# We start processing the different hashes

$time='['. timestamp(). ']'."\n";
print "Start PROCESSING:$time\n";

	foreach my $ENSG_tok(sort keys%hash1)
	{
	foreach my $SYMBOL_tok(sort keys%{$hash1{$ENSG_tok}})
	{
	foreach my $strand_tok(sort keys%{$hash1{$ENSG_tok}{$SYMBOL_tok}})
	{
		foreach my $CHROM_tok(sort keys%{$hash2{$ENSG_tok}})
		{
		foreach my $ENST_tok(sort keys%{$hash2{$ENSG_tok}{$CHROM_tok}})
		{
		foreach my $BIOTYPE_tok(sort keys%{$hash2{$ENSG_tok}{$CHROM_tok}{$ENST_tok}})
		{
			my $string="NaN";
			
			# order_hash is a hash created to order the UTRS and classsify them between 5 prime and 3 prime UTRs.
			
			my %order_hash=();
			foreach my $CCDS_tok(sort keys%{$hash3{$ENST_tok}})
			{
			foreach my $ENSP_tok(sort keys%{$hash3{$ENST_tok}{$CCDS_tok}})
			{
				$string=join("\t",$ENST_tok,$CHROM_tok,$SYMBOL_tok,$strand_tok,$BIOTYPE_tok,$CCDS_tok,$ENSP_tok);
				
				# print "$ENST_tok\t$CHROM_tok\t$SYMBOL_tok\t$strand_tok\t$BIOTYPE_tok\t$CCDS_tok\t$ENSP_tok\n";
			
			foreach my $fields_tok(sort keys%{$hash3{$ENST_tok}{$CCDS_tok}{$ENSP_tok}})
			{	
				#~ # print "Fields_CDS:$fields_tok\n";
				my @CDS_tmp=split(/\t/,$fields_tok);
				foreach my $CDS_tmp_tok(@CDS_tmp)
				{
					if($CDS_tmp_tok=~/CDS:\[(.+)\]\-\[(.+)\]/)
					{
						my $begin=$1;
						my $end=$2;
						$order_hash{'CDS'}{$begin}{$end}=1;
						$ENST_hash{$strand_tok}{$string}{'CDS'}{$begin}{$end}=1;
						# print "'CDS'\t$begin\t$end**\n";
					}
				}
			}
			}
			}
			
			# If the transcript has UTRS we store them in order_hash
			
			if(exists($hash4{$ENST_tok}))
			{
				foreach my $UTR_fields_tok(sort keys%{$hash4{$ENST_tok}})
				{
					my @UTR_tmp=split(/\t/,$UTR_fields_tok);
					foreach my $UTR_tmp_tok(@UTR_tmp)
					{
						if($UTR_tmp_tok=~/UTR:\[(.+)\]\-\[(.+)\]/)
						{
							my $begin=$1;
							my $end=$2;
							$order_hash{'UTR'}{$begin}{$end}=1;
							# print "'UTR'\t$begin\t$end**\n";
						}
					}
				}
			}
			
			# If the trasncript has not UTRS their coordinates will be noted as NaN and stored in order hash
			
			else
			{
				my $begin="NaN";
				my $end="NaN";
				$order_hash{'UTR'}{$begin}{$end}=1;
				# print "'UTR'\t$begin\t$end**\n";
			}
			
			# If the transcript has a START codon we stor the coordinates in ENST hash
			
			if(exists($hash5{$ENST_tok}))
			{
				foreach my $START_fields_tok(sort keys%{$hash5{$ENST_tok}})
				{
					my @START_tmp=split(/\t/,$START_fields_tok);
					foreach my $START_tmp_tok(@START_tmp)
					{
						if($START_tmp_tok=~/START:\d+:\[(.+)\]\-\[(.+)\]/)
						{
							my $begin=$1;
							my $end=$2;
							$ENST_hash{$strand_tok}{$string}{'START'}{$begin}{$end}=1;
							# print "'START'\t$begin\t$end**\n";
						}
					}
				}
			}
			
			# If the transcript has not a START codon their coordinates will be noted as NaN
			
			else
			{
				my $begin="NaN";
				my $end="NaN";
				$ENST_hash{$strand_tok}{$string}{'START'}{$begin}{$end}=1;
				# print "'START'\t$begin\t$end**\n";
			}
			
			# If the transcript has a STOP codon we store the coordinates in ENST hash

			if(exists($hash6{$ENST_tok}))
			{
				foreach my $STOP_fields_tok(sort keys%{$hash6{$ENST_tok}})
				{
					my @STOP_tmp=split(/\t/,$STOP_fields_tok);
					foreach my $STOP_tmp_tok(@STOP_tmp)
					{
						if($STOP_tmp_tok=~/STOP:\d+:\[(.+)\]\-\[(.+)\]/)
						{
							my $begin=$1;
							my $end=$2;
							$ENST_hash{$strand_tok}{$string}{'STOP'}{$begin}{$end}=1;
							# print "'STOP'\t$begin\t$end**\n";
						}
					}
				}
			}
			
			# If the transcript has not a START codon their coordinates will be noted as NaN

			else
			{
				my $begin="NaN";
				my $end="NaN";
				$ENST_hash{$strand_tok}{$string}{'STOP'}{$begin}{$end}=1;
				# print "'STOP'\t$begin\t$end**\n";
			}
			
			# We recover the first coding position of the CDS. We don't use the first START position as not all the protein coding transcripts
			# have a canonical start but all of them have by definition CDS. When a transcript has a START and CDS the first coordinate of them both is the same.
			
			my @begin_CDS_tmp=sort{$a<=>$b} keys%{$order_hash{'CDS'}};
			my $counter_5_prime=0;
			my $counter_3_prime=0;
			
			# Ordering of the UTRs in 5 prime and 3 prime is strand dependent
			
			if($strand_tok eq '+')
			{
				my $START_CDS=$begin_CDS_tmp[0];
				foreach my $begin_UTR_tok(sort{$a<=>$b} keys%{$order_hash{'UTR'}})
				{
				foreach my $end_UTR_tok(sort{$a<=>$b} keys%{$order_hash{'UTR'}{$begin_UTR_tok}})
				{	
					if($begin_UTR_tok !~ 'NaN' && $end_UTR_tok !~ 'NaN')
					{
						if($begin_UTR_tok < $START_CDS)
						{
							# When strand is +; if the begin UTR coordinate is smaller than the START position of the CDS, then the UTR segment (there might be more than one) is UTR 5 prime
							
							$counter_5_prime++;
							$ENST_hash{$strand_tok}{$string}{'UTR_5_prime'}{$begin_UTR_tok}{$end_UTR_tok}=1;
							# print "AA'UTR_5_prime'\t$counter_5_prime\t$begin_UTR_tok\t$end_UTR_tok**\n";
						}
						elsif($begin_UTR_tok > $START_CDS)
						{
							# When strand is +; if the begin UTR coordinate is bigger than the START position of the CDS, then the UTR segment (there might be more than one) is UTR 3 prime
							
							$counter_3_prime++;
							$ENST_hash{$strand_tok}{$string}{'UTR_3_prime'}{$begin_UTR_tok}{$end_UTR_tok}=1;
							# print "AA'UTR_3_prime'\t$counter_3_prime\t$begin_UTR_tok\t$end_UTR_tok**\n";
						}
					}
					else
					{
						$counter_5_prime++;
						$counter_3_prime++;
						$ENST_hash{$strand_tok}{$string}{'UTR_5_prime'}{$begin_UTR_tok}{$end_UTR_tok}=1;
						$ENST_hash{$strand_tok}{$string}{'UTR_3_prime'}{$begin_UTR_tok}{$end_UTR_tok}=1;
						# print "BB'UTR_5_prime'\t$counter_5_prime\t$begin_UTR_tok\t$end_UTR_tok**\n";
						# print "BB'UTR_3_prime'\t$counter_3_prime\t$begin_UTR_tok\t$end_UTR_tok**\n";
					}
				}
				}
				if($counter_5_prime ==0)
				{
					my $begin_5_prime="NaN";
					my $end_5_prime="NaN";
					$ENST_hash{$strand_tok}{$string}{'UTR_5_prime'}{$begin_5_prime}{$end_5_prime}=1;
					# print "CC'UTR_5_prime'\t$counter_5_prime\t$begin_5_prime\t$end_5_prime**\n";
				}
				elsif($counter_3_prime ==0)
				{
					my $begin_3_prime="NaN";
					my $end_3_prime="NaN";
					$ENST_hash{$strand_tok}{$string}{'UTR_3_prime'}{$begin_3_prime}{$end_3_prime}=1;
					# print "CC'UTR_3_prime'\t$counter_3_prime\t$begin_3_prime\t$end_3_prime**\n";
				}
			}
			if($strand_tok eq '-')
			{
				my @begin_CDS_tmp_reversed=reverse(@begin_CDS_tmp);
				my $START_CDS=$begin_CDS_tmp_reversed[0];
				foreach my $begin_UTR_tok(sort{$a<=>$b} keys%{$order_hash{'UTR'}})
				{
				foreach my $end_UTR_tok(sort{$a<=>$b} keys%{$order_hash{'UTR'}{$begin_UTR_tok}})
				{
					if($begin_UTR_tok !~ 'NaN' && $end_UTR_tok !~ 'NaN')
					{
						# When strand is -; if the begin UTR coordinate is bigger than the START position of the CDS, then the UTR segment (there might be more than one) is UTR 5 prime

						if($begin_UTR_tok > $START_CDS)
						{
							$counter_5_prime++;
							$ENST_hash{$strand_tok}{$string}{'UTR_5_prime'}{$begin_UTR_tok}{$end_UTR_tok}=1;
							# print "AA'UTR_5_prime'\t$counter_5_prime\t$begin_UTR_tok\t$end_UTR_tok**\n";
						}
						
						# When strand is -; if the begin UTR coordinate is smaller than the START position of the CDS, then the UTR segment (there might be more than one) is UTR 3 prime
						
						elsif($begin_UTR_tok < $START_CDS)
						{
							$counter_3_prime++;
							$ENST_hash{$strand_tok}{$string}{'UTR_3_prime'}{$begin_UTR_tok}{$end_UTR_tok}=1;
							# print "AA'UTR_3_prime'\t$counter_3_prime\t$begin_UTR_tok\t$end_UTR_tok**\n";
						}
					}
					else
					{
						$counter_5_prime++;
						$counter_3_prime++;
						$ENST_hash{$strand_tok}{$string}{'UTR_5_prime'}{$begin_UTR_tok}{$end_UTR_tok}=1;
						$ENST_hash{$strand_tok}{$string}{'UTR_3_prime'}{$begin_UTR_tok}{$end_UTR_tok}=1;
						# print "BB'UTR_5_prime'\t$counter_5_prime\t$begin_UTR_tok\t$end_UTR_tok**\n";
						# print "BB'UTR_3_prime'\t$counter_3_prime\t$begin_UTR_tok\t$end_UTR_tok**\n";
					}
				}
				}
				if($counter_5_prime ==0)
				{
					my $begin_5_prime="NaN";
					my $end_5_prime="NaN";
					$ENST_hash{$strand_tok}{$string}{'UTR_5_prime'}{$begin_5_prime}{$end_5_prime}=1;
					# print "CC'UTR_5_prime'\t$counter_5_prime\t$begin_5_prime\t$end_5_prime**\n";
				}
				elsif($counter_3_prime ==0)
				{
					my $begin_3_prime="NaN";
					my $end_3_prime="NaN";
					$ENST_hash{$strand_tok}{$string}{'UTR_3_prime'}{$begin_3_prime}{$end_3_prime}=1;
					# print "CC'UTR_3_prime'\t$counter_3_prime\t$begin_3_prime\t$end_3_prime**\n";
				}
			}
		}#
		}
		}
	}
	}	
	}

$time='['. timestamp(). ']'."\n";
print "Start PRINTING:$time\n";

# Now we print the results of the hash of hashes in a strand specific way; + is printed in an increasing way (sort); - is printed in a decreasing way (reverse sort).

if (open(OUTPUT, '>'.$output))
{
	foreach my $strand_tok(sort keys %ENST_hash)
	{
		if($strand_tok eq '+')
		{
			foreach my $string_tok(sort keys %{$ENST_hash{$strand_tok}})
			{
				print OUTPUT "$string_tok\t";
				
				foreach my $begin_5_prime_tok(sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'UTR_5_prime'}})
				{
				foreach my $end_5_prime_tok(sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'UTR_5_prime'}{$begin_5_prime_tok}})
				{
					print OUTPUT "UTR_5_prime:[$begin_5_prime_tok]-[$end_5_prime_tok]\t";
				}
				}
				foreach my $begin_START_tok(sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'START'}})
				{
				foreach my $end_START_tok(sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'START'}{$begin_START_tok}})
				{
					print OUTPUT "START:[$begin_START_tok]-[$end_START_tok]\t";
				}
				}
				foreach my $begin_CDS_tok(sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'CDS'}})
				{
				foreach my $end_CDS_tok(sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'CDS'}{$begin_CDS_tok}})
				{
					print OUTPUT "CDS:[$begin_CDS_tok]-[$end_CDS_tok]\t";
				}
				}
				foreach my $begin_STOP_tok(sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'STOP'}})
				{
				foreach my $end_STOP_tok(sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'STOP'}{$begin_STOP_tok}})
				{
					print OUTPUT "STOP:[$begin_STOP_tok]-[$end_STOP_tok]\t";
				}
				}
				foreach my $begin_3_prime_tok(sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'UTR_3_prime'}})
				{
				foreach my $end_3_prime_tok(sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'UTR_3_prime'}{$begin_3_prime_tok}})
				{
					print OUTPUT "UTR_3_prime:[$begin_3_prime_tok]-[$end_3_prime_tok]\t";
				}
				}
				print OUTPUT "\n";
			}
		}
		elsif($strand_tok eq '-')
		{
			foreach my $string_tok(sort keys %{$ENST_hash{$strand_tok}})
			{
				print OUTPUT "$string_tok\t";
				
				foreach my $begin_5_prime_tok(reverse sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'UTR_5_prime'}})
				{
				foreach my $end_5_prime_tok(reverse sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'UTR_5_prime'}{$begin_5_prime_tok}})
				{
					print OUTPUT "UTR_5_prime:[$begin_5_prime_tok]-[$end_5_prime_tok]\t";
				}
				}
				foreach my $begin_START_tok(reverse sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'START'}})
				{
				foreach my $end_START_tok(reverse sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'START'}{$begin_START_tok}})
				{
					print OUTPUT "START:[$begin_START_tok]-[$end_START_tok]\t";
				}
				}
				foreach my $begin_CDS_tok(reverse sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'CDS'}})
				{
				foreach my $end_CDS_tok(reverse sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'CDS'}{$begin_CDS_tok}})
				{
					print OUTPUT "CDS:[$begin_CDS_tok]-[$end_CDS_tok]\t";
				}
				}
				foreach my $begin_STOP_tok(reverse sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'STOP'}})
				{
				foreach my $end_STOP_tok(reverse sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'STOP'}{$begin_STOP_tok}})
				{
					print OUTPUT "STOP:[$begin_STOP_tok]-[$end_STOP_tok]\t";
				}
				}
				foreach my $begin_3_prime_tok(reverse sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'UTR_3_prime'}})
				{
				foreach my $end_3_prime_tok(reverse sort{$a<=>$b} keys %{$ENST_hash{$strand_tok}{$string_tok}{'UTR_3_prime'}{$begin_3_prime_tok}})
				{
					print OUTPUT "UTR_3_prime:[$begin_3_prime_tok]-[$end_3_prime_tok]\t";
				}
				}
				print OUTPUT "\n";
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
