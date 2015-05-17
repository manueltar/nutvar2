##	Script to convert .vcf into .vcf with joint calls expressed as minimal representation and undo cryptic alleles. 2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne) 

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;

# We use the Registry of ENSEMBL to resolve cryptic alternative alleles (<DEL> or < . >)

use Bio::EnsEMBL::Registry;


# Input/output declaration

my $input=$ARGV[0];
my $output_minimal_representation= $ARGV[1];

## Here we declare the Heather-to-index (HI) and Index-to-Heather hashes and the global variables for alleles that will be trimmed
## to yield minimal representation.

my %HI0=();
my %IH0=();
my $REF_post_trimming="NaN";
my $ALT_post_trimming="NaN";
my $POS_trimmed="NaN";

my $time='['. timestamp(). ']'."\n";
print "Start loading the Registry:$time\n";

my $registry = "Bio::EnsEMBL::Registry";
$registry->load_registry_from_url('mysql://anonymous@ensembldb.ensembl.org');


$time='['. timestamp(). ']'."\n";
print "Start OPENING FILES:$time\n";

# We create a counter to print the number of iterations we conduct with ENSEMBL API

my $l=0;

if(open(INPUT, $input) && open(OUTPUT, '>'.$output_minimal_representation))
{
	# We print the vcf standard header line in the output file 
	
	print OUTPUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";	
	
	# We open first input file
while (my $line = <INPUT>)
	{
		
		# Input: original vcf of the user. In the test provided: example.vcf
		
		#~ #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
		#~ 1       10583   rs58108140      G       A       100     PASS    AVGPOST=0.7707;RSQ=0.4319;LDAF=0.2327;ERATE=0.0161;AN=2184;VT=SNP;AA=.;THETA=0.0046;AC=314;SNPSOURCE=LOWCOV;AF=0.14;ASN_AF=0.13;AMR_AF=0.17;AFR_AF=0.04;EUR_AF=0.21
		#~ 1       10611   rs189107123     C       G       100     PASS    AN=2184;THETA=0.0077;VT=SNP;AA=.;AC=41;ERATE=0.0048;SNPSOURCE=LOWCOV;AVGPOST=0.9330;LDAF=0.0479;RSQ=0.3475;AF=0.02;ASN_AF=0.01;AMR_AF=0.03;AFR_AF=0.01;EUR_AF=0.02
		#~ 1       13302   rs180734498     C       T       100     PASS    THETA=0.0048;AN=2184;AC=249;VT=SNP;AA=.;RSQ=0.6281;LDAF=0.1573;SNPSOURCE=LOWCOV;AVGPOST=0.8895;ERATE=0.0058;AF=0.11;ASN_AF=0.02;AMR_AF=0.08;AFR_AF=0.21;EUR_AF=0.14


##	From input file we load HI & IH hashes with the heather fields.

		chomp ($line);
		if (index($line,'#CHROM')==0) 
				
			{
				$line=~/^#(.+)/;
				my @tmp_HI0 = split ("\t",$1);
				#~ my $string=join("*",@tmp_HI0);
				#~ print OUTPUT "The array is $string\n";
				for(my $i=0;$i<scalar(@tmp_HI0);$i++)
										
				{
					$HI0{$tmp_HI0[$i]}=$i;
					$IH0{$i}=$tmp_HI0[$i];
					#print OUTPUT "$HI0{$tmp_HI0[$i]} \n";
					#print OUTPUT "$IH0{$i} \n";
				}
			}

##	From lines carrying data we store every column in an array
						
		elsif ($line !~ /^#/)
			{	
				##	Here we define and re-initialize the arrays on a per-line basis.
				
				my %hash_REF_ALT=();
				my $trim="NaN";
				my $pos_trimmed="NaN";

				#print OUTPUT "$line\n";
				my @tmp= split(/\t/,$line);
				#~ my $string=join("*",@tmp);
				#~ print OUTPUT "The array is $string\n";
				
				my $CHROM="NaN";if (exists($HI0{"CHROM"})){$CHROM=$tmp[$HI0{"CHROM"}];}unless(defined($CHROM)){$CHROM="NaNein"; print "ERROR in CHROM\n";}
				my $POS="NaN";if (exists($HI0{"POS"})){$POS=$tmp[$HI0{"POS"}];}unless(defined($POS)){$POS="NaNein"; print "ERROR in POS\n";}
				my $ID="NaN";if (exists($HI0{"ID"})){$ID=$tmp[$HI0{"ID"}];}unless(defined($ID)){$ID="NaNein"; print "ERROR in ID\n";}
				my $REF="NaN";if (exists($HI0{"REF"})){$REF=$tmp[$HI0{"REF"}];}unless(defined($REF)){$REF="NaNein"; print "ERROR in REF\n";}
				my $ALT="NaN";if (exists($HI0{"ALT"})){$ALT=$tmp[$HI0{"ALT"}];}unless(defined($ALT)){$ALT="NaNein"; print "ERROR in ALT\n";}
				my $QUAL="NaN";if (exists($HI0{"QUAL"})){$QUAL=$tmp[$HI0{"QUAL"}];}unless(defined($QUAL)){$QUAL="NaNein"; print "ERROR in QUAL\n";}
				my $FILTER="NaN";if (exists($HI0{"FILTER"})){$FILTER=$tmp[$HI0{"FILTER"}];}unless(defined($FILTER)){$FILTER="NaNein"; print "ERROR in FILTER\n";}
				my $INFO="NaN";if (exists($HI0{"INFO"})){$INFO=$tmp[$HI0{"INFO"}];}unless(defined($INFO)){$INFO="NaNein"; print "ERROR in INFO\n";}
				my $FORMAT="NaN";if (exists($HI0{"FORMAT"})){$FORMAT=$tmp[$HI0{"FORMAT"}];}unless(defined($FORMAT)){$FORMAT="NaNein"; print "ERROR in FORMAT\n";}
				my $NA12878="NaN";if (exists($HI0{"NA12878"})){$NA12878=$tmp[$HI0{"NA12878"}];}unless(defined($NA12878)){$NA12878="NaNein"; print "ERROR in NA12878\n";}
				
				#~ print OUTPUT "$CHROM\t$POS\t$ID\t$REF\t$ALT\n";
				#~ print OUTPUT "La otra línea es:$QUAL\t$FILTER\t$INFO\t$FORMAT\t$NA12878\n";
				#~ exit;
				
				# We declare an array to store possible alternative allele joint calls
				# p.e. ATG,A,ATTTGGG
				
				my @ALT_tmp=();
				
				# Should there be joint calls we split them
				
				if($ALT=~/\,/)
				{
					@ALT_tmp=split(',',$ALT);
				}
				
				# Else we store the simple allele call
				else {push (@ALT_tmp,$ALT);}
				
				foreach my $ALT_tmp_tok(@ALT_tmp)
				{
					$ALT_tmp_tok=uc($ALT_tmp_tok);
					
					# Should there be any cryptic allele (<DEL> or < . >)
					
					if ($ALT_tmp_tok=~/DEL/ || $ALT_tmp_tok=~/\./)
					{
						$l++;
						
						# We initialize the API iteration
						
						print "Starting API $l\n";
						
						# The position we want to retrieve is the previuos one to the alternative allele
						
						my $POS_transferred=$POS-1;
						print "Positions:$POS\t$POS_transferred\n";

						# Get the human *core* Adaptor for Slices
						my $slice_adaptor =Bio::EnsEMBL::Registry->get_adaptor("human", "core", "Slice");

						# Get the slice corresponding to the region of interest
						my $slice = $slice_adaptor->fetch_by_region("chromosome", $CHROM, $POS_transferred, $POS_transferred);
						
						# We recover the nucleotide in the position of interest
						my $initial_nucleotide=$slice->seq();
						print "The nucleotide is:$initial_nucleotide\n";
						
						# We define the new ref and alternative alleles in the minimal representation (see M&M); A->DEL turns into TA->T
						my $REF_reconstructed=join('',$initial_nucleotide,$REF);
						my $ALT_reconstructed=$initial_nucleotide;
						
						# We store values in a hash of hashes
						$hash_REF_ALT{$CHROM}{$POS_transferred}{$ID}{$REF_reconstructed}{$ALT_reconstructed}{$QUAL}{$FILTER}{$INFO}{$FORMAT}{$NA12878}=1;
					}
					else
					{
						# We store values in a hash of hashes
						$hash_REF_ALT{$CHROM}{$POS}{$ID}{$REF}{$ALT_tmp_tok}{$QUAL}{$FILTER}{$INFO}{$FORMAT}{$NA12878}=1;
					}
				}
				
				# We undo the hash of hashes on a per line basis
				foreach my $CHROM_tok(sort{$a<=>$b} keys %hash_REF_ALT)
				{
					foreach my $POS_tok(sort{$a<=>$b} keys %{$hash_REF_ALT{$CHROM_tok}})
					{
					foreach my $ID_tok(sort keys %{$hash_REF_ALT{$CHROM_tok}{$POS_tok}})
					{
					foreach my $REF_tok(sort keys %{$hash_REF_ALT{$CHROM_tok}{$POS_tok}{$ID_tok}})
					{
					foreach my $ALT_tok(sort keys %{$hash_REF_ALT{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}})
					{
					foreach my $QUAL_tok(sort keys %{$hash_REF_ALT{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}})
					{	
					foreach my $FILTER_tok(sort keys %{$hash_REF_ALT{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$QUAL_tok}})
					{	
					foreach my $INFO_tok(sort keys %{$hash_REF_ALT{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$QUAL_tok}{$FILTER_tok}})
					{
					foreach my $FORMAT_tok(sort keys %{$hash_REF_ALT{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$QUAL_tok}{$FILTER_tok}{$INFO_tok}})
					{	
					foreach my $NA12878_tok(sort keys %{$hash_REF_ALT{$CHROM_tok}{$POS_tok}{$ID_tok}{$REF_tok}{$ALT_tok}{$QUAL_tok}{$FILTER_tok}{$INFO_tok}{$FORMAT_tok}})
					{	
						# If REF and ALT alleles are already in their minimal representation
						
						if (length($REF_tok)==1 && length($ALT_tok)==1)
								{
										#~ print OUTPUT "El alelo partido por comas es y que cumple las NOMODIF es:$ALT_tok\n";
										# We print OUTPUT a label of it being left untouched and add it to info field
										$trim= "No_modification_to_get_minimal_representation";
										my $string=join(";",$INFO_tok,$trim);
										# We print OUTPUT the line in vcf format in the output file
										print OUTPUT "$CHROM_tok\t$POS_tok\t$ID_tok\t$REF_tok\t$ALT_tok\t$QUAL_tok\t$FILTER_tok\t$string\n";	
								}
						elsif (length($REF_tok) == 1 && length($ALT_tok) != 1)
								{
										#~ print OUTPUT "El alelo partido por comas es y que cumple las NOMODIF es:$ALT_tok\n";
										# We print OUTPUT a label of it being left untouched and add it to info field
										$trim= "No_modification_to_get_minimal_representation";
										my $string=join(";",$INFO_tok,$trim);
										# We print OUTPUT the line in vcf format in the output file
										print OUTPUT "$CHROM_tok\t$POS_tok\t$ID_tok\t$REF_tok\t$ALT_tok\t$QUAL_tok\t$FILTER_tok\t$string\n";
								}
						elsif (length($REF_tok) != 1 && length($ALT_tok) == 1)
								{
										#~ print OUTPUT "El alelo partido por comas es y que cumple las NOMODIF es:$ALT_tok\n";
										# We print OUTPUT a label of it being left untouched and add it to info field
										$trim= "No_modification_to_get_minimal_representation";
										my $string=join(";",$INFO_tok,$trim);
										# We print OUTPUT the line in vcf format in the output file
										print OUTPUT "$CHROM_tok\t$POS_tok\t$ID_tok\t$REF_tok\t$ALT_tok\t$QUAL_tok\t$FILTER_tok\t$string\n";
								}
								
						# If REF and ALT alleles are not already in their minimal representation
						elsif(length($REF_tok) > 1 && length($ALT_tok) > 1)
								{
										##	We split by letters the components of REF_tok and ALT

										#~ print OUTPUT "$REF_tok\t$ALT_tok\n";
										my @REF_tok_tmp=split('',$REF_tok);
										#~ print OUTPUT "@REF_tok_tmp\n";
										my @ALT_tmp=split('',$ALT_tok);
										#~ print OUTPUT "@ALT_tmp\n";
										$POS_trimmed=$POS_tok;
										
										# Conditions of trimming downstream; same element in last position of REF_tok and ALT and REF_tok and ALT more than one in length. 
										
										while (($REF_tok_tmp[scalar(@REF_tok_tmp)-1] eq $ALT_tmp[scalar(@ALT_tmp)-1] && scalar(@REF_tok_tmp)>1 && scalar(@ALT_tmp)>1))
											{
												## Extract last element of the array.
												#print OUTPUT "Hello_world_1\n";
												pop (@REF_tok_tmp);
												pop (@ALT_tmp);
											}
										
										# Conditions of trimming upstream; same element in last position of REF_tok and ALT and REF_tok and ALT more than one in length.
										
										while ($REF_tok_tmp[0] eq $ALT_tmp[0] && scalar(@REF_tok_tmp)>1 && scalar(@ALT_tmp)>1)
											{
												#print OUTPUT "Hello_world_2\n";
												shift (@REF_tok_tmp);
												shift (@ALT_tmp);
												$POS_trimmed++;
												#print OUTPUT "$POS_tok\n";
												#print OUTPUT "$POS_tok\n";
											}
										
										# Once the trimming has been completed we recompose whats been left in each array
										
										$REF_post_trimming=join('',@REF_tok_tmp);
										$ALT_post_trimming=join('',@ALT_tmp);
										
										# We print OUTPUT a label of it having been modified and add it to info field
										$trim= "This_line_has_been_modified_to_get_minimal_representation";
										my $string=join(";",$INFO_tok,$trim);
										# We print OUTPUT the line in vcf format in the output file
										print OUTPUT "$CHROM_tok\t$POS_tok\t$ID_tok\t$REF_tok\t$ALT_tok\t$QUAL_tok\t$FILTER_tok\t$string\n";
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
		
}else {print "unable to open $output_minimal_representation \n";}
	
##or $output_minimal_representation\n";}

$time='['. timestamp(). ']'."\n";
print "THE END:$time\n";

sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
