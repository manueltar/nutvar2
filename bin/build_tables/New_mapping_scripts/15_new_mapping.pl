##	Script to ascertain to divide ENSEMBL transcripts between those for whom the UniProt displayed isoform matches excatly
## or with an N-ter offset in the ENSEMBL polypetide and those that do no have a match at all. 2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne)

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;

my %hash1=();
my %hash_no_aligned=();

my $input=$ARGV[0];
my $output1=$ARGV[1];
my $output2=$ARGV[2];
my $output3=$ARGV[3];
my $Flag=0;

my $time='['. timestamp(). ']'."\n";
print "Tiempo de incio de carga_hash1:$time\n";
	
my $AC="NaN";
my $ENSG="NaN";
my $ENST="NaN";
my $coordinates_up="NaN";
my $coordinates_down="NaN";

# We open the file with all the correspondences
		
if(open(INPUT,$input))
{
	# Input= Equiv_ENSG_seq.txt
	
	#~ UNIPROT>O43657  ENSG00000000003
	#~ MASPSRRLQTKPVITCFKSVLLIYTFIFWITGVILLAVGIWGKVSLENYFSLLNEKATNVPFVLIATGTVIILLGTFGCFATCRASAWMLKLYAMFLTLVFLVELVAAIVGFVFRHEIKNSFKNNYEKALKQYNSTGDYRSHA
	#~ UPSTREAM=0;DOWNSTREAM=0>O43657  ENSG00000000003 ENST00000373020
	#~ MASPSRRLQTKPVITCFKSVLLIYTFIFWITGVILLAVGIWGKVSLENYFSLLNEKATNVPFVLIATGTVIILLGTFGCFATCRASAWMLKLYAMFLTLVFLVELVAAIVGFVFRHEIKNSFKNNYEKALKQYNSTGDYRSHA
	#~ 
	#~ UNIPROT>Q9H2S6  ENSG00000000005
	#~ MAKNPPENCEDCHILNAEAFKSKKICKSLKICGLVFGILALTLIVLFWGSKHFWPEVPKKAYDMEHTFYSNGEKKKIYMEIDPVTRTEIFRSGNGTDETLEVHDFKNGYTGIYFVGLQKCFIKTQIKVIPEFSEPEEEIDENE
	#~ UPSTREAM=0;DOWNSTREAM=0>Q9H2S6  ENSG00000000005 ENST00000373031
	#~ MAKNPPENCEDCHILNAEAFKSKKICKSLKICGLVFGILALTLIVLFWGSKHFWPEVPKKAYDMEHTFYSNGEKKKIYMEIDPVTRTEIFRSGNGTDETLEVHDFKNGYTGIYFVGLQKCFIKTQIKVIPEFSEPEEEIDENE
	#~ 
	#~ UNIPROT>O60762  ENSG00000000419
	#~ MASLEVSRSPRRSRRELEVRSPRQNKYSVLLPTYNERENLPLIVWLLVKSFSESGINYEIIIIDDGSPDGTRDVAEQLEKIYGSDRILLRPREKKLGLGTAYIHGMKHATGNYIIIMDADLSHHPKFIPEFIRKQKEGNFDIV
	#~ UPSTREAM=0;DOWNSTREAM=0>O60762  ENSG00000000419 ENST00000371588
	#~ MASLEVSRSPRRSRRELEVRSPRQNKYSVLLPTYNERENLPLIVWLLVKSFSESGINYEIIIIDDGSPDGTRDVAEQLEKIYGSDRILLRPREKKLGLGTAYIHGMKHATGNYIIIMDADLSHHPKFIPEFIRKQKEGNFDIV
	#~ UPSTREAM=NO_ALIGN;DOWNSTREAM=NO_ALIGN>O60762    ENSG00000000419 ENST00000371582
	#~ MASLEVSRSPRRSRRELEVRSPRQNKYSVLLPTYNERENLPLIVWLLVKSFSESGINYEIIIIDDGSPDGTRDVAEQLEKIYGSDRILLRPREKKLGLGTAYIHGMKHATGNYIIIMDADLSHHPKFIPEFIRKQKEGNFDIV
	#~ UPSTREAM=NO_ALIGN;DOWNSTREAM=NO_ALIGN>O60762    ENSG00000000419 ENST00000371584
	#~ XASLEVSRSPRRSRRELEVRSPRQNKYSVLLPTYNERENLPLIVWLLVKSFSESGINYEIIIIDDGSPDGTRDVAEQLEKIYGSDRILLRPREKKLGLGTAYIHGMKHATGNYIIIMDADLSHHPKFIPEFIRKQKEGNFDIV
	#~ UPSTREAM=NO_ALIGN;DOWNSTREAM=NO_ALIGN>O60762    ENSG00000000419 ENST00000413082
	#~ MASLEVSRSPRRSRRELEVRSPRQNKYSVLLPTYNERENLPLIVWLLVKSFSESGINYEIIIIDDGSPDGTRDVAEQLEKIYGSDRILLRPREKKLGLGTAYIHGMKHATGNYIIIMDADLSHHPKFIPEFISDGVLPCCPGW
	#~ UPSTREAM=NO_ALIGN;DOWNSTREAM=NO_ALIGN>O60762    ENSG00000000419 ENST00000371583
	#~ MASLEVSRSPRRSRRELEVRSPRQNKYSVLLPTYNERENLPLIVWLLVKSFSESGINYEIIIIDDGSPDGTRDVAEQLEKIYGSDRILLRPREKKLGLGTAYIHGMKHATGNYIIIMDADLSHHPKFIPEFISDGVLPCCPGW
	#~ 
	#~ UNIPROT>Q9Y6X5  ENSG00000001561
	#~ MKLLVILLFSGLITGFRSDSSSSLPPKLLLVSFDGFRADYLKNYEFPHLQNFIKEGVLVEHVKNVFITKTFPNHYSIVTGLYEESHGIVANSMYDAVTKKHFSDSNDKDPFWWNEAVPIWVTNQLQENRSSAAAMWPGTDVPI
	#~ UPSTREAM=0;DOWNSTREAM=0>Q9Y6X5  ENSG00000001561 ENST00000321037
	#~ MKLLVILLFSGLITGFRSDSSSSLPPKLLLVSFDGFRADYLKNYEFPHLQNFIKEGVLVEHVKNVFITKTFPNHYSIVTGLYEESHGIVANSMYDAVTKKHFSDSNDKDPFWWNEAVPIWVTNQLQENRSSAAAMWPGTDVPI
	
	while(my $line=<INPUT>)
	{
		chomp $line;
		#~ print "$line\n";
		if ($line=~/^(.+)>(.+)/)
		{
			#~ print "Hello_world:$line\n";
			my $label=$1;
			my $AC_ENSG_ENST=$2;
			
			# If the line does not correspond to a UniProt identifier
			
			if ($label ne 'UNIPROT')
			{
				my @label_tmp=split(";",$label);
				my @tmp=split("\t",$AC_ENSG_ENST);
				
				# Obtain the corresponding UniProt identifier and the ENSG and ENST id's
				
				$AC=$tmp[0];
				$ENSG=$tmp[1];
				$ENST=$tmp[2];
				
				# Obtain the offsets and separate them if they are No Align
				
				foreach my $label_tmp_tok(@label_tmp)
				{
					#~ print "**$label_tmp_tok**\n";
					if($label_tmp_tok =~ /UPSTREAM=(.+)/)
					{
						my $upstream=$1;
						#~ print "**$upstream**\n";
						if($upstream ne 'NO_ALIGN')
						{
							$coordinates_up=$upstream;
							$Flag=1;
							#~ print "**$AC\t$ENSG\t'ALIGN'\t$coordinates_up\t$Flag**\n";
						}
						elsif($upstream eq 'NO_ALIGN')
						{
							$Flag=2;
							#~ print "$AC\t$ENSG\t'NO_ALIGN'\t$ENST\t$Flag\n";
						}
					}
					elsif($label_tmp_tok =~ /DOWNSTREAM=(.+)/)
					{
						my $downstream=$1;
						if($downstream =~/\d+/)
						{
							$coordinates_down=$downstream;
							#~ print "**$AC\t$ENSG\t'ALIGN'\t$coordinates_down\t$Flag**\n";
						}
					}
				}
				
			}
			
			# If the sequence belongs to UniProt displayed isoform
			 
			elsif($label eq 'UNIPROT')
			{
				my @tmp=split("\t",$AC_ENSG_ENST);
				$AC=$tmp[0];
				$ENSG=$tmp[1];
				$Flag=3;
				#~ print "$AC\t$ENSG\t'UNIPROT'\t$Flag\n";
			}
		}
		
		# Extract all the peptide sequences
		
		elsif($line !~ /^(.+)>(.+)/)
		{
			# Store each sequence in a hash o hashes identifying ENSEMBL isoforms that match the UniProt (Align) and ENSEMBL isoforms that do not match the UniProt displayed
			# isoform (No Align). Of note for the same genes there are usually isoforms matching the UniProt and shorter isoforms that do not.
			# Store also the UniProt displayed isoform and its sequence.
			if($Flag==1)
			{
				unless($coordinates_up eq 'NaN' && $coordinates_down eq 'NaN')
				{
					$hash1{$AC}{$ENSG}{'ALIGN'}{$ENST}{$coordinates_up}{$coordinates_down}=1;
					$Flag=0;
					#~ print "RESUME:$AC\t$ENSG\t'ALIGN'\t$ENST\t$coordinates_up\t$coordinates_down\t$Flag\n";
				}	
			}
			elsif($Flag==2)
			{
					$hash1{$AC}{$ENSG}{'NO_ALIGN'}{$ENST}{$line}=1;
					$Flag=0;
					#~ print "RESUME:$AC\t$ENSG\t'NO_ALIGN'\t$ENST\t$line\t$Flag\n";	
			}
			elsif($Flag==3)
			{
					$hash1{$AC}{$ENSG}{'UNIPROT'}{$line}=1;
					$Flag=0;
					#~ print "$AC\t$ENSG\t'UNIPROT'\t$line\t$Flag\n";
			}
		}
	}
}else{print "Unable to open $input\n";}

$time='['. timestamp(). ']'."\n";
print "Tiempo de IMPRESIÓN_ALIGNED:$time\n";

if(open(OUTPUT, '>'.$output1))
{
foreach my $AC_tok(sort keys %hash1)
{
foreach my $ENSG_tok(sort keys %{$hash1{$AC_tok}})
{
	my @tmp_aligned=sort keys %{$hash1{$AC_tok}{$ENSG_tok}{'ALIGN'}};
	
	# If there are none ENSEMBL polypetides matching the UniProt displayed isoform for a given gene then store those ids in a hash
	
	if (scalar(@tmp_aligned) == 0)
	{
		$hash_no_aligned{$AC_tok}{$ENSG_tok}=1;
	}
	
	# If there are ENSEMBL polypetides matching the UniProt isoform sequence
	
	elsif(scalar(@tmp_aligned) >= 0)
	{
		# Define an order hash on a key basis
		
		my %order_hash=();
		#~ print "El array_aligned es:$AC_tok\t$ENSG_tok\t@tmp_aligned\n";
		foreach my $ENST_tok(@tmp_aligned)
		{
		foreach my $upstream_coordinate_tok(sort keys %{$hash1{$AC_tok}{$ENSG_tok}{'ALIGN'}{$ENST_tok}})
		{
		foreach my $downstream_coordinate_tok(sort keys %{$hash1{$AC_tok}{$ENSG_tok}{'ALIGN'}{$ENST_tok}{$upstream_coordinate_tok}})
		{
			# Store in the order hash all the offsets and the corresponding ENST id
			
			$order_hash{$upstream_coordinate_tok}{$downstream_coordinate_tok}{$ENST_tok}=1;
			#~ print "ORDER_HASH:$upstream_coordinate_tok\t$downstream_coordinate_tok\t$ENST_tok\n";
			
		}
		}
		}
		my @upstream_tmp=sort{$a<=>$b}keys%order_hash;
		my @downstream_tmp=sort{$a<=>$b}keys%{$order_hash{$upstream_tmp[0]}};
		my @ENST_tmp=sort keys %{$order_hash{$upstream_tmp[0]}{$downstream_tmp[0]}};
		
		# Select from the order hash the first trasncript with the lowest N-ter and C-ter offsets. This is the fists trasncript with 0,0 offsets for 18.074 transcripts
		# while in 132 the offsets are no 0 so the minimum offset is selected 
		
		print  OUTPUT "$AC_tok\t$ENSG_tok\t$ENST_tmp[0]\tUPSTREAM=$upstream_tmp[0]\tDOWNSTREAM=$downstream_tmp[0]\n";
	}

}
}
}

$time='['. timestamp(). ']'."\n";
print "Tiempo de Reannealing:$time\n";

if(open(OUTPUT2, '>'.$output2) && open(OUTPUT3, '>'.$output3))
{
	# For those ENSEMBL genes whose different polypetide sequences do not match their cognate UniProt displayed isoform
	# that we have stored in the hash %hash_no_aligned
	
foreach my $AC_tok(sort keys %hash_no_aligned)
{
foreach my $ENSG_tok(sort keys %{$hash_no_aligned{$AC_tok}})
{
foreach  my $UNIPROT_seq_tok(sort keys %{$hash1{$AC_tok}{$ENSG_tok}{'UNIPROT'}})	
{
	# As we are going to align them using BLAST we need the UniProt displayed isoform to align them against.
	# We need to print the sequence in FASTA format (70 residues per line)
	
	print OUTPUT2 ">sp|$AC_tok|$ENSG_tok"."_UNIPROT\n";
	while(length($UNIPROT_seq_tok)>70) 
	{
		print OUTPUT2 substr($UNIPROT_seq_tok,0,70),"\n";
		$UNIPROT_seq_tok=substr($UNIPROT_seq_tok,70);
	}
	print OUTPUT2 "$UNIPROT_seq_tok\n";
	
	# We print all the corresponding ENSEMBL polypetides for the same gene with therir polypetide sequence in FASTA format.
	
foreach  my $ENST_tok(sort keys %{$hash1{$AC_tok}{$ENSG_tok}{'NO_ALIGN'}})	
{
foreach  my $ENST_seq_tok(sort keys %{$hash1{$AC_tok}{$ENSG_tok}{'NO_ALIGN'}{$ENST_tok}})	
{
	print OUTPUT3 ">sp|$AC_tok|$ENSG_tok"."_"."$ENST_tok\n";
	while(length($ENST_seq_tok)>70) 
	{
		print OUTPUT3 substr($ENST_seq_tok,0,70),"\n";
		$ENST_seq_tok=substr($ENST_seq_tok,70);
	}
	print OUTPUT3 "$ENST_seq_tok\n";
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



