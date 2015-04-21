use strict;
use warnings;
use Time::localtime;

my %hash1=();
my %hash2=();



my $time="NaN";

$time='['. timestamp(). ']'."\n";
print "Tiempo de inicio de carga_hash1:$time\n";

my $input1=$ARGV[0];
my $output1=$ARGV[1];
my $input2=$ARGV[2];
my $output2=$ARGV[3];


if (open(INPUT1, $input1))
{
	## Input1= UNIPROT_GLOBAL_PRE_UNIQUE.txt

#	ENSG00000000003 O43657  ENST00000373020 ENSP00000362111 Unique  Displayed       FEATURE:Tetraspanin-6.__0__735__CHAIN   IPR:IPR000301__Tetraspa
#	ENSG00000000419 O60762  ENST00000371588 ENSP00000360644 Unique  Displayed       FEATURE:Dolichol-phosphate mannosyltransferase__3__777__CHAIN
#	ENSG00000000457 Q8IZE3  ENST00000367770 ENSP00000356744 Q8IZE3-1        Displayed       FEATURE:HEAT 1.__594__120__REPEAT       FEATURE:HEAT 2.

	while(my $line=<INPUT1>)
	{
		chomp($line);
		my @tmp=split("\t",$line);
		my $AC=$tmp[1];
		my $ENSG=$tmp[0];
		my $ENST=$tmp[2];
		my $ENSP=$tmp[3];
		my $IPR="NaN";
		my $PFAM="NaN";
		my $FEATURE="NaN";
		foreach my $tmp_tok(@tmp)
		{
			if($tmp_tok=~/^IPR:/)
			{
				$IPR=$tmp_tok;
				$hash1{$AC}{$ENSG}{$ENST}{$ENSP}{'IPR'}{$IPR}=1;
				#print OUTPUT"$AC\t$ENSG\t$ENST\t$ENSP"."IPR\t"."$IPR\n";
			}
			#elsif($tmp_tok=~/^PFAM:/)
			#{
			#	$PFAM=$tmp_tok;
			#	$hash1{$AC}{$ENSG}{$ENST}{$ENSP}{'PFAM'}{$PFAM}=1;
			#}
			elsif($tmp_tok=~/^FEATURE:/)
			{
				$FEATURE=$tmp_tok;
				$hash1{$AC}{$ENSG}{$ENST}{$ENSP}{'FEATURE'}{$FEATURE}=1;
			}
			elsif($tmp_tok=~/^CCDS:/)
			{
				$FEATURE=$tmp_tok;
				$hash1{$AC}{$ENSG}{$ENST}{$ENSP}{'CCDS'}{$FEATURE}=1;
			}
		}
	}
}else {print OUTPUT "impossible to open INPUT1\n";die;}


$time='['. timestamp(). ']'."\n";
print "Tiempo de inicio de impresión1:$time\n";

if (open(OUTPUT, '>'.$output1))
{
	foreach my $AC_tok(sort keys %hash1)
	{
		#print OUTPUT "**************************$AC_tok\n";
		foreach my $ENSG_tok(sort keys %{$hash1{$AC_tok}})
		{
			#print OUTPUT "**************************$ENSG_tok\n";
			my @ENST_tmp=sort keys %{$hash1{$AC_tok}{$ENSG_tok}};
			#print OUTPUT "**************************$ENST_tmp[0]\n";
			foreach my $ENSP_tok(sort keys %{$hash1{$AC_tok}{$ENSG_tok}{$ENST_tmp[0]}})
			{
				print OUTPUT "$AC_tok\t$ENSG_tok\t$ENST_tmp[0]\t$ENSP_tok\t";
				foreach my $FEATURE_tok(sort keys %{$hash1{$AC_tok}{$ENSG_tok}{$ENST_tmp[0]}{$ENSP_tok}{'FEATURE'}})
				{
					print OUTPUT "$FEATURE_tok\t";
				}
				foreach my $IPR_tok(sort keys %{$hash1{$AC_tok}{$ENSG_tok}{$ENST_tmp[0]}{$ENSP_tok}{'IPR'}})
				{
					print OUTPUT "$IPR_tok\t";
				}
				#foreach my $PFAM_tok(sort keys %{$hash1{$AC_tok}{$ENSG_tok}{$ENST_tmp[0]}{$ENSP_tok}{'PFAM'}})
				#{
				#	print OUTPUT "$PFAM_tok\t";
				#}
				foreach my $CCDS_tok(sort keys %{$hash1{$AC_tok}{$ENSG_tok}{$ENST_tmp[0]}{$ENSP_tok}{'CCDS'}})
				{
					print OUTPUT "$CCDS_tok\t";
				}
				print OUTPUT "\n";
			}
		}	
	}
}else{print "Impossible to open OUTPUT\n";die;}

my %hash_unique_post=();

$time='['. timestamp(). ']'."\n";
print "Tiempo de inicio de carga_hash2:$time\n";

if(open(OUTPUT, $output1))
{
	#	OUTPUT=UNIPROT_GLOBAL_POST_UNIQUE_2.txt
	#	A0A183  ENSG00000235942 ENST00000431011 ENSP00000411070 FEATURE:Late cornified envelope protein 6A.__0__240__CHAIN      
	#	A0AV96  ENSG00000163694 ENST00000295971 ENSP00000295971 FEATURE:RNA-binding protein 47.__0__1779__CHAIN IPR:IPR000504__RRM_dom  IPR:IPR006535__HnRNP_R/Q_splicing_fac   IPR:IPR012677__Nucleotide-bd_a/b_plait  PFAM:PF00076__RRM_1__3  
	#	A0AVF1  ENSG00000105948 ENST00000464848 ENSP00000419279 FEATURE:TPR 1.__168__102__REPEAT        FEATURE:TPR 2.__273__102__REPEAT        FEATURE:TPR 3.__450__102__REPEAT        FEATURE:TPR 4.__1401__102__REPEAT       FEATURE:Tetratricopeptide repeat protein 26.__0__1662__CHAIN    IPR:IPR011990__TPR-like_helical IPR:IPR013026__TPR-contain_dom  
	#	A0AVI4  ENSG00000168936 ENST00000382936 ENSP00000372394 FEATURE:Transmembrane protein 129.__84__1002__CHAIN     IPR:IPR018801__Tmpp129  PFAM:PF10272__Tmpp129__1        
	#	A0AVK6  ENSG00000129173 ENST00000250024 ENSP00000250024 FEATURE:Transcription factor E2F8.__0__2601__CHAIN      IPR:IPR003316__E2F_TDP  IPR:IPR011991__WHTH_DNA-bd_dom  IPR:IPR015633__E2F      PFAM:PF02319__E2F_TDP__2        

	while(my $line=<OUTPUT>)
	{
		chomp $line;
		my @output_tmp=split("\t", $line);
		my $AC=$output_tmp[0];
		my $ENSG=$output_tmp[1];
		my $ENST=$output_tmp[2];
		my $ENSP=$output_tmp[3];
		foreach my $tmp_tok(@output_tmp)
		{
			if($tmp_tok=~/^IPR:/){my $IPR=$tmp_tok; $hash_unique_post{$AC}{$ENSG}{$ENST}{$ENSP}{'IPR'}{$IPR}=1;}
			elsif($tmp_tok=~/^PFAM:/){my $PFAM=$tmp_tok; $hash_unique_post{$AC}{$ENSG}{$ENST}{$ENSP}{'PFAM'}{$PFAM}=1;}
			elsif($tmp_tok=~/^FEATURE:/){my $FEATURE=$tmp_tok; $hash_unique_post{$AC}{$ENSG}{$ENST}{$ENSP}{'FEATURE'}{$FEATURE}=1;}
			elsif($tmp_tok=~/^CCDS:/){my $CCDS=$tmp_tok; $hash_unique_post{$AC}{$ENSG}{$ENST}{$ENSP}{'CCDS'}{$CCDS}=1;}
		}	
	}
}else{print "Impossible to open OUTPUT second time\n";die;}

my %INTERPRO_coordinates_hash=();

$time='['. timestamp(). ']'."\n";
print "Tiempo de inicio de carga_hash3:$time\n";

if(open(INPUT2,$input2))
{
	# Input 2=~/Escritorio/Proyecto_clasificador/Documentos/INTERPRO/protein2ipr.dat
	
	# Q8NH21	IPR000276	G protein-coupled receptor, rhodopsin-like	PR00237	19	43
	# Q8NH21	IPR000276	G protein-coupled receptor, rhodopsin-like	PR00237	52	73
	# Q8NH21	IPR000276	G protein-coupled receptor, rhodopsin-like	PR00237	97	119
	# Q8NH21	IPR000276	G protein-coupled receptor, rhodopsin-like	PR00237	262	288
	# Q8NH21	IPR000276	G protein-coupled receptor, rhodopsin-like	PS00237	103	119
	# Q8NH21	IPR000725	Olfactory receptor	PR00245	85	96
	# Q8NH21	IPR000725	Olfactory receptor	PR00245	122	134
	# Q8NH21	IPR000725	Olfactory receptor	PR00245	169	185
	# Q8NH21	IPR000725	Olfactory receptor	PR00245	228	237
	# Q8NH21	IPR000725	Olfactory receptor	PR00245	273	284
	# Q8NH21	IPR017452	GPCR, rhodopsin-like, 7TM	G3DSA:1.20.1070.10	3	297
	# Q8NH21	IPR017452	GPCR, rhodopsin-like, 7TM	PS50262	34	280
	
	while(my $line=<INPUT2>)
	{
		chomp $line;
		#print "$line\n";
		$line=~/([^\t]+)\t/;
		my $AC_INTERPRO=$1;
		#print "$AC_INTERPRO\n";
		if(exists($hash_unique_post{$AC_INTERPRO}))
		{
			#print "HUMAN:$line\n";
			my @tmp_INTERPRO=split(/\t/,$line);
			#print "El array es:@tmp_INTERPRO\n";
			my $IPR_INTERPRO=$tmp_INTERPRO[1];
			my $coordinate_begin=$tmp_INTERPRO[4];
			my $coordinate_end=$tmp_INTERPRO[5];
			for (my$i=$coordinate_begin;$i<=$coordinate_end;$i++)
			{
				$INTERPRO_coordinates_hash{$AC_INTERPRO}{$IPR_INTERPRO}{$i}=1;
			}
		}else{next;}
	}# while	
}else{print "Unable to open $input2\n";}

$time='['. timestamp(). ']'."\n";
print "Tiempo de inicio de impresión 2:$time\n";

if (open(OUTPUT2, '>'.$output2))
{
	foreach my $AC_tok(sort keys %hash_unique_post)
	{
		foreach my $ENSG_tok(sort keys %{$hash_unique_post{$AC_tok}})
		{
			foreach my $ENST_tok(sort keys%{$hash_unique_post{$AC_tok}{$ENSG_tok}})
			{
				foreach my $ENSP_tok(sort keys %{$hash_unique_post{$AC_tok}{$ENSG_tok}{$ENST_tok}})
				{	
					print OUTPUT2 "$AC_tok\t$ENSG_tok\t$ENST_tok\t$ENSP_tok\t";
					
					foreach my $FEATURE_tok(sort keys %{$hash_unique_post{$AC_tok}{$ENSG_tok}{$ENST_tok}{$ENSP_tok}{'FEATURE'}})
					{
						print OUTPUT2 "$FEATURE_tok\t";
					}
					foreach my $IPR_tok(sort keys %{$hash_unique_post{$AC_tok}{$ENSG_tok}{$ENST_tok}{$ENSP_tok}{'IPR'}})
					{
						#print OUTPUT2 "$IPR_tok\n";
						if($IPR_tok=~/IPR:([^\_\_]+)/)
						{
							my $IPR_ID=$1;
							my $counter=1;
							my $START="NaN";
							#print OUTPUT2 "Hello_world2:$IPR_ID\n";
							my @POS_tmp=sort{ $a <=> $b }keys %{$INTERPRO_coordinates_hash{$AC_tok}{$IPR_ID}};
							for (my $i=0;$i<scalar(@POS_tmp);$i++)
							{
								if ($i==0)
									{
										$START=$POS_tmp[$i];
										my $distance_begin=($START-1)*3;
										print OUTPUT2 "$IPR_tok"."__"."$counter"."__"."$distance_begin"."__";
									}
									else
									{
										if($POS_tmp[$i]-$POS_tmp[$i-1] ==1 && $i<scalar(@POS_tmp)-1)
										{
										}
										elsif($POS_tmp[$i]-$POS_tmp[$i-1] !=1)
										{
											my $distance_feature=($POS_tmp[$i-1]-$START +1)*3;
											print OUTPUT2 "$distance_feature\t";
											$counter++;
											$START=$POS_tmp[$i];
											my $distance_begin=($START-1)*3;
											print OUTPUT2 "$IPR_tok"."__"."$counter"."__"."$distance_begin"."__";
										}
										elsif($i==scalar(@POS_tmp)-1)
										{
											my $distance_feature=($POS_tmp[$i]-$START +1)*3;
											print OUTPUT2 "$distance_feature\t"
										}
									}
										
										
								}		
									
						}	
					}
					#foreach my $PFAM_tok(sort keys %{$hash_unique_post{$AC_tok}{$ENSG_tok}{$ENST_tok}{$ENSP_tok}{'PFAM'}})
					#{
					#	print OUTPUT2 "$PFAM_tok\t";
					#}
					foreach my $CCDS_tok(sort keys %{$hash_unique_post{$AC_tok}{$ENSG_tok}{$ENST_tok}{$ENSP_tok}{'CCDS'}})
					{
						print OUTPUT2 "$CCDS_tok\t";
					}
					print OUTPUT2 "\n";
				}
			}	
		}
	}
}else{print "Impossible to open OUTPUT2\n";die;}



sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
