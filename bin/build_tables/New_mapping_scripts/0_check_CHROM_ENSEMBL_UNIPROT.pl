###############
###############

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;


my $input=$ARGV[0];
my $input2=$ARGV[1];
my $input3=$ARGV[2];
my $output1=$ARGV[3];
my %hash2=();

my $time='['. timestamp(). ']'."\n";
print "Tiempo de incio de carga_hash1:$time\n";

if(open(INPUT,$input))
{
	# Input=  gtf_output_ENSG.txt
	
	#~ ENSG00000186092 OR4F5   +       69091   70008
	#~ ENSG00000237683 AL627309.1      -       134901  139379
	#~ ENSG00000235249 OR4F29  +       367640  368634
	#~ ENSG00000185097 OR4F16  -       621059  622053
	
	while(my $line=<INPUT>)
	{
		chomp $line;
		#~ print "$line\n";
		my @tmp=split(/\t/,$line);
		#~ print "El array es:@tmp\n";
		my $ENSG=$tmp[0];
		#~ print "************$ENSG\n";
		my $SYMBOL=$tmp[1];
		
		$hash2{$ENSG}{$SYMBOL}=1;
	}	
}else{print "Unable to open $input\n";}


$time='['. timestamp(). ']'."\n";
print "Tiempo de incio de carga_hash2:$time\n";

if(open(INPUT2,$input2) && open(OUTPUT,'>'.$output1))
{
	# Input=  ENSEMBL_includes_UNIPROT.txt
	
	#~ A0A183  ENSG00000235942 ENST00000431011 UPSTREAM=0      DOWNSTREAM=0
	#~ A0AUZ9  ENSG00000144445 ENST00000281772 UPSTREAM=0      DOWNSTREAM=0
	#~ A0AV02  ENSG00000221955 ENST00000393469 UPSTREAM=0      DOWNSTREAM=0
	#~ A0AV96  ENSG00000163694 ENST00000295971 UPSTREAM=0      DOWNSTREAM=0
	#~ A0AVF1  ENSG00000105948 ENST00000464848 UPSTREAM=0      DOWNSTREAM=0
	#~ A0AVI4  ENSG00000168936 ENST00000382936 UPSTREAM=0      DOWNSTREAM=0
	#~ A0AVK6  ENSG00000129173 ENST00000250024 UPSTREAM=0      DOWNSTREAM=0
	#~ A0AVT1  ENSG00000033178 ENST00000322244 UPSTREAM=0      DOWNSTREAM=0
	#~ A0FGR9  ENSG00000158220 ENST00000389567 UPSTREAM=0      DOWNSTREAM=0

	
	while(my $line=<INPUT2>)
	{
		chomp $line;
		#~ print "$line\n";
		if($line =~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\tUPSTREAM=(.+)\tDOWNSTREAM=(.+)/)
		{
			#~ print "Hello_world:$line\n";
			my $AC=$1;
			my $ENSG=$2;
			my $ENST=$3;
			my $upstream=$4;
			my $downstream=$5;
			if(exists($hash2{$ENSG}))
			{
				print OUTPUT ">$ENSG\t$ENST\n";
			}
		}
		
	}	
}else{print "Unable to open $input2\n";}
$time='['. timestamp(). ']'."\n";
print "Tiempo de incio de carga_hash3:$time\n";

if(open(INPUT3,$input3)&& open(OUTPUT,'>>'.$output1))
{
	# Input=  aligned_no_display_offset.txt
	
	#~ O60763  ENSG00000138768 UPSTREAM=LIMIT:-72      DOWNSTREAM=0    NO_DISPLAY=0    OFFSET=0        
	#~ O60774  ENSG00000117507 UPSTREAM=0      DOWNSTREAM=LIMIT:418    NO_DISPLAY=0    OFFSET=0        
	#~ O75140  ENSG00000100150 UPSTREAM=0      DOWNSTREAM=0    NO_DISPLAY=725__733     OFFSET=734__-9  
	#~ O75949  ENSG00000130054 UPSTREAM=0      DOWNSTREAM=0    NO_DISPLAY=369__369     OFFSET=370__-1  
	#~ O95072  ENSG00000100918 UPSTREAM=0      DOWNSTREAM=0    NO_DISPLAY=236__236     OFFSET=237__-1  

	
	while(my $line=<INPUT3>)
	{
		chomp $line;
		#~ print "$line\n";
		if($line =~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\tUPSTREAM=(.+)\tDOWNSTREAM=(.+)\tNO_DISPLAY=(.+)\tOFFSET=(.+)/)
		{
			#~ print "Hello_world:$line\n";
			my $AC=$1;
			my $ENSG=$2;
			my $ENST=$3;
			my $upstream=$4;
			my $downstream=$5;
			my $no_display=$6;
			my $offset=$7;
			if(exists($hash2{$ENSG}) && $upstream ne 'XX')
			{
				print OUTPUT "$ENSG\t$ENST\n";
			}
		}
	}	
}else{print "Unable to open $input3\n";}



sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
