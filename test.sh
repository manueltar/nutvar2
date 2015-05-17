# Define time and stamp it in every result directory

echo `ulimit -a`

Softwaredir=$1
vcfinput=$2
output=$3
workingdir=$4

snpEFFdir=SOFTWARE/snpEff
VEPdir=ensembl-tools-release-75/scripts/variant_effect_predictor/

bindir1=bin/shared
bindir2=bin/snpEff
bindir3=bin/VEP
datadir1=data/intermediate
datadir2=data/build_tables
datadir3=data/external
datadir4=data/final

mypwd=$(pwd)

	# Open file and eliminate header lines

#~ cat $2 | perl -ne 'chomp;unless($_=~/^##/){$_=~s/^[Cc]hr//;print "$_\n";}' > $4/vcfinput.vcf 

	# Run the minimal representation script

	# perl $1/${bindir1}/2_Script_minimal_representation_vcf_7.0.pl $4/vcfinput.vcf ${datadir1}/vcfinput_mr.vcf
#~ perl $1/${bindir1}/2_Script_minimal_representation_vcf_7.0.pl $4/vcfinput.vcf $1/${datadir1}/vcfinput_mr.vcf

	# Ask the user whether she wants to run SnpEff, VEP or both

	# NOTICE GRCh37.75 database for snpEff and the ENSEMBL version of genome 37.75 are installed through the install script

	# Run SnpEff

echo "Runing Snpeff"

java -Xmx4g -jar $1/${snpEFFdir}/snpEff.jar eff -c $1/${snpEFFdir}/snpEff.config -v GRCh37.75 -lof -csvStats -nextProt -sequenceOntology $1/${datadir1}/vcfinput_mr.vcf > $1/${datadir1}/vcfinput_mr_eff.vcf

mv snpEff_genes.txt $1/${datadir1}/snpEff_genes.txt
mv snpEff_summary.csv $1/${datadir1}/snpEff_summary.csv
