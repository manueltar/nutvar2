echo `ulimit -a`

Softwaredir=$1
vcfinput=$2
output=$3
workingdir=$4

snpEFFdir=snpEff

bindir=bin
datadir=data

mypwd=$(pwd)


# NOTICE
# The GRCh37.75 database for snpEff and ENSEMBL human 37.75 for VEP are installed by default with the install.sh script
# java -jar snpEff.jar download -v GRCh37.75
# cd ${Softwaredir}/${snpEFFdir}

echo "Runing Snpeff"

## cat ${vcfinput} | perl -ne 'chomp;unless($_=~/^##/){$_=~s/^[Cc]hr//;print "$_\n";}' >${workingdir}/vcfinput.vcf
## java -Xmx4g -jar ${Softwaredir}/${snpEFFdir}/snpEff.jar eff -c ${Softwaredir}/${snpEFFdir}/snpEff.config -v GRCh37.71 -noLog ${workingdir}/vcfinput.vcf >${workingdir}/vcfinput_snpeffout.vcf

# java -Xmx4g -jar ${Softwaredir}/${snpEFFdir}/snpEff.jar eff -c ${Softwaredir}/${snpEFFdir}/snpEff.config -v GRCh37.71 -noLog ${vcfinput} >${workingdir}/vcfinput_snpeffout.vcf

# perl ${Softwaredir}/${bindir}/Snpeff_parser.pl ${workingdir}/vcfinput_snpeffout.vcf >${workingdir}/vcfinput_snpeffout_Snpeff_parsed.txt

### HERE TO BYPASS THE SNPEFF ASSESSEMENT AND PARSE AN EXTERNAL SNPEFF.out
perl ${Softwaredir}/${bindir}/Snpeff_parser.pl ${vcfinput} >${workingdir}/vcfinput_snpeffout_Snpeff_parsed.txt

echo "Snpeff done. Runing NUTVAR perl scripts"

# perl ${Softwaredir}/${bindir}/Mapping2CCDS_plusSeqFeatures.pl ${Softwaredir}/${datadir} ${vcfinput} >${workingdir}/mapped2CCDS_plusSeqFeatures.txt

perl ${Softwaredir}/${bindir}/Mapping2CCDS_plusSeqFeatures.pl ${Softwaredir}/${datadir} ${workingdir}/vcfinput_snpeffout_Snpeff_parsed.txt >${workingdir}/mapped2CCDS_plusSeqFeatures.txt

perl ${Softwaredir}/${bindir}/InterProParsing_Continuous.pl ${Softwaredir}/${datadir} ${workingdir}/mapped2CCDS_plusSeqFeatures.txt >${workingdir}/mapped2CCDS_plusSeqFeatures_InterPro.txt

less -S ${workingdir}/mapped2CCDS_plusSeqFeatures_InterPro.txt | head -1 | cut -f2-47 >${workingdir}/mapped2CCDS_plusSeqFeatures_InterPro_uniq.txt

less -S ${workingdir}/mapped2CCDS_plusSeqFeatures_InterPro.txt | egrep -v "^dummyIDforR" | cut -f2-47 | sort -u  >>${workingdir}/mapped2CCDS_plusSeqFeatures_InterPro_uniq.txt

perl ${Softwaredir}/${bindir}/NMDannotation.pl ${Softwaredir}/${datadir} ${workingdir}/mapped2CCDS_plusSeqFeatures_InterPro_uniq.txt >${workingdir}/mapped2CCDS_plusSeqFeatures_InterPro_uniq_NMD.txt

perl ${Softwaredir}/${bindir}/LongestCCDSannotation.pl ${Softwaredir}/${datadir} ${workingdir}/mapped2CCDS_plusSeqFeatures_InterPro_uniq_NMD.txt >${workingdir}/mapped2CCDS_plusSeqFeatures_InterPro_uniq_NMD_LongestCCDS.txt

perl ${Softwaredir}/${bindir}/AddingGeuvadisData.pl ${Softwaredir}/${datadir} ${workingdir} ${workingdir}/mapped2CCDS_plusSeqFeatures_InterPro_uniq_NMD_LongestCCDS.txt >${workingdir}/mapped2CCDS_plusSeqFeatures_InterPro_uniq_NMD_LongestCCDS_plusGeuvadis.txt

echo "Before"

perl ${Softwaredir}/${bindir}/Collapsing_per_gene.pl ${Softwaredir}/${datadir} ${workingdir}/mapped2CCDS_plusSeqFeatures_InterPro_uniq_NMD_LongestCCDS_plusGeuvadis.txt >${workingdir}/OUTPUTVCF_CollapsedPerGene.txt

echo "After"

less ${workingdir}/OUTPUTVCF_CollapsedPerGene.txt | perl -ne 'chomp;@tr=split/\t/;if($_=~/^baseNCBI37/){print "$_\n";};if(($tr[5]>0)){print "$_\n";}' >${workingdir}/OUTPUTVCF_CollapsedPerGene_plusGeuvadis_STOPS.txt
less ${workingdir}/OUTPUTVCF_CollapsedPerGene.txt | perl -ne 'chomp;@tr=split/\t/;if($_=~/^baseNCBI37/){print "$_\n";};if(($tr[6]>0)){print "$_\n";}' >${workingdir}/OUTPUTVCF_CollapsedPerGene_plusGeuvadis_FRAMESHIFTS.txt

perl ${Softwaredir}/${bindir}/BayesianProbCalculator_plusRANKPercentileESP1000G.pl ${Softwaredir}/${datadir} ${workingdir}/OUTPUTVCF_CollapsedPerGene_plusGeuvadis_STOPS.txt stop-gains NewPpal.txt_Stop-gains_Scores_.txt  >${workingdir}/OUTPUTVCF_CollapsedPerGene_plusGeuvadis_STOPS_scores.txt
perl ${Softwaredir}/${bindir}/BayesianProbCalculator_plusRANKPercentileESP1000G.pl ${Softwaredir}/${datadir} ${workingdir}/OUTPUTVCF_CollapsedPerGene_plusGeuvadis_FRAMESHIFTS.txt frameshifts NewPpal.txt_Frame-Shifts_Scores_.txt  >${workingdir}/OUTPUTVCF_CollapsedPerGene_plusGeuvadis_FRAMESHIFTS_scores.txt

mkdir ${workingdir}/${5}_NUTvar_results
mv ${workingdir}/OUTPUTVCF_CollapsedPerGene_plusGeuvadis_STOPS_scores.txt ${workingdir}/${5}_NUTvar_results/${5}_NUTvar_Stop-gains_output.txt
mv ${workingdir}/OUTPUTVCF_CollapsedPerGene_plusGeuvadis_FRAMESHIFTS_scores.txt ${workingdir}/${5}_NUTvar_results/${5}_NUTvar_Frameshifts_output.txt

cd ${workingdir}/
tar -zcvf ${5}_NUTvar_results.tar.gz ${5}_NUTvar_results

## rm -r ${5}_NUTvar_results
## rm mapped2CCDS_plusSeqFeatures_InterPro_uniq_NMD_LongestCCDS_plusGeuvadis.txt
## rm mapped2CCDS_plusSeqFeatures_InterPro_uniq_NMD_LongestCCDS.txt
## rm mapped2CCDS_plusSeqFeatures_InterPro_uniq_NMD.txt
## rm mapped2CCDS_plusSeqFeatures_InterPro_uniq.txt
## rm mapped2CCDS_plusSeqFeatures_InterPro.txt
## rm mapped2CCDS_plusSeqFeatures.txt
## rm Merging_GEUVADIS_plus_VariationDatasets_ErrorFile.txt
## rm OUTPUTVCF_CollapsedPerGene_plusGeuvadis_FRAMESHIFTS.txt
## rm OUTPUTVCF_CollapsedPerGene_plusGeuvadis_STOPS.txt
## rm OUTPUTVCF_CollapsedPerGene.txt
## rm vcfinput_snpeffout_Snpeff_parsed.txt
## rm vcfinput_snpeffout.vcf
## rm vcfinput.vcf

cd ${mypwd}

mv ${workingdir}/${5}_NUTvar_results.tar.gz $3

echo "Done!"
