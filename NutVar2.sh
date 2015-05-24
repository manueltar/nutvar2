# Define time and stamp it in every result directory

echo `ulimit -a`

Softwaredir=$1
# ~/Downloads/nutvar2-master
vcfinput=$2
# USER INPUT p.e. vep_example.vcf
output=$3
# data/final

snpEFFdir=SOFTWARE/snpEff
VEPdir=SOFTWARE/ensembl-tools-release-75/scripts/variant_effect_predictor

bindir1=bin/shared
bindir2=bin/snpEff
bindir3=bin/VEP
datadir1=data/intermediate
datadir2=data/build_tables
datadir3=data/external
datadir4=data/final

mypwd=$(pwd)

	# Open file and eliminate header lines

cat $2 | perl -ne 'chomp;unless($_=~/^##/){$_=~s/^[Cc]hr//;print "$_\n";}' > ${datadir1}/vcfinput.vcf

	# Run the minimal representation script

	perl $1/${bindir1}/2_Script_minimal_representation_vcf_7.0.pl ${datadir1}/vcfinput.vcf $1/${datadir1}/vcfinput_mr.vcf

	# Ask the user whether she wants to run SnpEff, VEP or both #ISSUE

	# NOTICE GRCh37.75 database for snpEff and the ENSEMBL version of genome 37.75 are installed through the install script

	# Run SnpEff

echo "Runing Snpeff"

java -Xmx4g -jar $1/${snpEFFdir}/snpEff.jar eff -c $1/${snpEFFdir}/snpEff.config -v GRCh37.75 -lof -csvStats -nextProt -sequenceOntology $1/${datadir1}/vcfinput_mr.vcf > $1/${datadir1}/vcfinput_mr_eff.vcf

mv snpEff_genes.txt $1/${datadir1}/snpEff_genes.txt
mv snpEff_summary.csv $1/${datadir1}/snpEff_summary.csv

#Parsing the results of snpEff

# ISSUE There are two more scripts in this folder /bin/snpEff/ 24 and 25 ---> Erase them?

echo "Parsing Snpeff results"

perl $1/${bindir2}/24_snpEff_parser_def_minus_heather_2.0.pl ${datadir1}/vcfinput_mr_eff.vcf ${datadir1}/out_snpeff_parsed.txt

# Run VEP

# VEP is very slow the user must know that if he wants to use VEP his files must be chopped

echo "Runing VEP"

perl $1/${VEPdir}/variant_effect_predictor.pl -i ${datadir1}/vcfinput_mr.vcf --offline --output_file ${datadir1}/vcfinput_mr_vep.vcf --everything --vcf --cache --dir $1/SOFTWARE/.vep/

#Parsing the results of VEP

# ISSUE There are two more scripts in this folder /bin/VEP/ versions of 24 ---> Erase them?

echo "Parsing VEP results"

perl $1/${bindir3}/24_VEP_parser_def_minus_heather_variante_def.pl ${datadir1}/vcfinput_mr_vep.vcf ${datadir1}/out_vep_parsed.txt 

echo "snpEff and VEP done. Runing NUTVAR perl scripts for snpEff results"

# Here there is a chance to MPI as script 25 is the longest the others can take place untill 27 while 25 is executing

perl $1/${bindir1}/25_Downstream_frameshift_API_independent_5.0.pl ${datadir1}/out_snpeff_parsed.txt ${datadir2}/CDS_genomic_coordinates_full_compresed.txt ${datadir1}/snpeff_derived_PTCS_API_independent.txt
#~ 
perl $1/${bindir1}/26_NEW_EXTRA_key_%_sequence_2.0.pl ${datadir1}/out_snpeff_parsed.txt ${datadir2}/ENST_table_full_condensed.txt ${datadir1}/snpeff_percentage.txt
#~ 
perl $1/${bindir1}/27_key_NMD_5.0_3.0.pl ${datadir1}/out_snpeff_parsed.txt ${datadir2}/NMD_table.txt ${datadir2}/ENST_table_full_condensed.txt ${datadir1}/snpeff_NMD.txt
#~ 
perl $1/${bindir1}/32_key_PROTEINS_8.0_GLOBAL_3.0.pl ${datadir1}/out_snpeff_parsed.txt ${datadir2}/ENST_table_full_condensed.txt ${datadir2}/ALL_ISOFORMS_PROTEIN_table_full.txt ${datadir1}/snpeff_detailed_ProtAndSite_Pre_step.txt

sort -k1,1 -k2,2 -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 ${datadir1}/snpeff_detailed_ProtAndSite_Pre_step.txt > ${datadir1}/snpeff_detailed_ProtAndSite_Pre_step_ordered.txt

perl $1/${bindir1}/32_key_PROTEINS_8.0_GLOBAL_ParteII.pl ${datadir1}/snpeff_detailed_ProtAndSite_Pre_step_ordered.txt ${datadir1}/snpeff_detailed_ProtAndSite_Post_step.txt

perl $1/${bindir1}/32_key_PROTEINS_8.0_GLOBAL_3.0.pl ${datadir1}/out_snpeff_parsed.txt ${datadir2}/ENST_table_full_condensed.txt ${datadir2}/ALL_ISOFORMS_DOMAIN_table_full.txt ${datadir1}/snpeff_DOMAINS_Pre_step.txt

sort -k1,1 -k2,2 -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 ${datadir1}/snpeff_DOMAINS_Pre_step.txt > ${datadir1}/snpeff_DOMAINS_Pre_step_ordered.txt

perl $1/${bindir1}/32_key_PROTEINS_8.0_GLOBAL_ParteII.pl ${datadir1}/snpeff_DOMAINS_Pre_step_ordered.txt ${datadir1}/snpeff_DOMAINS_Post_step.txt
#~ 
perl $1/${bindir1}/27_key_NMD_5.0_DERIVED_STOPS_2.0.pl ${datadir1}/snpeff_derived_PTCS_API_independent.txt ${datadir1}/out_snpeff_parsed.txt ${datadir2}/NMD_table.txt ${datadir2}/ENST_table_full_condensed.txt ${datadir1}/snpeff_derived_NMD.txt
#~ 
perl $1/${bindir1}/38_global_feature_table_1_4_paralell.pl ${datadir1}/snpeff_NMD.txt ${datadir1}/snpeff_derived_NMD.txt ${datadir1}/snpeff_detailed_ProtAndSite_Post_step.txt ${datadir1}/snpeff_DOMAINS_Post_step.txt ${datadir1}/snpeff_percentage.txt ${datadir1}/snpeff_first_table.txt

perl $1/${bindir1}/40_tabla_PEJMAN_16_def_2.0.pl ${datadir1}/snpeff_first_table.txt ${datadir1}/snpeff_NMD.txt ${datadir1}/snpeff_derived_NMD.txt ${datadir2}/gtf_output_ENST.txt ${datadir2}/gtf_output_ENSG.txt ${datadir2}/ENST_table_full_condensed.txt ${datadir3}/appris_principal_isoform_gencode_19_15_10_2014.txt ${datadir3}/Pervasive.txt ${datadir1}/Matrix_snpeff.txt
#~ 
perl $1/${bindir1}/41_CCDS_collapser_3.0.pl ${datadir2}/gtf_tabladef_sorted_by_SYMBOL.txt ${datadir1}/snpeff_NMD.txt ${datadir1}/snpeff_NMD_CCDS.txt ${datadir1}/snpeff_derived_NMD.txt ${datadir1}/snpeff_derived_NMD_CCDS.txt ${datadir2}/gtf_output_ENSG.txt ${datadir2}/gtf_output_ENSG_CCDS.txt ${datadir2}/gtf_output_ENST.txt ${datadir2}/gtf_output_ENST_CCDS.txt ${datadir2}/ENST_table_full_condensed.txt ${datadir2}/ENST_table_full_condensed_CCDS.txt ${datadir3}/appris_principal_isoform_gencode_19_15_10_2014.txt ${datadir3}/appris_principal_isoform_gencode_19_15_10_2014_CCDS.txt ${datadir3}/Pervasive.txt ${datadir3}/Pervasive_CCDS.txt ${datadir1}/snpeff_first_table.txt ${datadir1}/snpeff_first_table_CCDS.txt

perl $1/${bindir1}/42_tabla_PEJMAN_15.0_version_paralel_4.0.pl ${datadir1}/snpeff_first_table_CCDS.txt ${datadir1}/snpeff_NMD_CCDS.txt ${datadir1}/snpeff_derived_NMD_CCDS.txt ${datadir2}/gtf_output_ENST_CCDS.txt ${datadir2}/gtf_output_ENSG_CCDS.txt ${datadir2}/ENST_table_full_condensed_CCDS.txt ${datadir3}/appris_principal_isoform_gencode_19_15_10_2014_CCDS.txt ${datadir3}/Pervasive_CCDS.txt ${datadir1}/Matrix_snpeff_CCDS.txt

mkdir ${datadir4}/

perl $1/${bindir1}/53BIS_Fuse_Matrix\&Gene_based.pl ${datadir1}/Matrix_snpeff.txt ${datadir3}/pRDG2.txt ${datadir3}/Genes_AllInnateImmunity.txt ${datadir3}/Genes_Antiviral.txt ${datadir3}/Genes_ISGs.txt ${datadir3}/Genes_OMIMrecessive.txt ${datadir3}/RVIS2.txt ${datadir4}/Matrix_snpeff_added_gene_based_scores.txt

perl $1/${bindir1}/53BIS_Fuse_Matrix\&Gene_based.pl ${datadir1}/Matrix_snpeff_CCDS.txt ${datadir3}/pRDG2.txt ${datadir3}/Genes_AllInnateImmunity.txt ${datadir3}/Genes_Antiviral.txt ${datadir3}/Genes_ISGs.txt ${datadir3}/Genes_OMIMrecessive.txt ${datadir3}/RVIS2.txt ${datadir4}/Matrix_snpeff_CCDS_added_gene_based_scores.txt

echo "snpEff data matrix generated"


