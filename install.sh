
# First install and decompress files

# This line decompresses a big file needed to work

echo Decompressing CDS_genomic_coordinates_full_compresed.txt.tar.gz
tar -xzvf data/build_tables/CDS_genomic_coordinates_full_compresed.txt.tar.gz -C data/build_tables/
rm data/build_tables/CDS_genomic_coordinates_full_compresed.txt.tar.gz

echo Decompressing snpEff v3_6
mkdir SOFTWARE/snpEff
unzip SOFTWARE/snpEff_v3_6_core.zip -d ./SOFTWARE/snpEff/
rm SOFTWARE/snpEff_v3_6_core.zip

# The genome files are obtained from ELSEWHERE (#ISSUE 1)

echo obtaining the human genome version GRCh37.75
scp -r /home/manueltar/Desktop/Proyecto_NutVar2/snpEff_v3_6_GRCh37.75.zip SOFTWARE/snpEff
unzip SOFTWARE/snpEff/snpEff_v3_6_GRCh37.75.zip -d ./SOFTWARE/snpEff/
rm SOFTWARE/snpEff/snpEff_v3_6_GRCh37.75.zip



# snpEff and the human genome (release GRCh37.75)


