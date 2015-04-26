
# First install and decompress files

# This line decompresses a big file needed to work

echo Decompressing CDS_genomic_coordinates_full_compresed.txt.tar.gz
tar -xzvf data/build_tables/CDS_genomic_coordinates_full_compresed.txt.tar.gz -C data/build_tables/
rm data/build_tables/CDS_genomic_coordinates_full_compresed.txt.tar.gz

echo Decompressing snpEff v3_6
mkdir SOFTWARE/snpEff
unzip SOFTWARE/snpEff_v3_6_core.zip -d ./SOFTWARE/
rm SOFTWARE/snpEff_v3_6_core.zip

# The genome files are obtained from ELSEWHERE (#ISSUE 1.1)

echo obtaining the human genome version GRCh37.75
scp -r /home/manueltar/Desktop/Proyecto_NutVar2/snpEff_v3_6_GRCh37.75.zip SOFTWARE/snpEff
unzip SOFTWARE/snpEff/snpEff_v3_6_GRCh37.75.zip -d ./SOFTWARE/snpEff/
rm SOFTWARE/snpEff/snpEff_v3_6_GRCh37.75.zip



# VARIANT EFFECT PREDICTOR


# First, install the module DBI

# A zip file of DBI-1.633 is included in /SOFTWARE
echo THE PERL-DBI MODULE
tar -xvzf SOFTWARE/DBI-1.633.tar.gz -C SOFTWARE/
cd SOFTWARE/DBI-1.633/
perl Makefile.PL
make
make test
sudo make install
cd ../..
rm SOFTWARE/DBI-1.633.tar.gz

# Second, install the module DBD-mysql

# A zip file of DBI-1.633 is included in /SOFTWARE
#~ echo THE DBD-mysql MODULE
#~ tar -xvzf SOFTWARE/DBD-mysql-4.031.tar.gz -C SOFTWARE/
#~ cd SOFTWARE/DBD-mysql-4.031/
#~ perl Makefile.PL
#~ make
#~ make test
#~ sudo make install
#~ cd ../..
#~ rm SOFTWARE/DBD-mysql-4.031.tar.gz

# Decompress the VEP.zip file

unzip SOFTWARE/ensembl-tools-release-75.zip -d ./SOFTWARE/
echo INSTALLING VEP
perl SOFTWARE/ensembl-tools-release-75/scripts/variant_effect_predictor/INSTALL.pl

# printf n/n (#ISSUE 2)

# The genome files are obtained from ELSEWHERE (#ISSUE 1.2)

echo obtaining the human genome version GRCh37.75 for VEP

scp -r /home/manueltar/Desktop/Proyecto_NutVar2/homo_sapiens_vep_75.tar.gz SOFTWARE/

# mkdir -p SOFTWARE/{vep/{tmp/{homo_sapiens/{75,},},},}

mkdir SOFTWARE/vep/
mv SOFTWARE/vep/ SOFTWARE/.vep
tar -xvzf SOFTWARE/homo_sapiens_vep_75.tar.gz -C SOFTWARE/.vep/
cd SOFTWARE/.vep/
tar -xvf homo_sapiens_vep_75.tar

rm homo_sapiens_vep_75.tar
cd ../..

echo erasing homo_sapiens_vep_75.tar.gz
rm SOFTWARE/homo_sapiens_vep_75.tar.gz 


# Creating directories needed to run nutvar2

mkdir data/input/
mkdir data/intermediate/


