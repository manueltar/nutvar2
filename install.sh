
# First install and decompress files

# This line decompresses a big file needed to work

echo Decompressing CDS_genomic_coordinates_full_compresed.txt.tar.gz
tar -xzvf data/build_tables/CDS_genomic_coordinates_full_compresed.txt.tar.gz -C data/build_tables/
rm data/build_tables/CDS_genomic_coordinates_full_compresed.txt.tar.gz

echo Decompressing snpEff v3_6
mkdir SOFTWARE/snpEff
unzip SOFTWARE/snpEff_v3_6_core.zip -d ./SOFTWARE/
rm SOFTWARE/snpEff_v3_6_core.zip

# The genome files are obtained from ELSEWHERE (#ISSUE 1.1) ALL OF THESE FILES ARE OBTAINED FROM LOCAL THIS HAS TO CHANGE

echo obtaining the human genome version GRCh37.75
scp -r /media/manueltar/Data/Dropbox/Dropbox/0_p_NutVar2/snpEff_v3_6_GRCh37.75.zip SOFTWARE/snpEff
unzip SOFTWARE/snpEff/snpEff_v3_6_GRCh37.75.zip -d ./SOFTWARE/snpEff/
rm SOFTWARE/snpEff/snpEff_v3_6_GRCh37.75.zip


# VARIANT EFFECT PREDICTOR


# First, install the module DBI

# A zip file of DBI-1.633 is included in /SOFTWARE
#~ echo THE PERL-DBI MODULE
#~ tar -xvzf SOFTWARE/DBI-1.633.tar.gz -C SOFTWARE/
#~ cd SOFTWARE/DBI-1.633/
#~ perl Makefile.PL
#~ make
#~ make test
#~ sudo make install
#~ cd ../..
#~ rm SOFTWARE/DBI-1.633.tar.gz

# Second, install and configure mysql in case you don't have it

#~ sudo apt-get install mysql-server
#~ 
#~ sudo apt-get install libmysqlclient-dev

#~ sudo service mysql restart

# check if mysql is running

#~ sudo netstat -tap | grep mysql

# Set password

#~ sudo dpkg-reconfigure mysql-server-5.5

# To test the pass

# mysql -u root -p
# exit

# Third, install the module DBD-mysql


#~ echo THE DBD-mysql MODULE
#~ tar -xvzf SOFTWARE/DBD-mysql-4.031.tar.gz -C SOFTWARE/
#~ cd SOFTWARE/DBD-mysql-4.031/
#~ perl Makefile.PL
#~ make
#~ make test
#~ sudo make install
#~ cd ../..
#~ rm SOFTWARE/DBD-mysql-4.031.tar.gz

# Perl API


# Check first 

#~ $ perl -v
#~ 
#~ This is perl 5, version 18, subversion 2 (v5.18.2) built for x86_64-linux-gnu-thread-multi
#~ (with 41 registered patches, see perl -V for more detail)
#~ 
#~ $ perl -MDBI -e 'warn $DBI::VERSION'
#~ 1.633 at -e line 1.
#~ $ perl -MDBD::mysql -e 'warn $DBD::mysql::VERSION'
#~ 4.025 at -e line 1.

#scp -r

# The release of BioPerl that works with VEP75 is  BioPerl-1.6.1. not 1.2.3 as explained in the video of the Ensembl API installation Tutorial

#~ mkdir ~/src
#~ tar zxf SOFTWARE/BioPerl-1.6.1.tar.gz -C ~/src/
#~ tar zxf SOFTWARE/ensembl-api.tar.gz -C ~/src/

# Add the following lines to .profile and refresh .profile

#~ $ nano ~/.profile
#~ export PERL5LIB=$HOME/src/ensembl/modules:$PERL5LIB
#~ export PERL5LIB=$HOME/src/ensembl-variation/modules:$PERL5LIB
#~ export PERL5LIB=$HOME/src/ensembl-compara/modules:$PERL5LIB
#~ export PERL5LIB=$HOME/src/ensembl-funcgen/modules:$PERL5LIB
#~ export PERL5LIB=$HOME/src/ensembl-tools/modules:$PERL5LIB
#~ export PERL5LIB=$HOME/src/bioperl-1.6.1/:$PERL5LIB



#~ $ . ~/.profile

# Test

#~ perl ~/src/ensembl/misc-scripts/ping_ensembl.pl 
#~ 
#~ Installation is good. Connection to Ensembl works and you can query the human core database


# Decompress the VEP.zip file

unzip SOFTWARE/ensembl-tools-release-75.zip -d ./SOFTWARE/
echo INSTALLING VEP
perl SOFTWARE/ensembl-tools-release-75/scripts/variant_effect_predictor/INSTALL.pl

# printf n/n/n (#ISSUE 2)

# The genome files are obtained from ELSEWHERE (#ISSUE 1.2)

echo obtaining the human genome version GRCh37.75 (VEP)

scp -r /media/manueltar/Data/Dropbox/Dropbox/0_p_NutVar2/homo_sapiens_vep_75.tar.gz SOFTWARE/

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

