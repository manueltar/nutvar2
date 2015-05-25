# nutvar2
Classifier of the potential pathogenicity of human genomic truncations

# INSTALLATION

  # Dependencies
  
    #  Perl DBI module
        
        If you have not already installed this module a zip file of DBI-1.633 is included in /SOFTWARE module. to            install it follow the ensuing instructions:

     tar -xvzf SOFTWARE/DBI-1.633.tar.gz -C SOFTWARE/
     cd SOFTWARE/DBI-1.633/
     perl Makefile.PL
     make
     make test
     sudo make install
     cd ../..
     
      # To check the installation
 
      perl -MDBI -e 'warn $DBI::VERSION'
    
    #   MySQL
    
        If you have not installed and configured MySQL in your computer follow the ensuing instructions
      
    sudo apt-get install mysql-server
    
    sudo apt-get install libmysqlclient-dev
    
    sudo service mysql restart

    # check if mysql is running
    
     sudo netstat -tap | grep mysql

    # Set password
    
     sudo dpkg-reconfigure mysql-server-5.5

    # To test the pass
    
    mysql -u root -p
    mysql>exit
    
    #    DBD-mysql module

         If you have not already installed this module a zip file of DBI-1.633 is included in /SOFTWARE module. to          install it follow the ensuing instructions:

     tar -xvzf SOFTWARE/DBD-mysql-4.031.tar.gz -C SOFTWARE/
     cd SOFTWARE/DBD-mysql-4.031/
     perl Makefile.PL
     make
     make test
     sudo make install
     cd ../..
     
     # To check the installation
     
      perl -MDBD::mysql -e 'warn $DBD::mysql::VERSION'
      
    #   ggplot package in R
    
      RStudio> install.packages("ggplot2")
    
    # ROCR package in R
    
      RStudio> install.packages("ROCR")

 
# NUTVAR2 INSTALLATION
  
Download the .zip files with  git clone https://github.com/manueltar/NutVar2.git

unzip  nutvar2-master.zip

cd nutvar2-master

  #   Setting the path variables for ENSEMBL API.
    
    Export the following variables or add them to your .profile
    
     export PERL5LIB=$HOME/src/ensembl/modules:$PERL5LIB
     export PERL5LIB=$HOME/src/ensembl-variation/modules:$PERL5LIB
     export PERL5LIB=$HOME/src/ensembl-compara/modules:$PERL5LIB
     export PERL5LIB=$HOME/src/ensembl-funcgen/modules:$PERL5LIB
     export PERL5LIB=$HOME/src/ensembl-tools/modules:$PERL5LIB
     export PERL5LIB=$HOME/src/bioperl-1.6.1/:$PERL5LIB

  # Install NutVar2
  
  ./install.sh

  During the installation of ENSEMBL VEP for release 75 the user will be asked in he/she wants to install a cache version of the genome or Fasta files. In both cases the answer is NO, as the GRCh37.75 genome is provided within NuTVar2.
  

# RUNNING


# There are three different options to run NutVar2: using SnpEff as the only predictor for variant outcome, using VEP as the only predictor for variant outcome, and using both SnpEff and VEP. 

**Note: VEP is a very comprehensive and informative predictor of variant outcome, but its time of execution is way larger than that of SnpEff. If VEP is going to be used in large input vcf files, cut files in smaller subfiles prior to running NutVar2.


cd nutvar2-master

# 1- SnpEff

nutvar2-master$ ./NutVar2_snpEff.sh ~/Downloads/nutvar2-master user.vcf data/final

Test: nutvar2-master$ ./NutVar2_snpEff.sh ~/Downloads/nutvar2-master example.vcf data/final

# 2- VEP

nutvar2-master$ ./NutVar2_VEP.sh ~/Downloads/nutvar2-master user.vcf data/final

Test:nutvar2-master$ ./NutVar2_VEP.sh ~/Downloads/nutvar2-master example.vcf data/final

# 3- VEP and SnpEff

nutvar2-maste$ ./NutVar2_snpEff_and_VEP.sh ~/Downloads/nutvar2-master example.vcf data/final

ISSUE!! Create a relative path to setwd in R

# Issues

  # Installation ISSUE!! The genomes are huge where are we going to allocate them online for the user to download?. Right now I retrieve them from my local disk.
  
  # Running Issue 1 !!-> combine VEP and SnpEf results somehow in the third option.
  
  # Running Issue 2 !!-> Allow some MPI option specially for VEP and also intra bash script, run script 25 while running the rest of scripts.
  
  # Training issue: Cannot introduce NaN in the code. Right now I am imputing NaN values to different proxys. See script Processing_data.R.
  
  # Training issue: 4 non Pathogenic splice variants of the training set have NaN as Percentage value, why?
  
  # Training issue: Calssifysplice variants according to the length of the intron in which they lie in and the type of splice site: canonical versus non-canonical.
  
  # Output issue: Set an R top function to analyse the matrix of results of the user input.
  
  # Output issue: Integrate TAN-Bayes in this top function.




