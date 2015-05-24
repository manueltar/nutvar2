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

    #   Setting the path variables for ENSEMBL API.
    
    Export the following variables or add them to your .profile
    
     export PERL5LIB=$HOME/src/ensembl/modules:$PERL5LIB
     export PERL5LIB=$HOME/src/ensembl-variation/modules:$PERL5LIB
     export PERL5LIB=$HOME/src/ensembl-compara/modules:$PERL5LIB
     export PERL5LIB=$HOME/src/ensembl-funcgen/modules:$PERL5LIB
     export PERL5LIB=$HOME/src/ensembl-tools/modules:$PERL5LIB
     export PERL5LIB=$HOME/src/bioperl-1.6.1/:$PERL5LIB

 
  # NUTVAR2 INSTALLATION
  
Download the .zip files with  git clone https://github.com/manueltar/NutVar2.git

unzip  nutvar2-master.zip

cd nutvar2-master

./install.sh

  During the installation of ENSEMBL VEP for release 75 the user will be asked in he/she wants to install a cache version of the genome or Fasta files. In both cases the answer is NO, as the GRCh37.75 genome is provided within NuTVar2.

# RUNNING

cd nutvar2-master

nutvar2-master>




