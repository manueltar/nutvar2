#!/bin/bash


#BSUB -L /bin/bash
#BSUB -o snpeff_01-%J-output.txt		# Standard output [-%J-] includes jobname [-%I-] includes iteartion
#BSUB -e snpeff_01-%J-error.txt		# Standard error [-%J-] includes jobname [-%I-] includes iteration
#BSUB -J smpeff_01				# Name your job
#BSUB -u manueltar@hotmail.com	# Set email address
#BSUB -N				# Send an email when finished
#BSUB -n 4				# Select number of cores 
#BSUB –R "span[ptile=4]"		# All cores on the same host
#BSUB -R "rusage[mem=8192]" # 2 GB RAM
#BSUB -M 8388608

module add Development/java_jre/1.8.0_05;

java -Xmx4g -jar /home/mtardagu/Programs/snpEff/snpEff.jar eff -c /home/mtardagu/Programs/snpEff/snpEff.config -v GRCh37.75 -lof -csvStats -nextProt -sequenceOntology  ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf > OUT_1000gk.eff.vcf

