#!/bin/bash


#BSUB -L /bin/bash
#BSUB -o 2_vep_01-%J-output.txt		# Standard output [-%J-] includes jobname [-%I-] includes iteartion
#BSUB -e 2_vep_01-%J-error.txt		# Standard error [-%J-] includes jobname [-%I-] includes iteration
#BSUB -J 2_vep_01				# Name your job
#BSUB -u manueltar@hotmail.com	# Set email address
#BSUB -N				# Send an email when finished
#BSUB -n 4				# Select number of cores 
#BSUB –R "span[ptile=4]"		# All cores on the same host
#BSUB -R "rusage[mem=8192]" # 8 GB RAM
#BSUB -M 8388608

module add Development/Ensembl_API/75;

perl /scratch/cluster/monthly/mtardagu/Programs/ensembl-tools-release-75/scripts/variant_effect_predictor/variant_effect_predictor.pl -i 2GKaa --offline --output_file variant_effec_output_second_round_01.txt --everything --vcf --cache --dir /scratch/cluster/monthly/mtardagu/.vep/tmp/
