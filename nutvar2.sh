
# This line creates the subdirectory data/intermediate

mkdir data/intermediate/

# This line runs snpEff with example file

java -Xmx4g -jar SOFTWARE/snpEff/snpEff.jar eff -c SOFTWARE/snpEff/snpEff.config -v GRCh37.75 -lof -csvStats -nextProt -sequenceOntology test/example.vcf > data/intermediate/example_eff.vcf

mv snpEff_genes.txt data/intermediate/
mv snpEff_summary.csv data/intermediate/
