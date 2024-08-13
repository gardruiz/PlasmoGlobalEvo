 bcftools norm -a -f Genome/Pfalciparum.genome.fasta  -o cameroon_sera_splt.vcf annotated_Cameroon_SERA2.vcf 
awk '/^#/ || (length($4)==1 && length($5)==1) {print $0}' input.vcf > filtered_output.vcf

