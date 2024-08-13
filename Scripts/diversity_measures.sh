
#Nucleotide diversiy pi
vcftools --vcf input_file --window-pi 10000 --out output_file

#Tajima's D
vcftools --vcf input_file --TajimaD 10000 --out output_file

#Fst
vcftools --vcf all_samples.vcf --weir-fst-pop population1 --weir-fst-pop population2 --fst-window-size 10000 --out pop1_vs_pop2_FST_10kb
