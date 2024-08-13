#!/usr/bin/env python

#This script calculates the average pairwise nucleotide diversity (Dxy) between 2 populations. 
#Usage: python ./Dxy.py population1.vcf population2.vcf  genomic_region_start genomic_region_end.
#The VCF contains only the variants information of the genomic region of interest.

import numpy as np
import allel
import sys

def Dxy(pob1, pob2, start, end):
    callset1 = allel.read_vcf(pob1)
    callset2 = allel.read_vcf(pob2)
    genotype1 = callset1['calldata/GT']
    genotype2 = callset2['calldata/GT']
    gt1 = allel.GenotypeArray(genotype1)
    gt2 = allel.GenotypeArray(genotype2)
    ac1 = gt1.count_alleles()
    ac2 = gt2.count_alleles()
    dxy = allel.mean_pairwise_difference_between(ac1, ac2)
    dxy_sum = np.sum(dxy)
    region_length = int(end) - int(start)
    dxy_avg = dxy_sum / region_length

    return dxy_avg 

if len(sys.argv) != 5:
    print("This script calculates the average pairwise nucleotide diversity (Dxy) between 2 populations." + "\n\n" + "Usage: python ./Dxy.py population1.vcf population2.vcf  genomic_region_start genomic_region_end." + "\n\n" +"The VCF contains only the variants information of the genomic region of interest.Usage: python Dxy.py population1.vcf population2.vcf genomic_region_start genomic_region_end")

else:
    pob1 = sys.argv[1]
    pob2 = sys.argv[2]
    start = int(sys.argv[3])
    end = int(sys.argv[4])
    result = Dxy(pob1, pob2, start, end)
    print(result)

