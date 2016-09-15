#!/bin/bash

/cap1/tyler.kent/Recombination/Fijiensis/Scripts/vcftools/bin/vcftools --vcf /cap1/tyler.kent/Recombination/Fijiensis/Data/All_indiv_n.vcf --minDP 10 --maxDP 50 --chr scaffold_1 --out /cap1/tyler.kent/Recombination/Fijiensis/Data/all_n_quick --ldhelmet 
