#!/bin/bash

/cap1/tyler.kent/Software/LDhelmet_v1.7/ldhelmet rjmcmc --num_threads 20 -l $1 -p $2 -m /cap1/tyler.kent/Recombination/Fijiensis/Scripts/mut_mat.txt -s $3 -b 50 -w 50 --burn_in 10000 -n 100000 -o $4
