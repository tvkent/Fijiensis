#!/bin/bash

scripts=/cap1/tyler.kent/Recombination/Fijiensis/Scripts
data=/cap1/tyler.kent/Recombination/Fijiensis/Data
results=/cap1/tyler.kent/Recombination/Fijiensis/Results
stdout=/cap1/tyler.kent/Recombination/Fijiensis
pops=(oldworld newworld)

for pop in '$pops[@]'
do

${scripts}/find_confs.bash

${scripts}/table_gen.bash

${scripts/pade.bash

${scripts}/rjmcmc.bash

${scripts}/post_to_text.bash

${scripts}/max_lk.bash

done
