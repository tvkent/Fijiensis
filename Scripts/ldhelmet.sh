#!/bin/bash

scripts=/cap1/tyler.kent/Recombination/Fijiensis/Scripts
data=/cap1/tyler.kent/Recombination/Fijiensis/Data
results=/cap1/tyler.kent/Recombination/Fijiensis/Results
stdout=/cap1/tyler.kent/Recombination/Fijiensis
pops=(oldworld newworld)

for pop in "${pops[@]}"
do

seq=/cap1/tyler.kent/Recombination/Fijiensis/Data/${pop}_indiv_n.filtered.scaffold_1.ldhelmet.snps
echo "Beginning find_confs" >> ${stdout}/${pop}.out
echo date >> ${stdout}/${pop}.out
bash ${scripts}/find_confs.bash ${results}/${pop}.conf ${seq} >> ${stdout}/${pop}.out 2>> ${stdout}/${pop}.err
echo "find_confs finished" >> ${stdout}/${pop}.out
echo date >> ${stdout}/${pop}.out

echo "Beginning table_gen" >> ${stdout}/${pop}.out
echo date >> ${stdout}/${pop}.out
bash ${scripts}/table_gen.bash ${results}/${pop}.conf ${results}/${pop}.lk >> ${stdout}/${pop}.out 2>> ${stdout}/${pop}.err
echo "table_gen finished" >> ${stdout}/${pop}.out
echo date >> ${stdout}/${pop}.out

echo "Beginning pade" >> ${stdout}/${pop}.out
echo date >> ${stdout}/${pop}.out
bash ${scripts}/pade.bash ${results}/${pop}.conf ${results}/${pop}.pade >> ${stdout}/${pop}.out 2>> ${stdout}/${pop}.err
echo "pade finished" >> ${stdout}/${pop}.out
echo date >> ${stdout}/${pop}.out

echo "Beginning rjmcmc" >> ${stdout}/${pop}.out
echo date >> ${stdout}/${pop}.out
bash ${scripts}/rjmcmc.bash ${results}/${pop}.lk ${results}/${pop}.pade ${seq} ${results}/${pop}.post >> ${stdout}/${pop}.out 2>> ${stdout}/${pop}.err
echo "rjmcmc finished" >> ${stdout}/${pop}.out
echo date >> ${stdout}/${pop}.out

echo "Beginning post_to_text" >> ${stdout}/${pop}.out
echo date >> ${stdout}/${pop}.out
bash ${scripts}/post_to_text.bash ${results}/${pop}.out.txt ${results}/${pop}.post >> ${stdout}/${pop}.out 2>> ${stdout}/${pop}.err
echo "post_to_text finished" >> ${stdout}/${pop}.out
echo date >> ${stdout}/${pop}.out

#echo "Beginning max_lk" >> ${stdout}/${pop}.out
#echo date >> ${stdout}/${pop}.out
#bash ${scripts}/max_lk.bash ${results}/${pop}.lk ${pop}.post ${seq} >> ${stdout}/${pop}.out 2>> ${stdout}/${pop}.err
#echo "max_lk finished" >> ${stdout}/${pop}.out
#echo date >> ${stdout}/${pop}.out

done
