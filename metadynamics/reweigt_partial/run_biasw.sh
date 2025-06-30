#!/bin/bash

#clear
rm  weights.dat weights_norm.dat  eRMSDoutWeights.dat

#plumed
  rm bias.dat eRMSDout.dat
plumed driver --plumed plumed.dat --mf_xtc TRAJ.xtc --kt 2.494339 --pdb ../partial_compl/md0.pdb

#get max value of bias:
bmax=`awk 'BEGIN{max=0.}{if($1!="#!" && $2>max)max=$2}END{print max}' bias.dat`
echo $bmax
#calculate weights
awk '{if($1!="#!") print $1,(exp(($2-bmax)/kbt))}' kbt=2.494339 bmax=$bmax bias.dat > weights.dat
#get sum of weights to normalize them to 1
norm=`awk -F' ' '{sum+=$2;} END{print sum;}' weights.dat`
echo $norm
#calculate normalized weights (not needed to get K anyway)
awk '{print $1,$2/norm}' norm=$norm  weights.dat > weights_norm.dat
awk -F' ' '{sum+=$2;} END{print sum;}' weights_norm.dat
#add header to the file
sed -i '1s/^/#time weights\n/' weights_norm.dat
#megre the eRMSD and weight files
paste eRMSDout.dat weights_norm.dat >> eRMSDoutWeights.dat
#caclulate dG fold based on eRMSD<1 criteria 
# dG = kbTlnK where kbT is in kcal/mol and log in bash is natural logarithm
# check if in eRMSDoutWeights.dat the first line has weight 1 then remove that line
echo "energy difference between the two target states (kcal/mol):"
awk '{if($2<0.7 && $3>0.7) q1+=$5; else if($3<0.7 && $2>0.7) q2+=$5; else if($2<0.7 && $3 < 0.7 && $5) q3+=$5; else q4+=5} END{print q1, q2, q3, q4, -0.289*8.3144598*(log(q1/q2))/4.184}' eRMSDoutWeights.dat


