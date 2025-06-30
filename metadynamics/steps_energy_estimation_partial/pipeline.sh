#!/bin/bash

################################################################################
# first do to the reweigt_partial dir and caclulate weights for each frame there
################################################################################

####################################################################
# then calculate clusters corresponding to our step (based on eRMSD)
####################################################################
rm bck*
plumed driver --plumed plumed_print_ensembles.dat --mf_xtc TRAJ.xtc --kt 2.494339 --pdb ../partial_compl/md0.pdb
rm clusters.dat clusters_ss.dat


#only steps without slip-stranded states
awk '{if($6<0.7 && $10>0.7) cluster=3; else if($5<0.7 && $4>0.8 && $9<0.7) cluster=2; else if($2<0.7 && $3 > 0.8 && $4>0.8 && $8<0.7) cluster=1; else if($10<0.7 && $6>0.7) cluster=0; else cluster=9} {print cluster}' bps.dat  >> clusters.dat

#only steps with slip-stranded states
awk '{if($6<0.7 && $10>0.7) cluster=3; else if($5<0.7 && $4>0.8 && ($11<0.7 || $12<0.7)) cluster=2; else if($2<0.7 && $3 > 0.8 && $4>0.8 && $12<0.7) cluster=1; else if($13<0.7 && $6>0.7) cluster=0; else if(($10<0.7 && $6>0.7)) cluster=5; else cluster=9} {print cluster}' bps.dat  >> clusters_ss.dat
# cluster 5 here is BSL

rm eRMSDoutWeights.dat eRMSDoutWeights_ss.dat
paste clusters.dat ../weight/weights_norm.dat >> eRMSDoutWeights.dat
paste clusters_ss.dat ../weight/weights_norm.dat >> eRMSDoutWeights_ss.dat

# check if eRMSDoutWeights.dat  the first line has weight 1 then remove that line

echo "without slip-stranded states"
echo "step  dG_overall  dG_wrtBSL"
start=`awk '{if($1==0) q1+=$3; else q2+=$3} END{print q1}' eRMSDoutWeights.dat`
#awk '{if($1==0) q1+=$3; else q2+=$3} END{print "step 0: ", -0.289*8.3144598*(log(q1/q2))/4.184, -0.289*8.3144598*(log(q1/start))/4.184}' start=$start eRMSDoutWeights.dat
awk '{if($1==1) q1+=$3; else q2+=$3} END{print "step 1: ", -0.289*8.3144598*(log(q1/q2))/4.184, -0.289*8.3144598*(log(q1/start))/4.184}' start=$start eRMSDoutWeights.dat
awk '{if($1==2) q1+=$3; else q2+=$3} END{print "step 2: ", -0.289*8.3144598*(log(q1/q2))/4.184, -0.289*8.3144598*(log(q1/start))/4.184}' start=$start eRMSDoutWeights.dat
awk '{if($1==3) q1+=$3; else q2+=$3} END{print "step 3: ", -0.289*8.3144598*(log(q1/q2))/4.184, -0.289*8.3144598*(log(q1/start))/4.184}' start=$start eRMSDoutWeights.dat


echo ""
echo "with slip-stranded states"
echo "step  dG_overall  dG_wrtBSL"
awk '{if($1==0) q1+=$3; else q2+=$3} END{print "step 0: ", -0.289*8.3144598*(log(q1/q2))/4.184, -0.289*8.3144598*(log(q1/start))/4.184}' start=$start eRMSDoutWeights_ss.dat
awk '{if($1==1) q1+=$3; else q2+=$3} END{print "step 1: ", -0.289*8.3144598*(log(q1/q2))/4.184, -0.289*8.3144598*(log(q1/start))/4.184}' start=$start eRMSDoutWeights_ss.dat
awk '{if($1==2) q1+=$3; else q2+=$3} END{print "step 2: ", -0.289*8.3144598*(log(q1/q2))/4.184, -0.289*8.3144598*(log(q1/start))/4.184}' start=$start eRMSDoutWeights_ss.dat
awk '{if($1==3) q1+=$3; else q2+=$3} END{print "step 3: ", -0.289*8.3144598*(log(q1/q2))/4.184, -0.289*8.3144598*(log(q1/start))/4.184}' start=$start eRMSDoutWeights_ss.dat
echo "step 0 is going from canonical BSL to slip-stranded BSL"


################################################################
# now estimate the error and energy difference to a custom step:
################################################################

rm all_bias.dat all_bias2.dat all_bias_ss2.dat all_bias_ss.dat
paste ../weight/weights_norm.dat clusters.dat ../reweigt_partial/bias.dat  ../reweigt_partial/weights_norm.dat >> all_bias.dat
awk '{print $1, $3, $5, $7}' all_bias.dat >> all_bias2.dat

# modify header to: #time cluster t1.bias weight
paste ../weight/weights_norm.dat clusters_ss.dat ../reweigt_partial/bias.dat  ../reweigt_partial/weights_norm.dat >> all_bias_ss.dat
awk '{print $1, $3, $5, $7}' all_bias_ss.dat >> all_bias_ss2.dat

# example calculation for the pathway without slip-strand (all_bias.dat file)
# this will give energy differenrce (and error estimates) beween the selected step and a reference state
# update in the blocks_bs3.py code the reference state (bsl) manually !!
# the for loop below runs the bootstrapping code with different block sizes/numbers
step = 3
C="all_bias2.dat"
for i in 1 2 4 40 400 800 4000 6250 78125 10000 12500 15625 20000 25000 31250 50000 62500 100000 125000 250000; do python blocks_bs3.py $i cluster $step $C ; done
