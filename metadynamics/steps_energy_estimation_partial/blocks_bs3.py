#!/local/usr/bin/python3

'''
Used for Bootstrapping analysis with divided blocks 
1. Divide the COLVAR into n blocks
2. Calculate the free energies of different substates in each block
3. Random picking up x blocks for t times
4. Calculate the standard deviations of the free energies in the four substates

======
taken from Zhengue Zhang at https://github.com/sponerlab/WT-MetaD-HREX_HJ/blob/main/blocks_bs.py
 and modified by Pavlina Pokorna

Colvar file should have:
time cv1 metad.bias weight(1/0)

see also: https://gtribello.github.io/mathNET/block_averaging_video.html
====
'''

import os
import sys
import numpy as np
import math
import warnings
warnings.filterwarnings("ignore")

def is_number(string):
    # Check if string is a valid number (float)
    try:
        float(string)
    except ValueError:
        return False
    return True

def make_rand_traj(blocks, trajsize):
    # Generate trajectory indexes for one iteration, returns array of indexes of blocks that are to be concaternated to the traj
    blockpick = np.random.choice(blocks, size=trajsize, replace=True) #.tolist()
    return blockpick #traj

def cal_dG(trajblock):
    # Calculate free energies in different substates
    boundary_bsl = 2.0 # cluster of BSL 
    boundary_bph  = float(cluster) #cluster of my state
    kbt = 2.478957

    tot_w = 0.0
    bsl_w = 0.0
    bph_w = 0.0

    replicaweight = 1.0 # this is a leftover from multi-replica simulations I analysed before where the weight of non-reference ones was set to 0
#    print(boundary_bsl)
#    print(boundary_bph)
    for frame in trajblock:
#        print(frame)
        weight=(math.exp(frame[2]/kbt))*replicaweight
        tot_w += weight
        if frame[0] == boundary_bsl:
            bsl_w += weight
        elif frame[0] == boundary_bph:
            bph_w += weight

    try:
        diff = -kbt * np.log(bsl_w/bph_w)
    except:
        diff = float('inf') 

    return diff

if __name__ == "__main__":
    
    # First import arguments and check
    if len(sys.argv) != 5:
        print("Usage: python3 ./blocks_bs.py <block num> cluster <cluster> <COLVAR file>")
        quit()
    else:
        if not is_number(sys.argv[1]):
            print("Block or replica number should be number")
            quit()
        if not os.path.isfile(sys.argv[4]):
            print("Bias file not exists!")
            quit()
        bnum  = int(sys.argv[1])
        cv1   = sys.argv[2] #cluster
        cluster   = sys.argv[3] #cluster
        fname = sys.argv[4]

    # Get values from bias file
    colvar = np.loadtxt(fname, comments = '#')
    # Get the header
    with open(fname, 'r') as inp:
        header = inp.readline()
        if header.startswith('#'):
            header = header.replace("#! FIELDS ","").split()
            if not cv1 in header:
                print("Input CVs not in the file!")
                quit()
            ind_cv1 = header.index(cv1)
            ind_bias = header.index("t1.bias")
    
    # Only keep our CVs, bias and replica number
    colvar = colvar[:,[ind_cv1,ind_cv1,ind_bias]]
    # split this into blocks:
    timescale = colvar.shape[0]
    blocks=[]
    blocklength=round(timescale/bnum)
 #   print("block lenght: ")
 #   print(blocklength)
    i=0
    for blocknum in range(0,bnum):
        blocks.append(colvar[i:i+blocklength])
        i+=blocklength
    #calculate dG for each block
    block_dGs=[]
    for block in blocks:
        block_dGs.append(cal_dG(block))
 #   print("dG arr for the individual blocks: ")
#    print(block_dGs)
    #only considesr blocks where folding was sampled
    block_idxs=[]
    j=0
    for block in block_dGs:
        if not np.isfinite(block):
                pass
        else:
            block_idxs.append(j)
        j+=1
 #   print("indexes for dG arr for the individual blocks with both target states sampled: ")
 #   print(block_idxs)

    traj = make_rand_traj(block_idxs, 200)
    all_folded=[]

    for trajblock_idx in traj:
        all_folded.append(block_dGs[trajblock_idx])


    kjtokcal=0.238846
 #   print(all_folded)


    all_folded_array = np.array(all_folded)
    dg_fold=kjtokcal*np.average(all_folded)
    std_fold=kjtokcal*np.std(all_folded)
    error=std_fold/(math.sqrt(len(block_idxs)))
    #print(dg_fold)

    print(str(dg_fold) + "\t" + str(error) + "\t" + str(blocklength) )

