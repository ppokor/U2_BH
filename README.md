# U2-intron strand invasion modeling

Preprint is here: https://www.biorxiv.org/content/10.1101/2025.10.31.685824v1 
 
This directory contains the simulation setup and some analysis files

initial coordinates (md0.rst7 files) are those after equilibration

\*\_stripWAT.parm7 files are water-stripped topologies for the deposited trajectories

ZN parameters are in the U2SF3b\_complex/, see the example PDB files there for residue naming 
PSU (pseudouracil) parameters are in the force\_field\_PSU, see the example PDB file there for residue naming

Run the calculations with GROMACS:

    gmx_mpi mdrun -s MTD.tpr -deffnm md -ntomp 8 -v -nb gpu -pme gpu -npme 1 -pin on -cpi md.cpt -nsteps 200000000 -noappend
    gmx_mpi mdrun -s MTD.tpr -deffnm md -ntomp 8 -v -nb gpu -pme gpu -npme 1 -pin on -plumed plumed.dat -cpi md.cpt -nsteps 200000000 -noappend

trajectories and HILLS files are deposited at Zenodo (link will be included later)

Coarse-grained structures and trajectories can be visualized with the oxView web server.

