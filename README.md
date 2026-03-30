this directory contains simulation setup and some analysis files

initial coordinates (md0.rst7 files) are those after equilibration

\*\_stripWAT.parm7 files are water-stripped topologies for the deposited trajectories

ZN parameters are in the U2SF3b\_complex/, see the exmaple PDB files there for residue naming 
PSU (pseudouracil) parameters are in the force\_field\_PSU, see the exmaple PDB file there for residue naming

run the calculations with gromacs:
    gmx_mpi mdrun -s MTD.tpr -deffnm md -ntomp 8 -v -nb gpu -pme gpu -npme 1 -pin on -cpi md.cpt -nsteps 200000000 -noappend
    gmx_mpi mdrun -s MTD.tpr -deffnm md -ntomp 8 -v -nb gpu -pme gpu -npme 1 -pin on -plumed plumed.dat -cpi md.cpt -nsteps 200000000 -noappend

trajectories and HILLS files are deposited at Zenodo

coarse-grained structures and trajectories can be visualized with oxView webserver

