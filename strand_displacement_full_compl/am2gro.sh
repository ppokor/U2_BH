#!/bin/bash

TOP=u2BSL_intron_full.parm7
parmed -e <<EOF
!!
import parmed as pmd
amber = pmd.load_file('$TOP', 'md0.rst7')

# Save GROMACS topology and GRO file
amber.save('gromacs.top')
amber.save('gromacs.gro')
!!
EOF

#normal topology
gmx_mpi grompp -f ../MTD.mdp -c gromacs.gro -p gromacs.top -o MTD.tpr -maxwarn 1

#stafix-scaled topology
python ../scaleSTAFIX_gromacs.py gromacs.top 0.9
gmx_mpi grompp -f ../MTD.mdp -c gromacs.gro -p gromacsSTAFIX0.9.top -o MTD_stafix.tpr -maxwarn 1

