#!/bin/bash

TOP=u2BSL_intron_partial_PSU.parm7

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
gmx_mpi grompp -f ../MD.mdp -c gromacs.gro -p gromacs.top -o MD.tpr -maxwarn 1


