#!/bin/bash

TOP=7evo_noTATSF1.parm7
parmed -e <<EOF
!!
import parmed as pmd
amber = pmd.load_file('$TOP', 'md0.rst7')

# Save GROMACS topology and GRO file
amber.save('gromacs.top')
amber.save('gromacs.gro')
!!
EOF

gmx_mpi grompp -f ../MTD.mdp -c gromacs.gro -p gromacs.top -o MTD.tpr -maxwarn 1

