#!/bin/sh

vaspprog=/home/ivan/VASP/vasp.5.4.4/bin/vasp_std

## test if bash_profile is sourced
##showPaths.sh $LD_LIBRARY_PATH
echo "mpirun -np $1 $vaspprog > out_1" >> RUN
mpirun -np $1 $vaspprog > out_1
