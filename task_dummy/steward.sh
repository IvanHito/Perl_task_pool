#!/bin/bash

rm WAIT
touch RUN

##
## Some vasp work will be done here
##
cd VASP
chmod +x run.sh
./run.sh
cd ../
##

rm RUN
mkdir result
./extract_dyn.pl $1
touch DONE
