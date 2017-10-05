#!/bin/bash

rm WAIT
touch RUN

##
## Some vasp work will be done here
##
cd VASP
chmod +x run.sh
./run.sh $2
##chmod +x run_test.sh
##./run_test.sh $2
cd ../
##

rm RUN
mkdir result
./extract_dyn.pl $1 > extrct_log
touch DONE
