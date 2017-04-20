#!/bin/bash

#This bash script expects 2 arguments
#start is the power of 2 to start nx and nSteps with
#end is the power of 2 to end nx and nSteps with
#filename is the name of the file to store the output from the runs

filename=log_HG_Dirichlet_forward_simulation.txt

start=$1
let end=$2+1

rm $filename

for ((p=$start; p<$end; p++)); do
    echo Starting run for nx = $[ 2**p ], nt = $[ 2**p ]
    #	./uRO.exe --nx=$[ 2**nx ] --nSteps=$[ 2**nt ] --verbosity=0 >> $filename
    ./forward_problem_HG_Dirichlet.exe --nx=$[ 2**p ] --nSteps=$[ 2**p ] --verbosity=1 2>&1 | tee -a $filename
done


