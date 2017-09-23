#!/bin/bash

#This bash script expects 1 argument, which is the list of integer values
# to execute POD_Generator.exe with for nx and nSteps
# The list must be enclosed in {} and comma-separated.
#filename is the name of the file to store the output from the runs

filename=log_POD_Generation_tol0.999.txt
executable=POD_Generator.exe

rm $filename

for i in "$@"
do
    echo "Running for nx,nSteps = $i"
    ./$executable --nx=$i --nSteps=$i 2>&1 | tee -a $filename
done
