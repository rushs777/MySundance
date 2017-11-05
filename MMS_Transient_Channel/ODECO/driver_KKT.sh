#!/bin/bash

#This bash script expects 2 arguments
# The first argument is a list of the Reynolds numbers for which to generate the POD basis, with
# the first element in the list being the number of Reynolds numbers
# The second argument is a list of integers for nx,nSteps
# These lists can be enclosed in {} and comma-separated, or simply as a series of space
# separated values; but the {} notation is easier for readibility

#filename is the name of the file to store the output from the runs
filename=log_test.txt
executable=ODECO.exe

rm $filename

ReArray=( "${@:2:$1}" ); shift "$(( $1 + 1))"
nxnt=( "$@" )

for i in "${!ReArray[@]}"
do
#    echo "Value of Re = ${ReArray[$i]}"
    for j in "${!nxnt[@]}"
    do
#	echo "Value of nxnt = ${nxnt[$j]}"
	echo "Starting run for Re = ${ReArray[$i]}, nx=nSteps = ${nxnt[$j]}"
	./$executable --nx=${nxnt[$j]} --nSteps=${nxnt[$j]} --Re=${ReArray[$i]} --verbosity=0 2>&1 | tee -a $filename
    done
done
