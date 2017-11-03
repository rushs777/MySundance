#!/bin/bash

#This bash script expects 2 arguments
#start is the power of 2 to start nx and nSteps with
#end is the power of 2 to end nx and nSteps with
#filename is the name of the file to store the output from the runs

#filename=log_forward_problem_TransientChannel_tol8.txt
filename=log_test.txt
filenameShort=log_test_short.txt
executable=forward_problem_TransientChannel.exe

#start=$1
#let end=$2+1

rm $filename

#for ((p=$start; p<$end; p++)); do
#    echo Starting run for nx = $[ 2**p ], nt = $[ 2**p ]
    #	./uRO.exe --nx=$[ 2**nx ] --nSteps=$[ 2**nt ] --verbosity=0 >> $filename
#    ./$executable --nx=$[ 2**p ] --nSteps=$[ 2**p ] --verbosity=1 2>&1 | tee -a $filename
#done


# Trying to pass in the Reynolds Number and a list of integers for nx,nSteps
# The list can be enclosed in {} and comma-separated, or it seems that spaces work
# as well; however, the {} approach is easier for readibility

ReArray=( "${@:2:$1}" ); shift "$(( $1 + 1))"
nxnt=( "$@" )

for i in "${!ReArray[@]}"
do
#    echo "Value of Re = ${ReArray[$i]}"
    for j in "${!nxnt[@]}"
    do
#	echo "Value of nxnt = ${nxnt[$j]}"
	echo "Starting run for Re = ${ReArray[$i]}, nx=nSteps = ${nxnt[$j]}"
	./$executable --nx=${nxnt[$j]} --nSteps=${nxnt[$j]} --Re=${ReArray[$i]} --verbosity=1 2>&1 | tee -a $filename
    done
done

grep -i Re= $filename > $filenameShort

# Prints out everything in the array
#declare -p ReArray
#declare -p nxnt

# How get the indices for an array and then use them to access elements
#for i in "${!foo[@]}"; do 
#  printf "%s\t%s\n" "$i" "${foo[$i]}"
#done

# How to do multiple arrays where the first index of each array is the length
#array1=( "${@:2:$1}" ); shift "$(( $1 + 1 ))"
#array2=( "${@:2:$1}" ); shift "$(( $1 + 1 ))"
#array3=( "${@:2:$1}" ); shift "$(( $1 + 1 ))"

