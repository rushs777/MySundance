#!/bin/bash

#This bash script expects 2 arguments
# The first argument is a list of the Reynolds numbers for which to generate the POD basis, with
# the first element in the list being the number of Reynolds numbers
# The second argument is a list of integers for nx,nSteps
# These lists can be enclosed in {} and comma-separated, or simply as a series of space
# separated values; but the {} notation is easier for readibility

# UPDATE: After Christmas 2017, the {} seems to have broken inexplicitly. Thus have to use the
# space separated version now

#filename is the name of the file to store the output from the runs
filename=log_test.txt
filenameShort=log_test_short.txt
executable=forward_problem_TransientChannel.exe

echo "removing old files"
rm -f $filename

ReArray=( "${@:2:$1}" ); shift "$(( $1 + 1))"
nxnt=( "$@" )
# in the current iteration, we are letting nt be the number of timesteps per second
# Thus the real number of time steps is the time of the simulation (tf) times nt
tf=1

echo "start loop"
echo "nxnt=" ${nxnt}
echo "Re array=" "${ReArray[@]}"

for i in "${!ReArray[@]}"
do
    echo "Value of Re = ${ReArray[$i]}"
    for j in "${!nxnt[@]}"
    do
	echo "Value of nxnt = ${nxnt[$j]}"
	nx=${nxnt[$j]}
	nt=$((${tf} * ${nxnt[$j]}))
	echo "Starting run for Re = ${ReArray[$i]}, nx=$nx, nt=${nt}, tf=${tf}"
	./$executable --nx=${nx} --nSteps=${nt} --tFinal=${tf} --Re=${ReArray[$i]} --verbosity=1 2>&1 | tee -a $filename
    done
done

grep -i tFinal= $filename > $filenameShort

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

