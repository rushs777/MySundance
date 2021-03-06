#!/bin/bash

#This bash script expects 2 arguments
# The first argument is a list of the Reynolds numbers for which to generate the POD basis, with
# the first element in the list being the number of Reynolds numbers
# The second argument is a list of integers for nx,nSteps
# These lists can be enclosed in {} and comma-separated, or simply as a series of space
# separated values; but the {} notation is easier for readibility

# UPDATE: After Christmas 2017, the {} seems to have broken inexplicitly. Thus have to use the
# space separated version now

# Parse the passed array into the relevant pieces
ReArray=( "${@:2:$1}" ); shift "$(( $1 + 1))"
nxnt=( "$@" )

# filename is the prefix of the name of the file to store the output from the runs
filename="log_Re{"
# Name excutable the script will run
executable=POD_Generator.exe
# Define the directory where to store results
directory="Log_Files/"

# Build out the portion of the filename concering Re
for i in "${ReArray[@]:0:${#ReArray[@]}-1}"
do
    echo $i
    filename=$filename$i","
done
# Append the last Re value along with the next prefix in the filename
filename=$filename${ReArray[-1]}"}_nxnt{"
# Build out the portion of the filename concerning nx,nt
for j in "${nxnt[@]:0:${#nxnt[@]}-1}"
do
    filename=$filename$j","
done
# Append the closing brace
filename=$filename${nxnt[-1]}"}"
# Build the short filename
filenameShort=$filename"_short.txt"
# Add the file extension
filename=$filename".txt"

echo $filename
echo $filenameShort


# in the current iteration, we are letting nt be the number of timesteps per second
# Thus the real number of time steps is the time of the simulation (tf) times nt
tf=1

for i in "${!ReArray[@]}"
do
#    echo "Value of Re = ${ReArray[$i]}"
    for j in "${!nxnt[@]}"
    do
#	echo "Value of nxnt = ${nxnt[$j]}"
	nx=${nxnt[$j]}
	nt=$((${tf} * ${nxnt[$j]}))
	echo "Starting run for Re = ${ReArray[$i]}, nx=$nx, nt=${nt}, tf=${tf}"
	./$executable --nx=${nx} --nSteps=${nt} --tFinal=${tf} --Re=${ReArray[$i]} --verbosity=0 2>&1 | tee -a $filename
    done
done


# Only keep the main results
grep -i tFinal= $filename > $filenameShort

# Move the files nto the Log_Files directory
mv $filename $directory
mv $filenameShort $directory
