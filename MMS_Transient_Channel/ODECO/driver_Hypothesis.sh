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
executable=Hypothesis.exe

echo "removing old files"
rm $filename
rm $filenameShort

ReArray=( "${@:2:$1}" ); shift "$(( $1 + 1))"
nxnt=( "$@" )
# in the current iteration, we are letting nt be the number of timesteps per second
# Thus the real number of time steps is the time of the simulation (tf) times nt
tf=1

space='multiple'
for i in "${!ReArray[@]}"
do
#    echo "Value of Re = ${ReArray[$i]}"
    for j in "${!nxnt[@]}"
    do
#	echo "Value of nxnt = ${nxnt[$j]}"
	nx=${nxnt[$j]}
	nt=$((${tf} * ${nxnt[$j]}))
	for sens in 1
	do
	    echo "Starting run for Re = ${ReArray[$i]}, nx=$nx, nt=${nt}, tf=${tf} nSens=${sens}"
	    if [ "$space" == "single" ]; then
		./$executable --nx=${nx} --nSteps=${nt} --tFinal=${tf} --Re=${ReArray[$i]} --numSens=${sens} --verbosity=0 --writeVTK --y0=.75 2>&1 | tee -a $filename
	    elif [ "$space" == 'multiple' ]; then
		./$executable --nx=${nx} --nSteps=${nt} --tFinal=${tf} --Re=${ReArray[$i]} --numSens=${sens} --verbosity=0 --writeVTK --parameterSpace="multiple" --y0=.75 2>&1 | tee -a $filename
	    fi
	done
    done
done

grep -A 10 --group-separator="" "tFinal=" $filename > $filenameShort

# Here's how to format the rest
# Copy over an old naming of the headers
# Use ctrl+H to get rid of all the labels. Remember to copy new line characters too
# Once you have all the labels removed, to get rid of the extra new lines,
#           Find: \n\n1
#           Replace with: \n1


