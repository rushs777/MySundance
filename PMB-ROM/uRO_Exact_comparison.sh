#!/bin/bash

#This bash script expects 4 arguments
#nxStart is the power of 2 to start nx with
#nxEnd is the power of 2 to end nx with
#ntStart is the power of 2 to start nSteps with
#ntEnd is the power of 2 to end nSteps with
#filename is the name of the file to store the output from the runs

filename=uRO_comparisonLog.txt

nxStart=$1
let nxEnd=$2+1
ntStart=$3
let ntEnd=$4+1

echo $nxStart $nxEnd $ntStart $ntEnd
for ((nt=$ntStart; nt<$ntEnd; nt++)); do
    for ((nx=$nxStart; nx<$nxEnd; nx++)); do
	echo Starting run for nx = $[ 2**nx ], nt = $[ 2**nt ]
	./uRO.exe --nx=$[ 2**nx ] --nSteps=$[ 2**nt ] --verbosity=0 >> $filename
    done
done


