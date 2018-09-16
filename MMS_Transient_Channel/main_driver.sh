#!/bin/bash

# args is the parameters to pass
# First number is the number of Re
# Then the Re values to use
# Then the nx,nSteps values to use
# HAVE TO CHANGE tFinal IN EACH ONE MANUALLY

# args for working with the sensor/ODECO
args=(6 1 20 40 60 80 100 24 36 48 60 72)

# args for testing multiple KKT solution against Re not used to generate
#args=(4 10 25 50 75 24 36 48 60 72)

# args for the Viento verification process
#args=(8 1 2 4 8 16 32 64 128 8 16 32 64 128)

#echo "================== RUNNING SIMULATIONS ================= "
#cd ForwardProblem
#./driver_simulation.sh "${args[@]}"
#echo "================== FORMING PODS ================= "
#cd ../POD/SingleParameterSpace
#./driver_POD.sh "${args[@]}"
#echo "================== STAGING PODS ================= "
#./driver_stager.sh "${args[@]}"
#echo "================== RUNNING INVERSE PROBLEM ================= "
#cd InverseProblem
#./driver_InvProb.sh "${args[@]}"
echo "================== RUNNING KKT ================= "
cd ODECO
#cd ../../ODECO
./driver_KKT.sh "${args[@]}"
