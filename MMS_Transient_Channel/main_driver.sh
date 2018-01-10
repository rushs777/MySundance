#!/bin/bash

# args is the parameters to pass
# First number is the number of Re
# Then the Re values to use
# Then the nx,nSteps values to use
# HAVE TO CHANGE tFinal IN EACH ONE MANUALLY
args=(2 10 100 24 36)

#echo "================== RUNNING SIMULATIONS ================= "
#cd ForwardProblem
#./driver_simulation.sh "${args[@]}"
echo "================== FORMING PODS ================= "
#cd ../POD/SingleParameterSpace
#./driver_POD.sh "${args[@]}"
echo "================== STAGING PODS ================= "
#./driver_stager.sh "${args[@]}"
echo "================== RUNNING KKT ================= "
cd ODECO
#cd ../../ODECO
./driver_KKT.sh "${args[@]}"
