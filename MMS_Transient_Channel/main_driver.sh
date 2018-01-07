#!/bin/bash

# args is the parameters to pass
# First number is the number of Re
# Then the Re values to use
# Then the nx,nSteps values to use
# HAVE TO CHANGE tFinal IN EACH ONE MANUALLY
args=(2 10 100 24 36)

cd ForwardProblem
echo "================== RUNNING SIMULATIONS ================= "
./driver_simulation.sh "${args[@]}"
cd ../POD/SingleParameterSpace
echo "================== FORMING PODS ================= "
./driver_POD.sh "${args[@]}"
echo "================== STAGING PODS ================= "
./driver_stager.sh "${args[@]}"
echo "================== RUNNING KKT ================= "
cd ODECO #../../ODECO
./driver_KKT.sh "${args[@]}"
