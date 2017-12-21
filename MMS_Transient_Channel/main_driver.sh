#!/bin/bash

#cd ForwardProblem
#echo "================== RUNNING SIMULATIONS ================= "
#./driver_simulation.sh {2,10,100} {24,36,48,60,72}
#cd ../POD/SingleParameterSpace
#echo "================== FORMING PODS ================= "
#./driver_POD.sh {2,10,100} {24,36,48,60,72}
#echo "================== STAGING PODS ================= "
#./driver_stager.sh {2,10,100} {24,36,48,60,72}
echo "================== RUNNING KKT ================= "
cd ODECO #../../ODECO
./driver_KKT.sh {2,10,100} {24,36,48,60,72}
