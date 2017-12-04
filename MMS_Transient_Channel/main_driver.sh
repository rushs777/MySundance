#!/bin/bash

cd ForwardProblem
./driver_simulation.sh {4,1,10,20,40} {25,50}
cd ../POD/SingleParameterSpace
./driver_POD.sh {4,1,10,20,40} {25,50}
./driver_stager.sh {4,1,10,20,40} {25,50}
cd ../../ODECO
./driver_KKT.sh {4,1,10,20,40} {25,50}
