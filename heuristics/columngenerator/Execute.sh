#!/bin/bash

rm $1
rm $1.o


echo "Compiling.."
g++ -O0 -c -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -I/opt/ibm/ILOG/CPLEX_Studio128/cplex/include -I/opt/ibm/ILOG/CPLEX_Studio128/concert/include $1.cpp -o $1.o 

echo "Linking..."
g++ -O0 -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -I/opt/ibm/ILOG/CPLEX_Studio128/cplex/include -I/opt/ibm/ILOG/CPLEX_Studio128/concert/include -L/opt/ibm/ILOG/CPLEX_Studio128/cplex/lib/x86-64_linux/static_pic -L/opt/ibm/ILOG/CPLEX_Studio128/concert/lib/x86-64_linux/static_pic -o $1 $1.o -lconcert -lilocplex -lcplex -lm -lpthread -ldl

echo "Starting Execution!"
./$1 $2 $3 $4 $5 $6

