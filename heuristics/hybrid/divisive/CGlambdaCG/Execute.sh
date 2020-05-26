#!/bin/bash

rm $1
rm $1.o

echo "Compiling.."
g++ $1.cpp -o $1.o -c -m64 -O5 -fPIC -fexceptions -DNDEBUG -DIL_STD -DILOSTRICTPOD -I/opt/ibm/ILOG/CPLEX_Studio1261/cplex/include/ilcplex/ -I/opt/ibm/ILOG/CPLEX_Studio1261/cplex/include/ -I/opt/ibm/ILOG/CPLEX_Studio1261/concert/include/ -std=c++0x

echo "Linking..."
g++ $1.o -o $1 -m64 -O5 -fPIC -fexceptions -DNDEBUG -DIL_STD -DILOSTRICTPOD -I/opt/ibm/ILOG/CPLEX_Studio1261/cplex/include/ilcplex/ -I/opt/ibm/ILOG/CPLEX_Studio1261/cplex/include/ -L/opt/ibm/ILOG/CPLEX_Studio1261/cplex/lib/x86-64_linux/static_pic/ -lilocplex -lcplex -L /opt/ibm/ILOG/CPLEX_Studio1261/concert/lib/x86-64_linux/static_pic/ -lconcert -lm -lpthread -std=c++0x

echo "Starting Execution!"
./$1 $2 $3 $4
