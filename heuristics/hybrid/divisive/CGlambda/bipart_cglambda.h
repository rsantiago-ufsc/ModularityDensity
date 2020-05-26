#ifndef BIPART_CGI_ILS_H
#define BIPART_CGI_ILS_H


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <cstdlib>
#include <time.h>
#include <ilcplex/ilocplex.h>
#include <ilconcert/iloenv.h>
#include <chrono>
#include <list>

#include "../../../../graph/largegraph.h"
#include "../../../../graph/solution.h"
#include "../../../../utils/modularitylg.h"





ILOSTLBEGIN

#define RC_EPS 1.0e-9

/**
 * @brief The Bipart_CGI_ILS class
 *
 * Column generation to find the best cut of a given cluster.
 * This version uses the model CGI for auxiliar problem, and
 * the heuristic ILS.
 */

class Bipart_CGLambda
{
public:

    //graph instance
    LargeGraph * graph;

    //the first and the second half of the cluster optimized by
    //... 'execute' procedure.
    vector<unsigned> firstHalf;
    vector<unsigned> secondHalf;

    //auxiliar variables
    IloNum **Wij;   //adjacency matrix
    IloNum *Di;     //degree of nodes

    //maximal modularity obtained in the execution procedure
    IloNum maxDensity;

    //parameter of ILS heuristic
    double alpha;

    Bipart_CGLambda(LargeGraph * lg, double alpha);

    //executes the column generator
    void execute(vector<unsigned> &nodes);


    //auxiliar functions
    void coefOF(vector<unsigned> &cluster, IloNum &coef, IloNumArray &A);
    void prepareWij(IloNum **&wij);
    void prepareDi (IloNum *&di);

};

#include "bipart_cglambda.cpp"


#endif // BIPART_CGI_ILS_H
