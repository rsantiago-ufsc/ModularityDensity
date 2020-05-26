#ifndef DIVISIVE_CGI_ILS_H
#define DIVISIVE_CGI_ILS_H

#include "../../../../graph/largegraph.h"
#include "bipart_cglambdacg.h"

/**
 * @brief The Divisive_CGI_ILS class
 *
 * Divisive heuristic used to solve Modularity Density Maximization.
 * This heuristic uses auxiliar problem CGI and heutistic ILS.
 */

class Divisive_CGLambdaCG
{
public:
    //graph instance
    LargeGraph * graph;

    //parameter used for heuristic ILS
    float alpha;

    //maximal density found in the execution process
    long double maxDensity;

    Divisive_CGLambdaCG(LargeGraph * graph, float alpha);

    //executes the divisive heuristic CGI+ILS
    void execute();

    //calculate the D value of a cluster
    long double calculateD( vector<unsigned> * cluster );
};


#include "divisive_cglambdacg.cpp"
#endif // DIVISIVE_CGI_ILS_H
