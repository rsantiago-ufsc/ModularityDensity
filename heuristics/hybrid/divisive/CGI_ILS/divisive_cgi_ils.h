#ifndef DIVISIVE_CGI_ILS_H
#define DIVISIVE_CGI_ILS_H

#include "../../../../graph/largegraph.h"
#include "bipart_cgi_ils.h"


#define GROUND_TRUTH
#include "../../../../graph/largegraph.h"
#include "../../../../graph/solution.h"
#include "../../../../utils/modularitylg.h"
#include "../../mcnlnm/mcnlnm.h"
#include <list>

#include "../../../../graph/solution.h"

/**
 * @brief The Divisive_CGI_ILS class
 *
 * Divisive heuristic used to solve Modularity Density Maximization.
 * This heuristic uses auxiliar problem CGI and heutistic ILS.
 */

class Divisive_CGII_ILS
{
public:
    //graph instance
    LargeGraph * graph;

    //parameter used for heuristic ILS
    float alpha;

    //list of solutions
    list<Solution *> solutions;

    //maximal density found in the execution process
    long double maxDensity;

    Divisive_CGII_ILS(LargeGraph * graph, float alpha);

    //executes the divisive heuristic CGI+ILS
    void execute();

    //calculate the D value of a cluster
    long double calculateD( vector<unsigned> * cluster );
};


#include "divisive_cgi_ils.cpp"
#endif // DIVISIVE_CGI_ILS_H
