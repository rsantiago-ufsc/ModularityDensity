/// Hybrid heuristic which first applies a MCN (LVDENS), then LNM.
/// We will try to scape from the failure of the MCN for some instances.

#ifndef MCNLNM_H
#define MCNLNM_H



#include "../../../graph/largegraph.h"
#include "../../../heuristics/hybrid/LAMBDAlouvain/lmbdlouvain.h"
#include "../../../heuristics/constructive/lmbdfastgreedylg.h"
class LMBDMcnLnm
{
public:
    LargeGraph * graph;

    float lambda = 0.5;

    LMBDLouvain * hlv;
    LMBDFastGreedyLG *fglg;
    bool FastGreedyImproved = false;

    LMBDMcnLnm(LargeGraph *graph, float lambda);



    long double bestDensity;

    void execute();

};

#include "lmbdmcnlnm.cpp"

#endif // MCNLNM_H
