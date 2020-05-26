/// Hybrid heuristic which first applies a MCN (LVDENS), then LNM.
/// We will try to scape from the failure of the MCN for some instances.

#ifndef MCNLNM_H
#define MCNLNM_H



#include "../../../graph/largegraph.h"
#include "../../../heuristics/hybrid/louvain/louvain.h"
#include "../../../heuristics/constructive/fastgreedylg.h"
class McnLnm
{
public:
    LargeGraph * graph;
    Louvain * hlv;
    FastGreedyLG *fglg;
    bool FastGreedyImproved = false;

    McnLnm(LargeGraph *graph);



    long double bestDensity;

    void execute();

};

#include "mcnlnm.cpp"

#endif // MCNLNM_H
