#ifndef FASTGREEDYLG_H
#define FASTGREEDYLG_H


#include <vector>

#include "../../graph/graphcoarsenersinglesteplg.h"
#include "../../graph/solution.h"




struct FGCommunityLG{
    long double internalEdges;
    long double totalDegree;
    long double nnodes;
    long double density;

    FGCommunityLG(){}
    FGCommunityLG(long double ie, long double rtd, long double nn, long double dens):
        internalEdges(ie), totalDegree(rtd), nnodes(nn), density(dens){}

};

class FastGreedyLG
{
public:
    LargeGraph *graph;
    long double bestDensity;

    unsigned maxIt;
    unsigned bestIt;
    long long int totalTime;
    long long int bestTime;



    Solution * basedSolution=NULL;
    GraphCoarsenerSingleStepLG * gcss = NULL; //graph coarsened, which includes the original graph pointer
    vector<unsigned> nodeCommunity; //community which each node is
    vector<FGCommunityLG> communities; //communities information
    long double iter;

    Solution * executeAndConstruct(long double minimalDensity=MINUS_INFINITY);

    FastGreedyLG(LargeGraph * graph);
    FastGreedyLG(Solution *base);
    FastGreedyLG(GraphCoarsenerSingleStepLG * gcss);

    void execute(long double minimalDensity=MINUS_INFINITY); //execute the heuristic and returns a solution
    void update(GraphCoarsenerSingleStepLG *gcss); //just update: do not create new elements

    long long int t1 = 0;
    long long int t2 = 0;
    long long int t3 = 0;
    long long int t4 = 0;
    long long int t5 = 0;



    Solution * toSolution();
};

#include "fastgreedylg.cpp"

#endif // FASTGREEDYLG_H
