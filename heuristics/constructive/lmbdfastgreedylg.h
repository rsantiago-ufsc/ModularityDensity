#ifndef FASTGREEDYLG_H
#define FASTGREEDYLG_H


#include <vector>

#include "../../graph/lmbdgraphcoarsenersinglesteplg.h"
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

class LMBDFastGreedyLG
{
public:

    float lambda = 0.5;

    LargeGraph *graph;
    long double bestDensity= std::numeric_limits<long double>::min() ;
    string bestSolution;
    unsigned bestNComm = 0;


    unsigned maxIt;
    unsigned bestIt;
    long long int totalTime;
    long long int bestTime;



    Solution * basedSolution=NULL;
    LMBDGraphCoarsenerSingleStepLG * gcss = NULL; //graph coarsened, which includes the original graph pointer
    vector<unsigned> nodeCommunity; //community which each node is
    vector<FGCommunityLG> communities; //communities information
    long double iter;

    Solution * executeAndConstruct(long double minimalDensity=MINUS_INFINITY);

    LMBDFastGreedyLG(LargeGraph * graph, float lambda = 0.5);
    LMBDFastGreedyLG(Solution *base, float lambda = 0.5);
    LMBDFastGreedyLG(LMBDGraphCoarsenerSingleStepLG * gcss, float lambda = 0.5);

    void execute(long double minimalDensity=MINUS_INFINITY); //execute the heuristic and returns a solution
    void update(LMBDGraphCoarsenerSingleStepLG *gcss); //just update: do not create new elements

    long long int t1 = 0;
    long long int t2 = 0;
    long long int t3 = 0;
    long long int t4 = 0;
    long long int t5 = 0;



    Solution * toSolution();
};

#include "lmbdfastgreedylg.cpp"

#endif // FASTGREEDYLG_H
