#ifndef SINGLESTEPGREEDYLG_H
#define SINGLESTEPGREEDYLG_H

#include<vector>
#include<chrono>



#include "../../graph/largegraph.h"
#include "../../graph/solution.h"
#include "../../graph/lmbdgraphcoarsenersinglesteplg.h"
#include "../../utils/utils.h"

using namespace std;
using namespace std::chrono;



class LMBDSingleStepGreedyLG
{
public:

    float lambda = 0.5;

    long double bestDensity=MINUS_INFINITY;
    long double bestMod=MINUS_INFINITY;
    long double singletonMod=MINUS_INFINITY;
    unsigned bestCommSize;

#ifdef GROUND_TRUTH
    string bestCommunity = "";
#endif

    unsigned maxIt;
    unsigned bestIt;
    long long int totalTime;
    long long int bestTime;

    vector< std::pair< std::pair<unsigned, unsigned>, long double> > levels;

    LMBDGraphCoarsenerSingleStepLG * gcss=NULL;


    LMBDSingleStepGreedyLG(LargeGraph * graph, float lambda=0.5);


    void execute(); //execute the heuristic
    Solution * bestToSolution();
    string bestToString();


    //for multilevel heuristics
    void start();
    void nextLevel();
    void stop();


#ifdef COLLECT_DATA_CONSTRUCTIVE
    //collecting data to other metamethods, like multi-level heuristics
    vector< std::pair<NodeDetail, NodeDetail> > levelDetails;
#endif

private:
    LargeGraph * graph;



};

#include "lmbdsinglestepgreedylg.cpp"


#endif // SINGLESTEPGREEDYLG_H
