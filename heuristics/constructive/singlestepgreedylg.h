#ifndef SINGLESTEPGREEDYLG_H
#define SINGLESTEPGREEDYLG_H

#include<vector>
#include<chrono>



#include "../../graph/largegraph.h"
#include "../../graph/solution.h"
#include "../../graph/graphcoarsenersinglesteplg.h"
#include "../../utils/utils.h"

using namespace std;
using namespace std::chrono;



class SingleStepGreedyLG
{
public:

    long double bestDensity=MINUS_INFINITY;
    long double bestMod=MINUS_INFINITY;
    long double singletonMod=MINUS_INFINITY;
    unsigned bestCommSize;

    unsigned maxIt;
    unsigned bestIt;
    long long int totalTime;
    long long int bestTime;

    vector< std::pair< std::pair<unsigned, unsigned>, long double> > levels;

    GraphCoarsenerSingleStepLG * gcss=NULL;


    SingleStepGreedyLG(LargeGraph * graph);


    void execute(); //execute the heuristic
    Solution * bestToSolution();

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

#include "singlestepgreedylg.cpp"


#endif // SINGLESTEPGREEDYLG_H
