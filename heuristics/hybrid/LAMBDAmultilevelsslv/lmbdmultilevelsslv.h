#ifndef SSENDLV_H
#define SSENDLV_H

//SSEndLV makes a SingleSteep search, then the Louvain+Density search




#include "../../../graph/largegraph.h"
#include "../../../graph/solution.h"
#include "../../../utils/modularitylg.h"

//storing the nodes coarsened by each level with
// ... important data (internal edges, total degree, and number of nodes)
// ... by defining COLLECT_DATA_CONSTRUCTIVE
#define COLLECT_DATA_CONSTRUCTIVE


#include "../../constructive/lmbdsinglestepgreedylg.h"
#include "lmbdlouvainsselv.h"



class LMBDMultiLevelSSLV
{
public:

    LargeGraph * graph;

    unsigned typeLouvain;

    float lambda = 0.5;

    //predefined type=delta_density
    LMBDMultiLevelSSLV(LargeGraph * graph, float lambda=0.5, unsigned typeLouvain=2);


    //starting uncoarsening from last level partition
    void execute();
    //starting uncoarsening from level of the best partition
    void executeFromBest();
    //uses coarsened nodes from SingleStep on singleton partitions
    void executeAlwaysSingleton();

    void reorganization(LMBDGraphCoarsenerSingleStepLG *&gcss, vector <unsigned> & commByVertice);

    unsigned it;

    long double bestDensity;
    long double bestMod;
    unsigned bestCommSize;
    long long int bestTime;
    long long int totalTime;
    unsigned bestIt;
    string bestCommStr;

};

#include "lmbdmultilevelsslv.cpp"

#endif // SSENDLV_H
