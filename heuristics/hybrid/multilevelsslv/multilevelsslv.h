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


#include "../../constructive/singlestepgreedylg.h"
#include "louvainsselv.h"



class MultiLevelSSLV
{
public:

    LargeGraph * graph;

    unsigned typeLouvain;

    //predefined type=delta_density
    MultiLevelSSLV(LargeGraph * graph, unsigned typeLouvain=2);


    //starting uncoarsening from last level partition
    void execute();
    //starting uncoarsening from level of the best partition
    void executeFromBest();
    //uses coarsened nodes from SingleStep on singleton partitions
    void executeAlwaysSingleton();

    void reorganization(GraphCoarsenerSingleStepLG *&gcss, vector <unsigned> & commByVertice);

    unsigned it;

    long double bestDensity;
    long double bestMod;
    unsigned bestCommSize;
    long long int bestTime;
    long long int totalTime;
    unsigned bestIt;
    string bestCommStr;

};

#include "multilevelsslv.cpp"

#endif // SSENDLV_H
