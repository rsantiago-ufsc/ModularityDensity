
#ifndef SSENDLVMOD_H
#define SSENDLVMOD_H

//SSEndLV makes a SingleSteep search, then the Louvain+Density search




#include "../../../graph/largegraph.h"
#include "../../../graph/solution.h"
#include "../../../utils/modularitylg.h"

//storing the nodes coarsened by each level with
// ... important data (internal edges, total degree, and number of nodes)
// ... by defining COLLECT_DATA_CONSTRUCTIVE
#define COLLECT_DATA_CONSTRUCTIVE


#include "../../constructive/singlestepgreedylg.h"
#include "../multilevelsslv/louvainsselv.h"



class MultiLevelSSLVMod
{
public:

    LargeGraph * graph;

    unsigned typeLouvain;

    //predefined type=delta_density
    MultiLevelSSLVMod(LargeGraph * graph, unsigned typeLouvain=1);


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

#include "multilevelsslvmod.cpp"

#endif // SSENDLV_H
