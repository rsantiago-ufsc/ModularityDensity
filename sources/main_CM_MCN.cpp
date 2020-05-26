#include <iostream>
#include <chrono>



#include "../../utils/modularitylg.h"
#include "../../graph/largegraph.h"
#include "../../graph/solution.h"
#include "../../heuristics/hybrid/louvain/louvain.h"
#include "../../heuristics/constructive/singlestepgreedylg.h"
#include "../../heuristics/constructive/fastgreedylg.h"

using namespace std;
using namespace std::chrono;


#define TYPE_PRI_NOW TYPE_PRI

int main(int argc, char *argv[])
{
    chrono::system_clock::time_point before;
    before = chrono::system_clock::now();

    string fileInstance = argv[1];
    string filepath = "../../instances/"+fileInstance;

    string opt = argv[2];
    string exp = argv[3];

    srand(time(NULL)+atoi(exp.c_str()));

    LargeGraph lg(filepath);


    long long int totalTime;





    SingleStepGreedyLG sslg(&lg);
    sslg.execute();

    totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

    //construindo a solucao
    //ModularityLG mlg(&lg);
    Solution *solSS = sslg.bestToSolution();


    before = chrono::system_clock::now();
    Louvain hlv(solSS, &lg,TYPE_PRI_NOW);
    hlv.execute();

    //cout<<sslg.bestDensity<<", "<<hlv.bestDensity <<"\n";

    totalTime += chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

    vector<string> ix = Utils::splitString(fileInstance,'.');
    string instance= ix[ix.size()-2];
    ix = Utils::splitString(instance,'/');
    instance= ix[ix.size()-1];



    //storing the data
    string expFile="../data/CM_MCN"+instance+"_"+exp;
    ofstream f(expFile, std::fstream::app);
    f<<"CM+MCN,"<<instance<<","<<lg.numberOfNodes<<","<<lg.numberOfEdges<<","
           <<exp<<","
           <<hlv.it<<","<<hlv.it<<","
           <<opt<<","<< hlv.bestDensity << "," <<hlv.bestMod<<","
           <<hlv.bestNumberOfCommunities<<","
           <<hlv.it<<","<<hlv.it<<","
           <<totalTime<<","<<hlv.totalTime<<","<<hlv.totalTime
           <<"\n";
    f.close();


    return 0;
}


