#include <iostream>
#include <chrono>
#include "../../utils/utils.h"
#include "../../utils/modularitylg.h"
#include "../../graph/largegraph.h"
#include "../../graph/solution.h"
#include "../../heuristics/constructive/singlestepgreedylg.h"
#include "../../heuristics/constructive/fastgreedylg.h"

using namespace std;
using namespace std::chrono;




int main(int argc, char *argv[])
{


    string fileInstance = argv[1];
    string filepath = "../../instances/"+fileInstance;

    string opt = argv[2];
    string exp = argv[3];

    LargeGraph lg(filepath);

    ModularityLG mlg (&lg);

    chrono::system_clock::time_point before;
    before = chrono::system_clock::now();

    SingleStepGreedyLG sslg(&lg);
    sslg.execute();

    long long int totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();


    vector<string> ix = Utils::splitString(fileInstance,'.');
    string instance= ix[ix.size()-2];
    ix = Utils::splitString(instance,'/');
    instance= ix[ix.size()-1];

    Solution * sol = sslg.bestToSolution();

    //This following procedure calculates the number of communities
    unsigned numberCom = 0;
    for (unsigned com=0;com <sol->comunidades.size();com++){
        if (sol->comunidades[com]->qtd > 0){
            numberCom++;
        }
    }

    //storing the data
    string expFile="../data/"+instance+"_"+exp;
    ofstream f(expFile, std::fstream::app);
    f<<instance<<","<<lg.numberOfNodes<<","<<lg.numberOfEdges<<","
           <<exp<<","
           <<sslg.maxIt<<","<<sslg.bestIt<<","
           <<opt<<","<< sslg.bestDensity << "," <<mlg.calculateModularityND_NW(sol)<<","
           <<numberCom<<","
           <<sslg.maxIt<<","<<sslg.bestIt<<","
           <<totalTime<<","<<sslg.totalTime<<","<<sslg.bestTime
           <<"\n";
    f.close();
    return 0;
}

