#include <iostream>
#include <chrono>

#include "../../utils/modularitylg.h"
//#include "graph/graph.h"
#include "../../graph/largegraph.h"
#include "../../graph/solution.h"
#include "../../heuristics/constructive/singlestepgreedylg.h"
#include "../../heuristics/constructive/fastgreedylg.h"
//#include "heuristics/constructive/singlestepgreedy.h"
//#include "heuristics/constructive/fastgreedy.h"
#include "../../heuristics/constructive/singlestepmultilevellg.h"

using namespace std;
using namespace std::chrono;




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

    SingleStepGreedyLG sslg(&lg);
    sslg.execute();

    long long int totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

    //construindo a solucao
    Solution *solSS = sslg.bestToSolution();

    ModularityLG mlg(&lg);


    before = chrono::system_clock::now();
    FastGreedyLG fglg(solSS);

//    cout<<"\n"<<solSS->serialize();
//    cout<<"\n"<<mlg.calculateDensity(solSS);
    fglg.execute();
    totalTime += chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

    vector<string> ix = Utils::splitString(fileInstance,'.');
    string instance= ix[ix.size()-2];
    ix = Utils::splitString(instance,'/');
    instance= ix[ix.size()-1];


    //This following procedure calculates the number of communities
    Solution * sol= fglg.toSolution();
    unsigned numberCom = 0;
    for (unsigned com=0;com <sol->comunidades.size();com++){
        if (sol->comunidades[com]->qtd > 0){
            numberCom++;
        }
    }

//    cout<<"\n"<<sol->serialize();
//    cout<<"\n"<<mlg.calculateDensity(sol);

    //storing the data
    string expFile="../data/ssfg"+instance+"_"+exp;
    ofstream f(expFile, std::fstream::app);
    f<<instance<<","<<lg.numberOfNodes<<","<<lg.numberOfEdges<<","
           <<exp<<","
           <<fglg.maxIt<<","<<fglg.bestIt<<","
           <<opt<<","<< fglg.bestDensity << "," <<mlg.calculateModularityND_NW(sol)<<","
           <<numberCom<<","
           <<fglg.maxIt<<","<<fglg.bestIt<<","
           <<totalTime<<","<<fglg.totalTime<<","<<fglg.bestTime
           <<"\n";
    f.close();


    return 0;
}

