#include <iostream>
#include <chrono>

#include "../../utils/modularitylg.h"
#include "../../graph/largegraph.h"
#include "../../graph/solution.h"
#include "../../heuristics/hybrid/louvain/louvain.h"
using namespace std;
using namespace std::chrono;

#define TYPE_PRI_NOW TYPE_PRI


int main(int argc, char *argv[])
{
    chrono::system_clock::time_point before;
    before = chrono::system_clock::now();
    srand(time(NULL));
    string fileInstance = argv[1];
    string filepath = "../../instances/"+fileInstance;

    string opt = argv[2];
    string exp = argv[3];

    LargeGraph lg(filepath);


    long long int totalTime;

    //construindo a solucao
    ModularityLG mlg(&lg);


    //before = chrono::system_clock::now();
    Louvain hlv(&lg,TYPE_PRI_NOW);
    hlv.execute();

    totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

    vector<string> ix = Utils::splitString(fileInstance,'.');
    string instance= ix[ix.size()-2];
    ix = Utils::splitString(instance,'/');
    instance= ix[ix.size()-1];


    //This following procedure calculates the number of communities
    Solution * sol= hlv.bestToSolution();
    unsigned numberCom = 0;
    for (unsigned com=0;com <sol->comunidades.size();com++){
        if (sol->comunidades[com]->qtd > 0){
            numberCom++;
        }
    }

//    cout<<"\n"<<sol->serialize();
//    cout<<"\n"<<mlg.calculateDensity(sol);

    //storing the data
    string expFile="../data/lvmod"+instance+"_"+exp;
    ofstream f(expFile, std::fstream::app);
    f<<instance<<","<<lg.numberOfNodes<<","<<lg.numberOfEdges<<","
           <<exp<<","
           <<hlv.it<<","<<hlv.it<<","
           <<opt<<","<< mlg.calculateDensity(sol) << "," <<hlv.bestDensity<<","
           <<numberCom<<","
           <<hlv.it<<","<<hlv.it<<","
           <<totalTime<<","<<hlv.totalTime<<","<<hlv.totalTime
           <<"\n";
    f.close();


    return 0;
}

