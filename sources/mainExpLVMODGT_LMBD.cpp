#include <iostream>
#include <chrono>

//define the solution construction
#define GROUND_TRUTH


#include "../../utils/modularitylg.h"
#include "../../graph/largegraph.h"
#include "../../graph/solution.h"
#include "../../heuristics/hybrid/LAMBDAlouvain/lmbdlouvain.h"
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
    LMBDLouvain hlv(&lg,TYPE_PRI_NOW);
    hlv.lambda = 0.9;
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

    //cout<<"\n"<<sol->serialize();
    cout<<"\n"<<hlv.bestDensity;
    cout<<"\n"<<mlg.calculateDensity(sol,0.9)<<"\n";

    //storing the data
    string expFile="../data/GT.csv";
    //ofstream f(expFile, std::fstream::app);
    cout<<"Louvain;"<<instance<<";"<<lg.numberOfNodes<<";"<<lg.numberOfEdges<<";"
           <<exp<<";"
           <<hlv.it<<";"<<hlv.it<<";"
           <<opt<<";"<< mlg.calculateDensity(sol) << ";" <<hlv.bestDensity<<";"
           <<numberCom<<";"
           <<hlv.it<<";"<<hlv.it<<";"
           <<totalTime<<";"<<hlv.totalTime<<";"<<hlv.totalTime
           <<";"<<hlv.bestCommunityStr
           <<"\n";
    //f.close();


    return 0;
}

