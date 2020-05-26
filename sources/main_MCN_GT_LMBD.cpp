#include <iostream>
#include <chrono>

//define the solution construction
#define GROUND_TRUTH


#include "../../utils/modularitylg.h"
#include "../../graph/largegraph.h"
#include "../../graph/solution.h"
#include "../../heuristics/hybrid/LAMBDAlouvain/lmbdlouvain.h"
#include "../../heuristics/constructive/lmbdsinglestepgreedylg.h"
#include "../../heuristics/constructive/lmbdfastgreedylg.h"
using namespace std;
using namespace std::chrono;


#define TYPE_PRI_NOW TYPE_PRI

int main(int argc, char *argv[])
{
    chrono::system_clock::time_point before;
    before = chrono::system_clock::now();

    string fileInstance = argv[1];
    string filepath = "../../instances/"+fileInstance;

    string LAMBDA = argv[2];
    string exp = argv[3];

    srand(time(NULL)+atoi(exp.c_str()));

    LargeGraph lg(filepath);




    //construindo a solucao
    ModularityLG mlg(&lg);


    long long int totalTime;

    float lambda = atof(LAMBDA.c_str());

    totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

    //construindo a solucao
    //ModularityLG mlg(&lg);



    before = chrono::system_clock::now();
    LMBDLouvain hlv(&lg,TYPE_PRI_NOW);
    hlv.lambda = lambda;
    hlv.execute();


    totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

    vector<string> ix = Utils::splitString(fileInstance,'.');
    string instance= ix[ix.size()-2];
    ix = Utils::splitString(instance,'/');
    instance= ix[ix.size()-1];

/*Solution * sol = new Solution(&lg,"["+hlv.bestCommunityStr+"]");
cout<<"\n"<<hlv.bestDensity;//sol->serialize();
cout<<"\n"<<mlg.calculateDensity(sol,lambda);
cout<<"\n";*/

    //storing the data
    string expFile="../data/GT_lmbd.csv";
    ofstream f(expFile, std::fstream::app);
    f<<"MCN;"<<instance<<";"<<lg.numberOfNodes<<";"<<lg.numberOfEdges<<";"
           <<exp<<";"
           <<hlv.it<<";"<<hlv.it<<";"
           <<lambda<<";"<< hlv.bestDensity << ";" <<hlv.bestMod<<";"
           <<hlv.bestNumberOfCommunities<<";"
           <<hlv.it<<";"<<hlv.it<<";"
           <<totalTime<<";"<<hlv.totalTime<<";"<<hlv.totalTime
           <<";"<<hlv.bestCommunityStr
           <<"\n";
    f.close();


    return 0;
}



