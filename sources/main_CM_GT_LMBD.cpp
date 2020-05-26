#include <iostream>
#include <chrono>

#include "../../utils/modularitylg.h"
//#include "graph/graph.h"
#include "../../graph/largegraph.h"
#include "../../graph/solution.h"
#include "../../heuristics/constructive/lmbdsinglestepgreedylg.h"
//#include "../../heuristics/constructive/fastgreedylg.h"
//#include "heuristics/constructive/singlestepgreedy.h"
//#include "heuristics/constructive/fastgreedy.h"
//#include "../../heuristics/constructive/singlestepmultilevellg.h"

using namespace std;
using namespace std::chrono;




int main(int argc, char *argv[])
{
    chrono::system_clock::time_point before;
    before = chrono::system_clock::now();

    string fileInstance = argv[1];
    string filepath = "../../instances/"+fileInstance;



    string lambda = argv[2];
    string exp = argv[3];

    srand(time(NULL)+atoi(exp.c_str()));

    LargeGraph lg(filepath);

    const float LAMBDA = atof(lambda.c_str());

    ModularityLG mlg(&lg);

    LMBDSingleStepGreedyLG sslg(&lg, LAMBDA);
    sslg.execute();

    long long int totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

    //construindo a solucao
    Solution *sol = sslg.bestToSolution();



    vector<string> ix = Utils::splitString(fileInstance,'.');
    string instance= ix[ix.size()-2];
    ix = Utils::splitString(instance,'/');
    instance= ix[ix.size()-1];


    //This following procedure calculates the number of communities
    unsigned numberCom = 0;
    for (unsigned com=0;com <sol->comunidades.size();com++){
        if (sol->comunidades[com]->qtd > 0){
            numberCom++;
        }
    }
    string comm = sol->serialize();

//    cout<<"\n"<<sol->serialize();
//cout<<"\n"<<sslg.bestDensity;
//cout<<"\n"<<mlg.calculateDensity(sol,LAMBDA) <<"\n";

    //storing the data
    string expFile="../data/GT_lmbd.csv";
    ofstream f(expFile, std::fstream::app);
    f<<"CM;"<<instance<<";"<<lg.numberOfNodes<<";"<<lg.numberOfEdges<<";"
           <<exp<<";"
           <<sslg.maxIt<<";"<<sslg.bestIt<<";"
           <<lambda<<";"<< sslg.bestDensity << ";" <<mlg.calculateModularityND_NW(sol)<<";"
           <<numberCom<<";"
           <<sslg.maxIt<<";"<<sslg.bestIt<<";"
           <<totalTime<<";"<<sslg.totalTime<<";"<<sslg.bestTime
           <<";"<<comm
           <<"\n";
    f.close();


    return 0;
}




