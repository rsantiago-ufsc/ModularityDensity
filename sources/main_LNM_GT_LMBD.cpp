#include <iostream>
#include <chrono>

//define the solution construction
#define GROUND_TRUTH

#include "../../utils/modularitylg.h"
//#include "graph/graph.h"
#include "../../graph/largegraph.h"
#include "../../graph/solution.h"
#include "../../heuristics/constructive/lmbdsinglestepgreedylg.h"
#include "../../heuristics/constructive/lmbdfastgreedylg.h"
//#include "heuristics/constructive/lmbdsinglestepgreedy.h"
//#include "heuristics/constructive/lmbdfastgreedy.h"
//#include "../../heuristics/constructive/singlestepmultilevellg.h"

using namespace std;
using namespace std::chrono;




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



    ModularityLG mlg(&lg);

    float lambda = atof(LAMBDA.c_str());


    LMBDFastGreedyLG fglg(&lg, lambda);


//    cout<<"\n"<<solSS->serialize();
//    cout<<"\n"<<mlg.calculateDensity(solSS);

    fglg.execute();

    long double totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

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
/*
cout<<"\n"<<fglg.bestSolution;
cout<<"\n"<<fglg.bestDensity;//sol->serialize();
cout<<"\n"<<mlg.calculateDensity(new Solution(&lg, fglg.bestSolution),lambda);
cout<<"\n"<<mlg.calculateDensity(sol,lambda);
*/
    //storing the data
    string expFile="../data/GT_lmbd.csv";
    ofstream f(expFile, std::fstream::app);
    f<<"LNM;"<<instance<<";"<<lg.numberOfNodes<<";"<<lg.numberOfEdges<<";"
           <<exp<<";"
           <<fglg.maxIt<<";"<<fglg.bestIt<<";"
           <<lambda<<";"<< fglg.bestDensity << ";" <<mlg.calculateModularityND_NW(sol)<<";"
           <<numberCom<<";"
           <<fglg.maxIt<<";"<<fglg.bestIt<<";"
           <<totalTime<<";"<<fglg.totalTime<<";"<<fglg.bestTime
           <<";"<<fglg.bestSolution
           <<"\n";
    f.close();


    return 0;
}


