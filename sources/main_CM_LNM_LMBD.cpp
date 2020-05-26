#include <iostream>
#include <chrono>

#define GROUND_TRUTH

#include "../../utils/modularitylg.h"
//#include "graph/graph.h"
#include "../../graph/largegraph.h"
#include "../../graph/solution.h"
#include "../../heuristics/constructive/lmbdsinglestepgreedylg.h"

#undef GROUND_TRUTH


#include "../../heuristics/constructive/lmbdfastgreedylg.h"
//#include "heuristics/constructive/singlestepgreedy.h"
//#include "heuristics/constructive/fastgreedy.h"
//#include "../../heuristics/constructive/lmbdsinglestepgreedylg.h"

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

    float LAMBDA = atof(lambda.c_str());

    LargeGraph lg(filepath);
    ModularityLG mlg(&lg);
    /*Solution s1(&lg, "[4,5,20,14,18,6,13,2,25,12,21,3,0,23,22,19,7,10,15,1,24,11,17,28,16,8,27,9,26,30,29]");
    Solution s2(&lg, "[4,5,20,14,18,6,13,2,25,12,21,3,0,23,22,19,7,10,15,1,24,11,17,28,16,8,27,9,26,30,29]");
    cout<<"\n"<<mlg.calculateDensity(&s1, LAMBDA)<<"\n";
    cout<<"\n"<<mlg.calculateDensity(&s2, LAMBDA)<<"\n";*/





    LMBDSingleStepGreedyLG sslg(&lg, LAMBDA);
    sslg.execute();
    Solution *sol2 = new Solution(&lg, sslg.bestCommunity);
 //cout<<"\n"<<sslg.bestCommunity;
    //cout<<"\nSSG:"<<sslg.bestDensity;
    //cout<<"\nSSG:"<<mlg.calculateDensity(sol2, LAMBDA)<<"\n";


    //construindo a solucao
    Solution *solSS = sol2;




    before = chrono::system_clock::now();
    LMBDFastGreedyLG fglg(solSS, LAMBDA);
    fglg.execute();



    long double totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

    vector<string> ix = Utils::splitString(fileInstance,'.');
    string instance= ix[ix.size()-2];
    ix = Utils::splitString(instance,'/');
    instance= ix[ix.size()-1];


    //This following procedure calculates the number of communities
    Solution * sol= new Solution(&lg, fglg.bestSolution);
    //Solution * sol= fglg.toSolution();
    /*unsigned numberCom = 0;
    for (unsigned com=0;com < sol->comunidades.size();com++){
        if (sol->comunidades[com]->qtd > 0){
            numberCom++;
        }
    }*/
    //string comm = fglg.bestSolution;//sol->serialize();

//cout<<"\n"<<numberCom;
//cout<<"\n"<<fglg.bestNComm;
//cout<<"\n"<<fglg.bestDensity;
//cout<<"\n"<<mlg.calculateDensity(sol, LAMBDA)<<"\n";

    //storing the data
    string expFile="../data/data_lmbd.csv";
    ofstream f(expFile, std::fstream::app);
    f<<"CM+LNM;"<<instance<<";"<<lg.numberOfNodes<<";"<<lg.numberOfEdges<<";"
           <<exp<<";"
           <<fglg.maxIt<<";"<<fglg.bestIt<<";"
           <<lambda<<";"<< fglg.bestDensity << ";"
           <<mlg.calculateModularityND_NW(sol)<<";"
           <<fglg.bestNComm<<";"
           <<fglg.maxIt<<";"<<fglg.bestIt<<";"
           <<totalTime<<";"<<fglg.totalTime<<";"<<fglg.bestTime
           //<<";"<<comm
           <<"\n";
    f.close();


    return 0;
}



