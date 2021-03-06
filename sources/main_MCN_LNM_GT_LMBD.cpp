//define the solution construction
#define GROUND_TRUTH


#include <iostream>
#include <chrono>


#include "../../utils/modularitylg.h"
//#include "graph/graph.h"
#include "../../graph/largegraph.h"
#include "../../graph/solution.h"
//#include "../../heuristics/constructive/singlestepgreedylg.h"
//#include "../../heuristics/constructive/fastgreedylg.h"
//#include "heuristics/constructive/singlestepgreedy.h"
//#include "heuristics/constructive/fastgreedy.h"
//#include "../../heuristics/constructive/singlestepmultilevellg.h"
#include "../../heuristics/hybrid/LAMBDAmcnlnm/lmbdmcnlnm.h"
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

    float lambda = atof(LAMBDA.c_str());

    srand(time(NULL)+atoi(exp.c_str()));

    LargeGraph lg(filepath);

    string comm = "";
    long long int totalTime;

    ModularityLG mlg(&lg);


    LMBDMcnLnm mcnlnm(&lg, lambda);
    mcnlnm.execute();

    totalTime += chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

    vector<string> ix = Utils::splitString(fileInstance,'.');
    string instance= ix[ix.size()-2];
    ix = Utils::splitString(instance,'/');
    instance= ix[ix.size()-1];




    //This following procedure calculates the number of communities
    long double mod = mcnlnm.hlv->bestMod;
    unsigned numberCom = 0;
    if (mcnlnm.FastGreedyImproved == true){
        Solution * sol= mcnlnm.fglg->toSolution();
        for (unsigned com=0;com <sol->comunidades.size();com++){
            if (sol->comunidades[com]->qtd > 0){
                numberCom++;
            }
        }
        //cout<<"Dens."<< mlg.calculateDensity(sol)<<" - \n";
        comm = sol->serialize();
        mcnlnm.bestDensity = mlg.calculateDensity(sol, lambda);
        //cout<<sol->serialize()<<"\n";
        mod = mlg.calculateModularityND_NW(sol);
    }else{
        comm=mcnlnm.hlv->bestCommunityStr;
        Solution * sol = new Solution(&lg,"["+comm+"]");
        //cout<<sol->serialize()<<"\n";
        //cout<<"Dens."<< mlg.calculateDensity(sol)<<" + \n";;
        numberCom = mcnlnm.hlv->bestCommSize;
    }

/*Solution * sol = new Solution(&lg,"["+comm+"]");
cout<<"\n"<<mcnlnm.bestDensity;//sol->serialize();
cout<<"\n"<<mlg.calculateDensity(sol,lambda);
cout<<"\n";*/



    //storing the data
    string expFile="../data/GT_lmbd.csv";
    ofstream f(expFile, std::fstream::app);
    f<<"MCN+LNM;"<<instance<<";"<<lg.numberOfNodes<<";"<<lg.numberOfEdges<<";"
           <<exp<<";"
           <<mcnlnm.hlv->it+mcnlnm.fglg->maxIt<<";"<<mcnlnm.hlv->it+mcnlnm.fglg->maxIt<<";"
           <<lambda<<";"<< mcnlnm.bestDensity << ";" <<mod<<";"
           <<numberCom<<";"
           <<mcnlnm.hlv->it+mcnlnm.fglg->maxIt<<";"<<mcnlnm.hlv->it+mcnlnm.fglg->maxIt<<";"
           <<totalTime<<";"<<totalTime <<";"<<totalTime
           <<";"<<comm
           <<"\n";
    f.close();


    return 0;
}

