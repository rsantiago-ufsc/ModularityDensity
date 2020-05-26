#define GROUND_TRUTH
//MLSSLVLL: Multi Level Single Step LouVain Last Level
#include <iostream>
#include <chrono>



#include "../../utils/modularitylg.h"
#include "../../graph/largegraph.h"
#include "../../graph/solution.h"

#include "../../heuristics/hybrid/LAMBDAmultilevelsslvmod/lmbdmultilevelsslvmod.h"
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

    float lambda = atof(LAMBDA.c_str());

    long long int totalTime;


    LMBDMultiLevelSSLVMod mlsslv (&lg, lambda );
    mlsslv.execute();

    totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

    vector<string> ix = Utils::splitString(fileInstance,'.');
    string instance= ix[ix.size()-2];
    ix = Utils::splitString(instance,'/');
    instance= ix[ix.size()-1];

//Solution * sol = new Solution(&lg, mlsslv.bestCommStr);
//ModularityLG mlg(&lg);
//cout<<"\n"<<mlsslv.bestDensity;
//cout<<"\n"<<mlg.calculateDensity(sol, lambda);

    //storing the data
    string expFile="../data/GT_lmbd.csv";
    ofstream f(expFile, std::fstream::app);
    f<<"MDM;"<<instance<<";"<<lg.numberOfNodes<<";"<<lg.numberOfEdges<<";"
           <<exp<<";"
           <<mlsslv.it<<";"<<mlsslv.bestIt<<";"
           <<lambda<<";"<< mlsslv.bestDensity << ";" <<mlsslv.bestMod<<";"
           <<mlsslv.bestCommSize<<";"
           <<mlsslv.it<<";"<<mlsslv.bestIt<<";"
           <<totalTime<<";"<<mlsslv.totalTime<<";"<<mlsslv.bestTime
           <<";"<<mlsslv.bestCommStr
           <<"\n";
    f.close();


    return 0;
}




