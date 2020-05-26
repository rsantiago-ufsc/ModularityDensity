
//MLSSLVLL: Multi Level Single Step LouVain Last Level
#include <iostream>
#include <chrono>

#define GROUND_TRUTH

#include "../../utils/modularitylg.h"
#include "../../graph/largegraph.h"
#include "../../graph/solution.h"

#include "../../heuristics/hybrid/multilevelsslv/multilevelsslv.h"
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


    long long int totalTime;


    MultiLevelSSLV mlsslv (&lg);
    mlsslv.execute();

    totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

    vector<string> ix = Utils::splitString(fileInstance,'.');
    string instance= ix[ix.size()-2];
    ix = Utils::splitString(instance,'/');
    instance= ix[ix.size()-1];



    //storing the data
    string expFile="../data/GT.csv";
    ofstream f(expFile, std::fstream::app);
    f<<"MD;"<<instance<<";"<<lg.numberOfNodes<<";"<<lg.numberOfEdges<<";"
           <<exp<<";"
           <<mlsslv.it<<";"<<mlsslv.bestIt<<";"
           <<opt<<";"<< mlsslv.bestDensity << ";" <<mlsslv.bestMod<<";"
           <<mlsslv.bestCommSize<<";"
           <<mlsslv.it<<";"<<mlsslv.bestIt<<";"
           <<totalTime<<";"<<mlsslv.totalTime<<";"<<mlsslv.bestTime
           <<";"<<mlsslv.bestCommStr
           <<"\n";
    f.close();


    return 0;
}



