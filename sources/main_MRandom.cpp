#include <iostream>
#include <chrono>



#include "../../utils/modularitylg.h"
#include "../../graph/largegraph.h"
#include "../../graph/solution.h"

#include "../../heuristics/localsearch/monorandomsearch.h"
using namespace std;
using namespace std::chrono;


int main(int argc, char *argv[])
{
    srand(time(NULL));
    chrono::system_clock::time_point before;
    before = chrono::system_clock::now();

    string fileInstance = argv[1];
    string filepath = "../../instances/"+fileInstance;
    LargeGraph lg(filepath);


    string parameter = argv[2];
//    string exp = argv[3];

//    srand(time(NULL)+atoi(exp.c_str()));

//    LargeGraph lg(filepath);

//    float lambda = atof(LAMBDA.c_str());

//    long long int totalTime;

    ModularityLG modu(&lg);
//cout<<lg.numberOfEdges<<"aaa\n";

    MonoRandomSearch ts(&lg, &modu);
//cout<<lg.numberOfEdges<<"bbb\n";
    Solution sol(&lg);
    for (unsigned v = 0; v<lg.numberOfNodes; v++){
        sol.inserirVertice(v,v);
    }
    sol.modularidade = modu.calculateDensity(&sol,0.5);
    //cout<<"Antes: Modularity: "<<modu.calculateDensity(&sol)<<endl;
    Solution * res = ts.execute(&sol,atof(parameter.c_str()),0.0,10000.0,1000);


    //cout<<"Depois: Calculated Modularity: "<<ts.best->modularidade<<endl;
    //cout<<"Depois: Modularity: "<<modu.calculateDensity(res)<<endl;

    string result = "MonoRandom;"+parameter+";"
            +fileInstance+";"
            +to_string(lg.numberOfNodes)+";"
            +to_string(lg.numberOfEdges)+";"
            +to_string(sol.modularidade)+";"
            +to_string(ts.best->modularidade)+";"
            +to_string(ts.it)+";"
            +to_string(chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count())+""
            +"\n";
    //cout<<result;
    string expFile="../data/lsMD_results.csv";
    ofstream f(expFile, std::fstream::app);
    f<<result;
    f.close();



    return 0;
}




