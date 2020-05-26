#include <iostream>
#include <chrono>





#include "../../utils/modularitylg.h"
#include "../../graph/largegraph.h"
#include "../../graph/solution.h"
#include "../../heuristics/Qmodularity/louvainInternalDegree/louvainID.h"
using namespace std;
using namespace std::chrono;


int main(int argc, char *argv[])
{
    chrono::system_clock::time_point before;
    before = chrono::system_clock::now();
    srand(time(NULL));
    string fileInstance = argv[1];
    string filepath = "../../instances/"+fileInstance;

    string opt = argv[2];
    string exp = argv[3];
    bool MODE_INTERNAL_DEGREE = false;
    string op4 = "NOT";
    if (argc>4){
        string par4 = argv[4];
        if (par4 == "MODE_INTERNAL_DEGREE"){
            MODE_INTERNAL_DEGREE = true;
            op4 = "YES";
        }
    }

    LargeGraph lg(filepath);


    long long int totalTime;

    //construindo a solucao
    ModularityLG mlg(&lg);


    //before = chrono::system_clock::now();
    LouvainID hlv(&lg, MODE_INTERNAL_DEGREE);
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


    //storing the data
    string expFile=".GT_ID.csv";
    ofstream f(expFile, std::fstream::app);
    f<<"ImpLouvain;"<<instance<<";"<<lg.numberOfNodes<<";"<<lg.numberOfEdges<<";"
    //cout<<"Louvain;"<<instance<<";"<<lg.numberOfNodes<<";"<<lg.numberOfEdges<<";"
           <<exp<<";"
           <<hlv.it<<";"<<hlv.it<<";"
           <<opt<<";"<< mlg.calculateDensity(sol) << ";" <<hlv.bestDensity<<";"
           <<numberCom<<";"
           <<hlv.it<<";"<<hlv.it<<";"
           <<totalTime<<";"<<hlv.totalTime<<";"<<hlv.totalTime
           <<";"<<hlv.bestCommunityStr
           <<";"<<op4
           <<"\n";
    f.close();


    return 0;
}

