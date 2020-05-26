
#include <iostream>
#include <chrono>

#include "../../../../utils/modularitylg.h"
#include "../../../../graph/largegraph.h"
#include "../../../../graph/solution.h"

#include "divisive_cgi_ils.h"
using namespace std;
using namespace std::chrono;




int main(int argc, char *argv[])
{
    chrono::system_clock::time_point before;
    before = chrono::system_clock::now();

    string fileInstance = argv[1];
    string filepath = "../../../../instances/"+fileInstance;




    string opt = argv[2];

    string exp = argv[3];

    //srand(time(NULL)+atoi(exp.c_str()));
    srand(atoi(exp.c_str()));



    LargeGraph lg(filepath);

    Divisive_CGII_ILS divCGIILS(&lg,0.7);


    divCGIILS.execute();

    cout<<"MD: "<<divCGIILS.maxDensity;

    cout<<"\nTIME :"<<to_string(chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count()/1000.0);
    return 0;



    /*//storing the data
    string expFile="../data/mcnlnm"+instance+"_"+exp;
    ofstream f(expFile, std::fstream::app);
    f<<instance<<","<<lg.numberOfNodes<<","<<lg.numberOfEdges<<","
           <<exp<<","
           <<mcnlnm.hlv->it+mcnlnm.fglg->maxIt<<","<<mcnlnm.hlv->it+mcnlnm.fglg->maxIt<<","
           <<opt<<","<< mcnlnm.bestDensity << "," <<mod<<","
           <<numberCom<<","
           <<mcnlnm.hlv->it+mcnlnm.fglg->maxIt<<","<<mcnlnm.hlv->it+mcnlnm.fglg->maxIt<<","
           <<totalTime<<","<<totalTime <<","<<totalTime
           <<"\n";
    f.close();*/


    return 0;
}


