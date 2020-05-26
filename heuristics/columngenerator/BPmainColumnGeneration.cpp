#define GROUND_TRUTH

#include <boost/heap/fibonacci_heap.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <cstdlib>
#include <time.h>
#include <ilcplex/ilocplex.h>
#include <ilconcert/iloenv.h>
#include <chrono>
#include <list>





#include "../../graph/largegraph.h"
#include "../../graph/solution.h"
#include "../../utils/modularitylg.h"

//#include "../hybrid/mcnlnm/mcnlnm.h"
#include "../hybrid/LAMBDAmcnlnm/lmbdmcnlnm.h"
#include "minifastgreedybli.h"
using namespace std;
using namespace std::chrono;
using namespace boost::heap;

ILOSTLBEGIN

#define RC_EPS 1.0e-9


//number of columns
long long int ncInitial = 0;
long long int ncAPHeuristic = 0;
long long int ncAPExact = 0;

long long int totalTimeInitial      = 0;
long long int totalTimeLP           = 0;
long long int totalTimeAPExact      = 0;
long long int totalTimeAPHeuristic  = 0;
long double solveRMP = 0.0;
long double solveAux = 0.0;
string experimento = "0";

IloEnv env;

bool PAR_VERBOSE_ITER = false;

float alpha = 0.7;
IloNum *Di;
IloNum **Wij;
string parameter = "";

chrono::system_clock::time_point antes;

#define EPSILON 1E-6
bool DBL_EQL(double a, double b) {
    return fabs(a - b) < EPSILON;
}

static void coefOF(IloNum &coef, LargeGraph &graph, IloNumArray &A, float Plambda){


    vector<IloNum> nodes;

    for (IloInt v=0; v<A.getSize(); v++){
        if(A[v] >= 0.9){
           nodes.push_back(v);
        }

    }


    IloNum internalEdges = 0;
    IloNum totalDegree = 0;
    for (IloInt i=0; i<nodes.size(); i++){
        IloInt u = nodes[i];
        totalDegree += graph.degreeOfNode[u];
        for (IloInt j=0; j<nodes.size(); j++){
            IloInt v = nodes[j];
            internalEdges += graph.getAdj(u,v);

        }
    }
    coef = ( (internalEdges*4.0*Plambda)/2.0 -(2-2*Plambda)*(totalDegree-internalEdges))/nodes.size();
    //coef = ( (internalEdges*4.0)/2.0 - totalDegree)/nodes.size();
    //cout<<"\n Parte Cima: "<<(internalEdges*4.0)/2.0 - totalDegree;
    //cout<<"\n Parte deba: "<<nodes.size();

}



//"prepareWij" prepares modularity matrix from model
// --> Time complexity \theta( |V|^2 )
void prepareWij(LargeGraph &graph, IloNum **&wij){
    //memory allocation
    wij = new IloNum * [graph.numberOfNodes];
    for (IloInt i=0; i < graph.numberOfNodes; i++) {
        wij[i] = new IloNum[graph.numberOfNodes];
        for (IloInt j=0; j < graph.numberOfNodes; j++) {
            wij[i][j] =  graph.getAdj(i,j);
        }
    }
}

void prepareDi (LargeGraph &graph, IloNum *&di){
    di = new IloNum [ graph.numberOfNodes ];
    for (IloInt i=0; i < graph.numberOfNodes; i++) {
        di[i] = graph.degreeOfNode[i];
    }
}

using namespace std;

struct TSolution{
    pair<unsigned, unsigned> fixed;
    long long int parentId;
    long long int id;
    IloNum density;
    bool operator<(TSolution const & rhs) const{
        return density < rhs.density;
    }
};



IloNumArray * solverCGII_ILS(
            LargeGraph                  &graph,
            float                       &Plambda,
            fibonacci_heap<TSolution>   &heap,
            TSolution                   &bpnode,
            vector<TSolution *>         &bpnodes,
            IloNum                      &lowerbound,
            IloNumArray                 *&bestZs,
            IloArray<IloNumArray>       &A,
            IloNumVarArray              &Z,
            IloModel                    &masterOpt,
            IloCplex                    &masterCplex,
            IloModel                    &mauxOpt,
            IloNumVarArray              &AA,
            IloRangeArray               &masterConstraints,
            IloObjective                &masterObj,
            IloExpr                     &mauxOFExpConst
        ){




    //node and solution
    IloConstraintArray nodeConstraints(env);
cout<<"\n\n*******************\n";
    long long int id = bpnode.id;
    cout<<id<<" <- ";
    while (bpnodes[id]->parentId != -1){

        IloConstraint cons = Z[bpnodes[id]->fixed.first] == bpnodes[id]->fixed.second;
        masterOpt.add(cons);
        nodeConstraints.add(cons);
        id = bpnodes[id]->parentId;
        cout<<id<<" <- ";
    }
cout<<"\n*******************\n\n";
    IloInt it, i, t, j;

    //parameters and results (price and new column)
    IloNumArray lambda(env, graph.numberOfNodes); //dual variables from model
    IloNumArray newCommunity(env, graph.numberOfNodes);


    /*long long int firstTime = -1;
    IloNum firstValue = -1;
    long long int secondTime = -1;
    IloNum secondValue = -1;*/
    IloNum maxDensity = -IloInfinity;
    IloNumArray * Zs;

    //masterCplex.exportModel("antes.lp");

    for(it=1;;it++){
        chrono::system_clock::time_point antesRMP;
        antesRMP = chrono::system_clock::now();

        //solve the reduced master problem
        chrono::system_clock::time_point antesLP;
        antesLP = chrono::system_clock::now();
        masterCplex.solve();
        totalTimeLP += chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now()-antesLP).count();

        solveRMP += chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-antesRMP).count()/1000.0;
        //cout<<masterCplex.getObjValue();
        if (maxDensity < masterCplex.getObjValue() ){
            maxDensity = masterCplex.getObjValue();
        }
        //print the data of each iteration
        if (PAR_VERBOSE_ITER){
            //instance,method,dvalue,iteration
            string file = "evolutionByIter.csv";
            string texto = parameter+",CGI+ILS,"+to_string(maxDensity)+","+to_string(it)+","+experimento+"\n";
            ofstream ofs (file, std::fstream::out | std::fstream::app);
            ofs << texto;
            ofs.close();
        }

        /*if (it==1){
            firstTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-antes).count()/1000.0;
            firstValue = masterCplex.getObjValue();
        }
        if (it==2){
            secondTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-antes).count()/1000.0;
            secondValue = masterCplex.getObjValue();
        }*/

        //store the dual variables
        for (i=0;i<graph.numberOfNodes;i++){
            lambda[i] = masterCplex.getDual(masterConstraints[i]);
        }

        chrono::system_clock::time_point antesAux;
        antesAux = chrono::system_clock::now();

        MiniFastGreedyBLI mlcg(&graph, Wij, alpha, Plambda);

        chrono::system_clock::time_point antesAuxH = chrono::system_clock::now();
        mlcg.execute(lambda);
        totalTimeAPHeuristic += chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now()-antesAuxH).count();

        solveAux += chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-antesAux).count()/1000.0;

        long long int totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-antes).count()/1000.0;

        if (mlcg.maxValue < RC_EPS || totalTime >= 36000){
            //execute the exact algorithm
            IloExpr mauxOFExp(env);

            for (IloNum x=0;x<graph.numberOfNodes;x++){
                for (IloNum v=0;v<graph.numberOfNodes;v++){
                    mauxOFExp -= AA[x]*AA[v]*lambda[v];
                }
            }

            IloObjective obj = IloAdd(mauxOpt, IloMaximize(env, mauxOFExpConst + mauxOFExp));

            mauxOFExp.end();

            //solve column generator: and get the new column
            IloCplex mauxCplex(mauxOpt);
            mauxCplex.setParam(IloCplex::Param::Threads,1);
            mauxCplex.setParam(IloCplex::TiLim,36000);
            chrono::system_clock::time_point antesAux;
            antesAux = chrono::system_clock::now();

            chrono::system_clock::time_point antesAuxEx = chrono::system_clock::now();
            mauxCplex.solve();
            totalTimeAPExact += chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now()-antesAuxEx).count();

            solveAux += chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-antesAux).count()/1000.0;

            //copy the new column
            mauxCplex.getValues(newCommunity, AA);

            //if solution is optimal
            long long int totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-antes).count()/1000.0;

            if (mauxCplex.getObjValue() < RC_EPS || totalTime >= 36000){
                Zs =  new IloNumArray(env, AA.getSize());
                masterCplex.getValues((*Zs), Z);
                mauxOpt.remove(obj);
                obj.end();
                mauxCplex.end();
                break;
            }

            //pass AA to master model (inserting a new community in model.)
            Z.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT));
            A.add(IloNumArray(env,graph.numberOfNodes ));
            ncAPExact++;
            t=A.getSize()-1;
            for (j=0;j<graph.numberOfNodes;j++){
                    A[t][j] = newCommunity[j];
            }

            for (i=0;i<graph.numberOfNodes;i++){
                masterConstraints[i].setExpr(
                            masterConstraints[i].getExpr()+
                            A[t][i]*Z[t]
                            );
            }

            //the new contribution of Z_t to objetive function
            IloNum coef;
            coefOF(coef, graph, A[t], Plambda);
            masterObj.setExpr(masterObj.getExpr()+coef * Z[t]);
            mauxOpt.remove(obj);
            obj.end();
            mauxCplex.end();

            continue;
        }else{

            mlcg.constructSolutionArray(lambda);


            //pass AA to master model (inserting a new community in model.)
            Z.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT));
            A.add(IloNumArray(env,graph.numberOfNodes ));
            ncAPHeuristic++;
            t=A.getSize()-1;
            for (j=0;j<graph.numberOfNodes;j++){
                    A[t][j] = mlcg.maxNodes[j];
            }

            for (i=0;i<graph.numberOfNodes;i++){
                masterConstraints[i].setExpr(
                            masterConstraints[i].getExpr()+
                            A[t][i]*Z[t]
                            );
            }

            //the new contribution of Z_t to objetive function
            IloNum coef;
            coefOF(coef, graph, A[t], Plambda);
            masterObj.setExpr(masterObj.getExpr()+coef * Z[t]);

        }


    }

    //removing constraints
    for (IloInt x = 0; x < nodeConstraints.getSize(); x++){
         masterOpt.remove( nodeConstraints[x] );
    }

    //BRANCH
    bool isIntegral = true;
    IloNum bestFrac = 0.0;
    IloInt idFrac   = -1;
    for (IloInt c=0; c<Zs->getSize(); c++){
        if (!DBL_EQL((*Zs)[c], 0.0) && !DBL_EQL((*Zs)[c], 1.0)){
            isIntegral = false;
            IloInt intPart = (*Zs)[c];
            if (bestFrac < (*Zs)[c] - intPart){
                bestFrac = (*Zs)[c];
                idFrac = c;
            }
        }
    }
    if (isIntegral && maxDensity > lowerbound){
        lowerbound = maxDensity;
        if (bestZs != NULL){
            delete bestZs;
        }
        bestZs = Zs;
    }else{
        //branching
        TSolution *left = new TSolution, *right = new TSolution;

        IloNum coef;
        coefOF(coef,graph, A[idFrac], Plambda);
        left->density = maxDensity - coef *(*Zs)[idFrac] + 0.0 *coef;
        right->density = maxDensity - coef *(*Zs)[idFrac] + 1.0 *coef;
        right->parentId = left->parentId = bpnode.id;
        left->id = bpnodes.size();
        right->id = bpnodes.size()+1;
        left->fixed.first = idFrac;
        right->fixed.first = idFrac;
        left->fixed.second = 0;
        right->fixed.second = 1;
        bpnodes.push_back(left);
        bpnodes.push_back(right);
        heap.push(*left);
        heap.push(*right);
    }



}



int main(int argc, char **argv){
    srand (time(NULL));

    //starting chronometer
    chrono::system_clock::time_point tempoTudo;
    tempoTudo = chrono::system_clock::now();

    //reading instance from command line arguments


    if (argc > 1){
        parameter = argv[1];
        //alpha = atof(argv[2]);
    }
//


    float Plambda = atof(argv[2]);
    IloNumArray * bestZs = NULL;


    string filepath = "../../instances/"+parameter;
    LargeGraph graph(filepath);

    //constructing starting solution
    long double maxDensity = 0;
    //constructing starting solution
    Solution * startingComm=NULL;


    ModularityLG mlg(&graph);
    chrono::system_clock::time_point antesInitial;
    antesInitial = chrono::system_clock::now();
    for(unsigned i=0;i<30;i++){
        LMBDMcnLnm mcnlnm(&graph, Plambda);
        mcnlnm.execute();
        unsigned numberCom = 0;
        if (mcnlnm.FastGreedyImproved == true){
            if (startingComm!=NULL){
                delete startingComm;
            }
            Solution * sol= mcnlnm.fglg->toSolution();
            startingComm = sol;
            maxDensity = mlg.calculateDensity(sol, Plambda);
        }else{
            if (startingComm!=NULL){
                delete startingComm;
            }
            string comm=mcnlnm.hlv->bestCommunityStr;
            Solution * sol = new Solution(&graph,"["+comm+"]");
            maxDensity = mlg.calculateDensity(sol, Plambda);
            startingComm = sol;
        }
    }


    totalTimeInitial += chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now()-antesInitial).count();

    env.setOut(env.getNullStream() );
    try {
        IloInt  i, j, t, m=graph.numberOfEdges;


        prepareWij(graph, Wij);

        prepareDi (graph, Di);



        /// MASTER PROBLEM ///
        IloModel masterOpt (env);

        //creating the constraint coeficients
        IloArray<IloNumArray> A(env);
        unsigned commStCount = 0;
        for (unsigned c=0;c<startingComm->comunidades.size(); c++){
            LDE<unsigned> * comm = startingComm->comunidades[c];
            if(comm->getQtd() > 0){
                A.add(IloNumArray(env,graph.numberOfNodes ));
                ncInitial++;
                for (j=0;j<graph.numberOfNodes/2;j++){
                        A[commStCount][j]= 0;
                }
                itemLDE<unsigned> * nav = comm->getInicio();
                while(nav!=NULL){
                    A[commStCount][nav->id]= 1;
                    nav=nav->prox;
                }
                commStCount++;
            }
        }

        //creating vars from model
        //... Z
        IloNumVarArray Z(env, A.getSize(), 0.0, IloInfinity, ILOFLOAT);


        //define objective
        IloExpr masterOFExp(env);
        IloNum coef;
        for (t=0; t<A.getSize(); t++){
            coefOF(coef,graph, A[t], Plambda);
            masterOFExp +=  coef * Z[t];
        }


        IloObjective masterObj = IloAdd(masterOpt, IloMaximize(env, masterOFExp /*-C/(2.0*m)*/) );
        masterOFExp.end();



        //define master constraints
        IloRangeArray  masterConstraints(env);
        for (i=0;i<graph.numberOfNodes;i++){
            IloExpr  masterConsExp(env);
            for (t=0;t<A.getSize();t++){
                masterConsExp += A[t][i]*Z[t] ;

            }
            IloRange rng (env, 1.0, masterConsExp ,1.0);
            masterConsExp.end();
            masterConstraints.add(rng);
        }
        IloAdd(masterOpt ,masterConstraints);

        IloCplex masterCplex(masterOpt);
        masterCplex.setParam(IloCplex::Param::Threads,1);
        masterCplex.setParam(IloCplex::TiLim,36000);


        ///auxiliar problem



        //define AA_v \in {0,1}^|V|: our a_v of the auxiliar problem
        IloNumVarArray AA(env, graph.numberOfNodes, 0,1, ILOINT);//-IloInfinity,IloInfinity, ILOFLOAT);


        antes = chrono::system_clock::now();

        IloInt it;
        IloModel mauxOpt (env);

        //constant part of the OF
        IloExpr mauxOFExpConst(env);
        for (i=0;i<graph.numberOfNodes;i++){
            for (j=i+1;j<graph.numberOfNodes;j++){
                mauxOFExpConst += 4*Wij[i][j]*AA[i]*AA[j];
            }
        }

        for (i=0;i<graph.numberOfNodes;i++){
            mauxOFExpConst += (-2.0+2.0*Plambda)*Di[i]*AA[i]; //
        }



        //defining the starting node of the BP
        TSolution solInicial;
        solInicial.density = -IloInfinity;
        solInicial.parentId = -1;

        vector<TSolution *> bpnodes;
        solInicial.id = bpnodes.size();
        bpnodes.push_back(&solInicial);

        //the queue
        fibonacci_heap<TSolution> heap;
        heap.push(solInicial);

        IloNum lowerbound = -IloInfinity;

        bool terminate = false;

        while (!terminate){
            TSolution data = heap.top();
            heap.pop();
            if (data.density < lowerbound){
                break;
            }


            solverCGII_ILS(    graph,
                               Plambda,
                               heap,
                               data,
                               bpnodes,
                               lowerbound,
                               bestZs,
                               A,
                               Z,
                               masterOpt,
                               masterCplex,
                               mauxOpt,
                               AA,
                               masterConstraints,
                               masterObj,
                               mauxOFExpConst
                           );
            if (heap.empty()){
                break;
            }

        }







    string comms = "";
    for (IloInt i=0;i<bestZs->getSize(); i++){
        string comm = "";
        if ((*bestZs)[i] > 0.9){
            for (IloInt v=0;v<A[i].getSize(); v++){
                if (A[i][v] > 0.9){
                    if (comm == ""){
                        comm = to_string(v);
                    }else{
                        comm += ","+to_string(v);
                    }
                }
            }
            if (comm != ""){
                if (comms  == ""){
                    comms+= "["+comm+"]";
                }else{
                    comms+= ",["+comm+"]";
                }
            }
        }
    }

    Solution *sole = new Solution (&graph, comms);
    cout<<"\n*** True: "<<mlg.calculateDensity(sole, Plambda)<<endl;


//grafo,tipo,todaheuristica,todoprocesso,maxdensity,status,tempo1vez,valor1vez,tempo2vez,valor2vez,tempoRMP,tempoAux,it

    //masterCplex.exportModel("depois.lp");
    //cout<<"MaxDensity: "<<maxDensity;cout.flush();
    string file = "nossoCG_GT.csv";
    //string texto = parameter+",cgFG+BLI_hlsmd("+to_string(alpha)+")+G30,";
    string texto = parameter+";CGI+ILS;";
    texto += to_string(Plambda)+";";
    texto += to_string(chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-antes).count()/1000.0)+";";
    texto += to_string(chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-tempoTudo).count()/1000.0)+";";
    texto += to_string(maxDensity )+";";
    texto += to_string(masterCplex.getStatus())+";";
    /*texto += to_string(firstTime)+";";
    texto += to_string(firstValue)+";";
    texto += to_string(secondTime)+";";
    texto += to_string(secondValue)+";";*/
    texto += to_string(solveRMP)+";";
    texto += to_string(solveAux)+";";
    //texto += to_string(it)+"\n";
    texto += to_string(it)+";";

    texto += to_string(ncInitial)+";";
    texto += to_string(ncAPExact)+";";
    texto += to_string(ncAPHeuristic)+";";

    texto += to_string(totalTimeInitial)+";";
    texto += to_string(totalTimeLP)+";";
    texto += to_string(totalTimeAPExact)+";";
    texto += to_string(totalTimeAPHeuristic)+";";
    texto += comms+"\n";

    ofstream ofs (file, std::fstream::out | std::fstream::app);
    ofs << texto;
    ofs.close();




    env.end();
    //cout<<"\n\n";
    return 0;



}
catch (IloAlgorithm::CannotExtractException& e)
{ std::cerr <<
"CannoExtractException: " << e << std::endl; IloExtractableArray failed = e.getExtractables();

for (IloInt i = 0; i < failed.getSize(); ++i) std::cerr <<
"\t" << failed[i] << std::endl;
// Handle exception ...
}





catch (IloException& ex) {
   cerr << "Error: " << ex << endl;
}
catch (...) {
   cerr << "Error" << endl;
}


    return 0;
}
