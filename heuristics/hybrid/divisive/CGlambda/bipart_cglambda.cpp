#include <unordered_set>

#include "bipart_cglambda.h"

Bipart_CGLambda::Bipart_CGLambda(LargeGraph * lg, double alpha){
    this->graph = lg;
    this->alpha = alpha;

    this->prepareWij(Wij);
    this->prepareDi (Di);
}

//executes the column generator
void Bipart_CGLambda::execute(vector<unsigned> &nodes){
    this->firstHalf.clear();
    this->secondHalf.clear();

    //constructing starting solution
    this->maxDensity = 0;

    for (unsigned p=2; p<=nodes.size()-2;p++){
        IloEnv env;
        env.setOut(env.getNullStream() );
        try {
            IloInt  i, j, t, e, m=graph->numberOfEdges;

            /// MASTER PROBLEM ///
            IloModel masterOpt (env);

            IloNum P = p;//parameter//##

            //edges from the nodes
            vector< pair<unsigned, unsigned> >edges;
            for (i=0;i<nodes.size();i++){
                for (j=i+1;j<nodes.size();j++){
                    if (this->graph->getAdj(nodes[i], nodes[j]) > 0.9){
                        edges.push_back(pair<unsigned, unsigned>(i,j));
                    }
                }
            }


            //creating vars from model
            //... lambda vars
            IloNumVarArray L(env, edges.size()*3, 0, IloInfinity, ILOFLOAT);
            IloNumVarArray Y(env, nodes.size(), 0, 1, ILOINT);

            //define objective
            IloExpr masterOFExp(env);
            IloNum DIVISOR = p*(nodes.size() - P);
            //starts with Z_1 because Z_0 corresponds to an artificial cluster
            for (unsigned e=0; e<edges.size(); e++){
                masterOFExp +=  (4.0*nodes.size())/DIVISOR*L[edges.size()*2+e];
            }
            for (unsigned i=0; i<nodes.size(); i++){
                masterOFExp +=  (2.0*P  - nodes.size() )/DIVISOR*Y[i]*graph->degreeOfNode[nodes[i]];
            }
            for (unsigned e=0; e<edges.size(); e++){
                i = edges[e].first;
                j = edges[e].second;
                masterOFExp +=  (-4.0*P)/DIVISOR*(Y[i] + Y[j]);
            }
            for (unsigned i=0; i<nodes.size(); i++){
                masterOFExp +=  (- P  )/DIVISOR*graph->degreeOfNode[nodes[i]];
            }
            masterOFExp += 4.0*P*edges.size()/DIVISOR;
            IloObjective masterObj = IloAdd(masterOpt, IloMaximize(env, masterOFExp) );
            masterOFExp.end();


            //constraints


            //define master constraints
            IloRangeArray masterConstraints(env);
            IloExpr masterConsExp(env);
            for (i=0;i<nodes.size();i++){

                masterConsExp += Y[i];
            }
            masterConsExp -= P;
            IloRange rng (env, 0.0, masterConsExp ,0.0);
            masterConsExp.end();
            masterConstraints.add(rng);

            for (e=0;e<edges.size();e++){
                IloExpr masterConsExp(env);
                masterConsExp +=  L[edges.size()+e]+ L[edges.size()*2+e] - Y[edges[e].first];
                IloRange rng (env, 0.0, masterConsExp ,0.0);
                masterConsExp.end();
                masterConstraints.add(rng);
            }
            for (e=0;e<edges.size();e++){
                IloExpr masterConsExp(env);
                masterConsExp +=  L[edges.size()*0+e]+ L[edges.size()*2+e] - Y[edges[e].second];
                IloRange rng (env, 0.0, masterConsExp ,0.0);
                masterConsExp.end();
                masterConstraints.add(rng);
            }
            for (e=0;e<edges.size();e++){
                IloExpr masterConsExp(env);
                masterConsExp +=  L[edges.size()*0+e]+ L[edges.size()*1+e] +L[edges.size()*2+e];
                IloRange rng (env, 0.0, masterConsExp, 1.0);
                masterConsExp.end();
                masterConstraints.add(rng);
            }
            IloAdd(masterOpt, masterConstraints);

            IloCplex masterCplex(masterOpt);
            masterCplex.setParam(IloCplex::Param::Threads, 1);
            masterCplex.setParam(IloCplex::TiLim, 36000);

            masterCplex.solve();


            IloNum density =masterCplex.getObjValue();
            if (density > this->maxDensity){
                this-> maxDensity = density;
                //storing the two halves
                this->firstHalf.clear();
                this->secondHalf.clear();
                IloNumArray Ys(env, nodes.size());
                masterCplex.getValues(Ys, Y);
                for (j=0;j<nodes.size();j++){
                    if (Ys[j] >0.9){
                        this->firstHalf.push_back(nodes[j]);
                    }else{
                        this->secondHalf.push_back(nodes[j]);
                    }
                }
            }






            env.end();
        }
        catch (IloAlgorithm::CannotExtractException& e)
        {   std::cerr <<
            "CannoExtractException: " << e << std::endl;
            IloExtractableArray failed = e.getExtractables();
            for (IloInt i = 0; i < failed.getSize(); ++i)
                std::cerr << "\t" << failed[i] << std::endl;
                // Handle exception ...
        }
        catch (IloException& ex) {
           cerr << "Error: " << ex << endl;
        }
        catch (...) {
           cerr << "Error" << endl;
        }
    }

    unsigned j;
    cout<<"\n[";
    for (j=0;j<this->firstHalf.size();j++){
        cout<<", "<<this->firstHalf[j];
    }
    cout<<"][";
    for (j=0;j<this->secondHalf.size();j++){
        cout<<", "<<this->secondHalf[j];
    }
    cout<<"]";


}


void Bipart_CGLambda::coefOF(vector<unsigned> &cluster, IloNum &coef, IloNumArray &A){


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
        totalDegree += graph->degreeOfNode[cluster[u]];
        for (IloInt j=0; j<nodes.size(); j++){
            IloInt v = nodes[j];
            internalEdges += graph->getAdj(cluster[u],cluster[v]);

        }
    }
    if (nodes.size() > 0){
        coef = ( (internalEdges*4.0)/2.0 - totalDegree)/nodes.size();
    }else{
        coef = 0;
    }
//cout<<"\n Parte Cima: "<<(internalEdges*4.0)/2.0 - totalDegree;
//cout<<"\n Parte deba: "<<nodes.size();
//cout<<"\n Coef: "<<coef;

}



//"prepareWij" prepares modularity matrix from model
// --> Time complexity \theta( |V|^2 )
void Bipart_CGLambda::prepareWij(IloNum **&wij){
    //memory allocation
    wij = new IloNum * [graph->numberOfNodes];
    for (IloInt i=0; i < graph->numberOfNodes; i++) {
        wij[i] = new IloNum[graph->numberOfNodes];
        for (IloInt j=0; j < graph->numberOfNodes; j++) {
            wij[i][j] =  graph->getAdj(i,j);
        }
    }
}

void Bipart_CGLambda::prepareDi (IloNum *&di){
    di = new IloNum [ graph->numberOfNodes ];
    for (IloInt i=0; i < graph->numberOfNodes; i++) {
        di[i] = graph->degreeOfNode[i];
    }
}


