#include <unordered_set>

#include "bipart_cgii_ils.h"

Bipart_CGII_ILS::Bipart_CGII_ILS(LargeGraph * lg, double alpha){
    this->graph = lg;
    this->alpha = alpha;

    this->prepareWij(Wij);
    this->prepareDi (Di);
}

//executes the column generator
void Bipart_CGII_ILS::execute(vector<unsigned> &nodes){
    this->firstHalf.clear();
    this->secondHalf.clear();
//cout<<"\n --> 15";cout.flush();
    //reading instance from command line arguments
    float alpha = 0.7;

    //constructing starting solution
    this->maxDensity = 0;
    IloEnv env;
    env.setOut(env.getNullStream() );
    try {
        IloInt  i, j, t, m=graph->numberOfEdges;

        /// MASTER PROBLEM ///
        IloModel masterOpt (env);

        //creating the constraint coeficients
        IloArray<IloNumArray> A(env);

        //single cluster containing all 'nodes'
        A.add(IloNumArray(env,nodes.size()));
        for (j=0;j<nodes.size();j++){
            A[0][j]= 1;
        }

        //creating vars from model
        //... Z
        IloNumVarArray Z(env, 2, 0.0, IloInfinity, ILOFLOAT);

        //define objective
        IloExpr masterOFExp(env);
        IloNum coef;
        //starts with Z_1 because Z_0 corresponds to an artificial cluster
        for (t=0; t<A.getSize(); t++){
            coefOF(nodes, coef, A[t]);

            masterOFExp +=  coef * Z[t];
        }
        IloObjective masterObj = IloAdd(masterOpt, IloMaximize(env, masterOFExp) );
        masterOFExp.end();


        //constraints


        //define master constraints
        IloRangeArray masterConstraints(env);
        for (i=0;i<nodes.size();i++){
            IloExpr masterConsExp(env);
            for (t=0;t<A.getSize();t++){
                masterConsExp += A[t][i]*Z[t] ;
            }
            IloRange rng (env, 0.0, masterConsExp ,1.0);
            masterConsExp.end();
            masterConstraints.add(rng);
        }
        IloAdd(masterOpt, masterConstraints);


        IloRangeArray masterConstraintsY(env);
        IloExpr masterConsExpY(env);
        for (t=0;t<A.getSize();t++){
            masterConsExpY += Z[t] ;
        }
        IloRange rngY (env, 1.0, masterConsExpY ,1.0);
        masterConsExpY.end();
        masterConstraintsY.add(rngY);
        IloAdd(masterOpt, masterConstraintsY);


        IloCplex masterCplex(masterOpt);
        masterCplex.setParam(IloCplex::Param::Threads, 1);
        masterCplex.setParam(IloCplex::TiLim, 36000);


        ///heuristic for auxiliar problem
        Bipart_ILS bpils(nodes, graph, Wij, alpha);

        ///auxiliar problem
        //parameters and results (price and new column)
        IloNumArray lambda(env, nodes.size()); //dual variables from model
        IloNumArray newCommunity(env, nodes.size());

        unordered_set<unsigned> nodeFinder;
        for(i=0;i<nodes.size();i++){
            nodeFinder.insert(nodes[i]);
        }

        vector< pair<unsigned, unsigned> >edges;
        for (i=0;i<nodes.size();i++){
            for (j=i+1;j<nodes.size();j++){
                if (this->graph->getAdj(nodes[i], nodes[j]) > 0.9){
                    edges.push_back(pair<unsigned, unsigned>(i,j));
                }
            }
        }
        //define AA_v \in {0,1}^|V|: our a_v of the auxiliar problem
        IloNumVarArray AA(env, nodes.size(), 0,1, ILOINT);
        IloNumVarArray E(env, edges.size(), 0,1, ILOINT);
        IloNumVarArray E2(env, edges.size(), 0,1, ILOINT);
        IloNumVar Ea(env, 0,1, ILOINT);
        IloNumVar Eb(env, 0,1, ILOINT);
        IloNumVar Va(env, 0,IloInfinity, ILOFLOAT);
        IloNumVar Da(env, IloInfinity, ILOFLOAT);

        IloNumVar Wa(env,  0,IloInfinity, ILOFLOAT);
        IloNumVarArray G(env, nodes.size(), 0,IloInfinity, ILOFLOAT);


        IloInt D = 0.0;
        for (i=0;i<nodes.size();i++){
            D+= this->graph->degreeOfNode[nodes[i]];
        }
        IloInt V = nodes.size();


        IloInt it;
        IloModel mauxOpt (env);

        //constant part of the OF (firstHalf)
        IloExpr mauxOFExpConst(env);
        mauxOFExpConst += -4*Ea*V +Da*V ;
        mauxOFExpConst +=  4*Wa - 2*Da*Va;
        mauxOFExpConst += -4*Eb*Va +D*Va;

        mauxOpt.add(Wa <= Va);
        mauxOpt.add(Wa <= Ea*V);

        IloExpr mauxCntrVa(env);
        for (unsigned i=0; i<nodes.size();i++){
            mauxCntrVa += AA[i];
        }
        mauxOpt.add(Va == mauxCntrVa);

        for (unsigned i=0; i<nodes.size();i++){
            mauxOpt.add(G[i] <= Va);
            mauxOpt.add(G[i] <= AA[i]);
        }

        IloExpr mauxCntrEa(env);
        IloExpr mauxCntrEb(env);
        for (unsigned e=0; e<edges.size();e++){
            unsigned i = edges[e].first;
            unsigned j = edges[e].second;
            mauxOpt.add(E[e]<=AA[i]);
            mauxOpt.add(E[e]<=AA[j]);

            mauxOpt.add(E2[e]<=1-AA[i]);
            mauxOpt.add(E2[e]<=1-AA[j]);
            mauxCntrEa += E[e];
            mauxCntrEb+= E2[e];
        }
        mauxOpt.add(Ea == mauxCntrEa);
        mauxOpt.add(Eb == mauxCntrEb);

        long double solveRMP = 0.0;
        long double solveAux = 0.0;
        long long int firstTime = -1;
        IloNum firstValue = -1;
        long long int secondTime = -1;
        IloNum secondValue = -1;
        maxDensity = -IloInfinity;
        IloNumArray * Zs;

        chrono::system_clock::time_point antes;
        antes = chrono::system_clock::now();

//cout<<"\n --> 15 (antes do for)";cout.flush();
        for(it=1;;it++){

            //solve the reduced master problem
            masterCplex.solve();
//masterCplex.exportModel("antes.lp");


            if (maxDensity < masterCplex.getObjValue() ){
                maxDensity = masterCplex.getObjValue();
            }
//cout<<"\n\n ["<<it<<"]: "<<this->maxDensity;

            //store the dual variables
            for (i=0;i<nodes.size();i++){
                lambda[i] = masterCplex.getDual(masterConstraints[i]);
            }
            lambdaY = masterCplex.getDual(masterConstraintsY[0]);
//cout<<"\n --> 164 (antes do bpils.execute)";cout.flush();
            //bpils.execute(lambda, lambdaY);


            long long int totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-antes).count()/1000.0;

//            if (bpils.maxValue < RC_EPS || totalTime >= 36000){
//cout<<"\n --> 182 (tentanto EXATO)";cout.flush();
            if (true){
                //execute the exact algorithm
                //first part of non-constant OF
                IloExpr mauxOFExp(env);
                for (IloNum x=0;x<nodes.size();x++){
                    mauxOFExp -= Va*G[x]*lambda[x];
                }
                for (IloNum x=0;x<nodes.size();x++){
                    mauxOFExp += V*G[x]*lambda[x];
                }
                /*//second part of non-constant OF
                for (IloNum x=0;x<nodes.size();x++){
                    for (IloNum v=0;v<nodes.size();v++){
                        mauxOFExp -= (1.0-AA[x])*(1.0-AA[v])*lambda[v];
                    }
                }*/
                IloObjective obj = IloAdd(mauxOpt, IloMaximize(env, mauxOFExpConst + mauxOFExp));
                mauxOFExp.end();

                //solve column generator: and get the new column
                IloCplex mauxCplex(mauxOpt);
                mauxCplex.setParam(IloCplex::Param::Threads,1);
                mauxCplex.setParam(IloCplex::TiLim,36000);

//cout<<"\n   (e(x(aux)))...";cout.flush();
chrono::system_clock::time_point before;
before = chrono::system_clock::now();
                mauxCplex.solve();
cout<<"\n\tT(sec) :"<<to_string(chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count()/1000.0);
//cout<<"end!"; cout.flush();

                //copy the new column
                mauxCplex.getValues(newCommunity, AA);
cout<<"\n\tIt["<<it<<"]: "<<mauxCplex.getObjValue();
                //if solution is optimal
                if (mauxCplex.getObjValue() < RC_EPS || totalTime >= 36000){
                    Zs =  new IloNumArray(env, Z.getSize());
                    masterCplex.getValues((*Zs), Z);
                    break;
                }
                //pass AA to master model (inserting a new community in model.)
                //--> firstHalf
                Z.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT));
                A.add(IloNumArray(env,nodes.size() ));
                t=A.getSize()-1;
                for (j=0;j<nodes.size();j++){
                        A[t][j] = newCommunity[j];
                }
                for (i=0;i<nodes.size();i++){
                    masterConstraints[i].setExpr(
                                masterConstraints[i].getExpr()+
                                   A[t][i]*Z[t]
                                );
                }
                masterConstraintsY[0].setExpr(
                            masterConstraints[0].getExpr()+
                            Z[t]
                        );


                //the new contribution of Z_t to objetive function
                IloNum coef, coef2;
                coefOF(nodes, coef, A[t]);
                masterObj.setExpr(masterObj.getExpr()+coef * Z[t]);
                mauxOpt.remove(obj);
//cout<<"\n --> 182 (saindo EXATO)";cout.flush();
                continue;
            }else{
//cout<<"\n --> 264 (entrando EXATO)";cout.flush();
                //bpils.constructSolutionArray(lambda);


                //pass AA to master model (inserting a new community in model.)
                // --> firstHalf
                Z.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT));
                A.add(IloNumArray(env, nodes.size() ));
                t=A.getSize()-1;
                for (j=0;j<nodes.size();j++){
                        A[t][j] = bpils.maxNodes[j];
                }
                for (i=0;i<nodes.size();i++){
                    masterConstraints[i].setExpr(
                                masterConstraints[i].getExpr()+
                                A[t][i]*Z[t]
                                );
                }
                // --> secondHalf
                Z.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT));
                A.add(IloNumArray(env, nodes.size() ));
                for (j=0;j<nodes.size();j++){
                        A[t+1][j] = 1.0 - bpils.maxNodes[j];
                }
                for (i=0;i<nodes.size();i++){
                    masterConstraints[i].setExpr(
                                masterConstraints[i].getExpr()+
                                A[t+1][i]*Z[t+1]
                                );
                }


                //updating constraint \sum{z_t} = 2
                masterZt[0].setExpr(masterZt[0].getExpr()+Z[t]);
                masterZt[0].setExpr(masterZt[0].getExpr()+Z[t+1]);

                //the new contribution of Z_t to objetive function
                IloNum coef, coef2;
                coefOF(nodes, coef, A[t]);
                coefOF(nodes, coef2, A[t+1]);
                masterObj.setExpr(masterObj.getExpr()+coef * Z[t]+coef2 * Z[t+1]);
//cout<<"\n --> 264 (saindo EXATO)";cout.flush();
            }

        }

//cout<<"\n --> 310 (saindo)";cout.flush();
        //storing the two halves
        bool first = true;
        for (i=0;i<Zs->getSize(); i++ ){
            if ((*Zs)[i] >0.9){
                if (first){
                    first = false;
                    for (j=0;j<nodes.size();j++){
                        if (A[i][j] >0.9){
                            this->firstHalf.push_back(nodes[j]);
                        }
                    }
                }else{
                    for (j=0;j<nodes.size();j++){
                        if (A[i][j] >0.9){
                            this->secondHalf.push_back(nodes[j]);
                        }
                    }
                }


            }

        }
/*
cout<<"\n[";
for (j=0;j<this->firstHalf.size();j++){
    cout<<", "<<this->firstHalf[j];
}
cout<<"][";
for (j=0;j<this->secondHalf.size();j++){
    cout<<", "<<this->secondHalf[j];
}
cout<<"]";
*/
//cout<<"\n --> 310 (saindo final)";cout.flush();
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


void Bipart_CGII_ILS::coefOF(vector<unsigned> &cluster, IloNum &coef, IloNumArray &A){

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
    nodes.clear();
    for (IloInt v=0; v<A.getSize(); v++){
        if(A[v] < 0.9){
           nodes.push_back(v);
        }

    }

    internalEdges = 0;
    totalDegree = 0;
    for (IloInt i=0; i<nodes.size(); i++){
        IloInt u = nodes[i];
        totalDegree += graph->degreeOfNode[cluster[u]];
        for (IloInt j=0; j<nodes.size(); j++){
            IloInt v = nodes[j];
            internalEdges += graph->getAdj(cluster[u],cluster[v]);

        }
    }


    if (nodes.size() > 0){
        coef += ( (internalEdges*4.0)/2.0 - totalDegree)/nodes.size();
    }else{
        coef += 0;
    }


}



//"prepareWij" prepares modularity matrix from model
// --> Time complexity \theta( |V|^2 )
void Bipart_CGII_ILS::prepareWij(IloNum **&wij){
    //memory allocation
    wij = new IloNum * [graph->numberOfNodes];
    for (IloInt i=0; i < graph->numberOfNodes; i++) {
        wij[i] = new IloNum[graph->numberOfNodes];
        for (IloInt j=0; j < graph->numberOfNodes; j++) {
            wij[i][j] =  graph->getAdj(i,j);
        }
    }
}

void Bipart_CGII_ILS::prepareDi (IloNum *&di){
    di = new IloNum [ graph->numberOfNodes ];
    for (IloInt i=0; i < graph->numberOfNodes; i++) {
        di[i] = graph->degreeOfNode[i];
    }
}


