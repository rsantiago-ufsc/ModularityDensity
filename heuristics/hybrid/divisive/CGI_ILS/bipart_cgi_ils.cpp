#include "bipart_cgi_ils.h"

#include<unordered_set>
#include<unordered_map>


Bipart_CGII_ILS::Bipart_CGII_ILS(LargeGraph * lg, double alpha, list<Solution *> * solutions){
    this->graph = lg;
    this->alpha = alpha;

    this->prepareWij(Wij);
    this->prepareDi (Di);

    this->solutions=solutions;
}



//executes the column generator
void Bipart_CGII_ILS::execute(vector<unsigned> &nodes){
    this->firstHalf.clear();
    this->secondHalf.clear();
//cout<<"\n --> 15";cout.flush();
    //reading instance from command line arguments
    float alpha = 0.2;










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


        //empty cluster to assure \sum{z_t} = 2 (constraint)
        A.add(IloNumArray(env,nodes.size()));
        for (j=0;j<nodes.size();j++){
            A[0][j]= 0;
        }
        //single cluster containing all 'nodes'
        A.add(IloNumArray(env,nodes.size()));
        for (j=0;j<nodes.size();j++){
            A[1][j]= 1;
        }

        //new initial solutions
        /*unordered_set<unsigned> setNodes;
        setNodes.insert(nodes.begin(), nodes.end());
        unordered_map<unsigned, unsigned> mapNodes;
        for (i=0; i<nodes.size();i++){
            mapNodes[nodes[i]] = i;
        }
        for (Solution * sol: *this->solutions){
            for (LDE<unsigned>* comm: sol->comunidades){
                A.add(IloNumArray(env,nodes.size()));
                for (j=0;j<nodes.size();j++){
                    A[A.getSize()-1][j]= 0;
                }
//cout<<"\nAAAAA("<<comm->qtd<<")\n";cout.flush();
                itemLDE<unsigned> * aux = comm->getInicio();
                bool entrou = false;
                while(aux!=NULL){
                    if (setNodes.find(aux->id) != setNodes.end()){
                        entrou=true;
                        A[A.getSize()-1][mapNodes[aux->id]]=1;
                    }
                    aux=aux->prox;
                }
                if(! entrou){
                    A.remove(A.getSize()-1);
                }
            }

        }*/



        //creating vars from model
        //... Z
        IloNumVarArray Z(env, A.getSize(), 0.0, IloInfinity, ILOFLOAT);

        //define objective
        IloExpr masterOFExp(env);
        IloNum coef;
        //starts with Z_1 because Z_0 corresponds to an artificial cluster
        for (t=1; t<A.getSize(); t++){
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
            IloRange rng (env, 1.0, masterConsExp ,1.0);
            masterConsExp.end();
            masterConstraints.add(rng);
        }
        IloAdd(masterOpt, masterConstraints);

        //constraint \sum{z_t}=2
        IloRangeArray masterZt(env);
        IloExpr masterZtExp(env);
        for (t=0;t<A.getSize();t++){
            masterZtExp += Z[t];
        }
        IloRange rngZt (env, 2.0, masterZtExp ,2.0);
        masterZt.add(rngZt);
        IloAdd(masterOpt, masterZt);

        IloCplex masterCplex(masterOpt);
        //masterCplex.setParam(IloCplex::Param::Threads, 1);
        masterCplex.setParam(IloCplex::TiLim, 36000);


        ///heuristic for auxiliar problem
        Bipart_ILS bpils(nodes, graph, Wij, alpha);



        ///auxiliar problem
        //parameters and results (price and new column)
        IloNumArray lambda(env, nodes.size()); //dual variables from model
        IloNum lambdaY; //dual variable from model (specific bipartite)
        IloNumArray newCommunity(env, nodes.size());

        //define AA_v \in {0,1}^|V|: our a_v of the auxiliar problem
        IloNumVarArray AA(env, nodes.size(), 0,1, ILOINT);
        IloNumVar VV(env, 0,IloInfinity, ILOFLOAT);
vector< pair<unsigned, unsigned> >edges;
for (i=0;i<nodes.size();i++){
    for (j=i+1;j<nodes.size();j++){
        if (this->graph->getAdj(nodes[i], nodes[j]) > 0.9){
            edges.push_back(pair<unsigned, unsigned>(i,j));
        }
    }
}
IloNumVarArray E(env, edges.size(), 0,1, ILOINT);

        IloInt it;
        IloModel mauxOpt (env);

        //constant part of the OF (firstHalf)
        IloExpr mauxOFExpConst(env);
        /*for (i=0;i<nodes.size();i++){
            for (j=i+1;j<nodes.size();j++){
                mauxOFExpConst += 4*Wij[nodes[i]][nodes[j]]*AA[i]*AA[j];
            }
        }*/
for (i=0;i<edges.size();i++){
    mauxOFExpConst += 4*E[i];
}
        for (i=0;i<nodes.size();i++){
            mauxOFExpConst -= Di[nodes[i]]*AA[i];
        }
        //mauxOFExpConst -= Di[nodes[i]]*VV;
        /*//constant part of the OF (seconfHalf)
        for (i=0;i<nodes.size();i++){
            for (j=i+1;j<nodes.size();j++){
                mauxOFExpConst += 4*Wij[nodes[i]][nodes[j]]*(1.0-AA[i])*(1.0-AA[j]);
            }
        }
        for (i=0;i<nodes.size();i++){
            mauxOFExpConst -= Di[nodes[i]]*(1.0-AA[i]);
        }*/


        /*IloRangeArray mauxZeroA(env);
        IloExpr mauxExpZeroA(env);
        mauxExpZeroA+=AA[0];
        IloRange mauxRZeroA (env, 0.0, mauxExpZeroA, 0.0);
        mauxZeroA.add(mauxRZeroA);
        IloAdd(mauxOpt, mauxZeroA);*/

        IloExtractable zeroconst = mauxOpt.add(AA[0]==0);

        IloExpr mauxExpVV(env);
        mauxExpVV += VV;
        for (unsigned i=0;i<nodes.size();i++){
            mauxExpVV += -AA[i];
        }
        mauxOpt.add(mauxExpVV==0);
for (unsigned e=0; e<edges.size();e++){
    unsigned i = edges[e].first;
    unsigned j = edges[e].second;
    mauxOpt.add(E[e]<=AA[i]);
    mauxOpt.add(E[e]<=AA[j]);
}

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

        unsigned itOnlyHeuristic = 0;


        bool terminouZero=false;
//cout<<"\n --> 15 (antes do for)";cout.flush();
long long int tempoRMP = 0;
long long int tempoH = 0;
        for(it=1;;it++){
            itOnlyHeuristic ++;

            //solve the reduced master problem
chrono::system_clock::time_point beforeRMP;
beforeRMP = chrono::system_clock::now();
            masterCplex.solve();
tempoRMP += chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-beforeRMP).count();

//masterCplex.exportModel("antes.lp");


            if (maxDensity < masterCplex.getObjValue() ){
                maxDensity = masterCplex.getObjValue();
            }
//cout<<"\n\n ["<<it<<"]: "<<this->maxDensity;

            //store the dual variables
            for (i=0;i<nodes.size();i++){
                lambda[i] = masterCplex.getDual(masterConstraints[i]);
            }
            lambdaY = masterCplex.getDual(masterZt[0]);


//cout<<"\n --> 164 (antes do bpils.execute)";cout.flush();
            if (itOnlyHeuristic < graph->numberOfEdges){
                chrono::system_clock::time_point before;
                before = chrono::system_clock::now();
                bpils.execute(lambda, lambdaY);
                tempoH += chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
            }

//cout<<"\n --> 164 (depois do bpils.execute)";cout.flush();



            long long int totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-antes).count()/1000.0;

if(it%100 == 0){
cout<<"\n\t\t-->T(heur-sec - "<<nodes.size()<<" - "<<bpils.maxValue<<") :"<<tempoH/1000.0;cout.flush();
cout<<"\n\tT(RMP - "<<nodes.size()<<") :"<<tempoRMP/1000.0;cout.flush();
}

            if (itOnlyHeuristic >= graph->numberOfEdges || bpils.maxValue < RC_EPS || totalTime >= 36000){

cout<<"\n\tT(RMP - "<<nodes.size()<<") :"<<tempoRMP/1000.0;cout.flush();
cout<<"\n\tT(heur-sec - "<<nodes.size()<<" - "<<bpils.maxValue<<") :"<<totalTime/1000.0<<" -- "<<tempoH/1000.0;cout.flush();

//cout<<"\n --> 182 (tentanto EXATO)";cout.flush();
//            if (true){


                //execute the exact algorithm
                itOnlyHeuristic = 0;
                //first part of non-constant OF
                IloExpr mauxOFExp(env);
                IloExpr mauxOFExpFirst(env);
                for (IloNum v=0;v<nodes.size();v++){
                        mauxOFExpFirst -= AA[v]*lambda[v];
                }
                mauxOFExp += mauxOFExpFirst*VV;

                /*for (IloNum x=0;x<nodes.size();x++){
                    for (IloNum v=0;v<nodes.size();v++){
                        mauxOFExp -= AA[x]*AA[v]*lambda[v];
                    }
                }*/
                for (IloNum x=0;x<nodes.size();x++){
                        mauxOFExp -= AA[x]*lambdaY;
                }
                IloObjective obj = IloAdd(mauxOpt, IloMaximize(env, mauxOFExpConst + mauxOFExp));
                mauxOFExp.end();

                //solve column generator: and get the new column
                IloCplex mauxCplex(mauxOpt);
                //mauxCplex.setParam(IloCplex::Param::Threads,1);
                mauxCplex.setParam(IloCplex::TiLim,36000);

//cout<<"\n   (e(x(aux)))...";cout.flush();
chrono::system_clock::time_point before;
before = chrono::system_clock::now();
                mauxCplex.solve();
cout<<"\n\tT(aux-sec _ "<<mauxCplex.getObjValue()<<") :"<<to_string(chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count()/1000.0);
cout<<"end!"; cout.flush();

                //copy the new column
                mauxCplex.getValues(newCommunity, AA);

                //if solution is optimal
                if (mauxCplex.getObjValue() < RC_EPS || totalTime >= 36000){
                    if (terminouZero == false){
                        terminouZero = true;
                        mauxOpt.remove(zeroconst);
                        mauxOpt.remove(obj);
                        mauxOpt.add(AA[0]==1);//continua?
                        bpils.terminouZero = true;
                        continue;
                    }
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
                //--> secondHalf
                Z.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT));
                A.add(IloNumArray(env,nodes.size() ));
                for (j=0;j<nodes.size();j++){
                        A[t+1][j] = 1.0 - newCommunity[j];
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
                masterObj.setExpr(masterObj.getExpr()+coef * Z[t]+ coef2 * Z[t+1]);
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


