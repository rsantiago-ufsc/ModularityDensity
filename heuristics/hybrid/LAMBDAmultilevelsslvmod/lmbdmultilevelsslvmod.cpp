#include "lmbdmultilevelsslvmod.h"


LMBDMultiLevelSSLVMod::LMBDMultiLevelSSLVMod(LargeGraph * graph, float lambda, unsigned typeLouvain)
{
    this->lambda = lambda;
    this->graph = graph;
    this->typeLouvain = typeLouvain;
}



//starting uncoarsening from last level partition
void LMBDMultiLevelSSLVMod::execute(){
    chrono::system_clock::time_point before;
    before = chrono::system_clock::now();



    this->bestDensity = MINUS_INFINITY;
    long double m = this->graph->numberOfEdges;
    long double m2 = this->graph->numberOfEdges*2.0;

    //construction of starting solution
    LMBDSingleStepGreedyLG sslg(this->graph, lambda);
    sslg.execute();
#ifdef GROUND_TRUTH
    Solution *sol = sslg.bestToSolution();
    this->bestCommStr = sol->serialize();
#endif



    this->bestDensity  = sslg.bestDensity;
    this->bestCommSize = sslg.bestCommSize;


    this->bestMod      = sslg.bestMod;
    sslg.gcss->usesHeap = false;

    LMBDLouvainSSELV lvlg(this->graph, sslg.gcss, lambda, MLSSLV_MOD);

    //necessary preprocessing
    vector <unsigned> commByVertices (this->graph->numberOfNodes);

    unordered_set<unsigned>::iterator it;
    for (unsigned i=0;i<this->graph->numberOfNodes; i++){
        if (sslg.gcss->nodes[i].nodes.size()>0){
            //identifying where is each node
            it = sslg.gcss->nodes[i].nodes.begin();
            while (it != sslg.gcss->nodes[i].nodes.end()){
                commByVertices[*it] = i;
                it++;
            }
        }
    }


    sslg.gcss->size = this->graph->numberOfNodes;
    this->reorganization(sslg.gcss, commByVertices);
    lvlg.size = sslg.gcss->size = 1;

   //while
   this->it = 0;
   this->bestIt = 1;



   for (int l=sslg.levelDetails.size()-1; l>=0; l--){

        this->it++;
        //uncoarsening
        sslg.gcss->unmerge(sslg.levelDetails[l], lvlg.commByNode, commByVertices,this->it);

        //reorganizing the blank communities

        this->reorganization(sslg.gcss, commByVertices);


        //intensification

        lvlg.update(sslg.gcss);
        lvlg.execute(commByVertices);


        if (this->bestDensity < lvlg.bestDensity){
            this->bestDensity  = lvlg.bestDensity;
            this->bestMod      = lvlg.bestMod;
            this->bestIt       = this->it;
            this->bestCommSize = lvlg.bestNumberOfCommunities;
#ifdef GROUND_TRUTH
            this->bestCommStr  =  lvlg.bestCommunityStr;
//cout<<"\nDOI:"<<this->bestCommStr<<"\n";
//return;
#endif

//            cout<<bestCommStr<<"\n";
            this->bestTime     = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
        }

//if (this->it == 10)        return;
    }
   this->totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

}


void LMBDMultiLevelSSLVMod::reorganization(LMBDGraphCoarsenerSingleStepLG *&gcss, vector <unsigned> & commByVertice){
    /*vector<unsigned> translator(gcss->size, NO_COMMUNITY);
    for (unsigned x=0;x<translator.size();x++){
        translator[x] = x;
    }*/
    vector< pair<unsigned, unsigned> > trans;

    unsigned ileft=0;

    unsigned iright=gcss->size-1;

    //before = chrono::system_clock::now();

    while (ileft < iright){
        if (gcss->nodes[ileft].nodes.size() > 0 ){
//            translator[ileft] = ileft;
            ileft++;
            continue;
        }
        if (gcss->nodes[iright].nodes.size() == 0 ){
            iright--;
            continue;
        }

        swap(gcss->nodes[ileft], gcss->nodes[iright]);
        trans.push_back(pair<unsigned, unsigned>(ileft,iright));
//        translator[iright] = ileft;
        ileft++;
        iright--;
    }
    gcss->size =  iright+1;
//    timeReor1 += chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

    itemLDE<EdgeCoarsenerSingleStepLG> *aux;

    vector<pair<unsigned, unsigned>>::iterator it;

    unordered_set<unsigned>::iterator itNodes;

    it = trans.begin();

    while(it != trans.end()){
        aux = gcss->nodes[(*it).first].adjacencies.getInicio();
        while(aux!=NULL){
            aux->id.apointItem->id.anotherEndpoint = (*it).first;
            aux=aux->prox;
        }

        //updating commByVertice (vertices by community
        itNodes=gcss->nodes[(*it).first].nodes.begin();
        while(itNodes!=gcss->nodes[(*it).first].nodes.end()){
            commByVertice[*itNodes] = (*it).first;
            itNodes++;
        }
        it++;
    }
    /*
    for (unsigned in = 0;in<gcss->size; in++){
        aux = gcss->nodes[in].adjacencies.getInicio();
        while(aux!=NULL){
            aux->id.anotherEndpoint = translator[aux->id.anotherEndpoint];
            aux=aux->prox;
        }
    }*/
  //  timeReor2 += chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();


}




//starting uncoarsening from level of the best partition
void LMBDMultiLevelSSLVMod::executeFromBest(){

    this->bestDensity = MINUS_INFINITY;

    //construction of starting solution
    LMBDSingleStepGreedyLG sslg(this->graph, lambda);
    sslg.execute();

    this->bestDensity = sslg.bestDensity;
    //this->bestCommSize = sslg.

    LMBDGraphCoarsenerSingleStepLG * gcBest = new LMBDGraphCoarsenerSingleStepLG(this->graph, NULL,lambda, TYPE_PRI_DENSITY, false);
    for (unsigned l=0;l<=sslg.bestIt;l++){
        pair<unsigned,unsigned> p = sslg.levels[l].first;
        gcBest->merge(p.first, p.second, 0);
    }

    LMBDLouvainSSELV lvlg(this->graph, gcBest, lambda, MLSSLV_MOD);

    vector <unsigned> commByVertices (this->graph->numberOfNodes);
    //identifying where is each node
    unordered_set<unsigned>::iterator it;
    for (unsigned i=0;i<this->graph->numberOfNodes; i++){
//        cout<<"\n  ini.i:"<<i;cout.flush();
        if (gcBest->nodes[i].nodes.size()>0){
            it = gcBest->nodes[i].nodes.begin();
            while (it != gcBest->nodes[i].nodes.end()){
                commByVertices[*it] = i;
                it++;
            }
        }



    }


    gcBest->size = this->graph->numberOfNodes;
    this->reorganization(gcBest, commByVertices);
    lvlg.size = gcBest->size;// = 1;



    //while
    this->it = 0;

    unordered_set<unsigned>::iterator iit;

   for (int l=sslg.bestIt;l>=0;l--){
        this->it++;

        //uncoarsening
        gcBest->unmerge(sslg.levelDetails[l], lvlg.commByNode, commByVertices,this->it);

        //reorganizing the blank communities
        this->reorganization(gcBest, commByVertices);

        //intensification
        lvlg.update(gcBest);
        lvlg.execute(commByVertices);

        if (this->bestDensity < lvlg.bestDensity){
            this->bestDensity = lvlg.bestDensity;
        }
    }
}


//uses coarsened nodes from SingleStep on singleton partitions
void LMBDMultiLevelSSLVMod::executeAlwaysSingleton(){
    this->bestDensity = MINUS_INFINITY;

    //construction of starting solution
    LMBDSingleStepGreedyLG sslg(this->graph, lambda);
    sslg.execute();

    this->bestDensity = sslg.bestDensity;
    sslg.gcss->usesHeap = false;


    //identifying where is each node
    LMBDGraphCoarsenerSingleStepLG * gcSingleAux = new LMBDGraphCoarsenerSingleStepLG(this->graph, NULL, lambda, TYPE_PRI_DENSITY, false);
    LMBDGraphCoarsenerSingleStepLG * gcSingle = new LMBDGraphCoarsenerSingleStepLG(this->graph, NULL, lambda, TYPE_PRI_DENSITY, false);
    vector <unsigned> commByVertices (this->graph->numberOfNodes);
    unordered_set<unsigned>::iterator it;
    for (unsigned i=0;i<this->graph->numberOfNodes; i++){
        commByVertices[i] = i;
    }
    gcSingleAux->size = this->graph->numberOfNodes;
    gcSingle->size = this->graph->numberOfNodes;

    //while
    this->it = 0;
    unordered_set<unsigned>::iterator iit;
    for (int l=0;l<sslg.levelDetails.size()-1;l++){
        this->it++;
/*cout<<"\n 381 << iter["<<iter<<"]";cout.flush();
for (unsigned com = 0; com<gcSingle->nodes.size();com++){
    cout<<"\ncom["<<com<<"]"<<gcSingle->nodes[com].internalEdges<<", "
       <<gcSingle->nodes[com].totalDegree<<", "<<gcSingle->nodes[com].totalDegree
      <<", "<<gcSingle->nodes[com].nodes.size();
    cout<<"::: |"<<gcSingle->nodes[com].adjacencies.getQtd()<<"| ";
    itemLDE<EdgeCoarsenerSingleStepLG> * aux = gcSingle->nodes[com].adjacencies.getInicio();
    while(aux!=NULL){
        cout<<aux->id.anotherEndpoint<<", ";
        aux = aux->prox;
    }
}
cout.flush();*/


        //meging the nodes from level
        unordered_set<unsigned> all;
        all.insert(sslg.levelDetails[l].first.nodes->begin(), sslg.levelDetails[l].first.nodes->end());
        all.insert(sslg.levelDetails[l].second.nodes->begin(), sslg.levelDetails[l].second.nodes->end());

        unordered_set<unsigned>::iterator allIt = all.begin();
        unsigned master = *allIt;
        allIt++;
        while (allIt != all.end()){

            master = gcSingle->merge(master, *allIt, 0.0);
            allIt++;
        }

        allIt = all.begin();
        while (allIt != all.end()){

            commByVertices[*allIt] = master;
            allIt++;
        }


        gcSingle->size = this->graph->numberOfNodes - all.size() + 1;

        this->reorganization(gcSingle, commByVertices);

        //intensification
        LMBDLouvainSSELV lvlg(this->graph, gcSingle,lambda, MLSSLV_MOD);
        lvlg.execute(commByVertices);
        if (this->bestDensity < lvlg.bestDensity){
            this->bestDensity = lvlg.bestDensity;
        }

        /*I must to do unmerge to all nodes
         * (perhaps I must maintain a auxiliar
         * structure to some properties)
         */
        //Removing all edges from the gcSingle
        for(unsigned com=0; com<graph->numberOfNodes;com++){
            while (gcSingle->nodes[com].adjacencies.getQtd() > 0){
                gcSingle->nodes[com].adjacencies.remover(gcSingle->nodes[com].adjacencies.inicio);
            }
        }
        //Copying node metadatas from the gcSingleAux to the gcSingle
        for(unsigned com=0; com<gcSingleAux->size;com++){
            gcSingle->nodes[com].internalEdges = gcSingleAux->nodes[com].internalEdges;
            gcSingle->nodes[com].totalDegree = gcSingleAux->nodes[com].totalDegree;
            gcSingle->nodes[com].nodes.clear();
            gcSingle->nodes[com].nodes.insert(com);
            commByVertices[com]=com;
        }

        //Restating edges of the gcSingle
        for (unsigned e=0; e< graph->edges.size();e++){
            unsigned Ca = graph->edges[e].v1;
            unsigned Cb = graph->edges[e].v2;
            if (Ca == Cb){
                gcSingle->nodes[Ca].internalEdges = graph->getAdj(Ca,Cb);
                continue;
            }
            itemLDE<EdgeCoarsenerSingleStepLG> * itemA = new itemLDE<EdgeCoarsenerSingleStepLG>;
            itemLDE<EdgeCoarsenerSingleStepLG> * itemB = new itemLDE<EdgeCoarsenerSingleStepLG>;
            itemA->id.weight = itemB->id.weight = this->graph->getAdj(Ca,Cb);
            itemA->id.anotherEndpoint = Cb;
            itemB->id.anotherEndpoint = Ca;
            itemA->id.apointItem = itemB;
            itemB->id.apointItem = itemA;

            gcSingle->nodes[Ca].adjacencies.colocar(itemA);
            gcSingle->nodes[Cb].adjacencies.colocar(itemB);

        }
    }
}
