#include "lmbdgraphcoarsenersinglesteplg.h"



LMBDGraphCoarsenerSingleStepLG::LMBDGraphCoarsenerSingleStepLG(LargeGraph * graph, fibonacci_heap<DataGCSSLG> * heap, float lambda, unsigned type, bool usesHeap)
{
    this->lambda = lambda;
    this->usesHeap = usesHeap;
    this->heap = heap;
    this->density = 0.0;
    this->graph = graph;
    for (unsigned i=0; i< graph->numberOfNodes;i++){
        NodeCoarsenerSingleStepLG node;
        if (this->graph->getAdj(i,i)>0.0){
            //node.density = (4.0*this->graph->getAdj(i,i)) - graph->getDegree(i);
            node.density = (4.0*lambda*this->graph->getAdj(i,i)) -(2-2*lambda)*( graph->getDegree(i) - 2*this->graph->getAdj(i,i));
        }else{
            //node.density = graph->getDegree(i);
            node.density = -(2-2*lambda)*( graph->getDegree(i) );
            //node.density *= -1.0;
        }
        this->density += node.density;
        node.internalEdges = 0.0;
        node.totalDegree = graph->getDegree(i);
        node.nodes.insert(i);
        this->nodes.push_back(node);

    }

    for (unsigned e=0; e< graph->edges.size();e++){
        unsigned Ca = graph->edges[e].v1;
        unsigned Cb = graph->edges[e].v2;
/*if ((Ca == 998 || Ca == 349) && (Cb == 998 || Cb == 349)){
cout<<"\nCa["<<Ca<<"], Cb["<<Cb<<"]";
cout<<"\nAdj: "<<this->graph->getAdj(Ca,Cb);
cout<<"\ndeg("<<Ca<<") = "<<this->graph->degreeOfNode[Ca];
cout<<"\ndeg("<<Cb<<") = "<<this->graph->degreeOfNode[Cb];
}*/

        if (Ca == Cb){
            this->nodes[Ca].internalEdges = this->graph->getAdj(Ca,Cb);
            continue;
        }



        itemLDE<EdgeCoarsenerSingleStepLG> * itemA = new itemLDE<EdgeCoarsenerSingleStepLG>;
        itemLDE<EdgeCoarsenerSingleStepLG> * itemB = new itemLDE<EdgeCoarsenerSingleStepLG>;
        itemA->id.weight = itemB->id.weight = this->graph->getAdj(Ca,Cb);
        itemA->id.anotherEndpoint = Cb;
        itemB->id.anotherEndpoint = Ca;
        itemA->id.apointItem = itemB;
        itemB->id.apointItem = itemA;

        //creating references for the heap
        if (this->usesHeap){
            long double relDegree = graph->getDegree(Ca) + graph->getDegree(Cb);
            relDegree /= 2.0;
            DataGCSSLG data(this, Ca,Cb,1.0, relDegree, itemA->id.weight, type);
            itemA->id.heapReference = itemB->id.heapReference = this->heap->push(data);
        }
        this->nodes[Ca].adjacencies.colocar(itemA);
        this->nodes[Cb].adjacencies.colocar(itemB);
    }

}


//merge Ca with Cb and returns the enable master node
unsigned LMBDGraphCoarsenerSingleStepLG::merge(unsigned Ca, unsigned Cb, long double deltaDensity){
    //the basis will be the community with the greater number of edges
    if (this->nodes[Ca].adjacencies.getQtd() < this->nodes[Cb].adjacencies.getQtd() ){
        swap(Ca, Cb);
    }

    this->density+= deltaDensity;

    //metadata changing
    this->nodes[Ca].totalDegree += this->nodes[Cb].totalDegree;
    this->nodes[Ca].nodes.insert(this->nodes[Cb].nodes.begin(), this->nodes[Cb].nodes.end());
    this->nodes[Ca].density = deltaDensity + this->nodes[Ca].density+ this->nodes[Cb].density;
    this->nodes[Ca].internalEdges+=this->nodes[Cb].internalEdges;

    //passing edges
    itemLDE<EdgeCoarsenerSingleStepLG> * auxB = this->nodes[Cb].adjacencies.getInicio(), *erase, *next;
    while (auxB != NULL){
        //verifies if it is an internal edge
        if ( auxB->id.anotherEndpoint == Ca ){
            this->nodes[Ca].internalEdges += auxB->id.weight;
            this->nodes[Ca].adjacencies.remover(auxB->id.apointItem);
            erase = auxB;
            auxB = auxB->prox;
            this->nodes[Cb].adjacencies.remover(erase);
            continue; //continue because we pass to another adjacency.
        }

        //changing adjacent node of Cb
        unsigned Cc = auxB->id.anotherEndpoint;
        next = auxB->prox;
        itemLDE<EdgeCoarsenerSingleStepLG> * itemB = this->nodes[Cb].adjacencies.retirar(auxB);
        itemLDE<EdgeCoarsenerSingleStepLG> * anotherItemC = this->nodes[Cc].adjacencies.getInicio();
        while(anotherItemC != NULL){
            //it is not the same edge (that appints to Cb) and appoints to Ca
            if (anotherItemC != itemB->id.apointItem
                    && anotherItemC->id.anotherEndpoint == Ca){
                break;
            }
            anotherItemC=anotherItemC->prox;
        }

        if (anotherItemC == NULL){

            //the adjancent node do not be adjacent to Ca
            itemB->id.apointItem->id.anotherEndpoint = Ca;
            this->nodes[Ca].adjacencies.colocar(itemB);
        }else{
            //the adjacency already exists
            anotherItemC->id.anotherEndpoint = Ca;
            anotherItemC->id.weight += itemB->id.weight;
            anotherItemC->id.apointItem->id.weight = anotherItemC->id.weight;
            this->nodes[Cc].adjacencies.remover(itemB->id.apointItem);
            //removing from heap
            if (this->usesHeap){
                this->heap->erase(itemB->id.heapReference);
            }
            delete itemB;
        }

        auxB=next;
    }

    //disabling Cb
    this->nodes[Cb].density = MINUS_INFINITY;
    this->nodes[Cb].internalEdges = MINUS_INFINITY;
    this->nodes[Cb].totalDegree = MINUS_INFINITY;
    this->nodes[Cb].nodes.clear();
    //updating the adjacencies heap references
    itemLDE<EdgeCoarsenerSingleStepLG> * auxA = this->nodes[Ca].adjacencies.getInicio();
    while(auxA != NULL){

        long double nedges = this->nodes[Ca].internalEdges
                + this->nodes[auxA->id.anotherEndpoint].internalEdges
                + auxA->id.weight;
        long double nnodes = this->nodes[Ca].nodes.size()
                + this->nodes[auxA->id.anotherEndpoint].nodes.size();
        long double density = nedges/ (((nnodes*nnodes)-nnodes)/2.0);
        long double totalDegree =  this->nodes[Ca].totalDegree;
                + this->nodes[auxA->id.anotherEndpoint].totalDegree;
        if (this->usesHeap){
            (*auxA->id.heapReference).node1 = Ca;
            (*auxA->id.heapReference).node2 = auxA->id.anotherEndpoint;
            (*auxA->id.heapReference).weight = auxA->id.weight;
            (*auxA->id.heapReference).density = density;
            (*auxA->id.heapReference).totalRelativeDegree = totalDegree/nnodes;
            this->heap->update(auxA->id.heapReference);
        }
        auxA = auxA->prox;
    }

    return Ca;
}

Solution * LMBDGraphCoarsenerSingleStepLG::toSolution(){
    Solution * sol = new Solution(this->graph);
    unsigned com =0;
    sol->modularidade =0.0;
    for(unsigned i=0;i<this->graph->numberOfNodes;i++){
        if (this->nodes[i].density != MINUS_INFINITY){
            unordered_set<unsigned>::iterator it =   this->nodes[i].nodes.begin();
            while (it != this->nodes[i].nodes.end()){
                sol->inserirVertice(*it,com);
                it++;
            }
            sol->modularidade += this->nodes[i].density;
            com++;
        }

    }
    return sol;
}


#ifdef COLLECT_DATA_CONSTRUCTIVE
unsigned LMBDGraphCoarsenerSingleStepLG::unmerge(std::pair < NodeDetail, NodeDetail > & level,
                                             vector<unsigned> & commByNode, vector<unsigned> & commByVertice,
                                             unsigned &iter){
    //concentrate the nodes on the left --> I think that this is
    // ... already made by louvain method

    //new communities, generated
    unsigned ca = this->size;
    unsigned cb = this->size+1;
    this->size += 2;

    //updating the new two nodes splitted
    this->nodes[ca].internalEdges = level.first.edges;
    this->nodes[ca].totalDegree = level.first.degree;
    this->nodes[ca].nodes.insert(level.first.nodes->begin(), level.first.nodes->end());
    this->nodes[cb].internalEdges = level.second.edges;
    this->nodes[cb].totalDegree = level.second.degree;
    this->nodes[cb].nodes.insert(level.second.nodes->begin(), level.second.nodes->end());
    //restarting edges
    this->nodes[ca].adjacencies.qtd=0;
    this->nodes[ca].adjacencies.inicio=this->nodes[ca].adjacencies.fim=NULL;
    this->nodes[cb].adjacencies.qtd=0;
    this->nodes[cb].adjacencies.inicio=this->nodes[cb].adjacencies.fim=NULL;

//cout<<"\n 199";cout.flush();
    //computing all weights to remove
    unsigned node;
    long double weight;
    unordered_set<unsigned> all;
    all.insert(level.first.nodes->begin(), level.first.nodes->end());
    all.insert(level.second.nodes->begin(), level.second.nodes->end());
    unordered_map<unsigned, unordered_map<unsigned, long double> > toRemove;
    unordered_set<unsigned>::iterator allIt = all.begin();

    while (allIt != all.end()){
        node = *allIt;
        unsigned CNode1 = commByVertice[node];
        unordered_map<unsigned, LargeGraphEdge >::iterator adjIt = this->graph->adjNodes[node].begin();
        while (adjIt != this->graph->adjNodes[node].end()){
            unsigned CNode2 = commByVertice[adjIt->first];
            unordered_map<unsigned,long double>::const_iterator fit
                    = toRemove[CNode1].find(CNode2);
            if (fit == toRemove[CNode1].end()){
                weight = 0.0;
            }else{
                weight = toRemove[CNode1][CNode2];
            }
            if (CNode1 == CNode2 && all.find(adjIt->first)!=all.end()){
                toRemove[CNode1][CNode2] = weight + adjIt->second.weight/2.0;
            }else{
                toRemove[CNode1][CNode2] = weight + adjIt->second.weight;
            }
            adjIt++;
        }
        allIt ++;
    }
//cout<<"\n 238";cout.flush();

    //removing weights
    itemLDE<EdgeCoarsenerSingleStepLG> * rem;
    unordered_map<unsigned, unordered_map<unsigned, long double> >::iterator
          itCNode1  = toRemove.begin();
    while (itCNode1  != toRemove.end()){
        //if it is equal, then there is an internal edge to remove
        unordered_map<unsigned, long double>::iterator itCNode2;
        itCNode2 = toRemove[itCNode1->first].find(itCNode1->first);
        if (itCNode2 != toRemove[itCNode1->first].end()){
            this->nodes[itCNode1->first].internalEdges -= toRemove[itCNode1->first][itCNode2->first];

/*cout<<"\n\t["<<itCNode1->first<<"]";
unordered_set<unsigned>::iterator itdesp = this->nodes[itCNode1->first].nodes.begin();
float edg = 0;
while(itdesp != this->nodes[itCNode1->first].nodes.end()){
    unordered_set<unsigned>::iterator itdesp2 = itdesp;
    while(itdesp2 != this->nodes[itCNode1->first].nodes.end()){
        edg += this->graph->getAdj(*itdesp, *itdesp2);
        itdesp2++;
    }
    cout<<", "<<*itdesp;
    itdesp++;
}
cout<<"\nEdge: "<<edg;

itdesp = all.begin();
edg = 0;
cout<<"\n\t["<<itCNode1->first<<"]";
string det = "";
while(itdesp != all.end()){
    unordered_set<unsigned>::iterator itdesp2 = this->nodes[itCNode1->first].nodes.begin();
    while(itdesp2 != this->nodes[itCNode1->first].nodes.end()){
        if (commByVertice[*itdesp] == commByVertice[*itdesp2]){
            det+= "\n"+ to_string(*itdesp)+ " ,"+ to_string(*itdesp2) + " :"+to_string(this->graph->getAdj(*itdesp, *itdesp2));
            edg += this->graph->getAdj(*itdesp, *itdesp2);
            itdesp2++;
        }
    }
    cout<<", "<<*itdesp;
    itdesp++;
}
cout<<"\nEdge all: "<<edg;
cout<<det;
*/


//if (iter == 10) cout<<"\n"<<this->nodes[itCNode1->first].internalEdges;
        }
        //another adjacencies
        itemLDE<EdgeCoarsenerSingleStepLG> * nav = this->nodes[itCNode1->first].adjacencies.getInicio();
        while (nav!= NULL){
            rem = NULL;
            itCNode2= toRemove[itCNode1->first].find(nav->id.anotherEndpoint);
            if (itCNode2 != toRemove[itCNode1->first].end()){
                nav->id.weight-=toRemove[itCNode1->first][itCNode2->first];
                nav->id.apointItem->id.weight-=toRemove[itCNode1->first][itCNode2->first];
                if (nav->id.weight == 0.0){
                    rem = nav;
                    this->nodes[nav->id.anotherEndpoint].adjacencies.remover(nav->id.apointItem); //delete nav->id.apointItem;
                    nav = nav->prox;
                    this->nodes[itCNode1->first].adjacencies.remover(rem);//delete rem;
                }
            }
            if (rem == NULL){
                nav= nav->prox;
            }
        }
        itCNode1++;
    }
    //removing the nodes from the old containers
    allIt = all.begin();
    while (allIt!=all.end()){
        //cout<<"\n\t"<<*allIt<<" --> "<<commByVertice[*allIt];
        this->nodes[ commByVertice[*allIt] ].nodes.erase(*allIt);
        this->nodes[ commByVertice[*allIt] ].totalDegree -= this->graph->degreeOfNode[*allIt];
        allIt++;
    }
//cout<<"\n 278";cout.flush();
    //reorganizing
    allIt = this->nodes[ca].nodes.begin();
    while (allIt!=this->nodes[ca].nodes.end()){
        commByVertice[*allIt] = ca;
        allIt++;
    }
    allIt = this->nodes[cb].nodes.begin();
    while (allIt!=this->nodes[cb].nodes.end()){
        commByVertice[*allIt] = cb;
        allIt++;
    }
//cout<<"\n 290";cout.flush();
    //adding the new edges for ca
    unordered_map<unsigned, long double> toAddCa;
    unordered_map<unsigned, long double>::iterator finder;
    unordered_set<unsigned>::iterator itCa = this->nodes[ca].nodes.begin();

    while (itCa != this->nodes[ca].nodes.end()){
        node = *itCa;
        unordered_map<unsigned, LargeGraphEdge >::iterator adjIt = this->graph->adjNodes[node].begin();
        while(adjIt != this->graph->adjNodes[node].end()){
            unsigned anotherComm = commByVertice[adjIt->first];
            if (anotherComm != ca){
                finder = toAddCa.find(anotherComm);
                if (finder == toAddCa.end()){
                    toAddCa[anotherComm] = 0.0;
                }
                toAddCa[anotherComm] += adjIt->second.weight;
            }
            adjIt++;
        }
        itCa++;
    }
    unordered_map<unsigned, long double>::iterator itCaM = toAddCa.begin();
//cout<<"\n 328_";cout.flush();
    while (itCaM != toAddCa.end()){

        itemLDE<EdgeCoarsenerSingleStepLG> * itemA = new itemLDE<EdgeCoarsenerSingleStepLG>;
        itemLDE<EdgeCoarsenerSingleStepLG> * itemB = new itemLDE<EdgeCoarsenerSingleStepLG>;
        itemA->id.weight = itemB->id.weight = itCaM->second;
        itemA->id.anotherEndpoint = itCaM->first;
        itemB->id.anotherEndpoint = ca;
        itemA->id.apointItem = itemB;
        itemB->id.apointItem = itemA;

        this->nodes[ca].adjacencies.colocar(itemA);
        this->nodes[itCaM->first].adjacencies.colocar(itemB);

        itCaM++;
    }
//cout<<"\n 328a";cout.flush();

    //adding the new edges for cb
    unordered_map<unsigned, long double> toAddCb;
//cout<<"\n 331";cout.flush();
    //unordered_map<unsigned, long double>::iterator finder;
    unordered_set<unsigned>::iterator itCb = this->nodes[cb].nodes.begin();
//cout<<"\n 332";cout.flush();
    while (itCb != this->nodes[cb].nodes.end()){
//cout<<"\n 333";cout.flush();
        node = *itCb;
        unordered_map<unsigned, LargeGraphEdge >::iterator adjIt = this->graph->adjNodes[node].begin();
        while(adjIt != this->graph->adjNodes[node].end()){
//cout<<"\n 335";cout.flush();
            unsigned anotherComm = commByVertice[adjIt->first];
//cout<<"\n 336";cout.flush();
            if (anotherComm != cb && anotherComm != ca){
                finder = toAddCb.find(anotherComm);
                if (finder == toAddCb.end()){
                    toAddCb[anotherComm] = 0.0;
                }
                toAddCb[anotherComm] += adjIt->second.weight;
            }
            adjIt++;
        }
        itCb++;
    }
//cout<<"\n 350";cout.flush();
    unordered_map<unsigned, long double>::iterator itCbM = toAddCb.begin();
    while (itCbM != toAddCb.end()){
        itemLDE<EdgeCoarsenerSingleStepLG> * itemA = new itemLDE<EdgeCoarsenerSingleStepLG>;
        itemLDE<EdgeCoarsenerSingleStepLG> * itemB = new itemLDE<EdgeCoarsenerSingleStepLG>;
        itemA->id.weight = itemB->id.weight = itCbM->second;
        itemA->id.anotherEndpoint = itCbM->first;
        itemB->id.anotherEndpoint = cb;
        itemA->id.apointItem = itemB;
        itemB->id.apointItem = itemA;

        this->nodes[cb].adjacencies.colocar(itemA);
        this->nodes[itCbM->first].adjacencies.colocar(itemB);

        itCbM++;

    }
//cout<<"\n 367";cout.flush();
}
#endif


