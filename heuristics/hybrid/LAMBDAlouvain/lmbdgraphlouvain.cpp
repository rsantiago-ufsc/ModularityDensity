#include "lmbdgraphlouvain.h"



LMBDGraphLouvain::LMBDGraphLouvain(LargeGraph * graph)
{

    this->density = 0.0;
    this->graph = graph;
    for (unsigned i=0; i< graph->numberOfNodes;i++){
        NodeCoarsenerLouvain node;
        if (this->graph->getAdj(i,i)>0.0){
            node.density = (4.0*this->graph->getAdj(i,i)) - graph->getDegree(i);
        }else{
            node.density = graph->getDegree(i);
            node.density *= -1.0;
        }
        this->density += node.density;
        node.internalEdges = 0.0;
        node.totalDegree = graph->getDegree(i);
        node.nodes->insert(i);
        this->nodes.push_back(node);

    }
    for (unsigned e=0; e< graph->edges.size();e++){
        unsigned Ca = graph->edges[e].v1;
        unsigned Cb = graph->edges[e].v2;

        if (Ca == Cb){
            this->nodes[Ca].internalEdges = this->graph->getAdj(Ca,Cb);
            continue;
        }

        itemLDE<EdgeCoarsenerLouvain> * itemA = new itemLDE<EdgeCoarsenerLouvain>;
        itemLDE<EdgeCoarsenerLouvain> * itemB = new itemLDE<EdgeCoarsenerLouvain>;
        itemA->id.weight = itemB->id.weight = this->graph->getAdj(Ca,Cb);
        itemA->id.anotherEndpoint = Cb;
        itemB->id.anotherEndpoint = Ca;
        itemA->id.apointItem = itemB;
        itemB->id.apointItem = itemA;

        this->nodes[Ca].adjacencies.colocar(itemA);
        this->nodes[Cb].adjacencies.colocar(itemB);
    }
    this->size  = this->nodes.size();
}


//merge Ca with Cb and returns the enable master node
unsigned LMBDGraphLouvain::merge(unsigned Ca, unsigned Cb, long double deltaDensity){
    //the basis will be the community with the greater number of edges
    if (this->nodes[Ca].adjacencies.getQtd() < this->nodes[Cb].adjacencies.getQtd() ){
        swap(Ca, Cb);
    }

    this->density+= deltaDensity;

    //metadata changing
    this->nodes[Ca].totalDegree += this->nodes[Cb].totalDegree;
    this->nodes[Ca].nodes->insert(this->nodes[Cb].nodes->begin(), this->nodes[Cb].nodes->end());
    this->nodes[Ca].density = deltaDensity + this->nodes[Ca].density+ this->nodes[Cb].density;
    this->nodes[Ca].internalEdges+=this->nodes[Cb].internalEdges;

    //passing edges
    itemLDE<EdgeCoarsenerLouvain> * auxB = this->nodes[Cb].adjacencies.getInicio(), *erase, *next;
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
        itemLDE<EdgeCoarsenerLouvain> * itemB = this->nodes[Cb].adjacencies.retirar(auxB);
        itemLDE<EdgeCoarsenerLouvain> * anotherItemC = this->nodes[Cc].adjacencies.getInicio();
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
            delete itemB;
        }

        auxB=next;
    }


    //disabling Cb
    this->nodes[Cb].density = MINUS_INFINITY;
    this->nodes[Cb].internalEdges = MINUS_INFINITY;
    this->nodes[Cb].totalDegree = MINUS_INFINITY;
    this->nodes[Cb].nodes->clear();


    return Ca;
}

Solution * LMBDGraphLouvain::toSolution(){
    Solution * sol = new Solution(this->graph);
    unsigned com =0;
    sol->modularidade =0.0;
    for(unsigned i=0;i<this->graph->numberOfNodes;i++){
        if (this->nodes[i].density != MINUS_INFINITY){
            unordered_set<unsigned>::iterator it =   this->nodes[i].nodes->begin();
            while (it != this->nodes[i].nodes->end()){
                sol->inserirVertice(*it,com);
                it++;
            }
            sol->modularidade += this->nodes[i].density;
            com++;
        }

    }
    return sol;
}

