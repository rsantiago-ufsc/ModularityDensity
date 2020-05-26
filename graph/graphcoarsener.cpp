#include "graphcoarsener.h"

GraphCoarsener::GraphCoarsener(Graph * graph)
{
    this->density = 0.0;
    this->graph = graph;
    for (unsigned i=0; i< graph->qtd_vertices;i++){
        NodeCoarsener node;
        node.density = graph->listaVertices[i].getGrau();
        node.density *= -1.0;
        this->density += node.density;
        node.internalEdges = 0.0;
        node.totalDegree = graph->listaVertices[i].getGrau();
        node.nodes.insert(i);
        this->nodes.push_back(node);
    }
    for (unsigned e=0; e< graph->qtd_ligacoes;e++){
        unsigned Ca = graph->listaArestas[e].v1;
        unsigned Cb = graph->listaArestas[e].v2;
        if (Ca == Cb){
            this->nodes[Ca].internalEdges = this->graph->matrizAdj[Ca][Cb];
            continue;
        }

        itemLDE<EdgeCoarsener> * itemA = new itemLDE<EdgeCoarsener>;
        itemLDE<EdgeCoarsener> * itemB = new itemLDE<EdgeCoarsener>;
        itemA->id.weight = itemB->id.weight = this->graph->matrizAdj[Ca][Cb];
        itemA->id.anotherEndpoint = Cb;
        itemB->id.anotherEndpoint = Ca;
        itemA->id.apointItem = itemB;
        itemB->id.apointItem = itemA;

        this->nodes[Ca].adjacencies.colocar(itemA);
        this->nodes[Cb].adjacencies.colocar(itemB);
    }

}


//merge Ca with Cb and returns the enable master node
unsigned GraphCoarsener::merge(unsigned Ca, unsigned Cb, long double deltaDensity, unsigned it){
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
    itemLDE<EdgeCoarsener> * auxB = this->nodes[Cb].adjacencies.getInicio(), *erase, *next;
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
        itemLDE<EdgeCoarsener> * itemB = this->nodes[Cb].adjacencies.retirar(auxB);
        itemLDE<EdgeCoarsener> * anotherItemC = this->nodes[Cc].adjacencies.getInicio();
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

    return Ca;
}

Solution * GraphCoarsener::toSolution(){
    Solution * sol = new Solution(this->graph);
    unsigned com =0;
    sol->modularidade =0.0;
    for(unsigned i=0;i<this->graph->qtd_vertices;i++){
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
