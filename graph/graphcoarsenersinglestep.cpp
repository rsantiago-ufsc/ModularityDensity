#include "graphcoarsenersinglestep.h"

GraphCoarsenerSingleStep::GraphCoarsenerSingleStep(Graph * graph, fibonacci_heap<DataGCSS> * heap)
{
    this->heap = heap;
    this->density = 0.0;
    this->graph = graph;
    for (unsigned i=0; i< graph->qtd_vertices;i++){
        NodeCoarsenerSingleStep node;
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



        itemLDE<EdgeCoarsenerSingleStep> * itemA = new itemLDE<EdgeCoarsenerSingleStep>;
        itemLDE<EdgeCoarsenerSingleStep> * itemB = new itemLDE<EdgeCoarsenerSingleStep>;
        itemA->id.weight = itemB->id.weight = this->graph->matrizAdj[Ca][Cb];
        itemA->id.anotherEndpoint = Cb;
        itemB->id.anotherEndpoint = Ca;
        itemA->id.apointItem = itemB;
        itemB->id.apointItem = itemA;

        //creating references for the heap
        long double relDegree = graph->listaVertices[Ca].getGrau() + graph->listaVertices[Cb].getGrau();
        relDegree /= 2.0;
        DataGCSS data(Ca,Cb,1.0, relDegree, itemA->id.weight);
        itemA->id.heapReference = itemB->id.heapReference = this->heap->push(data);

        this->nodes[Ca].adjacencies.colocar(itemA);
        this->nodes[Cb].adjacencies.colocar(itemB);
    }

}


//merge Ca with Cb and returns the enable master node
unsigned GraphCoarsenerSingleStep::merge(unsigned Ca, unsigned Cb, long double deltaDensity){
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
    itemLDE<EdgeCoarsenerSingleStep> * auxB = this->nodes[Cb].adjacencies.getInicio(), *erase, *next;
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
        itemLDE<EdgeCoarsenerSingleStep> * itemB = this->nodes[Cb].adjacencies.retirar(auxB);
        itemLDE<EdgeCoarsenerSingleStep> * anotherItemC = this->nodes[Cc].adjacencies.getInicio();
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
            this->heap->erase(itemB->id.heapReference);

            delete itemB;
        }



        auxB=next;
    }


    //disabling Cb
    this->nodes[Cb].density = MINUS_INFINITY;
    this->nodes[Cb].internalEdges = MINUS_INFINITY;
    this->nodes[Cb].totalDegree = MINUS_INFINITY;

    //updating the adjacencies heap references
    itemLDE<EdgeCoarsenerSingleStep> * auxA = this->nodes[Ca].adjacencies.getInicio();
    while(auxA != NULL){

        long double nedges = this->nodes[Ca].internalEdges
                + this->nodes[auxA->id.anotherEndpoint].internalEdges
                + auxA->id.weight;
        long double nnodes = this->nodes[Ca].nodes.size()
                + this->nodes[auxA->id.anotherEndpoint].nodes.size();
        long double density = nedges/ (((nnodes*nnodes)-nnodes)/2.0);
        long double totalDegree =  this->nodes[Ca].totalDegree;
                + this->nodes[auxA->id.anotherEndpoint].totalDegree;

        (*auxA->id.heapReference).node1 = Ca;
        (*auxA->id.heapReference).node2 = auxA->id.anotherEndpoint;
        (*auxA->id.heapReference).weight = auxA->id.weight;
        (*auxA->id.heapReference).density = density;
        (*auxA->id.heapReference).totalRelativeDegree = totalDegree/nnodes;
//cout<<"\n  nedges["<<nedges<<"] nnodes["<<nnodes<<"] |"<<Ca<<"|  |"<<auxA->id.anotherEndpoint<<"|["<< this->nodes[auxA->id.anotherEndpoint].nodes.size()<<"]<< "<<density<<" | "<<totalDegree/nnodes<<" | "<<Ca<<"--"<<auxA->id.anotherEndpoint;
        this->heap->update(auxA->id.heapReference);

        auxA = auxA->prox;
    }


    return Ca;
}

Solution * GraphCoarsenerSingleStep::toSolution(){
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
