#include "modularitylg.h"


#include <math.h>
ModularityLG::ModularityLG(LargeGraph * grafo)
{
    this->grafo=grafo;
}

long double ModularityLG::calculateDensity(Solution*sol){
    LargeGraph * g = this->grafo;
    long double density = 0.0;
    for (unsigned com=0; com<sol->comunidades.size();com++){
        itemLDE <unsigned> * aux1 = sol->comunidades[com]->getInicio();
        long double edges       = 0.0;
        long double totalDegree = 0.0;
        while (aux1 != NULL){
            totalDegree += g->getDegree(aux1->id);
            itemLDE <unsigned> * aux2 = aux1;//aux1->prox;
            while (aux2 != NULL){
                edges += g->getAdj(aux1->id,aux2->id);
                aux2 = aux2->prox;
            }
            aux1 = aux1->prox;
        }
//        cout<<"\n:"<<( 4.0*edges - totalDegree ) / sol->comunidades[com]->qtd<<"  ==  ( 4.0*"<<edges<<" - "<<totalDegree<<" ) / "<<sol->comunidades[com]->qtd;
        long double nnodes = sol->comunidades[com]->qtd;
        if (nnodes > 0.0){
            density += ( 4.0*edges - totalDegree ) / nnodes ;
        }
    }
    return density;
}


long double ModularityLG::calculateDensity(Solution*sol, float lambda){
    LargeGraph * g = this->grafo;
    long double density = 0.0;
    for (unsigned com=0; com<sol->comunidades.size();com++){
        itemLDE <unsigned> * aux1 = sol->comunidades[com]->getInicio();
        long double edges       = 0.0;
        long double totalDegree = 0.0;
        while (aux1 != NULL){
            totalDegree += g->getDegree(aux1->id);
            itemLDE <unsigned> * aux2 = aux1;//aux1->prox;
            while (aux2 != NULL){
                edges += g->getAdj(aux1->id,aux2->id);
                aux2 = aux2->prox;
            }
            aux1 = aux1->prox;
        }
//        cout<<"\n:"<<( 4.0*edges - totalDegree ) / sol->comunidades[com]->qtd<<"  ==  ( 4.0*"<<edges<<" - "<<totalDegree<<" ) / "<<sol->comunidades[com]->qtd;
        long double nnodes = sol->comunidades[com]->qtd;
        if (nnodes > 0.0){
            //density += ( 4.0*edges - totalDegree ) / nnodes ;
            density += ( 4.0*lambda*edges - (2-2*lambda)*(totalDegree - 2*edges) ) / nnodes ;
        }
    }
    return density;
}





long double ModularityLG::calculateModularityND_NW(Solution * sol){
    long double mod = 0.0;
    for (unsigned ic=0;ic< sol->getComunidades()->size();ic++){
        itemLDE<unsigned> * auxV1 = (*sol->getComunidades())[ic]->getInicio();
        while (auxV1 != NULL){
            unsigned dg1 = this->grafo->getDegree(auxV1->id);
            itemLDE<unsigned> * auxV2 = (*sol->getComunidades())[ic]->getInicio();
            while (auxV2 != NULL){
                unsigned dg2 = this->grafo->getDegree(auxV2->id);
                long double aij = this->grafo->getAdj(auxV1->id,auxV2->id);

                mod += aij - ((long double)(dg1*dg2))/(2.0*grafo->numberOfEdges);
                auxV2 = auxV2->prox;
            }
            auxV1 = auxV1->prox;
        }
    }
    //cout<<"Qtd. de ligacoes: "<<(2*grafo->qtd_ligacoes)<<"\n";
    mod = mod /(2*grafo->numberOfEdges);
    return mod;
}


