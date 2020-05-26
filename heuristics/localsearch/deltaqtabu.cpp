#include "deltaqtabu.h"


long double DeltaQTabu::calculateFalseDeltaDOnceNeighborhood(Solution * sol,long double dens, unsigned node, unsigned comTo){

    unsigned cv = sol->getVerticeComunidade(node); //source community

    //number of edges inside of communities 'c' or 'cv'
    long double EC = 0.0;
    long double ECV = 0.0;


    //number of edges between 'node' and communities 'c' or 'cv'
    long double Ec = 0.0;
    long double Ecv = 0.0;

    //degree of communities 'c' or 'cv'
    long double Kc = 0.0;
    long double Kcv = 0.0;

    //number of nodes inside communities 'c' and 'cv'
    long double Nc = 0.0;
    long double Ncv = 0.0;

    long double dv = graph->degreeOfNode[node];

    //obtaining the values
    itemLDE<unsigned> * auxA = sol->comunidades[cv]->getInicio();
    while (auxA != NULL){
/*if ( ((comTo==0 && cv == 11) || (comTo==11 && cv == 0)) && node == 11 && auxA->id == 11){
cout<<"aaa "<<this->graph->getAdj(auxA->id, node)<<endl;
}*/
        Ecv += (double)this->graph->getAdj(auxA->id, node);
        Kcv += (double)this->graph->getDegree(auxA->id);
        Ncv++;
        itemLDE<unsigned> * auxB = auxA;
        while (auxB != NULL){
            ECV += (double)this->graph->getAdj(auxA->id, auxB->id);
            auxB = auxB->prox;
        }
        auxA = auxA->prox;
    }
    auxA = sol->comunidades[comTo]->getInicio();
    while (auxA != NULL){
        Ec += (double)this->graph->getAdj(auxA->id, node);
        Kc += (double)this->graph->getDegree(auxA->id);
        Nc++;
        itemLDE<unsigned> * auxB = auxA;
        while (auxB != NULL){
            EC += (double)this->graph->getAdj(auxA->id, auxB->id);
            auxB = auxB->prox;
        }
        auxA = auxA->prox;
    }

    long double result = dens;
    //past
    long double dens_oldc = 0;
    if  (Nc > 0.0){
        dens_oldc = - (4.0*EC -Kc)/Nc;
    }
    long double dens_oldcv = 0;
    if  (Ncv > 0.0){
        dens_oldcv = - (4.0*ECV -Kcv)/Ncv;
    }
    if  (Ncv > 0.0 && Nc > 0){
        result += dens_oldc + dens_oldcv ;
        result += (4.0*(EC+Ec) -Kc-dv)/(Nc+1);
        if (Ncv-1 > 0){
            result += (4.0*(ECV-Ecv) -Kcv+dv)/(Ncv-1);
        }
    }else{
        long double node = this->graph->numberOfNodes;
        result = -(node*node-node)/2;
    }
/*if ( ((comTo==0 && cv == 11) || (comTo==11 && cv == 0)) && node == 11){
    cout<<"\t -> "<<result<<endl;
    cout<<"\t    ["<<comTo<<"]: ";
    itemLDE<unsigned> * aux = sol->comunidades[comTo]->inicio;
    while (aux != NULL){
        cout<<aux->id<<",";
        aux = aux->prox;
    }
    cout<<endl;
    cout<<"\t    ["<<cv<<"]: ";
    aux = sol->comunidades[cv]->inicio;
    while (aux != NULL){
        cout<<aux->id<<",";
        aux = aux->prox;
    }
    cout<<endl;

    cout<<"EC("<<EC<<"); Kc("<<Kc<<"); Nc("<<Nc<<"); Ec("<<Ec<<");"<<endl;
    cout<<"ECV("<<ECV<<"); Kcv("<<Kcv<<"); Ncv("<<Ncv<<"); Ecv("<<Ecv<<");"<<endl;
    cout<<"dv("<<dv<<"); "<<endl;
    cout<<Ncv-1<<endl;


    itemLDE <unsigned> * aux1 = sol->comunidades[comTo]->getInicio();
    long double edges       = 0.0;
    long double eee       = 0.0;
    long double totalDegree = 0.0;
    while (aux1 != NULL){
        totalDegree += graph->getDegree(aux1->id);
        eee += graph->getAdj(aux1->id,node);
        itemLDE <unsigned> * aux2 = aux1;//aux1->prox;
        while (aux2 != NULL){
            edges += graph->getAdj(aux1->id,aux2->id);
            aux2 = aux2->prox;
        }
        aux1 = aux1->prox;
    }
//        cout<<"\n:"<<( 4.0*edges - totalDegree ) / sol->comunidades[com]->qtd<<"  ==  ( 4.0*"<<edges<<" - "<<totalDegree<<" ) / "<<sol->comunidades[com]->qtd;
    long double nnodes = sol->comunidades[comTo]->qtd;
    if (nnodes > 0.0){
        cout<<eee<<"( 4.0*"<<edges<<" - "<<totalDegree<<" ) / "<<nnodes<<endl;
    }
}*/

/*
if(comTo == 0 && cv == 1){
cout<<"\tdetlaqtabu.cpp#56: Nc("<<Nc<<") Ncv("<<Ncv<<")"<<endl;
cout<<"\tdetlaqtabu.cpp#56: EC("<<EC<<") ECV("<<ECV<<")"<<endl;
cout<<"\tdetlaqtabu.cpp#56: Ec("<<Ec<<") Ecv("<<Ecv<<")"<<endl;
cout<<"\tdetlaqtabu.cpp#56: Kc("<<Kc<<") Kcv("<<Kcv<<")"<<endl;
}

    if(Nc*Nc+Nc == 0 && Ncv*Ncv-Ncv == 0){
        long double nodes = this->graph->numberOfNodes;
        result = -(nodes *nodes  - nodes )/2;
    }else{

        long double divisor = Nc*Nc+Nc;
        if (divisor == 0){
            divisor = 1;
        }
        result += (Nc*(4.0*Ec-dv)-4.0*EC + Kc)/divisor;

if(comTo == 0 && cv == 1){
cout<<"\tdetlaqtabu.cpp#56: c "<<(Nc*(4.0*Ec-dv)-4.0*EC + Kc)/divisor<<endl;
}

        divisor = Ncv*Ncv-Ncv;
        if (divisor == 0){
            divisor = 1;
        }
        result += (Ncv*(-4.0*Ecv+dv)+4.0*ECV - Kcv)/divisor;
if(comTo == 0 && cv == 1){
cout<<"\tdetlaqtabu.cpp#56: cv"<<(Ncv*(-4.0*Ecv+dv)+4.0*ECV - Kcv)/divisor<<endl;}
    }
*/
    return result;
}






DeltaQTabu::DeltaQTabu(LargeGraph * graph, ModularityLG * modularity, Solution * startSol)
{
    unsigned i;

    this->graph = graph;
    this->modularity = modularity;

    this->deltaQGeral = new fibonacci_heap<heap_data>();


    //tabu value
    double nodes = graph->numberOfNodes;
    this->LowerBound = -std::numeric_limits<double>::infinity();//-(nodes*nodes-nodes)/2;


    //instanciando quem aponta
    for (i = 0; i< startSol->comunidades.size();i++){
        this->whoPointsCom.push_back(new fibonacci_heap<heap_data>::handle_type[this->graph->numberOfNodes]);
    }

    //delta data for each community
    for (i = 0; i< startSol->comunidades.size();i++){
        deltaCommunities.push_back( deltaDataValue() );
        deltaCommunities[i].nnodes = startSol->comunidades[i]->getQtd();
        deltaCommunities[i].internalEdges = 0;
        deltaCommunities[i].totalDegree = 0;
        if (deltaCommunities[i].nnodes  == 0){
            deltaCommunities[i].density = 0;
        }else{
            itemLDE<unsigned> * aux1 = startSol->comunidades[i]->getInicio();
            while (aux1 != NULL){
                deltaCommunities[i].totalDegree += this->graph->degreeOfNode[aux1->id];
                itemLDE<unsigned> * aux2 = aux1;
                while(aux2!=NULL){
                    deltaCommunities[i].internalEdges += this->graph->getAdj(aux1->id,aux2->id);
                    aux2 = aux2->prox;
                }
                aux1 = aux1->prox;
            }
            deltaCommunities[i].density = (4.0*deltaCommunities[i].internalEdges - deltaCommunities[i].totalDegree)/deltaCommunities[i].nnodes;

//cout<<"deltaqtabu.cpp#102: ["<<i<<"]"<<deltaCommunities[i].density<<endl;
        }
    }


    //alimentando os heaps com os deltaDs iniciais
    unsigned myCom, otherCom, v;
    long double deltaQ;
    for (v=0; v<this->graph->numberOfNodes;v++){

        myCom = startSol->verticeComunidade[v];
        //para cada outra comunidade
        for (otherCom = 0; otherCom < startSol->comunidades.size(); otherCom++){

            if (otherCom != myCom){
                deltaQ = this->calculateFalseDeltaDOnceNeighborhood(startSol,0.0, v, otherCom);
//cout<<"deltaqtabu.cpp#118: ["<<myCom<<"+"<<otherCom<<"]"<<deltaQ<<endl;
            }else{
                double nodes = graph->numberOfNodes;
                deltaQ = -(nodes*nodes-nodes)/2;
            }
            heap_data * dh= new  heap_data(deltaQ, v,  otherCom);
            fibonacci_heap<heap_data>::handle_type handle = this->deltaQGeral->push(*dh);
            (*handle).handle = handle;
            this->whoPointsCom[otherCom][v] = handle;
            delete dh;
        }
    }




}

//ao aceitar uma solução, ao mudar de comunidade (depois de mudar a solução),
//  deve-se atualizar o deltaQ: esta é a função para isto!
void DeltaQTabu::changeCommunity(Solution * solAfterUpdate, unsigned olderCom, unsigned newCom, unsigned it){
    unsigned node, com, nodeCom;
    fibonacci_heap<heap_data>::handle_type handle;
    long double deltaQ;


    //se aumentou uma comunidade
    if (solAfterUpdate->comunidades.size() > this->whoPointsCom.size()){
        this->whoPointsCom.push_back(new fibonacci_heap< heap_data >::handle_type[this->graph->numberOfNodes]);
        unsigned otherCom = this->whoPointsCom.size()-1;
        for (unsigned v=0; v<this->graph->numberOfNodes;v++){
            unsigned myCom = solAfterUpdate->verticeComunidade[v];
            //para cada outra comunidade
            if (otherCom != myCom){
                deltaQ = this->calculateFalseDeltaDOnceNeighborhood(solAfterUpdate,0.0, v, otherCom);
            }else{
                double nodes = graph->numberOfNodes;
                deltaQ = -(nodes*nodes-nodes)/2;
            }
            heap_data * dh= new  heap_data(deltaQ, v,  otherCom);
            fibonacci_heap<heap_data>::handle_type handle = this->deltaQGeral->push(*dh);
            (*handle).handle = handle;
            this->whoPointsCom[otherCom][v] = handle;
            delete dh;
        }
    }

    //comunidades alteradas
    // -- antiga
    itemLDE<unsigned> * nav = solAfterUpdate->comunidades[olderCom]->getInicio();
    while (nav!=NULL){
        node = nav->id;
        for (com = 0; com<solAfterUpdate->comunidades.size(); com++){
            if(com != olderCom){
                deltaQ = this->calculateFalseDeltaDOnceNeighborhood(solAfterUpdate, 0.0, node, com);
            }else{
                double nodes = graph->numberOfNodes;
                deltaQ = -(nodes*nodes-nodes)/2;
            }

            handle = this->whoPointsCom[com][node];
            //impedindo atualizacao: tabu
            if ((*handle).deltaq > this->LowerBound){//-2.0){
                (*handle).deltaq = deltaQ;
                //(*handle).node = node;
                (*handle).community = com;
                this->deltaQGeral->update(handle);
            }

        }
        nav=nav->prox;
    }

    //comunidades alteradas
    // -- nova
    nav = solAfterUpdate->comunidades[newCom]->getInicio();
    while (nav!=NULL){
        node = nav->id;
        for (com = 0; com<solAfterUpdate->comunidades.size(); com++){
            if(com != newCom){
                deltaQ = this->calculateFalseDeltaDOnceNeighborhood(solAfterUpdate, 0.0, node, com);
            }else{
                double nodes = graph->numberOfNodes;
                deltaQ = -(nodes*nodes-nodes)/2;
            }
            handle = this->whoPointsCom[com][node];
            //impedindo atualizacao: tabu
            if ((*handle).deltaq > this->LowerBound){//-2.0){
                (*handle).deltaq = deltaQ;
                (*handle).community = com;
                this->deltaQGeral->update(handle);
            }
        }
        nav=nav->prox;
    }

    //as que apontam
    // -- para a nova e para a antiga
    for (node=0; node<this->graph->numberOfNodes;node++){
        nodeCom = solAfterUpdate->verticeComunidade[node];
        if (nodeCom != newCom && nodeCom != olderCom){//somente nos não alterados
            //nova
            com = newCom;
            deltaQ = this->calculateFalseDeltaDOnceNeighborhood(solAfterUpdate, 0.0, node, com);
            handle = this->whoPointsCom[com][node];
            //impedindo atualizacao: tabu
            if ((*handle).deltaq > this->LowerBound){//-2.0){
                (*handle).deltaq = deltaQ;
                (*handle).community = com;
                this->deltaQGeral->update(handle);
            }
            //antiga
            com = olderCom;
            deltaQ = this->calculateFalseDeltaDOnceNeighborhood(solAfterUpdate, 0.0, node, com);
            handle = this->whoPointsCom[com][node];
            //impedindo atualizacao: tabu
            if ((*handle).deltaq > this->LowerBound){//-2.0){
                (*handle).deltaq = deltaQ;
                (*handle).community = com;
                this->deltaQGeral->update(handle);
            }
        }
    }
}



//desativa elemento para a Tabu
void DeltaQTabu::disableParTabu(unsigned &node, unsigned &community){
    fibonacci_heap<heap_data>::handle_type handle = this->whoPointsCom[community][node];
    (*handle).deltaq = this->LowerBound;//-2.0;
    (*handle).community = community;
    this->deltaQGeral->update(handle);
}

//ativa elemento para a Tabu
void DeltaQTabu::enableParTabu(Solution * current, unsigned &node, unsigned &community){
    fibonacci_heap<heap_data>::handle_type handle;
    handle = this->whoPointsCom[community][node];
    if (community == current->verticeComunidade[node]){
        double nodes = graph->numberOfNodes;
        (*handle).deltaq = -(nodes*nodes-nodes)/2;
    }else{
        (*handle).deltaq = this->calculateFalseDeltaDOnceNeighborhood(current, 0.0, node, community);
    }
    (*handle).community = community;
    this->deltaQGeral->update(handle);
}
