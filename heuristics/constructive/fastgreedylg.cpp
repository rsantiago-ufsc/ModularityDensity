#include "fastgreedylg.h"
#include "../../utils/modularitylg.h"


FastGreedyLG::FastGreedyLG(GraphCoarsenerSingleStepLG *gcss)
{


    this->gcss = gcss;
    this->graph = gcss->graph;

    //Storing the exact community where each node belongs.
    for (unsigned v=0;v < gcss->graph->numberOfNodes;v++){
        this->nodeCommunity.push_back(0);
    }

    unsigned cr=0;
    for (unsigned com=0;com < gcss->nodes.size();com++){
        if (gcss->nodes[com].density > MINUS_INFINITY){
            unordered_set<unsigned>::iterator aux = this->gcss->nodes[com].nodes.begin();
            while(aux != this->gcss->nodes[com].nodes.end()){
                this->nodeCommunity[*aux] = cr;
                aux++;
            }
        }
        cr++;
    }
//cout<<"28: ";cout.flush();
    //storing data about communities
    for (unsigned com=0;com < gcss->nodes.size();com++){
        unsigned nnodes = 0;
        if (gcss->nodes[com].density > MINUS_INFINITY){
            nnodes = gcss->nodes[com].nodes.size();
        }
        this->communities.push_back(FGCommunityLG(  gcss->nodes[com].internalEdges,
                                                     gcss->nodes[com].totalDegree,
                                                     nnodes,
                                                     gcss->nodes[com].density));
    }


//cout<<"42: ";cout.flush();
}



FastGreedyLG::FastGreedyLG(Solution *base)
{

    this->basedSolution = base;
    this->graph = base->lggrafo;

    //Storing the exact community where each node belongs.
    for (unsigned v=0;v < this->graph->numberOfNodes;v++){
        this->nodeCommunity.push_back(0);
    }

    for (unsigned v=0;v < this->graph->numberOfNodes; v++){
        this->nodeCommunity[v] = this->basedSolution->verticeComunidade[v];
    }



    //storing data about communities
    for (unsigned com=0;com < this->basedSolution->comunidades.size();com++){
        unsigned nnodes = this->basedSolution->comunidades[com]->qtd;
        itemLDE <unsigned> * aux1 = this->basedSolution->comunidades[com]->getInicio();
        long double edges       = 0.0;
        long double totalDegree = 0.0;
        while (aux1 != NULL){
            totalDegree += this->graph->getDegree(aux1->id);
            itemLDE <unsigned> * aux2 = aux1->prox;
            while (aux2 != NULL){
                edges += this->graph->getAdj(aux1->id,aux2->id);
                aux2 = aux2->prox;
            }
            aux1 = aux1->prox;
        }

        this->communities.push_back(FGCommunityLG(  edges,
                                                    totalDegree,
                                                    nnodes,
                                                    (4.0*edges-totalDegree)/nnodes ));
    }


}


FastGreedyLG::FastGreedyLG(LargeGraph * graph)
{

    //this->basedSolution = base;
    this->graph = graph;
    this->gcss = NULL;

    //Storing the exact community where each node belongs.
    for (unsigned v=0;v < this->graph->numberOfNodes;v++){
        this->nodeCommunity.push_back(0);
    }

    for (unsigned v=0;v < this->graph->numberOfNodes; v++){
        this->nodeCommunity[v] = v;
    }

    //storing data about communities
    this->bestDensity = 0;
    for (unsigned com=0;com < this->graph->numberOfNodes;com++){
        unsigned nnodes = 1;
        long double edges       = this->graph->getAdj(com,com);
        long double totalDegree = this->graph->getDegree(com);
        this->communities.push_back(FGCommunityLG(  edges,
                                                    totalDegree,
                                                    nnodes,
                                                    (4.0*edges-totalDegree)/nnodes ));
       this->bestDensity += (4.0*edges-totalDegree)/nnodes;
    }




}





//just update the structures for GraphCoarsenerSingleStepLG
void FastGreedyLG::update(GraphCoarsenerSingleStepLG *gcss){
    this->gcss = gcss;
    this->graph = gcss->graph;
    //storing the node and community information
    unsigned cr=0;
    for (unsigned com=0;com < gcss->nodes.size();com++){
        if (gcss->nodes[com].density > MINUS_INFINITY){
            unordered_set<unsigned>::iterator aux = this->gcss->nodes[com].nodes.begin();
            while(aux != this->gcss->nodes[com].nodes.end()){
                this->nodeCommunity[*aux] = com;
                aux++;
            }
            cr++;
        }
    }


    //storing data about communities
    cr = 0;
//long double x=0.0;
    for (unsigned com=0;com < gcss->nodes.size();com++){
        unsigned nnodes = 0;
        long double internalEdges   = 0.0;
        long double totalDegree     = 0.0;
        long double density         = 0.0;

        if (gcss->nodes[com].density > MINUS_INFINITY){
            nnodes = gcss->nodes[com].nodes.size();
            //internalEdges  = gcss->nodes[com].internalEdges;
            //totalDegree    = gcss->nodes[com].totalDegree;
            //density        = gcss->nodes[com].density;

        }
        if (cr == this->communities.size()){ this->communities.push_back(FGCommunityLG()); }
        this->communities[cr].internalEdges  = gcss->nodes[com].internalEdges;
        this->communities[cr].totalDegree    = gcss->nodes[com].totalDegree;
        this->communities[cr].nnodes         = nnodes;
        this->communities[cr].density        = gcss->nodes[com].density;
//x+=gcss->nodes[com].density;
        cr++;
    }
//cout<<"\nX: "<<x;

    //erasing the last elements
/*    if (gcss->nodes.size() < this->communities.size()){
        unsigned difference = this->communities.size() - gcss->nodes.size();
        vector<FGCommunityLG>::iterator del, it = this->communities.end();
        for (unsigned i=0;i<difference;i++){
            del = it;
            it--;
            this->communities.erase(del);

        }

    }*/

}

void FastGreedyLG::execute(long double minimalDensity){
    chrono::system_clock::time_point before;
    before = chrono::system_clock::now();

    unsigned randomNodes[this->graph->numberOfNodes];
    for (unsigned i=0;i<this->graph->numberOfNodes;i++){
        randomNodes[i]=i;
    }
    unsigned pos = this->graph->numberOfNodes;
    for (unsigned i=0;i<this->graph->numberOfNodes;i++){
        pos = rand()%this->graph->numberOfNodes;
        swap(randomNodes[i], randomNodes[pos]);
    }

    long double olderDensity;
    if (this->basedSolution!=NULL){
        olderDensity= this->basedSolution->modularidade;
    }else{
        if (this->gcss != NULL){
            olderDensity= this->gcss->density;
        }else{
            olderDensity= this->bestDensity;
        }
    }

    bool improvementFound=true;

    long double * tempEdgesInternal = new long double [this->communities.size()];
    //cleaning temporary
    for (unsigned com =0;com < this->communities.size(); com++){
        tempEdgesInternal[com] = 0.0;
    }

t1+=chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();


    unsigned iter = 0;
    while (improvementFound == true){
//cout<<"e["<<iter<<"]\n";cout.flush();

        this->bestTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
        iter++;
        improvementFound=false;
        for (unsigned x=0; x<this->nodeCommunity.size();x++){

            unsigned v = randomNodes[x];
            unordered_set<unsigned> adjCommunities;
            unsigned vcom = this->nodeCommunity[v];
before = chrono::system_clock::now();
            //what is the best community to stay?
            //answer: find out the best community for v
            unordered_map<unsigned, LargeGraphEdge >::const_iterator itAdj = this->graph->adjNodes[v].begin();
            while(itAdj != this->graph->adjNodes[v].end()){

                unsigned u = itAdj->first;
                unsigned ucom = this->nodeCommunity[u];
                tempEdgesInternal[ucom] += itAdj->second.weight;

                if (vcom != ucom){

                    adjCommunities.insert(ucom);

                }

                itAdj++;

            }
t2+=chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
before = chrono::system_clock::now();
            //best neighborhood //## há como melhorar (colocar em 'n' heaps e planejar atualização)
            long double maxGainDensity = MINUS_INFINITY;
            long double gain;
            unsigned maxCom;
            long double vdegree = this->graph->getDegree(v);
            unordered_set<unsigned>::iterator it = adjCommunities.begin();
            while (it != adjCommunities.end()){
                gain = (4.0*(tempEdgesInternal[*it] + communities[*it].internalEdges)
                        -communities[*it].totalDegree - vdegree)/(communities[*it].nnodes+1.0)
                        -communities[*it].density;
                if (gain > maxGainDensity){
                    maxGainDensity = gain;
                    maxCom = *it;

                }
                it++;
            }
t3+=chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
before = chrono::system_clock::now();
            //there are adjacencies for node 'v'?
            if (maxGainDensity > MINUS_INFINITY){
                //calculating the lost of old community
                long double lost = -communities[vcom].density;
                if (communities[vcom].nnodes-1.0 > 0.0){
                        lost =(4.0*(-tempEdgesInternal[vcom] + communities[vcom].internalEdges)
                        -communities[vcom].totalDegree + vdegree)/(communities[vcom].nnodes-1.0)
                        -communities[vcom].density;
                }
                long double maxDeltaDensity= lost+maxGainDensity;

                //changing community
                if (maxDeltaDensity > 0.0 ){


                    improvementFound = true;


                    //adapting the changes

                    this->nodeCommunity[v] = maxCom;

                    //old community
                    communities[vcom].internalEdges -= tempEdgesInternal[vcom];
                    communities[vcom].totalDegree   -= vdegree;
                    communities[vcom].nnodes        -= 1.0;
                    communities[vcom].density       += lost;

                    //new community
                    communities[maxCom].internalEdges  += tempEdgesInternal[maxCom];
                    communities[maxCom].totalDegree   += vdegree;
                    communities[maxCom].nnodes        += 1.0;
                    communities[maxCom].density       += maxGainDensity;




                }

            }
t4+=chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
before = chrono::system_clock::now();
            //cleaning temporary
            for (unsigned y=0;y<this->communities.size(); y++){
                tempEdgesInternal[y] = 0.0;
            }
t5+=chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
            /*it = adjCommunities.begin();
            while (it != adjCommunities.end()){
                tempEdgesInternal[*it] = 0.0;
                it++;
            }*/
        }


    }
    this->maxIt=iter;
    this->bestIt=iter-1;



    delete tempEdgesInternal;


    long double myDensity = 0.0;
    for (unsigned com=0; com<this->communities.size();com++){
        if (communities[com].nnodes > 0){
            myDensity += communities[com].density;
        }
    }
    if (olderDensity > myDensity){
        myDensity = olderDensity;
    }
/*cout<<"\nolderDensity:\t"<<olderDensity;
cout<<"\nmyDensity:\t"<<myDensity;
ModularityLG modul(this->graph);
Solution * sol = this->toSolution();
cout<<"\ntrue:\t"<<modul.calculateDensity(sol);
cout<<"\nsol:\t"<<sol->serialize();
*/
    if (myDensity>this->bestDensity){
        this->bestDensity = myDensity;
    }
    this->totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
}


Solution * FastGreedyLG::executeAndConstruct(long double minimalDensity){


    bool improvementFound=true;




    long double * tempEdgesInternal = new long double [this->communities.size()];
    //cleaning temporary
    for (unsigned com =0;com < this->communities.size(); com++){
        tempEdgesInternal[com] = 0.0;
    }

    unsigned iter = 0;
    while (improvementFound == true){
        iter++;

        improvementFound=false;

        for (unsigned v=0; v<this->nodeCommunity.size();v++){

            unordered_set<unsigned> adjCommunities;

            unsigned vcom = this->nodeCommunity[v];


            //what is the best community to stay?
            //answer: find out the best community for v
            unordered_map<unsigned, LargeGraphEdge >::const_iterator itAdj = this->graph->adjNodes[v].begin();

            while(itAdj != this->graph->adjNodes[v].end()){

                unsigned u = itAdj->first;
                unsigned ucom = this->nodeCommunity[u];
                tempEdgesInternal[ucom] += itAdj->second.weight;

                if (vcom != ucom){

                    adjCommunities.insert(ucom);

                }

                itAdj++;

            }

            //best neighborhood //## há como melhorar (colocar em 'n' heaps e planejar atualização)
            long double maxGainDensity = MINUS_INFINITY;
            long double gain;
            unsigned maxCom;
            long double vdegree = this->graph->getDegree(v);
            unordered_set<unsigned>::iterator it = adjCommunities.begin();
            while (it != adjCommunities.end()){
                gain = (4.0*(tempEdgesInternal[*it] + communities[*it].internalEdges)
                        -communities[*it].totalDegree - vdegree)/(communities[*it].nnodes+1.0)
                        -communities[*it].density;
                if (gain > maxGainDensity){
                    maxGainDensity = gain;
                    maxCom = *it;

                }
                it++;
            }

            //there are adjacencies for node 'v'?
            if (maxGainDensity > MINUS_INFINITY){
                //calculating the lost of old community
                long double lost = -communities[vcom].density;
                if (communities[vcom].nnodes-1.0 > 0.0){
                        lost =(4.0*(-tempEdgesInternal[vcom] + communities[vcom].internalEdges)
                        -communities[vcom].totalDegree + vdegree)/(communities[vcom].nnodes-1.0)
                        -communities[vcom].density;
                }
                long double maxDeltaDensity= lost+maxGainDensity;

                //changing community
                if (maxDeltaDensity > 0.0 ){


                    improvementFound = true;


                    //adapting the changes

                    this->nodeCommunity[v] = maxCom;

                    //old community
                    communities[vcom].internalEdges -= tempEdgesInternal[vcom];
                    communities[vcom].totalDegree   -= vdegree;
                    communities[vcom].nnodes        -= 1.0;
                    communities[vcom].density       += lost;

                    //new community
                    communities[maxCom].internalEdges  += tempEdgesInternal[maxCom];
                    communities[maxCom].totalDegree   += vdegree;
                    communities[maxCom].nnodes        += 1.0;
                    communities[maxCom].density       += maxGainDensity;




                }

            }

            //cleaning temporary
            for (unsigned y=0;y<this->communities.size(); y++){
                tempEdgesInternal[y] = 0.0;
            }
            /*it = adjCommunities.begin();
            while (it != adjCommunities.end()){
                tempEdgesInternal[*it] = 0.0;
                it++;
            }*/
        }


    }
    delete tempEdgesInternal;


    long double myDensity = 0.0;
    for (unsigned com=0; com<this->communities.size();com++){
        if (communities[com].nnodes > 0){
            myDensity += communities[com].density;
        }
    }

    if (minimalDensity > myDensity){
        return NULL;
    }


    Solution * sol = new Solution(this->graph);
    sol->modularidade = myDensity;
    unsigned * communitiesSol= new unsigned[this->communities.size()];
    for (unsigned com=0; com<this->communities.size();com++){
        communitiesSol[com] = MAX_UNSIGNED;
    }
    unsigned lastC = 0;
    for (unsigned v=0; v<this->nodeCommunity.size();v++){
        if (communitiesSol[this->nodeCommunity[v]]==MAX_UNSIGNED){
            communitiesSol[this->nodeCommunity[v]] = lastC;
            lastC++;
        }
        sol->inserirVertice(v,communitiesSol[this->nodeCommunity[v]]);
    }

    delete communitiesSol;
    this->iter = iter;

    return sol;

}




/*
//just update the structures for the merge between first and second solution
void FastGreedyLG::update(unsigned first, unsigned second){
//##é melhor partir da singlestepgreedylg, para evitar mais processamento inútil. Preciso criar um single step, passo a passo e permitir que este método melhore cada passo
    itemLDE<unsigned> * aux = this->basedSolution->comunidades[second]->inicio;
    while (aux != NULL){
        sol->verticeComunidade[aux->id] = first;
        this->nodeCommunity[*aux] = first;
        aux = aux->prox;
    }
    this->basedSolution->comunidades[first]->fim->prox = sol->comunidades[second]->inicio;
    this->basedSolution->comunidades[first]->fim = sol->comunidades[second]->fim;
    this->basedSolution->comunidades[first]->qtd += sol->comunidades[second]->qtd;
    this->basedSolution->comunidades[second]->qtd=0;
    this->basedSolution->comunidades[second]->inicio = sol->comunidades[second]->fim = NULL;

    //storing data about communities
    unsigned nnodes = this->basedSolution->comunidades[second]->qtd;
    this->communities[first].internalEdges  = gcss->nodes[com].internalEdges;
    this->communities[cr].totalDegree = gcss->nodes[com].totalDegree;
    this->communities[first].nnodes += nnodes;
    this->communities[cr].density   = gcss->nodes[com].density;


    //erasing the last elements
    if (gcss->nodes.size() < this->communities.size()){
        unsigned difference = this->communities.size() - gcss->nodes.size();
        vector<FGCommunityLG>::iterator del, it = this->communities.end();
        for (unsigned i=0;i<difference;i++){
            del = it;
            it--;
            this->communities.erase(del);

        }

    }

}
*/

Solution * FastGreedyLG::toSolution(){
    Solution * sol = new Solution(this->graph);
    sol->modularidade = this->bestDensity;
    unsigned * communitiesSol= new unsigned[this->communities.size()];
    for (unsigned com=0; com<this->communities.size();com++){
        communitiesSol[com] = MAX_UNSIGNED;
    }
    unsigned lastC = 0;
    for (unsigned v=0; v<this->nodeCommunity.size();v++){
        if (communitiesSol[this->nodeCommunity[v]]==MAX_UNSIGNED){
            communitiesSol[this->nodeCommunity[v]] = lastC;
            lastC++;
        }
        sol->inserirVertice(v,communitiesSol[this->nodeCommunity[v]]);
    }

    delete communitiesSol;
    return sol;
}
