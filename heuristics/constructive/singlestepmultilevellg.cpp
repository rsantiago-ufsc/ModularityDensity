#include "singlestepmultilevellg.h"

#include "fastgreedylg.h"

SingleStepMultiLevelLG::SingleStepMultiLevelLG(LargeGraph * graph)
{
    this->graph=graph;
}

void SingleStepMultiLevelLG::execute()
{
/*    SingleStepGreedyLG sslg(this->graph);
    sslg.start();
//cout<<"\naaa: ";cout.flush();
    FastGreedyLG fglg(sslg.gcss);
//cout<<"\nbbb: ";cout.flush();

//cout<<"\nccc: ";cout.flush();
    while (sslg.gcss->heap->empty() == false){
//        cout<<"\nIt: ";cout.flush();
        sslg.nextLevel();
        if (sslg.bestDensity > this->bestDensity){
            this->bestDensity = sslg.bestDensity;
        }
        fglg.update(sslg.gcss);
        fglg.execute();
        if (fglg.bestDensity > this->bestDensity){
            this->bestDensity = fglg.bestDensity;
        }
    }
    sslg.stop();

//    Solution * soll = sslg.bestToSolution();
//    cout<<"\nSSLG:\t "<<soll->modularidade;

*/


  SingleStepGreedyLG sslg(this->graph);
    sslg.execute();
    if (this->bestDensity < sslg.bestDensity){
        this->bestDensity = sslg.bestDensity;
    }

    if (sslg.levels.size() == 0){
        return;
    }

    //create solution by each level to apply fast greedy local search
    Solution *sol = new Solution(this->graph);
    //each node in the a singleton community
    FGCommunityLG * comms = new FGCommunityLG[this->graph->numberOfNodes];
    for (unsigned v=0; v < this->graph->numberOfNodes;v++){
        sol->inserirVertice(v,v);
        comms[v].totalDegree = this->graph->degreeOfNode[v];
        comms[v].internalEdges = this->graph->getAdj(v,v);
        comms[v].nnodes = 1;
        comms[v].density = comms[v].totalDegree * -1.0;
    }
    FastGreedyLG fglg(sol);
    //creating the new solution by the levels
    for (unsigned l=0; l<=min(sslg.bestIt+1,sslg.maxIt);l++){ //sslg.levels.size()-1; l++){

        pair<unsigned,unsigned> p = sslg.levels[l].first;
        long double betweenFirstAndSec = sslg.levels[l].second;

        //updating data structures
        itemLDE<unsigned> * aux = sol->comunidades[p.second]->inicio;
        while (aux != NULL){
            sol->verticeComunidade[aux->id] = p.first;
            aux = aux->prox;
        }

        //solution adaptation
        sol->comunidades[p.first]->fim->prox = sol->comunidades[p.second]->inicio;
        sol->comunidades[p.first]->fim = sol->comunidades[p.second]->fim;
        sol->comunidades[p.first]->qtd += sol->comunidades[p.second]->qtd;
        sol->comunidades[p.second]->qtd=0;
        sol->comunidades[p.second]->inicio = sol->comunidades[p.second]->fim = NULL;

        //start of update fglg FastGreedyLG fglg(sol);
        for (unsigned v=0;v < this->graph->numberOfNodes; v++){
            fglg.nodeCommunity[v] = sol->verticeComunidade[v];
        }
        // ... p.first

        //tenho que redefinir todos: não posso atualizar, pois é uma nova busca com novos valores.
        comms[p.first].internalEdges =  betweenFirstAndSec;
        comms[p.second].internalEdges = 0.0;
        comms[p.first].totalDegree += comms[p.second].totalDegree;
        comms[p.second].totalDegree = 0.0;
        comms[p.first].nnodes=sol->comunidades[p.first]->qtd;
        comms[p.second].nnodes= 0 ;
        comms[p.first].density = (4.0*comms[p.first].internalEdges-comms[p.first].totalDegree)/comms[p.first].nnodes;
        comms[p.second].density=   0.0;
        if (l >= sslg.bestIt-9){
            for (unsigned v=0; v < this->graph->numberOfNodes;v++){
                fglg.communities[v].internalEdges = comms[v].internalEdges;
                fglg.communities[v].nnodes = comms[v].nnodes;
                fglg.communities[v].totalDegree = comms[v].totalDegree;
                fglg.communities[v].density = comms[v].density;
            }
            //cout<<"\n"<<fglg.communities[p.first].density;
            //cout<<"\n"<<"(4.0*"<<fglg.communities[p.first].internalEdges<<" - "<<fglg.communities[p.first].totalDegree<<" ) / "<<fglg.communities[p.first].nnodes;
            //cout<<"\n"<<sol->serialize();
            //end of update fglg FastGreedyLG fglg(sol);

            fglg.execute(); //Solution * novel = fglg.executeAndConstruct();
            if (this->bestDensity < fglg.bestDensity){
                this->bestDensity = fglg.bestDensity;
            }
            //return;
         }
        //delete novel;
    }

}
