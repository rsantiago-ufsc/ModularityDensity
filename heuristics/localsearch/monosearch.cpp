#include "monosearch.h"

#include <chrono>

using namespace std;
using namespace std::chrono;


MonoSearch::MonoSearch(LargeGraph * graph, ModularityLG * modularity)
{
    //srand(time(NULL));
    this->graph = graph;
    this->modularity = modularity;
    this->best = NULL;
    this->current = NULL;
}


Solution * MonoSearch::execute( Solution * start)
{
    chrono::system_clock::time_point before;
    before = chrono::system_clock::now();

    //criando objeto com a variação da modularidade
    DeltaQTabu deltaQ(this->graph,this->modularity, start);


    //variaveis auxiliares
    unsigned e_node, e_community, older_community, timeLimit;
    long double e_deltaq;


    //recebendo solucao inicial
    if (this->best != NULL){
        delete this->best;
    }
    if (this->current != NULL){
        delete this->current;
    }
    //copia a solucao atual e a define como incumbente
    this->best = start->clone();
    this->current  = start->clone();
    this->best->garantirUmaVazia();
    this->current->garantirUmaVazia();

    this->it = 1;
    this->bestIt = 1;
    unsigned itSemNov = 1;
    long double lastMD = this->current->modularidade;
    while (true)
    {


        //encontrar melhor vizinho

        heap_data data  = deltaQ.deltaQGeral->top();
        e_node      = data.node;
        e_deltaq        = data.deltaq;
        e_community    = data.community;
        older_community = this->current->verticeComunidade[e_node];


        //atualiza solucao incumbente
        this->current->atualizaVertice(e_node,e_community);
        this->current->garantirUmaVazia();
        this->current->modularidade += e_deltaq;

        if (lastMD < this->current->modularidade){
            this->bestTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
        }else{
            this->time = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
            return this->current;
        }
        deltaQ.changeCommunity(this->current, older_community , e_community, it);
        lastMD = this->current->modularidade;


        this->it++;

    }
    this->time = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
    return this->current;
}

