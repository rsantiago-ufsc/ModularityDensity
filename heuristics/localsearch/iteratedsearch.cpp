#include "iteratedsearch.h"

#include <chrono>

using namespace std;
using namespace std::chrono;


IteratedSearch::IteratedSearch(LargeGraph * graph, ModularityLG * modularity)
{
    //srand(time(NULL));
    this->graph = graph;
    this->modularity = modularity;
    this->best = NULL;
    this->current = NULL;
}


Solution * IteratedSearch::execute( Solution * start,    //start solution
                        double strength,
                        long double duration,  //duration*|V|
                        unsigned int maxIterations, //[stop criteria] maximum number of iterations
                        unsigned int withoutbestIt, //[stop criteria] number of iterations without find best solutions
                        long double optimalValue //[stop criteria] best known value
                        )
{
    ModularityLG modul(graph);

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
    while (//critérios de parada
           it <= maxIterations
           && (itSemNov <= withoutbestIt || withoutbestIt == 0)
           )
    {

        //tempo de duracao
        if (duration == 0.0){
            double value = rand()%10;
            value /= 10.0;
            timeLimit = it + this->graph->numberOfNodes * value;
        }else{
            timeLimit = it + this->graph->numberOfNodes * duration;
        }


        //encontrar melhor vizinho

//cout<<"ts.cpp#72:"<<deltaQ.deltaQGeral->size()<<"\n";
        heap_data data  = deltaQ.deltaQGeral->top();
        e_node      = data.node;
        e_deltaq        = data.deltaq;
        e_community    = data.community;
        older_community = this->current->verticeComunidade[e_node];
//cout<<"ts.cpp#78: ("<<it<<") "<<e_deltaq<<"  node("<<e_node<<")[ "<< older_community <<"+ "<< e_community <<" ]\n";



        if (e_deltaq > 0.0){
            //atualiza solucao incumbente
            this->current->atualizaVertice(e_node,e_community);
            this->current->garantirUmaVazia();

            this->current->modularidade += e_deltaq;


            if (this->best->modularidade < this->current->modularidade){
                delete this->best; //O(n)
                this->best = this->current->clone();//O(n)
                this->bestIt = this->it;
                this->bestTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

                itSemNov = 0;
            }
            deltaQ.changeCommunity(this->current, older_community , e_community, it);
        }else{
            this->disturbMe(strength,this->modularity);
            deltaQ = DeltaQTabu(this->graph,this->modularity, current);
        }


        this->it++;
        itSemNov++;

    }
    this->time = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
    return this->best;
}


void IteratedSearch::disturbMe(double poder, ModularityLG * modul){

//cout<<"\nANTES DO DEMO ";
    unsigned qtd = poder*this->graph->numberOfNodes;
 //cout<<"\nqtd: " << qtd ;cout.flush();
    for (unsigned i=0; i<qtd; i++){
        //cout<<"\ni: " << i ;cout.flush();
        unsigned nodo = rand()%this->graph->numberOfNodes;
        //cout<<"\n" << nodo <<"  --  "<<this->atual->comunidades.size();cout.flush();
        unsigned comu = rand()%(this->current->comunidades.size());
        //cout<<"\n" << comu;cout.flush();
        while (comu == this->current->verticeComunidade[nodo]){
            comu = rand()%(this->current->comunidades.size());
        }
        this->current->atualizaVertice(nodo,comu);


    }
    this->current->modularidade = modul->calculateDensity(this->current, 0.5);
}
