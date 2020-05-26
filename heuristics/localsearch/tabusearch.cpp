#include "tabusearch.h"

#include <chrono>

using namespace std;
using namespace std::chrono;


TabuSearch::TabuSearch(LargeGraph * graph, ModularityLG * modularity)
{
    //srand(time(NULL));
    this->graph = graph;
    this->modularity = modularity;
    this->best = NULL;
    this->current = NULL;
}


Solution * TabuSearch::execute( Solution * start,    //start solution
                        long double duration,  //duration*|V|
                        unsigned int maxIterations, //[stop criteria] maximum number of iterations
                        unsigned int withoutbestIt, //[stop criteria] number of iterations without find best solutions
                        long double optimalValue //[stop criteria] best known value
                        )
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

        //colocando na lista tabu
        this->tabuInsert(deltaQ, timeLimit, e_node, e_community);

        //atualiza solucao incumbente
        this->current->atualizaVertice(e_node,e_community);
        this->current->garantirUmaVazia();
//cout<<"ts.cpp#82: "<<this->current->modularidade<<"\n";
        this->current->modularidade += e_deltaq;

/*ModularityLG modu(graph);
cout<<"\tts.cpp#87: ("<<it<<") MD = "<<modu.calculateDensity(current)<<" | calculated =  "<<this->current->modularidade<<endl;
if(!(this->current->modularidade <= modu.calculateDensity(current)+0.00001 && this->current->modularidade >= modu.calculateDensity(current)-0.00001)){
    cout<<"\tinfos: size of "<<older_community<<" ("<<current->comunidades[older_community]->qtd<<") and "<<e_community<<" ("<<current->comunidades[e_community]->qtd<<")";
    exit(0);
}*/


//cout<<"ts.cpp#87: if ("<<this->best->modularidade<<" < "<<this->current->modularidade<<"){...}\n";
//if (it == 5){return this->best;}
        if (this->best->modularidade < this->current->modularidade){
            delete this->best; //O(n)
            this->best = this->current->clone();//O(n)
            this->bestIt = this->it;
            this->bestTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

            /*if (this->best->modularidade >= optimalValue){
                this->time = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
                return this->best;
            }*/
            itSemNov = 0;
        }
        deltaQ.changeCommunity(this->current, older_community , e_community, it);
        this->tabuRemove(this->current, deltaQ, it);

        this->it++;
        itSemNov++;

    }
    this->time = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
    return this->best;
}


//colocando na lista tabu
void TabuSearch::tabuInsert(DeltaQTabu &deltaQ, unsigned & timeLimit, unsigned & e_node, unsigned & e_community){
    //insere na lista tabu
    fibonacci_heap<tabu_data>::handle_type handle;
    tabu_data dt(timeLimit, e_node, e_community);
    handle = this->tabuList.push(dt);
    (*handle).handle = handle;

    //desativa das vizinhancas
    deltaQ.disableParTabu(e_node, e_community);
}

//remove da lista tabu
void TabuSearch::tabuRemove(Solution * current, DeltaQTabu &deltaQ, unsigned & it){
    if (this->tabuList.size() > 0 ){
        tabu_data  data = this->tabuList.top();
        while(it  > data.duration){
            this->tabuList.pop();
            //desativa das vizinhancas
            deltaQ.enableParTabu(current, data.node, data.community);
            if (this->tabuList.size() > 0 ){
                data  = this->tabuList.top();
            }else{
                break;
            }

        }
    }
}
