#include "lmbdsinglestepgreedylg.h"
#include "../../utils/modularitylg.h"
using namespace std::chrono;
LMBDSingleStepGreedyLG::LMBDSingleStepGreedyLG(LargeGraph * graph, float lambda)
{
    this->lambda = lambda;
    this->graph = graph;

    //to calculate the modularity
    this->singletonMod = 0.0;
    long double dg;
    for(unsigned com = 0;com < graph->numberOfNodes;com++){
        dg = graph->degreeOfNode[com];
        singletonMod +=  - (dg*dg)/(4.0*graph->numberOfEdges*graph->numberOfEdges);
    }
}

//#include "../../utils/modularidade.h"
void LMBDSingleStepGreedyLG::execute(){



    chrono::system_clock::time_point before;
    before = chrono::system_clock::now();

    fibonacci_heap<DataGCSSLG> * heap = new fibonacci_heap<DataGCSSLG>;
    if (this->gcss != NULL){
        delete this->gcss;
    }
    this->gcss = new LMBDGraphCoarsenerSingleStepLG (this->graph,heap,this->lambda,TYPE_PRI_DENSITY);


    long double modularity = singletonMod;

    long double maxDeltaDensity=MINUS_INFINITY;
    this->bestCommSize = graph->numberOfNodes;
    long double size = graph->numberOfNodes;

    unsigned iter = 0;

    while (heap->empty() == false){


        DataGCSSLG data = heap->top();
        heap->pop();

        size--;
        unsigned Ca = data.node1;
        unsigned Cb = data.node2;
        //cout<<"\nIt: "<<iter<<"  ";
        //cout<<"\n   Ca:"<<Ca<<" / Cb: "<<Cb;
        //data.print();
        long double totalDegree = gcss->nodes[Ca].totalDegree + gcss->nodes[Cb].totalDegree;
        long double totalNodes = gcss->nodes[Ca].nodes.size() + gcss->nodes[Cb].nodes.size();
        long double edgesInside = data.weight + gcss->nodes[Ca].internalEdges + gcss->nodes[Cb].internalEdges;
        //long double density = (4.0*edgesInside - totalDegree)/totalNodes - gcss->nodes[Ca].density - gcss->nodes[Cb].density;
        long double density = (4.0*lambda*edgesInside -(2-2*lambda)*(totalDegree-2*edgesInside))/totalNodes - gcss->nodes[Ca].density - gcss->nodes[Cb].density;
        //cout<<"\ndensity>"<<density<<"\t"<<gcss->nodes[Ca].density<<" \t "<<gcss->nodes[Cb].density;
        modularity += data.weight/(this->graph->numberOfEdges) - (2.0*gcss->nodes[Ca].totalDegree*gcss->nodes[Cb].totalDegree)/(4.0*this->graph->numberOfEdges*this->graph->numberOfEdges);

#ifdef COLLECT_DATA_CONSTRUCTIVE
    //collecting data to other metamethods, like multi-level heuristics
    NodeDetail nca, ncb;

    nca.edges = gcss->nodes[Ca].internalEdges;
    nca.degree = gcss->nodes[Ca].totalDegree;
    nca.nodes = new unordered_set<unsigned>;
    nca.nodes->insert(gcss->nodes[Ca].nodes.begin(),gcss->nodes[Ca].nodes.end());

    ncb.edges = gcss->nodes[Cb].internalEdges;
    ncb.degree = gcss->nodes[Cb].totalDegree;
    ncb.nodes = new unordered_set<unsigned>;
    ncb.nodes->insert(gcss->nodes[Cb].nodes.begin(),gcss->nodes[Cb].nodes.end());

    this->levelDetails.push_back(std::pair<NodeDetail, NodeDetail>(nca, ncb));
#endif



        unsigned Cx = gcss->merge(Ca,Cb,density);
        unsigned Cy = Cb;
        if (Cx == Cb) {Cy = Ca;}
        //levels to construct a solution, or for another heuristic improvement
        this->levels.push_back(pair<pair<unsigned,unsigned>, long double> (pair<unsigned,unsigned> (Cx,Cy), edgesInside));

        //found a best partition
        if (gcss->density > maxDeltaDensity){
            this->bestIt = iter;
            this->bestTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
            maxDeltaDensity = gcss->density;
            bestDensity = maxDeltaDensity;


            this->bestCommSize = size;
            this->bestMod = modularity;

#ifdef GROUND_TRUTH

            this->bestCommunity = this->bestToString();

#endif
        }

        iter++;
    }
    this->maxIt = iter;
    this->totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

    delete heap;

}

Solution * LMBDSingleStepGreedyLG::bestToSolution(){
    if (this->levels.size() == 0){
        return NULL;
    }
    Solution *sol = new Solution(this->graph);

    //each node in the a singleton community
    for (unsigned v=0; v < this->graph->numberOfNodes;v++){
        sol->inserirVertice(v,v);
    }

    //creating the new solution by the levels

    for (unsigned l=0; l<=this->bestIt; l++){
        pair<unsigned,unsigned> p = this->levels[l].first;

        sol->comunidades[p.first]->fim->prox = sol->comunidades[p.second]->inicio;
        sol->comunidades[p.first]->fim = sol->comunidades[p.second]->fim;
        sol->comunidades[p.first]->qtd += sol->comunidades[p.second]->qtd;
        sol->comunidades[p.second]->qtd=0;
        sol->comunidades[p.second]->inicio = sol->comunidades[p.second]->fim = NULL;
    }

    //updating data structures
    for (unsigned com=0; com<sol->comunidades.size(); com++){
        itemLDE<unsigned> * aux = sol->comunidades[com]->inicio;
        while (aux != NULL){
            sol->verticeComunidade[aux->id] = com;
            aux = aux->prox;
        }
    }

    sol->modularidade = this->bestDensity;

    return sol;
}





string LMBDSingleStepGreedyLG::bestToString(){
    string solution = "";
    bool first = true;
    for (unsigned v=0;v<this->gcss->nodes.size();v++){
        if (this->gcss->nodes[v].nodes.size() == 0){continue;}
        if (!first){
            solution += ",";
        }
        first = false;
        string community = "[";
        unordered_set<unsigned>::iterator it = this->gcss->nodes[v].nodes.begin();
        unsigned i=0;
        while(it != this->gcss->nodes[v].nodes.end()){
            community += to_string(*it);
            it++;
            i++;
            if (i<this->gcss->nodes[v].nodes.size()){
                community += ",";
            }
        }
        community += "]";
        solution += community;

    }

    return solution;

}




long double maxDeltaDensity=MINUS_INFINITY;
chrono::system_clock::time_point before;
unsigned iter = 0;
void LMBDSingleStepGreedyLG::start(){

    before = chrono::system_clock::now();

    fibonacci_heap<DataGCSSLG> * heap = new fibonacci_heap<DataGCSSLG>;
    if (this->gcss != NULL){
        delete this->gcss;
    }
    this->gcss = new LMBDGraphCoarsenerSingleStepLG (this->graph,heap,lambda);



}
void LMBDSingleStepGreedyLG::nextLevel(){
//    while (heap->empty() == false){


        DataGCSSLG data = this->gcss->heap->top();
        this->gcss->heap->pop();

        unsigned Ca = data.node1;
        unsigned Cb = data.node2;

        long double totalDegree = gcss->nodes[Ca].totalDegree + gcss->nodes[Cb].totalDegree;
        long double totalNodes = gcss->nodes[Ca].nodes.size() + gcss->nodes[Cb].nodes.size();
        long double edgesInside = data.weight + gcss->nodes[Ca].internalEdges + gcss->nodes[Cb].internalEdges;
        //long double density = (4.0*edgesInside - totalDegree)/totalNodes - gcss->nodes[Ca].density - gcss->nodes[Cb].density;
        long double density = (4.0*lambda*edgesInside -(2-2*lambda)*( totalDegree-2*edgesInside))/totalNodes - gcss->nodes[Ca].density - gcss->nodes[Cb].density;

        unsigned Cx = gcss->merge(Ca,Cb,density);
        unsigned Cy = Cb;
        if (Cx == Cb) {Cy = Ca;}
        //levels to construct a solution, or for another heuristic improvement
        this->levels.push_back(pair<pair<unsigned,unsigned>, long double> (pair<unsigned,unsigned> (Cx,Cy), edgesInside));

        //found a best partition
        if (gcss->density > maxDeltaDensity){
            this->bestIt = iter;
            this->bestTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
            maxDeltaDensity = gcss->density;
            bestDensity = maxDeltaDensity;
        }
        iter++;
//    }
}

void LMBDSingleStepGreedyLG::stop(){
    this->maxIt = iter;
    this->totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
    delete this->gcss->heap;
}


