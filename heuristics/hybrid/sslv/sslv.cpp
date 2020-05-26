/***
 * DO NOT CONTINUE TO CODIFYING THIS
 * This code does not cost a penny, because we will coarsen the nodes two times:
 * one for single step, another for Louvain method. As the first will coarsen the nodes,
 * the Louvain method can only rearranges the already coarsened nodes.
 */
#include "sslv.h"

SSLV::SSLV(LargeGraph * graph)
{
    this->graph = graph;
}



void SSLV::execute(){


    //executing constructor
    SingleStepGreedyLG ss(this->graph);
    ss.execute();
    this->bestMod = ss.bestMod;
    this->bestCommSize = ss.bestCommSize;

    GraphCoarsenerSingleStepLG * gcBest = new GraphCoarsenerSingleStepLG(this->graph, NULL, TYPE_PRI_DENSITY, false);
    for (unsigned l=0;l<=sslg.bestIt;l++){
        pair<unsigned,unsigned> p = sslg.levels[l].first;
        gcBest->merge(p.first, p.second, 0);
    }

    LouvainSSELV lvlg(this->graph, gcBest);

    vector <unsigned> commByVertices (this->graph->numberOfNodes);
    //identifying where is each node
    unordered_set<unsigned>::iterator it;
    for (unsigned i=0;i<this->graph->numberOfNodes; i++){
//        cout<<"\n  ini.i:"<<i;cout.flush();
        if (gcBest->nodes[i].nodes.size()>0){
            it = gcBest->nodes[i].nodes.begin();
            while (it != gcBest->nodes[i].nodes.end()){
                commByVertices[*it] = i;
                it++;
            }
        }
    }

    gcBest->size = this->graph->numberOfNodes;
    this->reorganization(gcBest, commByVertices);
    lvlg.size = gcBest->size;// = 1;



}
