#include "mcnlnm.h"

#define TYPE_PRI_NOW TYPE_PRI

McnLnm::McnLnm(LargeGraph *graph)
{
    this->graph = graph;
}

void McnLnm::execute(){
    this->hlv = new Louvain(this->graph,TYPE_PRI_NOW);
    hlv->execute();
    this->bestDensity = hlv->bestDensity;

    Solution *solSS = hlv->bestToSolution();

    this->fglg = new FastGreedyLG(solSS);
    fglg->execute();

    if (this->bestDensity < fglg->bestDensity){
        this->bestDensity = fglg->bestDensity;
        this->FastGreedyImproved = true;
    }
}


