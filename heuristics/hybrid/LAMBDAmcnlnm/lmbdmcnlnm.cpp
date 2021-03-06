#include "lmbdmcnlnm.h"

#define TYPE_PRI_NOW TYPE_PRI

LMBDMcnLnm::LMBDMcnLnm(LargeGraph *graph, float lambda)
{
    this->lambda = lambda;
    this->graph = graph;
}

void LMBDMcnLnm::execute(){
    this->hlv = new LMBDLouvain(this->graph,TYPE_PRI_NOW);
    this->hlv->lambda = lambda;
    hlv->execute();
    this->bestDensity = hlv->bestDensity;

    Solution *solSS = hlv->bestToSolution();

    this->fglg = new LMBDFastGreedyLG(solSS, lambda);
    fglg->execute();

    if (this->bestDensity < fglg->bestDensity){
        this->bestDensity = fglg->bestDensity;
        this->FastGreedyImproved = true;
    }
}


