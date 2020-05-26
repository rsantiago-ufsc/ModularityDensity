#ifndef MINILOUVAINCG_H
#define MINILOUVAINCG_H

#include <unordered_set>
#include <list>
#include <vector>
#include "../../graph/largegraph.h"

#include <ilcplex/ilocplex.h>
#include <ilconcert/iloenv.h>

struct MLVCGNode{
    unsigned internalEdges = 0;
    unsigned degree = 0;
    unordered_set<unsigned> nodes;
    IloNum lambda;
    MLVCGNode(unsigned node, unsigned edges, unsigned deg, IloNum lbda):
        internalEdges(edges), degree(deg), lambda(lbda){
        nodes.insert(node);
    }

};

class MiniLouvainCG
{
public:
    LargeGraph * graph;
    vector<MLVCGNode> coarsedNodes;
    float lambda = 0.5;
    IloNum **Wij;
    float currentValue = 0.0;
    vector<unsigned> currentNodes; //binary array
    list<unsigned> currentLNodes; //list of nodes that belongs to the cluster

    float sumLambdaS = 0.0;

    vector<unsigned> a;//solution which is constructed by 'constructSolutionArray' method

    float maxValue = 0.0;
    vector<unsigned> maxNodes;



    MiniLouvainCG(LargeGraph * graph, IloNum **Wij, float lambda=0.5);

    void execute(IloNumArray &lambda, float strength=0.0);
    float gainIfIEnter(unsigned node, IloNumArray &lambda);
    float gainIfIExit(unsigned node, IloNumArray &lambda);
    float calculateAll(IloNumArray &lambda);//verify correctness
    void constructSolutionArray(IloNumArray &lambda);
    void shuffleCurrentSolution(IloNumArray &lambda, float strength);
    void copySolutionArray(vector<unsigned> &target, IloNumArray &lambda);
};

#include "minilouvaincg.cpp"

#endif // MINILOUVAINCG_H
