#ifndef GRAPHLOUVAIN_H
#define GRAPHLOUVAIN_H


#include <limits>
#include <vector>
#include <unordered_set>

#include "../../../graph/largegraph.h"
#include "../../../graph/lde.h"
#include "../../../graph/solution.h"


class EdgeCoarsenerLouvain{
public:
    long double weight;
    unsigned anotherEndpoint; //the another node from the edge
    itemLDE<EdgeCoarsenerLouvain> * apointItem; //the another endpoint item
};

class NodeCoarsenerLouvain{
public:
    long double internalEdges; //number of edges inside cluster
    long double totalDegree; // total degree
    unordered_set<unsigned> * nodes;// nodes inside
    long double density;
    LDE<EdgeCoarsenerLouvain> adjacencies; //adjacencies
    NodeCoarsenerLouvain(){this->nodes=new unordered_set<unsigned>;}
};

class LMBDGraphLouvain
{
public:
    LMBDGraphLouvain(LargeGraph * graph);

    //nodes
    vector<NodeCoarsenerLouvain> nodes;

    //merge Ca with Cb and returns the disabled master node
    unsigned merge(unsigned Ca, unsigned Cb, long double deltaDensity);

    Solution * toSolution();

    long double density;

    LargeGraph * graph;

    unsigned size;

};



#include "lmbdgraphlouvain.cpp"

#endif // GRAPHLOUVAIN_H
