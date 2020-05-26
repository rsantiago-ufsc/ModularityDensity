#ifndef GRAPHCOARSENER_H
#define GRAPHCOARSENER_H

#include <limits>
#include <vector>
#include <unordered_set>

#include "graph.h"
#include "lde.h"
#include "solution.h"

#define GC_IN_ARC 0
#define GC_OUT_ARC 1

const long double PLUS_INFINITY = std::numeric_limits<long double>::max();
const long double MINUS_INFINITY = - std::numeric_limits<long double>::max();


struct EdgeCoarsener{
    long double weight;
    unsigned anotherEndpoint; //the another node from the edge
    itemLDE<EdgeCoarsener> * apointItem; //the another endpoint item
};

struct NodeCoarsener{
    long double internalEdges; //number of edges inside cluster
    long double totalDegree; // total degree
    unordered_set<unsigned> nodes;// nodes inside
    long double density;
    LDE<EdgeCoarsener> adjacencies; //adjacencies
};

class GraphCoarsener
{
public:
    GraphCoarsener(Graph * graph);

    //nodes
    vector<NodeCoarsener> nodes;

    //merge Ca with Cb and returns the disabled master node
    unsigned merge(unsigned Ca, unsigned Cb, long double deltaDensity, unsigned it);



    Solution * toSolution();



    long double density;

private:
    Graph * graph;



};

#include "graphcoarsener.cpp"

#endif // GRAPHCOARSENER_H
