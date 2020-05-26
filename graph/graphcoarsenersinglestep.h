#ifndef GRAPHCOARSENERSINGLESTEP_H
#define GRAPHCOARSENERSINGLESTEP_H

#include <limits>
#include <vector>
#include <unordered_set>

#include <boost/heap/fibonacci_heap.hpp>
using namespace boost::heap;

#include "graph.h"
#include "lde.h"
#include "solution.h"

#define GC_IN_ARC 0
#define GC_OUT_ARC 1

const long double PLUS_INFINITY = std::numeric_limits<long double>::max();
const long double MINUS_INFINITY = - std::numeric_limits<long double>::max();
const long double MAX_UNSIGNED = std::numeric_limits<unsigned>::max();

struct DataGCSS{
    unsigned node1;
    unsigned node2;
    long double density;
    long double totalRelativeDegree;
    long double weight; //weight of edge between them

    DataGCSS(unsigned v1, unsigned v2, long double dens, long double relDegree, long double w ):
        node1(v1),node2(v2), density(dens), totalRelativeDegree(relDegree), weight(w)
    {}

    bool operator<(DataGCSS const & rhs) const
    {
        if (density < rhs.density){
            return true;
        }else{
            if (density == rhs.density && totalRelativeDegree > rhs.totalRelativeDegree){
                return true;
            }
        }
        return false;
    }

};

struct EdgeCoarsenerSingleStep{
    long double weight;
    unsigned anotherEndpoint; //the another node from the edge
    itemLDE<EdgeCoarsenerSingleStep> * apointItem; //the another endpoint item
    fibonacci_heap<DataGCSS>::handle_type heapReference; //reference from heap
};

struct NodeCoarsenerSingleStep{
    long double internalEdges; //number of edges inside cluster
    long double totalDegree; // total degree
    unordered_set<unsigned> nodes;// nodes inside
    long double density;
    LDE<EdgeCoarsenerSingleStep> adjacencies; //adjacencies
};

class GraphCoarsenerSingleStep
{
public:
    GraphCoarsenerSingleStep(Graph * graph, fibonacci_heap<DataGCSS> * heap);

    //nodes
    vector<NodeCoarsenerSingleStep> nodes;

    //merge Ca with Cb and returns the disabled master node
    unsigned merge(unsigned Ca, unsigned Cb, long double deltaDensity);

    Solution * toSolution();

    long double density;


    Graph * graph;
    fibonacci_heap<DataGCSS> * heap;


};


#include "graphcoarsenersinglestep.cpp"

#endif // GRAPHCOARSENERSINGLESTEP_H
