#ifndef GRAPHCOARSENERSINGLESTEPLG_H
#define GRAPHCOARSENERSINGLESTEPLG_H

#include <limits>
#include <vector>
#include <unordered_set>

#include <boost/heap/fibonacci_heap.hpp>
using namespace boost::heap;

#include "largegraph.h"
#include "lde.h"
#include "solution.h"

#define GC_IN_ARC 0
#define GC_OUT_ARC 1

const unsigned TYPE_PRI_DENSITY = 1;
const unsigned TYPE_PRI_MODGAIN = 2;

#ifdef COLLECT_DATA_CONSTRUCTIVE
struct NodeDetail{
    long double edges;
    long double degree;
    unordered_set<unsigned> * nodes;
};
#endif


//const long double PLUS_INFINITY = std::numeric_limits<long double>::max();
//const long double MINUS_INFINITY = - std::numeric_limits<long double>::max();
//const long double MAX_UNSIGNED = std::numeric_limits<unsigned>::max();

class GraphCoarsenerSingleStepLG;
class EdgeCoarsenerSingleStepLG;
class NodeCoarsenerSingleStepLG;
class DataGCSSLG;

class DataGCSSLG{
public:
    GraphCoarsenerSingleStepLG * gcss;

    unsigned type;
    unsigned node1;
    unsigned node2;
    long double density;
    long double totalRelativeDegree;
    long double weight; //weight of edge between them

    DataGCSSLG(GraphCoarsenerSingleStepLG * g, unsigned v1, unsigned v2, long double dens, long double relDegree, long double w, unsigned t = TYPE_PRI_DENSITY);

    //DataGCSSLG(GraphCoarsenerSingleStepLG * g, unsigned v1, unsigned v2, long double dens, long double relDegree, long double w, long double o1, long double o2, unsigned t = TYPE_PRI_DENSITY);
    bool operator<(DataGCSSLG const & rhs) const;

    void print();

};

class EdgeCoarsenerSingleStepLG{
public:
    long double weight;
    unsigned anotherEndpoint; //the another node from the edge
    itemLDE<EdgeCoarsenerSingleStepLG> * apointItem; //the another endpoint item
    fibonacci_heap<DataGCSSLG>::handle_type heapReference; //reference from heap
};

class NodeCoarsenerSingleStepLG{
public:
    long double internalEdges; //number of edges inside cluster
    long double totalDegree; // total degree
    unordered_set<unsigned> nodes;// nodes inside
    long double density;
    LDE<EdgeCoarsenerSingleStepLG> adjacencies; //adjacencies
};

class GraphCoarsenerSingleStepLG
{
public:
    GraphCoarsenerSingleStepLG(LargeGraph * graph, fibonacci_heap<DataGCSSLG> * heap, unsigned type = TYPE_PRI_DENSITY, bool usesHeap=true);
    unsigned type;

    //nodes
    vector<NodeCoarsenerSingleStepLG> nodes;

    //merge Ca with Cb and returns the disabled master node
    unsigned merge(unsigned Ca, unsigned Cb, long double deltaDensity);

    Solution * toSolution();

    long double density;

#ifdef COLLECT_DATA_CONSTRUCTIVE
    unsigned unmerge(std::pair < NodeDetail, NodeDetail > & level,
                     vector<unsigned> & commByNode, vector<unsigned> & commByVertice,
                     unsigned &iter);
#endif

    //this var is true if the heap is enabled
    bool usesHeap=true;

    LargeGraph * graph;
    fibonacci_heap<DataGCSSLG> * heap;

    unsigned size;
};




DataGCSSLG::DataGCSSLG(GraphCoarsenerSingleStepLG * g, unsigned v1, unsigned v2, long double dens, long double relDegree, long double w, unsigned t):
    gcss(g), node1(v1),node2(v2), density(dens), totalRelativeDegree(relDegree), weight(w), type(t)
{}

//DataGCSSLG::DataGCSSLG(GraphCoarsenerSingleStepLG * g, unsigned v1, unsigned v2, long double dens, long double relDegree, long double w, long double o1, long double o2, unsigned t):
//    gcss(g), node1(v1),node2(v2), density(dens), totalRelativeDegree(relDegree), weight(w), older1(o1), older2(o2), type(t)
//{}

bool DataGCSSLG::operator<(DataGCSSLG const & rhs) const
{
    switch(type){
        case TYPE_PRI_DENSITY:
            if (density < rhs.density){
                return true;
            }else{
                if (density == rhs.density && totalRelativeDegree > rhs.totalRelativeDegree){
                    return true;
                }
            }
            return false;

        case TYPE_PRI_MODGAIN:
            long double older1 = (gcss->nodes[node1].internalEdges)/(gcss->graph->numberOfEdges)-(gcss->nodes[node1].totalDegree) * (gcss->nodes[node1].totalDegree) / (4.0*gcss->graph->numberOfEdges*gcss->graph->numberOfEdges);
            long double older2 = (gcss->nodes[node2].internalEdges)/(gcss->graph->numberOfEdges)-(gcss->nodes[node2].totalDegree) * (gcss->nodes[node2].totalDegree) / (4.0*gcss->graph->numberOfEdges*gcss->graph->numberOfEdges);
            long double mine = (gcss->nodes[node1].internalEdges+gcss->nodes[node2].internalEdges +weight)/(gcss->graph->numberOfEdges)
                    -(gcss->nodes[node1].totalDegree+gcss->nodes[node2].totalDegree) * (gcss->nodes[node1].totalDegree+gcss->nodes[node2].totalDegree) / (4.0*gcss->graph->numberOfEdges*gcss->graph->numberOfEdges)
                    -older1-older2;
            long double other = (rhs.gcss->nodes[rhs.node1].internalEdges+rhs.gcss->nodes[rhs.node2].internalEdges +weight)/(rhs.gcss->graph->numberOfEdges)
                    -(rhs.gcss->nodes[node1].totalDegree+rhs.gcss->nodes[rhs.node2].totalDegree) * (rhs.gcss->nodes[rhs.node1].totalDegree+gcss->nodes[rhs.node2].totalDegree) / (4.0*rhs.gcss->graph->numberOfEdges*rhs.gcss->graph->numberOfEdges)
                    -older1-older2;

            if (mine < other){
                return true;
            }
            return false;
    }
}
void DataGCSSLG::print(){
long double older1 = (gcss->nodes[node1].internalEdges)/(gcss->graph->numberOfEdges)-(gcss->nodes[node1].totalDegree) * (gcss->nodes[node1].totalDegree) / (4.0*gcss->graph->numberOfEdges*gcss->graph->numberOfEdges);
long double older2 = (gcss->nodes[node2].internalEdges)/(gcss->graph->numberOfEdges)-(gcss->nodes[node2].totalDegree) * (gcss->nodes[node2].totalDegree) / (4.0*gcss->graph->numberOfEdges*gcss->graph->numberOfEdges);
long double mine = (gcss->nodes[node1].internalEdges+gcss->nodes[node2].internalEdges +weight)/(gcss->graph->numberOfEdges)
        -(gcss->nodes[node1].totalDegree+gcss->nodes[node2].totalDegree) * (gcss->nodes[node1].totalDegree+gcss->nodes[node2].totalDegree) / (4.0*gcss->graph->numberOfEdges*gcss->graph->numberOfEdges)
        -older1-older2;
    cout<<"\nOlder1: "<<older1;
    cout<<"\nOlder2: "<<older2;
    cout<<"\nDeltaQ: "<<mine;
    cout<<"\n DegreNode1: "<<gcss->nodes[node1].totalDegree;
    cout<<"\n Debaixo: "<<gcss->nodes[node1].internalEdges+gcss->nodes[node2].internalEdges +weight;
    unordered_set<unsigned> n;
    n.insert(gcss->nodes[node1].nodes.begin(), gcss->nodes[node1].nodes.end());
    n.insert(gcss->nodes[node2].nodes.begin(), gcss->nodes[node2].nodes.end());
    unordered_set<unsigned>::iterator x, y;

    long double mod=0.0;

    x= n.begin();
    while(x!= n.end()){
        unsigned dg1 = gcss->graph->getDegree(*x);
        y= n.begin();
        while(y!= n.end()){
            unsigned dg2 = gcss->graph->getDegree(*y);
            long double aij = gcss->graph->getAdj(*x,*y);
            mod += aij - ((long double)(dg1*dg2))/(2.0*gcss->graph->numberOfEdges);
            y++;
        }
        x++;
    }
    mod = mod /(2*gcss->graph->numberOfEdges);
    cout<<"\nMod(real):"<<mod;
    cout<<"\nMod(mine):"<<(gcss->nodes[node1].internalEdges+gcss->nodes[node2].internalEdges +weight)/(gcss->graph->numberOfEdges)
            -(gcss->nodes[node1].totalDegree+gcss->nodes[node2].totalDegree) * (gcss->nodes[node1].totalDegree+gcss->nodes[node2].totalDegree) / (4.0*gcss->graph->numberOfEdges*gcss->graph->numberOfEdges);
}

#include "graphcoarsenersinglesteplg.cpp"



#endif // GRAPHCOARSENERSINGLESTEPLG_H
