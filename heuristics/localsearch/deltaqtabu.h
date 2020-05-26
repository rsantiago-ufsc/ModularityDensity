#ifndef DELTAQ_H
#define DELTAQ_H
#include <algorithm>
#include <boost/heap/fibonacci_heap.hpp>
#include <map>

using namespace boost::heap;


#include"../../graph/largegraph.h"
#include"../../utils/modularitylg.h"

struct heap_data
{
    fibonacci_heap<heap_data>::handle_type handle;
    long double deltaq;
    unsigned community;
    unsigned node;

    heap_data(long double i, unsigned v, unsigned c):
        deltaq(i), node(v), community(c)
    {}

    bool operator<(heap_data const & rhs) const
    {
        return deltaq < rhs.deltaq;
    }
};

struct deltaDataValue{
    long double internalEdges;
    long double totalDegree;
    unsigned nnodes;
    long double density;
    deltaDataValue(){}
};

class DeltaQTabu
{

public:
    long double LowerBound;

    DeltaQTabu(LargeGraph * graph, ModularityLG * modularity, Solution * startSol);
    void changeCommunity(Solution * solAfterUpdate, unsigned olderCom, unsigned newCom, unsigned it=0);

    fibonacci_heap<heap_data> *    deltaQGeral;//http://stackoverflow.com/questions/17823495/how-do-i-get-elements-from-boostfibonacci-heap-by-handle
    vector< fibonacci_heap<heap_data>::handle_type * > whoPointsCom;

    //data store for delta calculations
    vector<deltaDataValue> deltaCommunities;

    //for Tabu search
    void enableParTabu(Solution * current, unsigned &node, unsigned &community);
    void disableParTabu(unsigned &node, unsigned &community);

    long double calculateFalseDeltaDOnceNeighborhood(Solution * sol,long double dens, unsigned node, unsigned comTo);


private:
    ModularityLG * modularity;
    LargeGraph * graph;


};


#include "deltaqtabu.cpp"

#endif // DELTAQ_H
