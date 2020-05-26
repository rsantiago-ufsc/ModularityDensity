#ifndef MonoSearch_H
#define MonoSearch_H



#include "../../graph/largegraph.h"
#include "../../graph/solution.h"
#include "../../utils/modularitylg.h"
#include "deltaqtabu.h"


struct tabu_data
{
    fibonacci_heap<tabu_data>::handle_type handle;
    unsigned duration;
    unsigned community;
    unsigned node;



    tabu_data(unsigned dur, unsigned v, unsigned c):
        duration(dur), node(v), community(c)
    {}

    bool operator<(tabu_data const & rhs) const
    {
        //inverti para inverter prioridade
        return duration > rhs.duration;
    }
};


//Busca Tabu
class MonoSearch
{
public:

    Solution * best;
    Solution * current;

    //report
    unsigned it;
    unsigned bestIt;
    long long int time;
    long long int bestTime;

    fibonacci_heap<tabu_data> tabuList;
    MonoSearch(LargeGraph * graph, ModularityLG * modularity);

    Solution * execute( Solution * start );



private:
    LargeGraph * graph;
    ModularityLG * modularity;
};

#include "monosearch.cpp"



#endif // MonoSearch_H
