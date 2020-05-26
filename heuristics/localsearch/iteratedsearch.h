#ifndef IteratedSearch_H
#define IteratedSearch_H



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
class IteratedSearch
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
    IteratedSearch(LargeGraph * graph, ModularityLG * modularity);

    Solution * execute( Solution * start,    //start solution
                            double strength,
                            long double duration,  //duration*|V|
                            unsigned int maxIterations, //[stop criteria] maximum number of iterations
                            unsigned int withoutbestIt = 0, //[stop criteria] number of iterations without find best solutions
                            long double optimalValue = 1.0 //[stop criteria] best known value
                            );
    void disturbMe(double poder, ModularityLG * modul);

private:
    LargeGraph * graph;
    ModularityLG * modularity;
};

#include "iteratedsearch.cpp"



#endif // IteratedSearch_H
