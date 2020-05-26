#ifndef SINGLESTEPMULTILEVELLG_H
#define SINGLESTEPMULTILEVELLG_H

#include "../../graph/largegraph.h"
#include "singlestepgreedylg.h"
class SingleStepMultiLevelLG
{
public:
    LargeGraph * graph;
    long double bestDensity = MINUS_INFINITY;

    SingleStepMultiLevelLG(LargeGraph * graph);

    void execute();
};

#include "singlestepmultilevellg.cpp"

#endif // SINGLESTEPMULTILEVELLG_H
