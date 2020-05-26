#ifndef SSLV_H
#define SSLV_H

#include "../../../graph/largegraph.h"
#include "../../constructive/singlestepgreedylg.h"
#include "../../hybrid/louvain/louvain.h"

/***
 * DO NOT CONTINUE TO CODIFYING THIS
 * This code does not cost a penny, because we will coarsen the nodes two times:
 * one for single step, another for Louvain method. As the first will coarsen the nodes,
 * the Louvain method can only rearranges the already coarsened nodes.
 */
class SSLV
{
public:
    LargeGraph * graph;
    long double bestMod=MINUS_INFINITY;
    unsigned bestCommSize;


    SSLV(LargeGraph * graph);

    void execute();
};

#endif // SSLV_H
