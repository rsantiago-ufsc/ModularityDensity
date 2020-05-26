#ifndef MINIFASTGREEDYBLI_H
#define MINIFASTGREEDYBLI_H

#include "minilouvaincg.h"



class MiniFastGreedyBLI: public MiniLouvainCG
{
public:
    MiniFastGreedyBLI(LargeGraph *graph, IloNum **Wij, float alpha, float lambda=0.5);
    void execute(IloNumArray &lambda);

    float alpha = 0.0;
};


#include "minifastgreedybli.cpp"
#endif // MINIFASTGREEDYBLI_H
