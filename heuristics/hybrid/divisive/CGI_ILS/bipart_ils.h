#ifndef BIPART_ILS_H
#define BIPART_ILS_H


#include <unordered_set>
#include <list>
#include <vector>
#include "../../../../graph/largegraph.h"

#include <ilcplex/ilocplex.h>
#include <ilconcert/iloenv.h>



/**
 * @brief The Bipart_ILS class
 *
 * Iterative Local Search to find the best cut for Modularity Density Maximization
 */

struct BipILSNode{
    unsigned internalEdges = 0;
    unsigned degree = 0;
    unordered_set<unsigned> nodes;
    IloNum lambda;
    BipILSNode(unsigned node, unsigned edges, unsigned deg, IloNum lbda):
        internalEdges(edges), degree(deg), lambda(lbda){
        nodes.insert(node);
    }

};


class Bipart_ILS
{
public:

        LargeGraph * graph;
        vector<BipILSNode> coarsedNodes;

        float alpha = 0.0;

        bool terminouZero;

        IloNum **Wij;
        float currentValue = 0.0;
        vector<unsigned> currentNodes; //binary array
        list<unsigned> currentLNodes; //list of nodes that belongs to the first half cluster
        list<unsigned> currentLNodes2; //list of nodes that belongs to the second half cluster


        vector<unsigned> * nodes; //binary array

        float sumLambdaS = 0.0;
        float sumLambdaS2 = 0.0;

        vector<unsigned> a;//solution which is constructed by 'constructSolutionArray' method

        float maxValue = 0.0;
        vector<unsigned> maxNodes;

        Bipart_ILS(vector<unsigned> &nodes, LargeGraph * graph, IloNum **Wij, float alpha);

        void execute(IloNumArray &lambda, IloNum &lambdaY);
        float gainIfIEnter(unsigned node, IloNumArray &lambda, IloNum &lambdaY);
        float gainIfIExit(unsigned node, IloNumArray &lambda, IloNum &lambdaY);
        float calculateAll(IloNumArray &lambda, IloNum &lambdaY);//verify correctness
        void constructSolutionArray(IloNumArray &lambda);
        void shuffleCurrentSolution(IloNumArray &lambda, float strength, IloNum &lambdaY);
        void copySolutionArray(vector<unsigned> &target, IloNumArray &lambda);

        long double calculateD( vector<unsigned> * cluster );
        float calculateAllSimulate(unsigned node, IloNumArray &lambda, IloNum &lambdaY);
        float calculateAllByVector(vector<unsigned> &a, IloNumArray &lambda, IloNum &lambdaY);

};


#include "bipart_ils.cpp"
#endif // BIPART_ILS_H
