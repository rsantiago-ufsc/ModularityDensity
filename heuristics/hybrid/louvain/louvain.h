#ifndef LOUVAIN_H
#define LOUVAIN_H

#include <vector>
#include <limits>

#include "../../../graph/largegraph.h"
#include "graphlouvain.h"

#define TYPE_PRI 2 //louvain priority


//type of prioritizers
#define LV_PRI_MOD_DELTA 1//delta mod.
#define LV_PRI_DNS_DELTA 2//delta density
#define LV_PRI_SAN_DELTA 3//delta made by me: real density delta
#define LV_PRI_SAN       4//made by me: real density //##parece ter erros

const long double NO_COMMUNITY = std::numeric_limits<unsigned>::max();


class CommunityLouvain{
public:

    long double degree;
    long double internalEdges;
    unsigned    nnodes;
    long double density;
    CommunityLouvain (long double d, long double ie, unsigned n, long double dens):
        degree (d),internalEdges(ie), nnodes(n), density(dens){}
};


class Louvain
{
public:
    unsigned typePrioritizer;

    Louvain(LargeGraph * graph, unsigned type = LV_PRI_MOD_DELTA);
    Louvain(Solution * sol, LargeGraph * graph, unsigned type = LV_PRI_MOD_DELTA);



    GraphLouvain * gl;

    //original graph
    LargeGraph * graph;

    //original graph
    vector<CommunityLouvain> communities;

    //current number of nodes (coarsener) from the solution graph
    unsigned size;

    //iterations
    unsigned it;


    //community by node
    vector<unsigned> commByNode;
    vector<unsigned> commByNodeAux; //auxiliar

    //temporary variables for neighborhood one level calculation
    vector<long double> neighWeight;
    vector<unsigned>    neighPos;
    unsigned            neighLast;

    //one level
    bool oneLevel();

    //calculate the density:: old modularity()
    long double calculateDensity();
    long double calculateModularity();
    //calculate the density gain for a node movement from its community to another
    long double calculateDensityGain(unsigned node, unsigned comm, long double dnodecomm, long double wdegree);
    long double calculateDensityGain2(unsigned node, unsigned comm, long double dnodecomm, long double wdegree);

    //return the string of the current solution
    string currentToString();

    //preparing the temporary variables
    void neighsPrepare(unsigned node);

    //removing the node from its community (one level)
    void removeFromCommunity(unsigned node, unsigned comm, long double dnodecomm);
    //inserting the node inside a community (one level)
    void insertInCommunity(unsigned node, unsigned comm, long double dnodecomm);
    //updating the main data structures of the Louvain heuristic
    void updateLouvainGraphAndCommunities();


    //this method executes the Louvain heuristic
    void execute();
    //this method returns the best solution
    Solution * bestToSolution();

    long double bestDensity = - std::numeric_limits<long double>::max();
    long double bestMod = - std::numeric_limits<long double>::max();
    unsigned bestCommSize = 0;
    unsigned bestNumberOfCommunities;
    long long int totalTime;
    string bestCommunityStr;

};




#include "louvain.cpp"

#endif // LOUVAIN_H
