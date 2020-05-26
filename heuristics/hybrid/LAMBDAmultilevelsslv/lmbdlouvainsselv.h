#ifndef LOUVAINSSELV_H
#define LOUVAINSSELV_H

#include "../../../graph/largegraph.h"
#include "../../../graph/lmbdgraphcoarsenersinglesteplg.h"

const long double NO_COMMUNITY = std::numeric_limits<unsigned>::max();

const unsigned MLSSLV_DENS = 1;
const unsigned MLSSLV_MOD = 2;



class CommunityLouvainSSELV{
public:

    long double degree;
    long double internalEdges;
    unsigned    nnodes;
    long double density;
    CommunityLouvainSSELV(){}
    CommunityLouvainSSELV (long double d, long double ie, unsigned n, long double dens):
        degree (d),internalEdges(ie), nnodes(n), density(dens){}

};


class LMBDLouvainSSELV
{
public:

    unsigned type = MLSSLV_DENS;
    float lambda = 0.5;

    LMBDLouvainSSELV(LargeGraph * graph, LMBDGraphCoarsenerSingleStepLG * gcss, float lambda=0.5, unsigned type=MLSSLV_DENS);
    void update(LMBDGraphCoarsenerSingleStepLG * gcss);


    LMBDGraphCoarsenerSingleStepLG * gl;

    //original graph
    LargeGraph * graph;

    //original graph
    vector<CommunityLouvainSSELV> communities;

    //current number of nodes (coarsener) from the solution graph
    unsigned size;

    //iterations
    unsigned it;


    //community by node
    vector<unsigned> commByNode;

    //temporary variables for neighborhood one level calculation
    vector<long double> neighWeight;
    vector<unsigned>    neighPos;
    unsigned            neighLast;

    //one level
    inline bool oneLevel();

    //calculate the density:: old modularity()
    long double calculateDensity();
    //calculate the density gain for a node movement from its community to another
    long double calculateDensityGain(unsigned node, unsigned comm, long double dnodecomm, long double wdegree);
    ///long double calculateDensityGain2(unsigned node, unsigned comm, long double dnodecomm, long double wdegree);
    long double calculateModularityGain(unsigned node, unsigned comm, long double dnodecomm, long double wdegree);
    long double calculateModularity();
    //preparing the temporary variables
    void neighsPrepare(unsigned node);

    //removing the node from its community (one level)
    void removeFromCommunity(unsigned node, unsigned comm, long double dnodecomm);
    //inserting the node inside a community (one level)
    void insertInCommunity(unsigned node, unsigned comm, long double dnodecomm);
    //updating the main data structures of the Louvain heuristic
    void updateLouvainGraphAndCommunities(vector<unsigned> & commByVertice);
    //return the string of the current solution
    string currentToString();


    //this method executes the Louvain heuristic
    void execute(vector<unsigned> & commByVertice);

    long double bestDensity = - std::numeric_limits<long double>::max();
    long double bestMod = - std::numeric_limits<long double>::max();
    string bestCommunityStr;
    unsigned bestCommSize = 0;
    long long int totalTime;
    unsigned bestNumberOfCommunities;

    long long int timeOneLevel = 0;
    long long int timeUpdate = 0;

    long long int timeOL1 = 0;
    long long int timeOL2 = 0;
    long long int timeOL3 = 0;
    long long int timeOL4 = 0;
    long long int timeOL5 = 0;
    long long int timeOL6 = 0;
    long long int timeOL7 = 0;


    long long int timeUpTranslate =0;
    long long int timeUpCopy =0;
    long long int timeUpSwap =0;
    long long int timeUpMerge =0;
    long long int timeUpRenumber =0;

};






#include "lmbdlouvainsselv.cpp"

#endif // LOUVAINSSELV_H
