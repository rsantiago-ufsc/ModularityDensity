#ifndef LARGEGRAPH_H
#define LARGEGRAPH_H


#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <string.h>
#include <unordered_map>
#include "../utils/utils.h"

using namespace std;

const long double PLUS_INFINITY = std::numeric_limits<long double>::max();
const long double MINUS_INFINITY = - std::numeric_limits<long double>::max();
const long double MAX_UNSIGNED = std::numeric_limits<unsigned>::max();


struct LargeGraphEdge{
    unsigned sourceNode;
    unsigned targetNode;
    float weight;
    LargeGraphEdge(unsigned s, unsigned t, float w):
        sourceNode(s), targetNode(t),weight (w){}

    bool operator== (LargeGraphEdge const & rhs) const {
        return targetNode == rhs.targetNode;}



};

struct SimpleEdgeLG{
    unsigned v1;
    unsigned v2;
    SimpleEdgeLG(unsigned n1, unsigned n2):v1(n1),v2(n2){}
};

class LargeGraph
{
public:
    LargeGraph(string path);

    void readFile(string path);

    bool        directed; //dirigido;
    bool        weighted; //ponderado;
    unsigned    numberOfNodes; //qtd_vertices;
    unsigned    numberOfEdges;//qtd_ligacoes;
    float       sumOfWeights;

    vector< unordered_map<unsigned, LargeGraphEdge > > adjNodes;
    vector< string > labelOfNode;
    vector< float > degreeOfNode;


    vector< SimpleEdgeLG > edges;

    float getAdj(unsigned s, unsigned t);
    float getDegree(unsigned v);

};

#include "largegraph.cpp"

#endif // LARGEGRAPH_H
