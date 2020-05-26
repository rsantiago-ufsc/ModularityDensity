#ifndef Solution_H
#define Solution_H

#include <vector>
#include <limits>
#include "lde.h"
#include "../utils/utils.h"
#include "largegraph.h"
#include "graph.h"



const unsigned COMUNIDADE_NULA =  std::numeric_limits<unsigned>::max();

class Solution
{
public:
    Solution(Graph * grafo);
    Solution(Graph * grafo, string texto);
    ~Solution();
    Solution(LargeGraph * grafo);
    Solution(LargeGraph * grafo, string texto);

    bool inserirVertice(unsigned id, unsigned comunidade);
    bool removerVertice(unsigned id);
    bool atualizaVertice(unsigned id, unsigned comunidade);
    bool garantirUmaVazia();//garante uma comunidade vazia, para as buscas locais
    string serialize();
    Solution * clone();
    vector<LDE<unsigned> *> * getComunidades(){return &(this->comunidades);}
    itemLDE<unsigned> * getVerticeEndereco(unsigned id){return this->verticeEndereco[id];}
    unsigned getVerticeComunidade(unsigned id){return this->verticeComunidade[id];}
    long double modularidade;

    //retorna em qual comunidade esta o vertice
    vector<unsigned> verticeComunidade;

    //o elemento da lista de vertices do vetor comunidades
    vector<itemLDE<unsigned> *>    verticeEndereco;

    //lista de vertices por comunidade: esta lista é composta por elementos de uma lista duplamente encadeada
    vector<LDE<unsigned> *>        comunidades;
    unsigned haComunidadeVazia = 0;

    Graph * grafo;
    LargeGraph * lggrafo;
};

#include "solution.cpp"

#endif // Solution_H
