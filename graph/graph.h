#ifndef Graph_H
#define Graph_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include "vertice.h"
#include <string.h>
#include "../utils/utils.h"

using namespace std;

struct Aresta{
  unsigned v1;
  unsigned v2;
  Aresta(unsigned x, unsigned y): v1(x), v2(y){}
};

class Graph
{
public:
    Graph(string caminho);

    void lerArquivo(string caminho);
    vector<Vertice>  getListaVertices(){return this->listaVertices;}
    int getQtdVertices(){return this->qtd_vertices;}
    int getQtdLigacoes(){return this->qtd_ligacoes;}

    bool dirigido;
    bool ponderado;
    int qtd_vertices;
    int qtd_ligacoes;
    float **matrizAdj; //matriz de adjacência
    float somatorio_pesos;
    vector<Vertice> listaVertices;    

    vector<Aresta> listaArestas;
};

#include "graph.cpp"

#endif // Graph_H
