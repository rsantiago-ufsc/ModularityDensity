#ifndef VERTICE_H
#define VERTICE_H

#include <vector>
#include <iostream>
using namespace std;

class Vertice
{
public:
    Vertice(unsigned id);
    void addVerticeAdjacente(Vertice *vertice);
    void addVerticeAdjacenteEntrada(Vertice *vertice);
    void addVerticeAdjacenteSaida(Vertice *vertice);
    void setId(int id){this->id = id;}
    void setRotulo(string rotulo){this->rotulo = rotulo;}
    void setForca(float forca){this->forca = forca;}
    void setForcaEntrada(float forca_entrada){this->forca_entrada = forca_entrada;}
    void setForcaSaida(float forca_saida){this->forca_saida = forca_saida;}
    unsigned getId() {return this->id;}
    unsigned getGrau(){return this->grau;}
    unsigned getGrauEntrada(){return this->grau_entrada;}
    unsigned getGrauSaida(){return this->grau_saida;}
    float getForca(){return this->forca;}
    float getForcaEntrada(){return this->forca_entrada;}
    float getForcaSaida(){return this->forca_saida;}
    string getRotulo(){return this->rotulo;}
    vector<Vertice*> getListaAdjacencia(){return this->lista_adjacencia;}
    vector<Vertice*> getListaAdjacenciaEntrada(){return this->lista_adjacencia_entrada;}
    vector<Vertice*> getListaAdjacenciaSaida(){return this->lista_adjacencia_saida;}    
    vector<Vertice*> lista_adjacencia;
    vector<Vertice*> lista_adjacencia_entrada; //contém os vértices adjacentes que apontam para este
    vector<Vertice*> lista_adjacencia_saida; //contém os vértices adjacentes que são apontados por este

private:
    unsigned id;
    unsigned grau;
    unsigned grau_entrada;
    unsigned grau_saida;
    float forca;
    float forca_entrada;
    float forca_saida;
    string rotulo;
};

#include "vertice.cpp"

#endif // VERTICE_H
