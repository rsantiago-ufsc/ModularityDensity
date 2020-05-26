#include "vertice.h"

Vertice::Vertice(unsigned id)
    :id(id), rotulo(""), grau(0), grau_entrada(0), grau_saida(0), forca(0), forca_entrada(0), forca_saida(0)
{
}

void Vertice::addVerticeAdjacente(Vertice *vertice)
{
    this->lista_adjacencia.push_back(vertice);
    vertice->lista_adjacencia.push_back(this);
    this->grau++;
    vertice->grau++;
}

void Vertice::addVerticeAdjacenteEntrada(Vertice *vertice)
{
    this->lista_adjacencia_entrada.push_back(vertice);
    this->grau_entrada++;
}

void Vertice::addVerticeAdjacenteSaida(Vertice *vertice)
{
    this->lista_adjacencia_saida.push_back(vertice);
    this->grau_saida++;
}
