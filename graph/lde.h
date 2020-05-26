#ifndef LDE_H
#define LDE_H

template<class T>
struct itemLDE{
    T id;
    itemLDE<T> * prox;
    itemLDE<T> * ant;
};


template<class T>
class LDE
{
public:

    LDE();
    itemLDE<T> * inserir(T id);
    void colocar(itemLDE<T> * novo);

    bool remover(itemLDE<T> * endereco);
    itemLDE<T> * retirar(itemLDE<T> * endereco);

    itemLDE<T> * getInicio();
    itemLDE<T> * getFim();
    unsigned getQtd();


    ~LDE();
    itemLDE<T> * inicio;
    itemLDE<T> * fim;
    int qtd;

};

#include "lde.cpp"

#endif // LDE_H
