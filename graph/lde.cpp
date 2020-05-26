#include "lde.h"

//O(1)
template<class T>
LDE<T>::LDE()
{
    this->inicio=NULL;
    this->fim=NULL;
    this->qtd = 0;
}

//O(1)
template<class T>
itemLDE<T> * LDE<T>::inserir(T id){
    itemLDE<T> * novo = new itemLDE<T>;
    novo->id = id;
    novo->prox = NULL;
    novo->ant = this->fim;

    if (this->inicio == NULL){
        this->inicio = this->fim = novo;
    }else{
        this->fim->prox = novo;
        this->fim = novo;
    }
    this->qtd++;
    return novo;
}

//O(1)
template<class T>
bool LDE<T>::remover(itemLDE<T> * endereco){
    itemLDE<T> * aremover = endereco;
    if (endereco  == this->inicio){
        this->inicio = this->inicio->prox;
    }
    if (endereco == this->fim){
        this->fim = this->fim->ant;
    }
    if (endereco->ant != NULL){
        endereco->ant->prox = endereco->prox;
    }
    if (endereco->prox != NULL){
        endereco->prox->ant = endereco->ant;
    }
    this->qtd--;
    delete aremover;
    return true;
}


template<class T>
itemLDE<T> * LDE<T>::retirar(itemLDE<T> * endereco){
    itemLDE<T> * aremover = endereco;
    if (endereco  == this->inicio){
        this->inicio = this->inicio->prox;
    }
    if (endereco == this->fim){
        this->fim = this->fim->ant;
    }
    if (endereco->ant != NULL){
        endereco->ant->prox = endereco->prox;
    }
    if (endereco->prox != NULL){
        endereco->prox->ant = endereco->ant;
    }
    this->qtd--;
    return aremover;
    //return true;
}

template<class T>
void LDE<T>::colocar(itemLDE<T> * novo){
    novo->prox = NULL;
    novo->ant = this->fim;

    if (this->inicio == NULL){
        this->inicio = this->fim = novo;
    }else{
        this->fim->prox = novo;
        this->fim = novo;
    }
    this->qtd++;
}

//O(n)
template<class T>
LDE<T>::~LDE()
{
    itemLDE<T> * aux = this->inicio, *aremover;
    while(aux != NULL){
        aremover = aux;
        aux = aux->prox;
        delete aremover;
    }
}

template<class T>
itemLDE<T> * LDE<T>::getInicio(){
    return this->inicio;
}

template<class T>
itemLDE<T> * LDE<T>::getFim(){
    return this->fim;
}

template<class T>
unsigned LDE<T>::getQtd(){
    return this->qtd;
}
