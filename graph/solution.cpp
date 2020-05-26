#include "solution.h"



//O(n)
Solution::Solution(Graph * grafo)
{
    this->grafo = grafo;
    this->lggrafo = NULL;
    for (unsigned i=0; i<grafo->getQtdVertices(); i++){
        this->verticeComunidade.push_back(COMUNIDADE_NULA );
        this->verticeEndereco.push_back(NULL);
    }
}


Solution::Solution(LargeGraph * grafo)
{
    this->grafo = NULL;
    this->lggrafo = grafo;
    for (unsigned i=0; i<grafo->numberOfNodes; i++){
        this->verticeComunidade.push_back(COMUNIDADE_NULA );
        this->verticeEndereco.push_back(NULL);
    }
}


//O(n)
Solution::Solution(Graph * grafo, string texto){
    this->grafo = grafo;
    this->lggrafo = NULL;
    for (unsigned i=0; i<grafo->getQtdVertices(); i++){
        this->verticeComunidade.push_back(COMUNIDADE_NULA );
        this->verticeEndereco.push_back(NULL);
    }
    vector<string> txtComunidades = Utils::splitString(texto,']');
    for (unsigned ic=0; ic<txtComunidades.size();ic++){
        this->comunidades.push_back(new LDE<unsigned>());
        vector<string> txtVertices = Utils::splitString(txtComunidades[ic],']');
        for(unsigned iv=0; iv<txtVertices.size();iv++ ){
            Utils::replace(txtVertices [iv], "]", "");
            Utils::replace(txtVertices [iv], "[", "");
            Utils::replace(txtVertices [iv], ",", "");
            if (txtVertices [iv].size() > 0){
                this->inserirVertice((unsigned) atoi(txtVertices[iv].c_str()),ic);
            }
        }
    }
}

Solution::Solution(LargeGraph * grafo, string texto){
    this->grafo = NULL;
    this->lggrafo = grafo;
    for (unsigned i=0; i<grafo->numberOfNodes; i++){
        this->verticeComunidade.push_back(COMUNIDADE_NULA );
        this->verticeEndereco.push_back(NULL);
    }
    vector<string> txtComunidades = Utils::splitString(texto,']');
    for (unsigned ic=0; ic<txtComunidades.size();ic++){
        this->comunidades.push_back(new LDE<unsigned>());
        vector<string> txtVertices = Utils::splitString(txtComunidades[ic],',');

        for(unsigned iv=0; iv<txtVertices.size();iv++ ){
            Utils::replaceAll(txtVertices [iv], "]", "");
            Utils::replaceAll(txtVertices [iv], "[", "");
            Utils::replaceAll(txtVertices [iv], ",", "");
            if (txtVertices [iv].size() > 0){
                this->inserirVertice((unsigned) atoi(txtVertices[iv].c_str()),ic);
            }
        }
    }
}




//O(1)
bool Solution::inserirVertice(unsigned id, unsigned comunidade){
    bool comunidadeValida = false;
    if (comunidade < this->comunidades.size()){
        //verifica se comunidade já existe
        comunidadeValida = true;
    }else{
        //se a comunidade não existir, então posso criar apenas uma nova
        if (comunidade == this->comunidades.size()){
            this->comunidades.push_back(new LDE<unsigned>());
            comunidadeValida = true;
        }
    }
    bool idValido = false;
    if (comunidadeValida){
        //soh pode inserir vertices nao inseridos
        if (this->verticeComunidade[id] == COMUNIDADE_NULA){
            idValido = true;
            this->verticeComunidade[id] = comunidade;
            this->verticeEndereco[id] = this->comunidades[comunidade]->inserir(id);

            if (this->haComunidadeVazia > 0
                && this->comunidades[comunidade]->getQtd() == 1){
                this->haComunidadeVazia--;
            }
        }
    }

    return idValido;
}

//O(1)
bool Solution::removerVertice(unsigned id){
    //para remover ele deve já ter sido inserido
    if (this->verticeComunidade[id] != COMUNIDADE_NULA){
        unsigned comunidade = this->verticeComunidade[id];
        this->verticeComunidade[id] = COMUNIDADE_NULA;
        this->comunidades[comunidade]->remover(this->verticeEndereco[id]);
        this->verticeEndereco[id] = NULL;

        if (this->comunidades[comunidade]->getQtd() == 0){
            this->haComunidadeVazia++;
        }

        return true;
    }else{
        return false;
    }
}


//O(1)
bool Solution::atualizaVertice(unsigned id, unsigned comunidade){
    bool remover = false, inserir = false;
    remover = this->removerVertice(id);
    if (remover){
        inserir = this->inserirVertice(id,comunidade);
    }
    return inserir && remover;
}



//O(n)
string Solution::serialize(){

    string res = "";
    string com;

    for (int ic=0;ic<this->comunidades.size();ic++){
        itemLDE<unsigned> * aux = this->comunidades[ic]->getInicio();
        com = "";
        while(aux!=NULL){
            com += std::to_string(aux->id)+",";
            aux = aux->prox;
        }
        com = com.substr(0, com.size()-1);
        res = res + "],[" + com;
    }
    res = res.substr(2, res.size()-2) + "]";
    return res;
}

//O(n)
Solution * Solution::clone(){
    Solution * cl;
    if (this->grafo != NULL){
        cl = new Solution(this->grafo);
        for (int ic=0;ic<this->comunidades.size();ic++){
            cl->comunidades.push_back(new LDE<unsigned>());
            itemLDE<unsigned> * aux = this->comunidades[ic]->getInicio();
            while(aux!=NULL){
                cl->inserirVertice(aux->id, ic);
                aux = aux->prox;
            }
        }
        cl->modularidade = this->modularidade;
        cl->haComunidadeVazia = this->haComunidadeVazia;
    }else{
        cl = new Solution(this->lggrafo);
        for (int ic=0;ic<this->comunidades.size();ic++){
            cl->comunidades.push_back(new LDE<unsigned>());
            itemLDE<unsigned> * aux = this->comunidades[ic]->getInicio();
            while(aux!=NULL){
                cl->inserirVertice(aux->id, ic);
                aux = aux->prox;
            }
        }
        cl->modularidade = this->modularidade;
        cl->haComunidadeVazia = this->haComunidadeVazia;
    }
    return cl;
}



//destrutor
Solution::~Solution()
{
    itemLDE<unsigned> * aux, *save;
    for (unsigned i=0; i<this->comunidades.size(); i++){
        aux = this->comunidades[i]->getInicio();
        while (aux != NULL){
            save = aux;
            aux = aux->prox;
            delete save;
        }
    }
}

//O(1)
bool Solution::garantirUmaVazia(){
    if (this->haComunidadeVazia == 0){
        //criar comunidade vazia para buscas locais
        this->comunidades.push_back(new LDE<unsigned>());
        this->haComunidadeVazia++;
        return true;
    }else{
        return false;
    }

}
