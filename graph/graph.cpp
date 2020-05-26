#include "graph.h"
#include<vector>

Graph::Graph(string caminho)
    :dirigido(false),ponderado(false),qtd_vertices(0),qtd_ligacoes(0),somatorio_pesos(0)
{

    this->lerArquivo(caminho);
}

//método responsável por ler o arquivo que contém o Graph
void Graph::lerArquivo(string caminho)
{

    ifstream arquivo(caminho.c_str());
    if(arquivo.is_open()){
        vector<string> toks = Utils::splitString(caminho,'.');
        string tipo_arquivo = toks[toks.size()-1];
        if (tipo_arquivo == "net" || tipo_arquivo == "paj"){
            string linha = "";
            string parte;
            bool primeira = false;
            do{

                getline(arquivo, linha);
                parte = linha.substr(0,8);
                if (Utils::strToUpper(parte ) == "*VERTICE"){
                    vector<string> qtd= Utils::splitString(linha, ' ');

                    this->qtd_vertices = atoi(qtd[1].c_str());
                    primeira = true;
                }
            }while(primeira == false);



            //lendo vertices
            bool semDetalhes=false;
            unsigned int comecaCom = this->qtd_vertices;
            unsigned int * ids = new unsigned int[this->qtd_vertices];
            string * rotulos = new string[this->qtd_vertices];
            for (unsigned int i=0; i < this->qtd_vertices; i++){
                getline(arquivo, linha);

                if (linha.size() >= 4){
                    parte = linha.substr(0,4);
                    if (Utils::strToUpper(parte) == "*ARC" or Utils::strToUpper(parte) == "*EDG"){
                        comecaCom = 1;
                        semDetalhes = true;
                        for (unsigned int j=0; j < this->qtd_vertices; j++){
                            ids[j] = 1;
                            rotulos[j] = to_string(ids[j]);
                        }
                        break;
                    }
                }

                vector<string> campos = Utils::splitString(linha,' ');
                ids[i] = atoi(campos[0].c_str());
                rotulos[i] = campos[1];
                if (ids[i]< comecaCom){
                    comecaCom = ids[i];
                }
            }

            //alocando vertices
            for (unsigned int i=0; i < this->qtd_vertices; i++){
                this->listaVertices.push_back(Vertice(0));
            }
            //associando valores
            for (unsigned int i=0; i < this->qtd_vertices; i++){
                this->listaVertices[ids[i]-comecaCom].setId(ids[i]-comecaCom);
                this->listaVertices[ids[i]-comecaCom].setRotulo(rotulos[i]);
            }
            //saber se eh dirigido ou nao
            if (semDetalhes == false){
                getline(arquivo, linha);
            }

            parte =linha.substr(0,4);
            if (Utils::strToUpper(parte) == "*ARC"){
                this->dirigido = true;
            }else{
                if (Utils::strToUpper(parte) == "*EDG"){
                    this->dirigido = false;
                }else{
                    cout<<linha<<"\n";
                    cout << "Arquivo não suportado!\n";
                    delete[] ids;
                    delete[] rotulos;
                    return;
                }
            }
            //leitura das arestas ou arcos
            getline(arquivo, linha);
            vector<string> campos = Utils::splitString(linha,' ');
            if (campos.size()<3){
                this->ponderado = false;
            }else{
                this->ponderado = true;
            }
            //declarando matrix de adj;
            //aloca matriz e a inicializa com zeros
            matrizAdj = new float *[qtd_vertices];
            for(int i = 0;i < qtd_vertices; ++i)
                    matrizAdj[i] = new float[qtd_vertices];
            for(int i = 0; i < qtd_vertices; i++){
                for(int j = 0; j < qtd_vertices; j++){
                    matrizAdj[i][j] = 0;
                }
            }
            float verificaPonderado = 0.0;//faz verificação com pesos exclusivamente 1: no final da rotina, o classifica como não ponderado
            this->qtd_ligacoes=0;
            unsigned int v1, v2;
            do{

                v1 = atoi(campos[0].c_str())-comecaCom;
                v2 = atoi(campos[1].c_str())-comecaCom;
                this->listaArestas.push_back(Aresta(v1,v2));

                if (this->ponderado){
                    this->somatorio_pesos += atof(campos[2].c_str());
                    this->matrizAdj[v1][v2] = atof(campos[2].c_str());
                    verificaPonderado += atof(campos[2].c_str());
                }else{
                    this->matrizAdj[v1][v2] = 1.0;
                }


                if (this->dirigido == false){
                    if (this->ponderado){
                        this->matrizAdj[v2][v1] = atof(campos[2].c_str());
                        this->listaVertices[v1].setForca(this->listaVertices[v1].getForca()+atof(campos[2].c_str()));
                    }else{
                        this->matrizAdj[v2][v1] = 1.0;
                    }
                    //graus
                    this->listaVertices[v1].addVerticeAdjacente(&this->listaVertices[v2]);
                }else{
                    if (this->ponderado){
                        this->listaVertices[v1].setForcaSaida(this->listaVertices[v1].getForcaSaida()+atof(campos[2].c_str()));
                        this->listaVertices[v2].setForcaEntrada(this->listaVertices[v1].getForcaEntrada()+atof(campos[2].c_str()));
                    }
                    this->listaVertices[v1].addVerticeAdjacenteSaida(&this->listaVertices[v2]);
                    this->listaVertices[v2].addVerticeAdjacenteEntrada(&this->listaVertices[v1]);
                }
                this->qtd_ligacoes++;
                getline(arquivo, linha);
                campos = Utils::splitString(linha,' ');

            }while(campos.size()>1);
            if (verificaPonderado == this->qtd_ligacoes){
                this->ponderado = false;
            }
            delete[] ids;
            delete[] rotulos;

        }else{
            cout << "Arquivo não suportado!";
        }


    }
}

