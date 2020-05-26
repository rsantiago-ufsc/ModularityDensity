#include "largegraph.h"

#include<vector>

LargeGraph::LargeGraph(string path)
    :directed(false),weighted(false),numberOfNodes(0),numberOfEdges(0),sumOfWeights(0)
{

    this->readFile(path);
}

//método responsável por ler o arquivo que contém o Graph
void LargeGraph::readFile(string path)
{

    ifstream file(path.c_str());
    if(file.is_open()){
        vector<string> toks = Utils::splitString(path,'.');
        string file_ext = toks[toks.size()-1];
        if (file_ext == "net" || file_ext == "paj"){
            string line = "";
            string part;
            bool first = false;
            do{

                getline(file, line);
                part = line.substr(0,8);
                if (Utils::strToUpper(part ) == "*VERTICE"){
                    vector<string> qtd= Utils::splitString(line, ' ');
                    this->numberOfNodes = atoi(qtd[1].c_str());
                    first = true;
                }
            }while(first == false);

            //reading nodes
            bool withoutDetails=false;
            unsigned int startsWith = this->numberOfNodes;
            unsigned int * ids = new unsigned int[this->numberOfNodes];
            string * labels = new string[this->numberOfNodes];
            for (unsigned int i=0; i < this->numberOfNodes; i++){

                getline(file, line);

                if (line.size() >= 4){
                    part = line.substr(0,4);
                    if (Utils::strToUpper(part) == "*ARC" or Utils::strToUpper(part) == "*EDG"){
                        startsWith = 1;
                        withoutDetails = true;
                        for (unsigned int j=0; j < this->numberOfNodes; j++){
                            ids[j] = 1;
                            labels[j] = to_string(ids[j]);
                        }
                        break;
                    }
                }

                vector<string> fields = Utils::splitString(line,' ');
                ids[i] = atoi(fields[0].c_str());
                if (fields.size() > 1){
                    labels[i] = fields[1];
                }else{
                    labels[i] = fields[0];
                }
                if (ids[i]< startsWith){
                    startsWith = ids[i];
                }
            }

            //alocando vertices
            for (unsigned int i=0; i < this->numberOfNodes; i++){
                this->adjNodes.push_back(unordered_map<unsigned, LargeGraphEdge > ());
                this->degreeOfNode.push_back(0.0);
                if (withoutDetails == false){
                    this->labelOfNode.push_back(labels[i]);
                }else{
                    this->labelOfNode.push_back(std::to_string(i+startsWith));
                }
            }

            //to know if the graph is directed
            if (withoutDetails == false){
                getline(file, line);
            }
            part =line.substr(0,4);
            if (Utils::strToUpper(part) == "*ARC"){
                this->directed = true;
            }else{
                if (Utils::strToUpper(part) == "*EDG"){
                    this->directed = false;
                }else{
                    cout<<line<<"\n";
                    cout << "Unsupported file!\n";
                    delete[] ids;
                    delete[] labels;
                    return;
                }
            }
            //reading edges and arcs
            getline(file, line);
            vector<string> fields = Utils::splitString(line,' ');
            if (fields.size()<3){
                this->weighted = false;
            }else{
                this->weighted = true;
            }

            float weightedVerification = 0.0;//faz verificação com pesos exclusivamente 1: no final da rotina, o classifica como não ponderado
            this->numberOfEdges=0;
            unsigned int v1, v2;
            do{

                v1 = atoi(fields[0].c_str())-startsWith;
                v2 = atoi(fields[1].c_str())-startsWith;

                edges.push_back(SimpleEdgeLG(v1,v2));

                LargeGraphEdge edge1(v1,v2,1.0);


                if (this->weighted){
                    this->sumOfWeights += atof(fields[2].c_str());
                    weightedVerification += atof(fields[2].c_str());
                    edge1.weight = atof(fields[2].c_str());
                    this->degreeOfNode[v1]+=atof(fields[2].c_str());
                }else{
                    this->sumOfWeights += 1.0;
                    this->degreeOfNode[v1]+=edge1.weight;
                }

                //if edges already exists, increases the weight
                unordered_map<unsigned,LargeGraphEdge>::iterator it = this->adjNodes[v1].find(v2);
                if (it!=this->adjNodes[v1].end()){
                    it->second.weight += edge1.weight;
                }else{
                    this->adjNodes[v1].insert( pair<unsigned, LargeGraphEdge>(v2, edge1));
                }
                if (this->directed == false){
                    LargeGraphEdge edge2(v2,v1,edge1.weight);
                    unordered_map<unsigned,LargeGraphEdge>::iterator it  = this->adjNodes[v2].find(v1);
                    if (it!=this->adjNodes[v2].end()){
                        it->second.weight += edge2.weight;
                    }else{
                        this->adjNodes[v2].insert(pair<unsigned, LargeGraphEdge>(v1, edge2));
                    }
                    if (v1 != v2){
                        this->degreeOfNode[v2]+=edge2.weight;
                    }
                }
                this->numberOfEdges++;
                getline(file, line);
                fields = Utils::splitString(line,' ');

            }while(fields.size()>1);
            if (weightedVerification == this->numberOfEdges){
                this->weighted = false;
            }
            delete[] ids;
            delete[] labels;

        }else{
            cout << "Unsupported extension!";
        }


    }
}

float LargeGraph::getDegree(unsigned v){
    return this->degreeOfNode[v];
}


float LargeGraph::getAdj(unsigned s, unsigned t){
    unordered_map<unsigned,LargeGraphEdge>::const_iterator it = this->adjNodes[s].find(t);
    if (it != this->adjNodes[s].end()){
        return it->second.weight;
    }else{
        return 0.0;
    }
}
