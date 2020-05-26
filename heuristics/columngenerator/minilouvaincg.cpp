#include "minilouvaincg.h"

#include <algorithm>

using namespace std;

MiniLouvainCG::MiniLouvainCG(LargeGraph *graph, IloNum **Wij, float lambda)
{
    this->lambda = lambda;
    this->graph = graph;
    this->Wij = Wij;
    for (unsigned i=0;i<graph->numberOfNodes;i++){
        this->maxNodes.push_back(0);
        this->currentNodes.push_back(0);
        a.push_back(0);
    }
}

void MiniLouvainCG::execute(IloNumArray &lambda, float strength){
    this->coarsedNodes.clear();
    this->currentLNodes.clear();
    ////starting a new search
    this->currentValue = 0.0;
    this->maxValue = 0.0;

    //... all nodes are not in the new cluster (new column)
    for (unsigned i=0;i<graph->numberOfNodes;i++){
        this->maxNodes[i]=0;
        this->currentNodes[i]=0;

    }
    //... start the coarsened nodes that are not coarsened yeat
    this->coarsedNodes.clear();
    for (unsigned i=0;i<graph->numberOfNodes;i++){
        unsigned deg = this->graph->degreeOfNode[i];
        this->coarsedNodes.push_back(MLVCGNode(i, this->Wij[i][i], deg, lambda[i]));
    }

    this->sumLambdaS = 0.0;

    //start random solution
    shuffleCurrentSolution(lambda, 0.5);

    ////mini louvain method
    vector<unsigned> randomOrder(this->coarsedNodes.size());

//vector<float> TESTE(this->coarsedNodes.size()); //##

    unsigned notImproved = 0;

    while( notImproved < 2){//graph->numberOfNodes){
//        cout<<"| ";//##
        ////level phase
        int notImprovedLvL = 0;
        while (notImprovedLvL < 2){//graph->numberOfNodes){
            for (int i=0 ; i<this->coarsedNodes.size() ; i++){
                randomOrder[i]=i;
            }
            for (int i=0 ; i<this->coarsedNodes.size()-1 ; i++) {
                int randpos = rand()%(this->coarsedNodes.size()-i)+i;
                swap(randomOrder[i], randomOrder[randpos]);
            }

            for(unsigned inode = 0; inode < this->coarsedNodes.size(); inode++){
                unsigned node = randomOrder[inode];
                float value;
                if (this->currentNodes[node] == 0){
                    value = this->gainIfIEnter(node, lambda);

/*if(notImproved == 0){
    TESTE[node] = value; //##
}*/

                    if(value >0.0000001){
                        this->currentNodes[node] = 1;
                        //an improvement happens
                        notImprovedLvL = -1;
                        notImproved = -1;
                        this->currentValue += value;
                        //each subnode belongs to the solution
                        unordered_set<unsigned>::iterator itnews = this->coarsedNodes[node].nodes.begin();
                        while (itnews != this->coarsedNodes[node].nodes.end()){
                            this->currentLNodes.push_back(*itnews);
                            itnews++;
                        }
                        this->sumLambdaS += this->coarsedNodes[node].lambda;
                    }
                }else{
                    value = this->gainIfIExit(node, lambda);
                    if(value >0.0000001){
                        this->currentNodes[node] = 0;
                        //an improvement happens
                        notImprovedLvL = -1;
                        notImproved = -1;
                        this->currentValue += value;
                        //each subnode belongs to the solution
                        unordered_set<unsigned>::iterator itnews = this->coarsedNodes[node].nodes.begin();
                        while (itnews != this->coarsedNodes[node].nodes.end()){
                            this->currentLNodes.remove(*itnews);
                            itnews++;
                        }
                        this->sumLambdaS -= this->coarsedNodes[node].lambda;
                    }
                }
                notImprovedLvL++;
                node++;
                if (this->currentValue > this->maxValue){
                    this->maxValue = this->currentValue;
                    this->copySolutionArray(this->maxNodes, lambda);
                }

            }
            /*if (this->currentValue > this->maxValue){
                this->maxValue =

            }*/

            notImproved++;


            ////coarsening phase

        }
    }

/*cout<<"\n";
for (unsigned i=0; i< TESTE.size();i++){
    cout<<"["<<i<<"] "<<TESTE[i]<<", ";
}*/

}




float MiniLouvainCG::gainIfIEnter(unsigned node, IloNumArray &lambda){
    //first part
    float value = this->sumLambdaS * this->coarsedNodes[node].nodes.size()
            +  (this->currentLNodes.size()+this->coarsedNodes[node].nodes.size())
            *this->coarsedNodes[node].lambda;


    //second part
    //value -= 4.0*this->coarsedNodes[node].internalEdges;
    unordered_set<unsigned>::iterator ita = this->coarsedNodes[node].nodes.begin();
    while (ita != this->coarsedNodes[node].nodes.end()){
        list<unsigned>::iterator itb = this->currentLNodes.begin();
        while(itb != this->currentLNodes.end()){
            value -= 4.0*this->Wij[*itb][*ita];
            itb++;
        }
        ita++;
    }

    //third part
    value -= this->coarsedNodes[node].degree*(-2.0+2.0*this->lambda);

    return /*max*/ value *-1.0;
}

float MiniLouvainCG::gainIfIExit(unsigned node, IloNumArray &lambda){
    //first part
    float value = 0.0;
    value -= (this->sumLambdaS-this->coarsedNodes[node].lambda) * this->coarsedNodes[node].nodes.size()
                +  (this->currentLNodes.size())*this->coarsedNodes[node].lambda;
                ;

        //second part
        //value -= 4.0*this->coarsedNodes[node].internalEdges;
        unordered_set<unsigned>::iterator ita = this->coarsedNodes[node].nodes.begin();
        while (ita != this->coarsedNodes[node].nodes.end()){
            list<unsigned>::iterator itb = this->currentLNodes.begin();
            while(itb != this->currentLNodes.end()){
                value += 4.0*this->Wij[*itb][*ita];
                itb++;
            }
            ita++;
        }

        //third part
        value += this->coarsedNodes[node].degree*(-2.0+2.0*this->lambda);

    return /*max*/ value *-1.0;
}


float MiniLouvainCG::calculateAll(IloNumArray &lambda){

/*vector<float> TESTE(this->coarsedNodes.size()); //##
for (unsigned i=0; i< TESTE.size();i++){
    TESTE[i] = 0.0;
}*/


    vector<unsigned> a(this->graph->numberOfNodes);
    for (unsigned i=0;i<this->coarsedNodes.size();i++){
        unordered_set<unsigned>::iterator ita = this->coarsedNodes[i].nodes.begin();
        while (ita != this->coarsedNodes[i].nodes.end()){
            a[*ita] = this->currentNodes[i];
            ita++;
        }
    }
    float value = 0.0;
    //first part
    for (unsigned u=0;u<this->graph->numberOfNodes; u++){
        for (unsigned v=0;v<this->graph->numberOfNodes; v++){
            value += a[u] * a[v] * lambda[v];
//TESTE[v] += a[u] * a[v] * lambda[v];
//TESTE[u] += a[u] * a[v] * lambda[v];
        }
    }
    //second part
    for (unsigned v=0;v<this->graph->numberOfNodes; v++){
        for (unsigned u=v+1;u<this->graph->numberOfNodes; u++){
            value -= 4.0*a[u] * a[v] * this->Wij[u][v];
//TESTE[v] -=  4.0*a[u] * a[v] * this->Wij[u][v];
        }
    }
    //third part
    for (unsigned v=0;v<this->graph->numberOfNodes; v++){
        value -= a[v]*this->graph->degreeOfNode[v]*(-2.0+2.0*this->lambda);
//TESTE[v]+= a[v]*this->graph->degreeOfNode[v];
    }

    cout<<"\n";
    for (unsigned v=0;v<this->graph->numberOfNodes; v++){
        cout<<a[v]<<", ";
    }

/*cout<<"\n";
for (unsigned i=0; i< TESTE.size();i++){
    cout<<"["<<i<<"] "<<TESTE[i]<<", ";
}*/

    return /*max*/ value * -1;

}

void MiniLouvainCG::constructSolutionArray(IloNumArray &lambda){
    for (unsigned i=0;i<this->coarsedNodes.size();i++){
        unordered_set<unsigned>::iterator ita = this->coarsedNodes[i].nodes.begin();
        while (ita != this->coarsedNodes[i].nodes.end()){
            this->a[*ita] = this->currentNodes[i];
            ita++;
        }
    }
}

void MiniLouvainCG::shuffleCurrentSolution(IloNumArray &lambda, float strength){
    unsigned nchanges = this->graph->numberOfNodes * strength;
    float value;
    for (unsigned ic=0;ic<nchanges;ic++){
        unsigned node = rand()%this->coarsedNodes.size();
        if (this->currentNodes[node] == 0){
            value = this->gainIfIEnter(node, lambda);
            this->currentNodes[node] = 1;
            //an improvement happens
            this->currentValue += value;
            //each subnode belongs to the solution
            unordered_set<unsigned>::iterator itnews = this->coarsedNodes[node].nodes.begin();
            while (itnews != this->coarsedNodes[node].nodes.end()){
                this->currentLNodes.push_back(*itnews);
                itnews++;
            }
            this->sumLambdaS += this->coarsedNodes[node].lambda;

        }else{
            value = this->gainIfIExit(node, lambda);
            this->currentNodes[node] = 0;
            //an improvement happens
            this->currentValue += value;
            //each subnode belongs to the solution
            unordered_set<unsigned>::iterator itnews = this->coarsedNodes[node].nodes.begin();
            while (itnews != this->coarsedNodes[node].nodes.end()){
                this->currentLNodes.remove(*itnews);
                itnews++;
            }
            this->sumLambdaS -= this->coarsedNodes[node].lambda;
        }

    }

}

void MiniLouvainCG::copySolutionArray(vector<unsigned> &target, IloNumArray &lambda){
    for (unsigned i=0;i<this->coarsedNodes.size();i++){
        unordered_set<unsigned>::iterator ita = this->coarsedNodes[i].nodes.begin();
        while (ita != this->coarsedNodes[i].nodes.end()){
            target[*ita] = this->currentNodes[i];
            ita++;
        }
    }
}
