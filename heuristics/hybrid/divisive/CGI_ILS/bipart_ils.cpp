#include "bipart_ils.h"

#include <math.h>

long double Bipart_ILS::calculateD( vector<unsigned> * cluster ){
    long double value = 0.0;
    long double size  = cluster->size();

    long double edges=0.0, degrees = 0.0;
    for (unsigned u=0; u < size; u++){
        degrees += this->graph->getDegree( (*cluster)[u] );
        for (unsigned v=u+1; v < size; v++){
            edges += this->graph->getAdj((*cluster)[u], (*cluster)[v]);
        }
    }
    value = (4.0*edges-degrees)/size;

    return value;
}



Bipart_ILS::Bipart_ILS(vector<unsigned> &nodes, LargeGraph * graph, IloNum **Wij, float alpha){
    this->graph = graph;
    this->Wij = Wij;
    this->alpha = alpha;
    this->nodes = &nodes;
    /*for (unsigned i=0;i<nodes.size();i++){
        this->maxNodes.push_back(0);
        this->currentNodes.push_back(0);
        a.push_back(0);
    }*/
    this->terminouZero = false;
}

void Bipart_ILS::execute(IloNumArray &lambda, IloNum &lambdaY){

    //this->terminouZero = false;
//cout<<"\n   (((( 1 ...";cout.flush();
    this->maxNodes.clear();
    this->currentNodes.clear();
    this->a.clear();

    for (unsigned i=0;i<nodes->size();i++){
        this->maxNodes.push_back(0);
        this->currentNodes.push_back(0);
        a.push_back(0);
    }
//cout<<"\n   (((( 2 ...";cout.flush();
    //choosing the alpha
    this->coarsedNodes.clear();
    this->currentLNodes.clear();
    this->currentLNodes2.clear();
    for (unsigned i=0; i<this->nodes->size(); i++){
        this->currentLNodes2.push_back(i);
    }

    ////starting a new search
    //this->currentValue = 0.0;
    this->currentValue = this->calculateAll(lambda, lambdaY);
/*if (isnan(currentValue)){
    cout<<"\nAAA:"<<lambda<<", "<<lambdaY;
    std::cin.get();
}*/
    //this->maxValue = 0.0;
    this->maxValue = this->currentValue ;


    //... all nodes are not in the new cluster (new column)
    //for (unsigned i=0;i<nodes->size();i++){
    //    this->maxNodes[i]=0;
    //    this->currentNodes[i]=0;
    //}
    //... start the coarsened nodes that are not coarsened yeat
    this->coarsedNodes.clear();
    for (unsigned i=0;i<nodes->size();i++){
        unsigned deg = this->graph->degreeOfNode[(*nodes)[i]];
        this->coarsedNodes.push_back(BipILSNode(i, this->Wij[(*nodes)[i]][(*nodes)[i]], deg, lambda[i]));
    }
//cout<<"\n   (((( 3 ...";cout.flush();
    this->sumLambdaS = 0.0;
    this->sumLambdaS2 = 0.0; //all nodes are in the second group
    for (unsigned i =0; i< lambda.getSize(); i++){
        this->sumLambdaS2 += lambda[i];
    }

//cout<<"\n   (((( 4 ...";cout.flush();
    float lastValue=0.0;

//cout<<"\nAtual(antesShuffle):  "<<this->currentValue;cout.flush();
//cout<<"\nReal (antesShuffle):   "<<this->calculateAll(lambda);


    //start random solution
    shuffleCurrentSolution(lambda, 0.5, lambdaY);

//cout<<"\n   (((( 5 ...";cout.flush();

//cout<<"\nAtual:  "<<this->currentValue;cout.flush();
//cout<<"\nReal:   "<<this->calculateAll(lambda);
//return;


    ////mini louvain method
    vector<unsigned> randomOrder(this->coarsedNodes.size());

    unsigned notImproved = 0;
    unsigned it=0;
    while( notImproved < nodes->size()){

//cout<<"\n   (((( 6("<<it<<") ...";cout.flush();

        it++;


        ////level phase
        int notImprovedLvL = 0;
        unsigned itInternal = 0;
        while (notImprovedLvL < 2 && itInternal<graph->numberOfEdges){//graph->numberOfNodes){
            itInternal++;
/*cout<<"\nAtual:  "<<this->currentValue;cout.flush();
cout<<"\nReal:   "<<this->calculateAll(lambda, lambdaY);
cout<<"\n size("<<currentLNodes.size()<<")";
*/
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

                //fixing the first node.

                if(node == 0){
                    if (this->terminouZero == false){
                         if(this->currentNodes[node] ==1){
                             value = this->gainIfIExit(node, lambda, lambdaY);
//if(it % 100 == 0){
//cout<<"\n \t\t\t------>insideH -- it["<<it<<"] value     ("<<value<<")";cout.flush();
//}
                            this->currentNodes[node] = 0;
                            this->currentValue += value;
                            //each subnode belongs to the solution
                            unordered_set<unsigned>::iterator itnews = this->coarsedNodes[node].nodes.begin();
                            while (itnews != this->coarsedNodes[node].nodes.end()){
                                this->currentLNodes.remove(*itnews);
                                this->currentLNodes2.push_back(*itnews);
                                itnews++;
                            }
                            this->sumLambdaS  -= this->coarsedNodes[node].lambda;
                            this->sumLambdaS2 += this->coarsedNodes[node].lambda;
                         }
                    }else{
                        if(this->currentNodes[node] ==0){
                            value = this->gainIfIEnter(node, lambda, lambdaY);
//if(it % 100 == 0){
//cout<<"\n \t\t\t------>insideH -- it["<<it<<"] value     ("<<value<<")";cout.flush();
//}
                            this->currentNodes[node] = 1;
                            this->currentValue += value;
                            //each subnode belongs to the solution
                            unordered_set<unsigned>::iterator itnews = this->coarsedNodes[node].nodes.begin();
                            while (itnews != this->coarsedNodes[node].nodes.end()){
                                this->currentLNodes.push_back(*itnews);
                                this->currentLNodes2.remove(*itnews);
                                itnews++;
                            }
                            this->sumLambdaS  += this->coarsedNodes[node].lambda;
                            this->sumLambdaS2 -= this->coarsedNodes[node].lambda;
                        }
                    }
                    continue;
                 }

/*cout<<"\nNode: ["<<node<<"] ==> "<<(*nodes)[node];
cout<<"\n";
this->constructSolutionArray(lambda);
vector<unsigned> x;
for (unsigned i=0; i< a.size();i++){
    cout<<", "<<this->a[i];
//    if (this->a[i] == 1)
//        x.push_back((*nodes)[i]);
}
//long double Da=calculateD(&x);
//cout<<"\nCORRECT (a):"<<Da;
//vector<unsigned> b;
//for (unsigned i=0; i< this->a.size();i++){

//    if (this->a[i] == 0)
//        b.push_back((*nodes)[i]);
//}
//long double Db=calculateD(&b);
//cout<<"\nCORRECT (b):"<<Db;
//cout<<"\nCORRECT (a+b):"<<Da+Db;
float cas = this->calculateAllSimulate(node,lambda);
cout<<"\nCalculateSimul: "<<this->calculateAllSimulate(node,lambda);
*/
                if (this->currentNodes[node] == 0){
                    value = this->gainIfIEnter(node, lambda, lambdaY);

//if(itInternal % 100 == 0){
//cout<<"\n \t\t\t------>insideH -- it["<<it<<"] value     ("<<value<<")";cout.flush();
//}

/*cout<<"\n*********BIen*** it["<<it<<"] cur.value ("<<this->currentValue<<")";cout.flush();
cout<<"\n*********BIen*** it["<<it<<"] value     ("<<value<<")";cout.flush();
cout<<"\nCalculateAll:   "<<this->calculateAll(lambda);
cout<<"\nRealValue:   "<<cas-this->calculateAll(lambda);
*/


                    if(value >0.0000001){
                        this->currentNodes[node] = 1;
                        //an improvement happens
                        notImprovedLvL = 0;

                        this->currentValue += value;
                        //each subnode belongs to the solution
                        unordered_set<unsigned>::iterator itnews = this->coarsedNodes[node].nodes.begin();
                        while (itnews != this->coarsedNodes[node].nodes.end()){
                            this->currentLNodes.push_back(*itnews);
                            this->currentLNodes2.remove(*itnews);
                            itnews++;
                        }
                        this->sumLambdaS  += this->coarsedNodes[node].lambda;
                        this->sumLambdaS2 -= this->coarsedNodes[node].lambda;
                    }
                }else{

                    value = this->gainIfIExit(node, lambda, lambdaY);
//if(itInternal % 100 == 0){
//cout<<"\n \t\t\t------>insideH -- it["<<it<<"] value     ("<<value<<")";cout.flush();
//}

/*cout<<"\n*********BIex*** it["<<it<<"] cur.value ("<<this->currentValue<<")";cout.flush();
cout<<"\n*********BIex*** it["<<it<<"] value     ("<<value<<")";cout.flush();
cout<<"\nCalculateAll: "<<this->calculateAll(lambda);
cout<<"\nRealValue:   "<<cas-this->calculateAll(lambda);
*/

                    if(value >0.0000001){
                        this->currentNodes[node] = 0;
                        //an improvement happens
                        notImprovedLvL = 0;

                        this->currentValue += value;
                        //each subnode belongs to the solution
                        unordered_set<unsigned>::iterator itnews = this->coarsedNodes[node].nodes.begin();
                        while (itnews != this->coarsedNodes[node].nodes.end()){
                            this->currentLNodes.remove(*itnews);
                            this->currentLNodes2.push_back(*itnews);
                            itnews++;
                        }
                        this->sumLambdaS  -= this->coarsedNodes[node].lambda;
                        this->sumLambdaS2 += this->coarsedNodes[node].lambda;
                    }
                }
                node++;
            }
            notImprovedLvL++;



//return;//##

        }
        if (lastValue >= this->currentValue){
            notImproved++;
        }
        lastValue = this->currentValue;

        if (this->currentValue > this->maxValue){
            this->maxValue = this->currentValue;
            this->copySolutionArray(this->maxNodes, lambda);
        }
        shuffleCurrentSolution(lambda, this->alpha, lambdaY);

        /*if( notImproved == nodes->size() && this->terminouZero == false){
            this->terminouZero = true;
            notImproved = 0;
        }*/
    }

}

float Bipart_ILS::gainIfIEnter(unsigned node, IloNumArray &lambda, IloNum &lambdaY){
///FIRST HALF
    //first part
    float value = this->sumLambdaS * this->coarsedNodes[node].nodes.size()
            +  (this->currentLNodes.size()+this->coarsedNodes[node].nodes.size())
            *this->coarsedNodes[node].lambda;

    //new part
    value   += lambdaY;


    //second part
    //value -= 4.0*this->coarsedNodes[node].internalEdges;
    unordered_set<unsigned>::iterator ita = this->coarsedNodes[node].nodes.begin();
    while (ita != this->coarsedNodes[node].nodes.end()){
        list<unsigned>::iterator itb = this->currentLNodes.begin();
        while(itb != this->currentLNodes.end()){
            value -= 4.0*this->Wij[(*nodes)[*itb]][(*nodes)[*ita]];
            itb++;
        }
        ita++;
    }

    //third part
    value += this->coarsedNodes[node].degree;


///SECOND HALF
/*    //first part
    value -= (this->sumLambdaS2-this->coarsedNodes[node].lambda) * this->coarsedNodes[node].nodes.size()
                +  (this->currentLNodes2.size())*this->coarsedNodes[node].lambda;
                ;

    //second part
    ita = this->coarsedNodes[node].nodes.begin();
    while (ita != this->coarsedNodes[node].nodes.end()){
        list<unsigned>::iterator itb = this->currentLNodes2.begin();
        while(itb != this->currentLNodes2.end()){
            value += 4.0*this->Wij[(*nodes)[*itb]][(*nodes)[*ita]];
            itb++;
        }
        ita++;
    }

    //third part
    value -= this->coarsedNodes[node].degree;
*/
    return /*max*/ value *-1.0;
}

float Bipart_ILS::gainIfIExit(unsigned node, IloNumArray &lambda, IloNum &lambdaY){
///FIRST HALF
    //first part
    float value = 0.0;
    value -= (this->sumLambdaS-this->coarsedNodes[node].lambda) * this->coarsedNodes[node].nodes.size()
                +  (this->currentLNodes.size())*this->coarsedNodes[node].lambda;
    //new part
    value   -= lambdaY;


        //second part
        //value -= 4.0*this->coarsedNodes[node].internalEdges;
        unordered_set<unsigned>::iterator ita = this->coarsedNodes[node].nodes.begin();
        while (ita != this->coarsedNodes[node].nodes.end()){
            list<unsigned>::iterator itb = this->currentLNodes.begin();
            while(itb != this->currentLNodes.end()){
                value += 4.0*this->Wij[(*nodes)[*itb]][(*nodes)[*ita]];
                itb++;
            }
            ita++;
        }

        //third part
        value -= this->coarsedNodes[node].degree;
///SECOND HALF
/*    //first part
    value += this->sumLambdaS2 * this->coarsedNodes[node].nodes.size()
            +  (this->currentLNodes2.size()+this->coarsedNodes[node].nodes.size())
            *this->coarsedNodes[node].lambda;


    //second part
    //value -= 4.0*this->coarsedNodes[node].internalEdges;
    ita = this->coarsedNodes[node].nodes.begin();
    while (ita != this->coarsedNodes[node].nodes.end()){
        list<unsigned>::iterator itb = this->currentLNodes2.begin();
        while(itb != this->currentLNodes2.end()){
            value -= 4.0*this->Wij[(*nodes)[*itb]][(*nodes)[*ita]];
            itb++;
        }
        ita++;
    }

    //third part
    value += this->coarsedNodes[node].degree;
*/
    return /*max*/ value *-1.0;


}

float Bipart_ILS::calculateAll(IloNumArray &lambda, IloNum &lambdaY){
    vector<unsigned> a(nodes->size());
    for (unsigned i=0;i<this->coarsedNodes.size();i++){
        unordered_set<unsigned>::iterator ita = this->coarsedNodes[i].nodes.begin();
        while (ita != this->coarsedNodes[i].nodes.end()){
            a[*ita] = this->currentNodes[i];
            ita++;
        }
    }
    float value = 0.0;
///FIRST HALF
    //first part
    for (unsigned u=0;u<nodes->size(); u++){
        for (unsigned v=0;v<nodes->size(); v++){
            value += a[u] * a[v] * lambda[v];
        }
    }
    //new part
    for (unsigned v=0;v<nodes->size(); v++){
        value += a[v] * lambdaY;
    }
    //second part
    for (unsigned v=0;v<this->nodes->size(); v++){
        for (unsigned u=v+1;u<this->nodes->size(); u++){
            value -= 4.0*a[u] * a[v] * this->Wij[(*nodes)[u]][(*nodes)[v]];
        }
    }
    //third part
    for (unsigned v=0;v<this->nodes->size(); v++){
        value += a[v]*this->graph->degreeOfNode[(*nodes)[v]];
    }
/*
///SECOND HALF
    //first part
    for (unsigned u=0;u<this->graph->numberOfNodes; u++){
        for (unsigned v=0;v<this->graph->numberOfNodes; v++){
            value += (1.0-a[u]) * (1.0-a[v]) * lambda[v];
        }
    }
    //second part
    for (unsigned v=0;v<this->graph->numberOfNodes; v++){
        for (unsigned u=v+1;u<this->graph->numberOfNodes; u++){
            value -= 4.0*(1.0-a[u]) * (1.0-a[v]) * this->Wij[u][v];
        }
    }
    //third part
    for (unsigned v=0;v<this->graph->numberOfNodes; v++){
        value += (1.0-a[v])*this->graph->degreeOfNode[v];
    }
*/

    return /*max*/ value * -1;


}

void Bipart_ILS::constructSolutionArray(IloNumArray &lambda){
    for (unsigned i=0;i<this->coarsedNodes.size();i++){
        unordered_set<unsigned>::iterator ita = this->coarsedNodes[i].nodes.begin();
        while (ita != this->coarsedNodes[i].nodes.end()){
            this->a[*ita] = this->currentNodes[i];
            ita++;
        }
    }

}

void Bipart_ILS::shuffleCurrentSolution(IloNumArray &lambda, float strength, IloNum &lambdaY){
    unsigned nchanges = this->graph->numberOfNodes * strength;

    float value;
    for (unsigned ic=0;ic<nchanges;ic++){
        unsigned node = rand()%this->coarsedNodes.size();
        if (this->currentNodes[node] == 0){
            value = this->gainIfIEnter(node, lambda, lambdaY);
            this->currentNodes[node] = 1;
            //an improvement happens
            this->currentValue += value;
            //each subnode belongs to the solution
            unordered_set<unsigned>::iterator itnews = this->coarsedNodes[node].nodes.begin();
            while (itnews != this->coarsedNodes[node].nodes.end()){
                this->currentLNodes.push_back(*itnews);
                this->currentLNodes2.remove(*itnews);
                itnews++;
            }
            this->sumLambdaS  += this->coarsedNodes[node].lambda;
            this->sumLambdaS2 -= this->coarsedNodes[node].lambda;
//cout<<"\nENTRA";
        }else{
            value = this->gainIfIExit(node, lambda, lambdaY);
            this->currentNodes[node] = 0;
            //an improvement happens
            this->currentValue += value;
            //each subnode belongs to the solution
            unordered_set<unsigned>::iterator itnews = this->coarsedNodes[node].nodes.begin();
            while (itnews != this->coarsedNodes[node].nodes.end()){
                this->currentLNodes.remove(*itnews);
                this->currentLNodes2.push_back(*itnews);
                itnews++;
            }
            this->sumLambdaS  -= this->coarsedNodes[node].lambda;
            this->sumLambdaS2 += this->coarsedNodes[node].lambda;
//cout<<"\nSAI";
        }
//cout<<"\n-> ic"<<ic<<", Atual:  "<<this->currentValue;cout.flush();
//cout<<"\n-> ic"<<ic<<", Real:   "<<this->calculateAll(lambda);


    }


}

void Bipart_ILS::copySolutionArray(vector<unsigned> &target, IloNumArray &lambda){
    for (unsigned i=0;i<this->coarsedNodes.size();i++){
        unordered_set<unsigned>::iterator ita = this->coarsedNodes[i].nodes.begin();
        while (ita != this->coarsedNodes[i].nodes.end()){
            target[*ita] = this->currentNodes[i];
            ita++;
        }
    }
}



float Bipart_ILS::calculateAllSimulate(unsigned posNode, IloNumArray &lambda, IloNum &lambdaY){
    vector<unsigned> a(this->graph->numberOfNodes);
    for (unsigned i=0;i<this->coarsedNodes.size();i++){
        unordered_set<unsigned>::iterator ita = this->coarsedNodes[i].nodes.begin();
        while (ita != this->coarsedNodes[i].nodes.end()){
            a[*ita] = this->currentNodes[i];
            ita++;
        }
    }
    float value = 0.0;
///FIRST HALF
    //first part
    unsigned au, av;
    for (unsigned u=0;u<this->graph->numberOfNodes; u++){
        for (unsigned v=0;v<this->graph->numberOfNodes; v++){
            au=a[u], av=a[v];
            if (u == posNode){
                au = !au;
            }
            if (v == posNode){
                av = !av;
            }
            value += au * av * lambda[v];
        }
    }
    //new part
    for (unsigned v=0;v<this->graph->numberOfNodes; v++){
        av=a[v];
        if (v == posNode){
            av = !av;
        }
        value += av * lambdaY;
    }
    //second part
    for (unsigned v=0;v<this->graph->numberOfNodes; v++){
        for (unsigned u=v+1;u<this->graph->numberOfNodes; u++){
            au=a[u], av=a[v];
            if (u == posNode){
                au = !au;
            }
            if (v == posNode){
                av = !av;
            }
            value -= 4.0*au * av * this->Wij[u][v];
        }
    }
    //third part
    for (unsigned v=0;v<this->graph->numberOfNodes; v++){
        av=a[v];

        if (v == posNode){
            av = !av;
        }
        value += av*this->graph->degreeOfNode[v];
    }
/*
///SECOND HALF
    //first part
    for (unsigned u=0;u<this->graph->numberOfNodes; u++){
        for (unsigned v=0;v<this->graph->numberOfNodes; v++){
            au=a[u], av=a[v];
            if (u == posNode){
                au = !au;
            }
            if (v == posNode){
                av = !av;
            }
            value += (1.0-au) * (1.0-av) * lambda[v];
        }
    }
    //second part
    for (unsigned v=0;v<this->graph->numberOfNodes; v++){
        for (unsigned u=v+1;u<this->graph->numberOfNodes; u++){
            au=a[u], av=a[v];
            if (u == posNode){
                au = !au;
            }
            if (v == posNode){
                av = !av;
            }
            value -= 4.0*(1.0-au) * (1.0-av) * this->Wij[u][v];
        }
    }
    //third part
    for (unsigned v=0;v<this->graph->numberOfNodes; v++){
        av=a[v];
        if (v == posNode){
            av = !av;
        }
        value += (1.0-av)*this->graph->degreeOfNode[v];
    }
*/
    return /*max*/ value * -1;


}



float Bipart_ILS::calculateAllByVector(vector<unsigned> &a, IloNumArray &lambda, IloNum &lambdaY){
cout<<"\n --> 1 ";cout.flush();
    float value = 0.0;

float primeiro, segundo, terceiro, quarto;

///FIRST HALF
    //first part
    for (unsigned u=0;u<this->nodes->size(); u++){
        for (unsigned v=0;v<this->nodes->size(); v++){
            value += a[u] * a[v] * lambda[v];
        }
    }
//primeiro = value;
//cout<<"\n --> 2 ";cout.flush();
    //new part
    for (unsigned v=0;v<this->nodes->size(); v++){
        value += a[v]*lambdaY;
    }
//segundo = value-primeiro;
//cout<<"\n --> 3 ";cout.flush();
    //second part
    for (unsigned v=0;v<this->nodes->size(); v++){
        for (unsigned u=v+1;u<this->nodes->size(); u++){

            value -= 4.0 * a[u] * a[v] * this->Wij[(*nodes)[u] ][(*nodes)[v] ];
        }
    }
//terceiro = value-segundo-primeiro;
//cout<<"\n --> 4 ";cout.flush();
    //third part
    for (unsigned v=0;v<this->nodes->size(); v++){
        unsigned x = this->graph->degreeOfNode[(*nodes)[v]];
        value += a[v]*x;
    }
//quarto = value-terceiro-segundo-primeiro;
//int test = this->nodes->size();
//cout<<"\ncccccREAL: ("<<primeiro<<", "<<segundo<<", "<<terceiro<<", "<<quarto<<")"<<test<<","<<lambdaY;
//cout.flush();


//cout<<"\n --> 5 ";cout.flush();
/*
///SECOND HALF
    //first part
    for (unsigned u=0;u<this->graph->numberOfNodes; u++){
        for (unsigned v=0;v<this->graph->numberOfNodes; v++){
            value += (1.0-a[u]) * (1.0-a[v]) * lambda[v];
        }
    }
    //second part
    for (unsigned v=0;v<this->graph->numberOfNodes; v++){
        for (unsigned u=v+1;u<this->graph->numberOfNodes; u++){
            value -= 4.0*(1.0-a[u]) * (1.0-a[v]) * this->Wij[u][v];
        }
    }
    //third part
    for (unsigned v=0;v<this->graph->numberOfNodes; v++){
        value += (1.0-a[v])*this->graph->degreeOfNode[v];
    }
*/
    return /*max*/ value * -1;


}
