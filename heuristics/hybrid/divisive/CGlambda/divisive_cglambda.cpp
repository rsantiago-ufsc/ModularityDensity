#include "divisive_cglambda.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>

Divisive_CGLambda::Divisive_CGLambda(LargeGraph * graph, float alpha){
    this->graph = graph;
    this->alpha = alpha;
}


void Divisive_CGLambda::execute(){

    Bipart_CGLambda bipCGLambda( this->graph, this->alpha );

    //creating partition and assigning the cluster with all nodes
    //vector< vector<unsigned> * > clusters;
    unordered_map< vector<unsigned> *, long double > clusters;
    vector<unsigned> *c1 = new vector<unsigned>;
    for (unsigned i=0; i<this->graph->numberOfNodes; i++){
        c1->push_back(i);

    }
    clusters.insert(pair< vector<unsigned> *, long double >(c1, calculateD(c1)));

    //visited clusters
    //vector<bool> visiteds;
    //visiteds.push_back(false);
    unordered_set< vector<unsigned> * > notVisiteds;
    notVisiteds.insert(c1);

unsigned level = 0;
    while (notVisiteds.size() > 0){
        level++;
//cout<<"37";cout.flush();

        unordered_set< vector<unsigned> * >::iterator current = notVisiteds.begin();
//cout<<"38";cout.flush();
        vector<unsigned> * currentCluster = *current;
//cout<<"39";cout.flush();
        long double Dvalue = clusters[currentCluster];

/*cout<<"\n \n\n::antes:: ("<<level<<") [";
for (unsigned j=0;j<currentCluster->size();j++){
    cout<<", "<<(*currentCluster)[j];
}
cout<<"]";
cout.flush();*/


        bipCGLambda.execute(*currentCluster);
//cout<<"40";cout.flush();

/*cout<<"\n - "<<level<<"[";
for (unsigned j=0;j<bipCGLambda.firstHalf.size();j++){
    cout<<", "<<bipCGLambda.firstHalf[j];
}
cout<<"][";
for (unsigned j=0;j<bipCGLambda.secondHalf.size();j++){
    cout<<", "<<bipCGLambda.secondHalf[j];
}
cout<<"]";*/
//if (level == 1) return;

//cout<<"56";cout.flush();

        notVisiteds.erase(currentCluster);
//cout<<"57";cout.flush();
        if (Dvalue < bipCGLambda.maxDensity){

            //if the cluster does not divide, then it is not visited again
            if (bipCGLambda.firstHalf.size() ==0 || bipCGLambda.secondHalf.size()==0){
                continue;
            }else{
                clusters.erase(currentCluster);
                delete currentCluster;
            }



            vector<unsigned> *ca = new vector<unsigned>;
            vector<unsigned> *cb = new vector<unsigned>;
            ca->insert(ca->begin(), bipCGLambda.firstHalf.begin(), bipCGLambda.firstHalf.end());
            cb->insert(cb->begin(), bipCGLambda.secondHalf.begin(), bipCGLambda.secondHalf.end());
            //clusters.push_back(ca);
            //clusters.push_back(cb);

            long double Dfirst, Dsecond;
            Dfirst = this->calculateD(& bipCGLambda.firstHalf);
            Dsecond= this->calculateD(& bipCGLambda.secondHalf);
            //D.push_back(Dfirst);
            //D.push_back(Dsecond);
            clusters.insert(pair< vector<unsigned> *, long double >(ca, Dfirst));
            clusters.insert(pair< vector<unsigned> *, long double >(cb, Dsecond));


            if (ca->size() > 3){
                //visiteds.push_back(false);
                notVisiteds.insert(ca);
            }
            if (cb->size() > 3){
                //visiteds.push_back(false);
                notVisiteds.insert(cb);
            }

        }

    }
    maxDensity = 0.0;
//cout<<"59";cout.flush();
    unordered_map< vector<unsigned> *, long double >::iterator current = clusters.begin();
    while (current != clusters.end()){
        maxDensity += current->second;
        current++;
    }
}

long double Divisive_CGLambda::calculateD( vector<unsigned> * cluster ){
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



