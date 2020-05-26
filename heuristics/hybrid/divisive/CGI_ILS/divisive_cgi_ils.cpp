#include "divisive_cgi_ils.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>

Divisive_CGII_ILS::Divisive_CGII_ILS(LargeGraph * graph, float alpha){
    this->graph = graph;
    this->alpha = alpha;




    ///////////////////////////////
    Solution * sol;
    for(unsigned i=0;i<30;i++){
        McnLnm mcnlnm(this->graph);
        mcnlnm.execute();
        if (mcnlnm.FastGreedyImproved == true){
            sol = mcnlnm.fglg->toSolution();
            sol->modularidade = mcnlnm.fglg->bestDensity;
            solutions.push_back(sol);
        }else{
            sol = new Solution(this->graph,"["+mcnlnm.hlv->bestCommunityStr+"]");
            sol->modularidade = mcnlnm.hlv->bestDensity;
            solutions.push_back(sol);
        }
    }
    ///////////////////////////////


}


void Divisive_CGII_ILS::execute(){

    Bipart_CGII_ILS bipCGIILS( this->graph, this->alpha, &this->solutions );

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


        bipCGIILS.execute(*currentCluster);
//cout<<"40";cout.flush();

/*cout<<"\n - "<<level<<"[";
for (unsigned j=0;j<bipCGIILS.firstHalf.size();j++){
    cout<<", "<<bipCGIILS.firstHalf[j];
}
cout<<"][";
for (unsigned j=0;j<bipCGIILS.secondHalf.size();j++){
    cout<<", "<<bipCGIILS.secondHalf[j];
}
cout<<"]";*/
//if (level == 1) return;

//cout<<"56";cout.flush();

        notVisiteds.erase(currentCluster);
//cout<<"57";cout.flush();
        if (Dvalue < bipCGIILS.maxDensity){

            //if the cluster does not divide, then it is not visited again
            if (bipCGIILS.firstHalf.size() ==0 || bipCGIILS.secondHalf.size()==0){
                continue;
            }else{
                clusters.erase(currentCluster);
                delete currentCluster;
            }



            vector<unsigned> *ca = new vector<unsigned>;
            vector<unsigned> *cb = new vector<unsigned>;
            ca->insert(ca->begin(), bipCGIILS.firstHalf.begin(), bipCGIILS.firstHalf.end());
            cb->insert(cb->begin(), bipCGIILS.secondHalf.begin(), bipCGIILS.secondHalf.end());
            //clusters.push_back(ca);
            //clusters.push_back(cb);

            long double Dfirst, Dsecond;
            Dfirst = this->calculateD(& bipCGIILS.firstHalf);
            Dsecond= this->calculateD(& bipCGIILS.secondHalf);
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

long double Divisive_CGII_ILS::calculateD( vector<unsigned> * cluster ){
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



