#include "louvain.h"

#include "../../../utils/modularitylg.h"



Louvain::Louvain(LargeGraph * graph, unsigned type)
{

    this->graph = graph;
    this->typePrioritizer = type;

    this->gl = new GraphLouvain(graph);


    unsigned i;
    unsigned m2 = this->graph->numberOfEdges*2;

    //preparing nodes
    this->size = this->graph->numberOfNodes;
    for (i=0;i<this->graph->numberOfNodes;i++){
        long double degree = this->graph->degreeOfNode[i];
        long double internalEdges = this->graph->getAdj(i,i);
        unsigned nnodes = 1;
        long double dens = internalEdges/m2 - (degree*degree)/(m2*m2);
        this->communities.push_back(CommunityLouvain(degree,internalEdges,nnodes,dens));

        this->commByNode.push_back(i);
        this->commByNodeAux.push_back(NO_COMMUNITY);
    }

    //preparing temporary
    this->neighWeight.resize(size,-1);
    this->neighPos.resize(size);
    this->neighLast=0;

}


Louvain::Louvain( Solution *sol, LargeGraph * graph, unsigned type)
{

    this->graph = graph;
    this->typePrioritizer = type;

    this->gl = new GraphLouvain(graph);


    unsigned i;
    unsigned m2 = this->graph->numberOfEdges*2;

    //preparing nodes

    for (i=0;i<this->graph->numberOfNodes;i++){
        this->commByNode.push_back(NO_COMMUNITY);
        this->commByNodeAux.push_back(NO_COMMUNITY);
    }

    this->size = sol->comunidades.size();
    for (i=0;i<this->size;i++){
        long double degree = 0;
        long double internalEdges = 0;
        unsigned nnodes = 0;
        itemLDE<unsigned> * aux = sol->comunidades[i]->getInicio();
        while (aux != NULL){
            nnodes++;
            this->commByNode[aux->id]=i;
            degree += this->graph->degreeOfNode[aux->id];
            itemLDE<unsigned> * aux2 = aux->prox;
            while (aux2 != NULL){
                internalEdges +=this->graph->getAdj(aux->id,aux2->id);
                aux2 = aux2->prox;
            }
            aux = aux->prox;
        }
        long double dens = internalEdges/m2 - (degree*degree)/(m2*m2);
        this->communities.push_back(CommunityLouvain(degree,internalEdges,nnodes,dens));
    }


    ///old
    /*this->size = this->graph->numberOfNodes;
    for (i=0;i<this->graph->numberOfNodes;i++){
        long double degree = this->graph->degreeOfNode[i];
        long double internalEdges = this->graph->getAdj(i,i);
        unsigned nnodes = 1;
        long double dens = internalEdges/m2 - (degree*degree)/(m2*m2);
        this->communities.push_back(CommunityLouvain(degree,internalEdges,nnodes,dens));

        this->commByNode.push_back(i);
        this->commByNodeAux.push_back(NO_COMMUNITY);
    }*/

    //preparing temporary
    this->neighWeight.resize(size,-1);
    this->neighPos.resize(size);
    this->neighLast=0;



}


long double Louvain::calculateDensity(){
    long double dens  = 0.0;
    if (this->typePrioritizer == LV_PRI_MOD_DELTA){
        long double m = this->graph->numberOfEdges;
        long double m2 = this->graph->numberOfEdges*2.0;
        //##we can avoid O(n) computation here by updating the vector communities at each pass
        for (unsigned i=0 ; i<this->gl->size ; i++) {
            if (this->communities[i].nnodes>0){
                dens += this->communities[i].internalEdges/m
                        - (this->communities[i].degree/m2)*(this->communities[i].degree/m2);

            }
        }

        return dens;
    }
    if (  this->typePrioritizer == LV_PRI_DNS_DELTA
          || this->typePrioritizer == LV_PRI_SAN_DELTA
          || this->typePrioritizer == LV_PRI_SAN){
        for (unsigned i=0 ; i<this->gl->size ; i++) {
            if (this->communities[i].nnodes>0){
                dens += (4.0*this->communities[i].internalEdges-this->communities[i].degree)
                        /this->communities[i].nnodes;


            }
        }
        return dens;

    }
}



long double Louvain::calculateModularity(){
    long double dens  = 0.0;
    long double m = this->graph->numberOfEdges;
    long double m2 = this->graph->numberOfEdges*2.0;
    //##we can avoid O(n) computation here by updating the vector communities at each pass
    for (unsigned i=0 ; i<this->gl->size ; i++) {
        if (this->communities[i].nnodes>0){
            dens += this->communities[i].internalEdges/m
                    - (this->communities[i].degree/m2)*(this->communities[i].degree/m2);

        }
    }

    return dens;

}


string Louvain::currentToString(){
    vector <string> comms;
    for (unsigned i=0; i<this->commByNode.size();i++){
        comms.push_back("");
    }
    for (unsigned i=0; i<this->commByNode.size();i++){

        unordered_set<unsigned>::iterator node = this->gl->nodes[i].nodes->begin();
        while(node != this->gl->nodes[i].nodes->end()){

            comms[ this->commByNode[i] ] += to_string(*node)+",";
            node++;
        }
    }
    string comm = "";
    bool first = true;
    for (unsigned i=0; i<this->commByNode.size();i++){
        if (comms[i].size() > 0){
            comms[i].erase(comms[i].size()-1,1);
            comms[i] += "]";
            if (first){
                comm += "["+comms[i];
                first=false;
            }else{
                comm += ",["+comms[i];
            }
        }
    }

    return comm;

}




bool Louvain::oneLevel(){

    unsigned i, randpos, inode, node;
    bool improvement=false;
    long double m = (this->graph->numberOfEdges);
    long double m2 = (this->graph->numberOfEdges*2.0);
    long double minModularity = 0.0000000001;
    int nmoves;
    int npass = 0;
    long double newDensity   = this->calculateDensity();
    double currDensity   = newDensity;
    long double delta = 0.0;

    unsigned numberOfCommunities = gl->size;

    //assuring random order
    vector<unsigned> randomOrder(this->gl->size);
    for (i=0 ; i<this->gl->size ; i++){
        randomOrder[i]=i;
    }

    for (int i=0 ; i<this->gl->size-1 ; i++) {
        randpos = rand()%(this->gl->size-i)+i;
        swap(randomOrder[i], randomOrder[randpos]);
    }

    do {
      currDensity = newDensity;
      nmoves = 0;
      npass++;
      delta = 0.0;

      //for each node (randomized)
      for (inode=0 ; inode<this->gl->size ; inode++) {

        node                = randomOrder[inode];


        unsigned nodeComm   = this->commByNode[node];
        long double wdegree     = this->gl->nodes[node].totalDegree;

        // computation of all neighboring communities of current node
        this->neighsPrepare(node);

        // remove node from its current community
        this->removeFromCommunity(node, nodeComm, neighWeight[nodeComm]);

        // compute the nearest community for node
        // default choice for future insertion is the former community
        unsigned bestComm  = nodeComm;
        long double bestNLinks  = 0.0;
        long double bestIncrease = 0.0;
#if TYPE_PRI == LV_PRI_SAN_DELTA || TYPE_PRI == LV_PRI_SAN
        long double bestIncrease2 = PLUS_INFINITY;
#endif
        for (i=0 ; i<this->neighLast ; i++) {
          long double increase = this->calculateDensityGain(node, neighPos[i], neighWeight[neighPos[i]], wdegree);



          //return false;

          if (increase>bestIncrease) {
            bestComm     = neighPos[i];
            bestNLinks  = neighWeight[neighPos[i]];
            bestIncrease = increase;
#if TYPE_PRI == LV_PRI_SAN_DELTA || TYPE_PRI == LV_PRI_SAN
            bestIncrease2 = this->calculateDensityGain2(node, neighPos[i], neighWeight[neighPos[i]], wdegree);
#endif
          }


#if TYPE_PRI == LV_PRI_SAN_DELTA || TYPE_PRI == LV_PRI_SAN //only apply the second facto if it is SAN_DELTA
          if (increase == bestIncrease){
              long double increase2 = this->calculateDensityGain2(node, neighPos[i], neighWeight[neighPos[i]], wdegree);
              //cout<<"\n"<<increase<<" == "<<bestIncrease;
              //cout<<"\n"<<increase2<<" < "<<bestIncrease2;
              if (increase2 < bestIncrease2){
                bestComm     = neighPos[i];
                bestNLinks  = neighWeight[neighPos[i]];
                bestIncrease = increase;
                bestIncrease2 = increase2;
              }
          }
#endif
        }


        // insert node in the nearest community
        this->insertInCommunity(node, bestComm, bestNLinks);

        if (bestComm!=nodeComm){
//if(this->it == 2)cout<<"\n\t\t node: "<<node<<" comm:"<<bestComm;cout.flush();
            if (this->communities[nodeComm].nnodes == 0){
                numberOfCommunities--;
            }
            if (this->communities[bestComm].nnodes == 1){
                numberOfCommunities++;
            }
            nmoves++;
        }

      }
      long double temp = this->calculateDensity();
      newDensity=temp;
      if (this->bestDensity < temp){
          this->bestDensity = temp;
          this->bestMod = this->calculateModularity();
#ifdef GROUND_TRUTH
          this->bestCommunityStr = this->currentToString();
#endif
          this->bestNumberOfCommunities = numberOfCommunities;
          improvement=true;
      }
      //if (this->typePrioritizer == LV_PRI_MOD_DELTA){
      //    this->bestDensity += delta;//deu problema com selloops... verificar depois
      //}

      if (nmoves>0){

      }

    } while (nmoves>0 && newDensity-currDensity>minModularity);
/*
 *
O ganho deve ser maior que zero para poder ser aceito (bestIncrease).
Se há ganhos maiores que zero, sempre a cada troca, significa que o delta está errado?

[  x  ] [  y  z ] = 10
[     ] [  y  z ] = 4
[     ] [  y x z ] = 6

CONSIDERAR ISTO: && newDensity-currDensity>min_modularity
 */
    return improvement;
}


void Louvain::neighsPrepare(unsigned node){

    unsigned neigh, i, neighComm;
    long double neighW;
    for (i=0 ; i<this->neighLast ; i++){
        this->neighWeight[this->neighPos[i]]=-1;
    }
    this->neighLast=0;

    neighPos[0]=this->commByNode[node];
    neighWeight[neighPos[0]]=0;
    neighLast=1;

    itemLDE<EdgeCoarsenerLouvain> * auxAdj = this->gl->nodes[node].adjacencies.getInicio();
    while (auxAdj != NULL){
        neigh = auxAdj->id.anotherEndpoint;
        neighComm   = this->commByNode[neigh];
        neighW = auxAdj->id.weight;

        if (neigh!=node) {
            if (this->neighWeight[neighComm]==-1) {
                neighWeight[neighComm]=0.;
                neighPos[neighLast]=neighComm;
                neighLast++;
            }
            neighWeight[neighComm]+=neighW;
        }

        auxAdj = auxAdj->prox;

    }
//cout<<"\nNL: "<<neighLast;

}

void Louvain::removeFromCommunity(unsigned node, unsigned comm, long double dnodecomm) {
    this->communities[comm].degree -= this->gl->nodes[node].totalDegree;
    this->communities[comm].internalEdges = this->communities[comm].internalEdges - (dnodecomm + this->gl->nodes[node].internalEdges);
    this->communities[comm].nnodes -= this->gl->nodes[node].nodes->size();
    this->commByNode[node]  = NO_COMMUNITY;
}

void Louvain::insertInCommunity(unsigned node, unsigned comm, long double dnodecomm) {
    this->communities[comm].degree += this->gl->nodes[node].totalDegree;
    this->communities[comm].nnodes += this->gl->nodes[node].nodes->size();
    this->communities[comm].internalEdges += dnodecomm + this->gl->nodes[node].internalEdges;
    this->commByNode[node]=comm;
}

long double Louvain::calculateDensityGain(unsigned node, unsigned comm, long double dnodecomm, long double wdegree) {

  if (this->typePrioritizer == LV_PRI_MOD_DELTA){
    long double m2   = this->graph->numberOfEdges*2.0;
    return (dnodecomm - this->communities[comm].degree * wdegree / m2);
  }
  if (this->typePrioritizer == LV_PRI_DNS_DELTA){
    /*return  (4.0*dnodecomm*this->communities[comm].nnodes
            -wdegree*this->communities[comm].nnodes
            -4.0*this->communities[comm].internalEdges
            +this->communities[comm].degree)/
            (this->communities[comm].nnodes*this->communities[comm].nnodes+this->communities[comm].nnodes);
    */
    return  (this->communities[comm].nnodes*4.0*dnodecomm
                  -this->communities[comm].nnodes*wdegree
                  -4.0*this->communities[comm].internalEdges*this->gl->nodes[node].nodes->size()
                  +this->communities[comm].degree*this->gl->nodes[node].nodes->size())/
                  (this->communities[comm].nnodes*this->communities[comm].nnodes+this->communities[comm].nnodes*this->gl->nodes[node].nodes->size());
  }
  if(this->typePrioritizer == LV_PRI_SAN_DELTA){
      //## ((verificar, acho que não precisa))preciso salvar o NNODES: preciso detectar quando está na mesma comunidade
      long double n = this->communities[comm].nnodes;
      long double x = n + this->gl->nodes[node].nodes->size();
      long double n2_n = (n*n-n);
      long double x2_x = (x*x-x);   

      if (n2_n>0.0){
          return  (2.0*(n2_n)*(this->communities[comm].internalEdges+dnodecomm)
                   - 2.0*x2_x*this->communities[comm].internalEdges
                   )/(n2_n*x2_x);
      }else{
          if (x2_x==0.0){
              if (wdegree == 0){
                  return 1.0; //elements with 0 degree should remain isolated
              }else{
                  return 0.0;
              }
          }else{
            return 2.0*(this->communities[comm].internalEdges+dnodecomm)/x2_x;
          }
      }
  }
  if(this->typePrioritizer == LV_PRI_SAN){
      long double n = this->communities[comm].nnodes;
      long double x = n + this->gl->nodes[node].nodes->size();
      long double x2_x = (x*x-x);
      if (x2_x>0.0){
        return 2.0*(this->communities[comm].internalEdges+dnodecomm)/x2_x;
      }else{
          if (wdegree == 0){
              return 1.0; //elements with 0 degree should remain isolated
          }else{
              return 0.0;
          }

      }
  }
}

long double Louvain::calculateDensityGain2(unsigned node, unsigned comm, long double dnodecomm, long double wdegree) {
    if(this->typePrioritizer == LV_PRI_SAN_DELTA){
        long double n = this->communities[comm].nnodes;
        long double x = this->gl->nodes[node].nodes->size();
        long double dv = wdegree;
        return  (n*dv - this->communities[comm].degree*x) / ((n+x)*n);
    }
    if(this->typePrioritizer == LV_PRI_SAN){
        long double n = this->communities[comm].nnodes;
        long double x = this->gl->nodes[node].nodes->size();
        return  (wdegree + this->communities[comm].degree) / (n+x);
    }
}


void Louvain::execute(){
    chrono::system_clock::time_point before;
    before = chrono::system_clock::now();
    if (this->typePrioritizer == LV_PRI_MOD_DELTA){
        this->bestDensity = 0.0;
    }

    this->it = 0;
    long double newDensity;
    bool improvement;

//ModularityLG mo(this->graph);

    do {
        this->it++;
//cout<<"\n Starting OneLevel "<<it<<chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count()/1000.0;cout.flush();
        improvement = this->oneLevel();
        newDensity = this->bestDensity;
//cout<<"\n Concluding OneLevel "<<it<<chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count()/1000.0;cout.flush();
//cout<<"\n Starting Update "<<it<<chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count()/1000.0;cout.flush();
        if(improvement){
            this->updateLouvainGraphAndCommunities();
        }
//cout<<"\n ConcludingUpdate "<<it<<chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count()/1000.0;cout.flush();
        if (it==1){ // do at least one more computation if partition is provided
           improvement=true;
        }
    } while(improvement);
    this->totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
}


void Louvain::updateLouvainGraphAndCommunities(){
//if(it==2){cout<<"\nNode:"<<node;cout.flush();}
    unsigned node, i, ic, in,master, ileft, iright;

    // Renumber communities
    vector<unsigned> renumber(gl->size, 0);
    for (node=0 ; node<gl->size; node++) {
      renumber[this->commByNode[node]]++;
    }
    unsigned final=0;
    for (i=0 ; i<gl->size ; i++){
      if (renumber[i]!=0){
        renumber[i]=final++;
      }
    }

    // Compute communities
    vector<vector<unsigned> > comm_nodes(final);
    for (node=0 ; node<gl->size ; node++) {
        comm_nodes[renumber[this->commByNode[node]]].push_back(node);
    }

    //first, the merge is made upon the nodes
    vector< pair<unsigned, unsigned> > trans;
//    vector<unsigned> translator(gl->size, NO_COMMUNITY);
    for (ic=0;ic<comm_nodes.size(); ic++){
        master = comm_nodes[ic][0];
        for (in=1;in<comm_nodes[ic].size(); in++){
            node = comm_nodes[ic][in];
            master = gl->merge(master, node, 0.0);
        }
        swap(this->communities[master],this->communities[this->commByNode[master]]);
//        translator[ic] = ic;

        //swap(gl->nodes[ic], gl->nodes[master]);

    }

    ileft=0;
    iright=gl->size-1;
    while (ileft < iright){

        if (gl->nodes[ileft].nodes->size() > 0 ){
//            translator[ileft] = ileft;
            ileft++;            
            continue;
        }
        if (gl->nodes[iright].nodes->size() == 0 ){
            iright--;
            continue;
        }
        //this->commByNodeAux[ileft] = this->commByNode[iright];
        swap(gl->nodes[ileft], gl->nodes[iright]);
        //swap(this->communities[ileft], this->communities[iright]);
//        translator[iright] = ileft;
        trans.push_back(pair<unsigned, unsigned>(ileft,iright));
        ileft++;
        iright--;
    }

   for (node=0;node<comm_nodes.size();node++){
       this->commByNode[node] = node;
       this->communities[node].internalEdges = gl->nodes[node].internalEdges;
       this->communities[node].degree = gl->nodes[node].totalDegree;
       this->communities[node].nnodes = gl->nodes[node].nodes->size();
       this->communities[node].density = gl->nodes[node].density;
   }
/*if (it == 2){
    cout<<"\n UP > GL4 edges["<<gl->nodes[4].internalEdges<<"] degree["<<gl->nodes[4].totalDegree<<"] nnodes["<<gl->nodes[4].nodes.size()<<"]";
    cout<<"\n UP > GL7 edges["<<gl->nodes[7].internalEdges<<"] degree["<<gl->nodes[7].totalDegree<<"] nnodes["<<gl->nodes[7].nodes.size()<<"]";
    cout<<"\n UP > COM4 edges["<<communities[4].internalEdges<<"] degree["<<communities[4].degree<<"] nnodes["<<communities[4].nnodes<<"]";
    cout<<"\n UP > COM7 edges["<<communities[7].internalEdges<<"] degree["<<communities[7].degree<<"] nnodes["<<communities[7].nnodes<<"]";
}*/


    //translating
   gl->size = comm_nodes.size();
   itemLDE<EdgeCoarsenerLouvain> *aux;
   vector<pair<unsigned, unsigned>>::iterator it;
   it = trans.begin();
   while(it != trans.end()){
       aux = gl->nodes[(*it).first].adjacencies.getInicio();
       while(aux!=NULL){
           aux->id.apointItem->id.anotherEndpoint = (*it).first;
           aux=aux->prox;
       }
       it++;
   }
/*
    for (in = 0;in<gl->size; in++){
        itemLDE<EdgeCoarsenerLouvain> *aux = gl->nodes[in].adjacencies.getInicio();
        while(aux!=NULL){
            aux->id.anotherEndpoint = translator[aux->id.anotherEndpoint];
            aux=aux->prox;
        }
    }*/

}

Solution * Louvain::bestToSolution(){
    return this->gl->toSolution();
}
