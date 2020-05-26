#include "louvainsselv.h"

LouvainSSELV::LouvainSSELV(LargeGraph * graph, GraphCoarsenerSingleStepLG * gcss, unsigned type)
{

    this->type = type;
    this->graph = graph;
    this->gl = gcss;

    //preparing nodes
    this->size = this->graph->numberOfNodes;

    for (unsigned i=0;i<this->graph->numberOfNodes;i++){
        this->communities.push_back(CommunityLouvainSSELV());
        this->commByNode.push_back(i);
    }

    //Define the nodes of each node.
    /*for (unsigned i=0;i<this->graph->numberOfNodes;i++){
        if(this->gl->nodes[i].nodes.size()>0){
            unordered_set<unsigned>::iterator itA = this->gl->nodes[i].nodes.begin();
            while (itA != this->gl->nodes[i].nodes.end()){
                this->commByNode[*itA] = i;
                itA ++;
            }
        }

    }*/



    //preparing temporary
    this->neighWeight.resize(size,-1);
    this->neighPos.resize(size);
    this->neighLast=0;

}

//stabilizes the data from GraphCoarsener used in SimpleStepGreedy
void LouvainSSELV::update(GraphCoarsenerSingleStepLG * gcss){

    this->gl = gcss;
    this->size = this->gl->size;
    for (unsigned node=0;node<this->gl->size; node++){
        this->commByNode[node] = node;
        this->communities[node].internalEdges = this->gl->nodes[node].internalEdges;
        this->communities[node].degree = this->gl->nodes[node].totalDegree;
        this->communities[node].nnodes = this->gl->nodes[node].nodes.size();
        this->communities[node].density = this->gl->nodes[node].density;
    }
    this->bestDensity = calculateDensity();

}



long double LouvainSSELV::calculateDensity(){
    long double dens  = 0.0;
    for (unsigned i=0 ; i<this->gl->size ; i++) {
        if (this->communities[i].nnodes>0){
            dens += (4.0*this->communities[i].internalEdges-this->communities[i].degree)
                    /this->communities[i].nnodes;


        }
    }
    return dens;
}

long double LouvainSSELV::calculateModularity(){
    long double mod  = 0.0;

    long double m = this->graph->numberOfEdges;
    long double m2 = this->graph->numberOfEdges*2.0;
    //##we can avoid O(n) computation here by updating the vector communities at each pass
//cout<<"\n";
//unsigned n=0, e=0;
    for (unsigned i=0 ; i<this->gl->size ; i++) {
        if (this->communities[i].nnodes>0){
//e+=this->communities[i].internalEdges;
//n+=this->communities[i].nnodes;
            mod += this->communities[i].internalEdges/m
                    - (this->communities[i].degree/m2)*(this->communities[i].degree/m2);
//cout<<":["<<i<<"] ie: "<<this->communities[i].internalEdges<<" / "<<m<<" - "<<this->communities[i].degree<<" / "<<m2;
        }
    }
/*    cout<<"\nnodes: "<<n;
    if (e> this->graph->numberOfEdges){
        cout<<"\n***************************************************************";
    }*/
    return mod;

}



bool LouvainSSELV::oneLevel(){


    unsigned i, randpos, inode, node;
    bool improvement=false;
    //long double minModularity = 23;//0.0000000001; --> 1% -->30 segundos/450
    //long double minModularity = 230;//0.0000000001;  --> 10% --> 60 segundos/435
    long double minModularity = 0.0000000001; // --> 10% -->
    int nmoves;
    int npass = 0;
    long double newDensity;

    newDensity = this->calculateDensity();

    double currDensity   = newDensity;
    unsigned numberOfCommunities = gl->size;

    long double delta = 0.0;

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

        for (i=0 ; i<this->neighLast ; i++) {

          long double increase;
          if (this->type == MLSSLV_DENS){
              increase = this->calculateDensityGain(node, neighPos[i], neighWeight[neighPos[i]], wdegree);
          }else{
              increase = this->calculateModularityGain(node, neighPos[i], neighWeight[neighPos[i]], wdegree);
          }

          if (increase>bestIncrease) {
            bestComm     = neighPos[i];
            bestNLinks  = neighWeight[neighPos[i]];
            bestIncrease = increase;

          }
        }
        // insert node in the nearest community
        this->insertInCommunity(node, bestComm, bestNLinks);

        if (bestComm!=nodeComm){
            if (this->communities[nodeComm].nnodes == 0){
                numberOfCommunities--;
            }
            if (this->communities[bestComm].nnodes == 1){
                numberOfCommunities++;
            }
            nmoves++;
        }

      }


      long double temp;

      temp = this->calculateDensity();

      newDensity=temp;


      if (this->bestDensity < temp){
          //if (this->type == MLSSLV_DENS){
          this->bestDensity = temp;
          this->bestMod = this->calculateModularity();
          //}else{
          //    this->bestDensity = temp;
          //    this->bestMod = this->calculateModularity();
          //}
#ifdef GROUND_TRUTH
          this->bestCommunityStr = this->currentToString();
//cout<<"\nONELEVEL{ "<<this->bestCommunityStr<<"}";cout.flush();
#endif

          this->bestNumberOfCommunities = numberOfCommunities;
          improvement=true;
      }
      /*if (nmoves>0){
        improvement=true;
      }*/

    } while (nmoves>0  && newDensity-currDensity>minModularity);

    return improvement;
}


string LouvainSSELV::currentToString(){
    vector <string> comms;
    for (unsigned i=0; i<this->commByNode.size();i++){
        comms.push_back("");
    }
    for (unsigned i=0; i<this->commByNode.size();i++){

        unordered_set<unsigned>::iterator node = this->gl->nodes[i].nodes.begin();
        while(node != this->gl->nodes[i].nodes.end()){

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



void LouvainSSELV::neighsPrepare(unsigned node){

    unsigned neigh, i, neighComm;
    long double neighW;
    for (i=0 ; i<this->neighLast ; i++){
        this->neighWeight[this->neighPos[i]]=-1;
    }
    this->neighLast=0;

    neighPos[0]=this->commByNode[node];
    neighWeight[neighPos[0]]=0;
    neighLast=1;

    itemLDE<EdgeCoarsenerSingleStepLG> * auxAdj = this->gl->nodes[node].adjacencies.getInicio();
    while (auxAdj != NULL){
//cout<<"\n\t\t ...> 154";cout.flush();
        neigh = auxAdj->id.anotherEndpoint;
//cout<<"\n\t\t ...> 155"<<neigh<<", "<<commByNode.size();cout.flush();
        neighComm   = this->commByNode[neigh];
//cout<<"\n\t\t ...> 156";cout.flush();
        neighW = auxAdj->id.weight;
//cout<<"\n\t\t ...> 158";cout.flush();
        if (neigh!=node) {
            if (this->neighWeight[neighComm]==-1) {
                neighWeight[neighComm]=0.;
                neighPos[neighLast]=neighComm;
                neighLast++;
            }
            neighWeight[neighComm]+=neighW;
        }
//cout<<"\n\t\t ...> 167";cout.flush();
        auxAdj = auxAdj->prox;
//cout<<"\n\t\t ...> 169";cout.flush();
    }
//cout<<"\nNL: "<<neighLast;

}

void LouvainSSELV::removeFromCommunity(unsigned node, unsigned comm, long double dnodecomm) {
    this->communities[comm].degree -= this->gl->nodes[node].totalDegree;
    this->communities[comm].internalEdges = this->communities[comm].internalEdges - (dnodecomm + this->gl->nodes[node].internalEdges);
    this->communities[comm].nnodes -= this->gl->nodes[node].nodes.size();
    this->commByNode[node]  = NO_COMMUNITY;
}

void LouvainSSELV::insertInCommunity(unsigned node, unsigned comm, long double dnodecomm) {
    this->communities[comm].degree += this->gl->nodes[node].totalDegree;
    this->communities[comm].nnodes += this->gl->nodes[node].nodes.size();
    this->communities[comm].internalEdges += dnodecomm + this->gl->nodes[node].internalEdges;
    this->commByNode[node]=comm;
}

long double LouvainSSELV::calculateDensityGain(unsigned node, unsigned comm, long double dnodecomm, long double wdegree) {
    return  (this->communities[comm].nnodes*4.0*dnodecomm
                  -this->communities[comm].nnodes*wdegree
                  -4.0*this->communities[comm].internalEdges*this->gl->nodes[node].nodes.size()
                  +this->communities[comm].degree*this->gl->nodes[node].nodes.size())/
                  (this->communities[comm].nnodes*this->communities[comm].nnodes+this->communities[comm].nnodes*this->gl->nodes[node].nodes.size());

}


long double LouvainSSELV::calculateModularityGain(unsigned node, unsigned comm, long double dnodecomm, long double wdegree) {
    long double m2   = this->graph->numberOfEdges*2.0;
    return (dnodecomm - this->communities[comm].degree * wdegree / m2);
}

void LouvainSSELV::execute(vector<unsigned> & commByVertice){
    chrono::system_clock::time_point before, before2;
    before = chrono::system_clock::now();



    this->it = 0;

    this->bestDensity = this->calculateDensity();
    this->bestMod = this->calculateModularity();
    this->bestCommunityStr = this->currentToString();

    this->bestNumberOfCommunities = this->gl->size;
    long double newDensity;
    bool improvement;

    do {
        this->it++;
        //Starting OneLevel..."
        improvement = this->oneLevel();
        newDensity = this->bestDensity;

        //Update
        if(improvement){
            this->updateLouvainGraphAndCommunities(commByVertice);
        }

        if (it==1){ // do at least one more computation if partition is provided
           improvement=true;
        }
    } while(improvement);
    this->totalTime = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();
}


void LouvainSSELV::updateLouvainGraphAndCommunities(vector<unsigned> &commByVertice){

    unsigned node, i, ic, in,master, ileft, iright;

    // Renumber communities
//    before = chrono::system_clock::now();
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
//    this->timeUpRenumber += chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();



    //first, the merge is made over the nodes
//    before = chrono::system_clock::now();
    vector< pair<unsigned, unsigned> > trans;
    unordered_set<unsigned>::iterator itNode;
    for (ic=0;ic<comm_nodes.size(); ic++){

        master = comm_nodes[ic][0];
        for (in=1;in<comm_nodes[ic].size(); in++){
            node = comm_nodes[ic][in];

            master = gl->merge(master, node, 0.0);

        }
        swap(this->communities[master],this->communities[this->commByNode[master]]);
//        translator[ic] = ic;

        //updating commByVertice
        itNode = gl->nodes[master].nodes.begin();
        while (itNode != gl->nodes[master].nodes.end()){
            commByVertice[*itNode] = master;
            itNode++;
        }
    }
//    this->timeUpMerge += chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

//    before = chrono::system_clock::now();
    ileft=0;
    iright=gl->size-1;
    while (ileft < iright){

        if (gl->nodes[ileft].nodes.size() > 0 ){
//            translator[ileft] = ileft;
            ileft++;
            continue;
        }
        if (gl->nodes[iright].nodes.size() == 0 ){
            iright--;
            continue;
        }

        swap(gl->nodes[ileft], gl->nodes[iright]);
        trans.push_back(pair<unsigned, unsigned>(ileft,iright));
//        translator[iright] = ileft;
        ileft++;
        iright--;
    }
//    this->timeUpSwap += chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

//   before = chrono::system_clock::now();
   for (node=0;node<comm_nodes.size();node++){
       this->commByNode[node] = node;
       this->communities[node].internalEdges = gl->nodes[node].internalEdges;
       this->communities[node].degree = gl->nodes[node].totalDegree;
       this->communities[node].nnodes = gl->nodes[node].nodes.size();
       this->communities[node].density = gl->nodes[node].density;
   }
//   this->timeUpCopy += chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

    //translating
//    before = chrono::system_clock::now();
    gl->size = comm_nodes.size();
    itemLDE<EdgeCoarsenerSingleStepLG> *aux;
    vector<pair<unsigned, unsigned>>::iterator it;
    it = trans.begin();
    while(it != trans.end()){
        aux = gl->nodes[(*it).first].adjacencies.getInicio();
        while(aux!=NULL){
            aux->id.apointItem->id.anotherEndpoint = (*it).first;
            aux=aux->prox;
        }

       //updating commByVertice
       itNode = gl->nodes[(*it).first].nodes.begin();
       while (itNode != gl->nodes[(*it).first].nodes.end()){
           commByVertice[*itNode] = (*it).first;
           itNode++;
       }
       it++;
    }
//    this->timeUpTranslate += chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-before).count();

    /*for (in = 0;in<gl->size; in++){
        itemLDE<EdgeCoarsenerSingleStepLG> *aux = gl->nodes[in].adjacencies.getInicio();
        while(aux!=NULL){
            aux->id.anotherEndpoint = translator[aux->id.anotherEndpoint];
            aux=aux->prox;
        }
    }*/

}

