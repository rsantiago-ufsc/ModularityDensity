#ifndef MODULARITYLG_H
#define MODULARITYLG_H

#include "../graph/largegraph.h"
#include "../graph/solution.h"

class ModularityLG
{
private:
    LargeGraph * grafo;
public:

    ModularityLG(LargeGraph * grafo);


    //DENSITY
    long double calculateDensity(Solution*sol);
    long double calculateDensity(Solution*sol, float lambda);

    //Modularity
    long double calculateModularityND_NW(Solution * sol);


};

#include "modularitylg.cpp"




#endif // MODULARITYLG_H
