#include <stdio.h>

#include "datatypes.h"
#include "density.h"
#include "kernel.h"

// Calculates density by simple sum (Equation 4.26 Liu 2003)
void sumDensity(particles_t *parts, int_pairs_t *pairs) {

    int k;
    int i, j; // Particles

    for(k = 0; k < parts->quant; k++) { // Sum the particle itself
        if(parts->particle[k].active) {
            parts->particle[k].rho = weight(0.0) * parts->particle[k].m;
        }
    }

    for(k = 0; k < pairs->quant; k++) { // Sums the other interacting particles

        i = pairs->int_pair[k].i;
        j = pairs->int_pair[k].j;
        
        parts->particle[i].rho += parts->particle[i].m * pairs->int_pair[k].w;
        parts->particle[j].rho += parts->particle[j].m * pairs->int_pair[k].w;
    }
}
