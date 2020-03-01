#include <stdio.h>
#include "datatypes.h"
#include "integration.h"
#include "vector.h"
#include "particles.h"
#include "kernel.h"
#include "density.h"
#include "pressureandsound.h"
#include "force.h"
#include "viscosity.h"
#include "temperature.h"

// Calculate the values for a time step
void singleStep(particles_t *particles, int_pairs_t *int_pairs, stars_t *stars) {

    // Interaction pairs
    directFind(particles, int_pairs);

    // Equations
    sumDensity(particles, int_pairs);
    pressureAndSoundSpeed(particles, 1.01); // gamma

    // Forces
    zerateForces(particles);
    artificialViscosity(particles, int_pairs);
    //artificialHeat(particles, int_pairs);
    internalForce(particles, int_pairs);
    externalForce(particles, stars);
    //betaCooling(particles, int_pairs, stars, 0.5);
}

// Find the particles that will make the interaction (that are within the SPH radius) and calculate the weight function for each pair
void directFind(particles_t *particles, int_pairs_t *int_pairs) {

    int i, j;
    int scale_k; // Scale factor for the SPH radius
    double dist[DIM]; // Distance vector
    double m_dist; // Distance module
    double mh; // Average particle interaction distance
    
    scale_k = 2; // Scale factor for: cubic spline kernel by W4 - Spline (Monaghan 1985)

    int_pairs->quant = 0; // Zerate number of interaction pairs

    for(i = 0; i < particles->quant; i++) { // All particles
        if(particles->particle[i].active) { // Active particles only
            for(j = i+1; j < particles->quant; j++) { // All pairs
                if(particles->particle[j].active) { // Active particles only
                    subVector(particles->particle[i].r, particles->particle[j].r, dist); // Find distance vector
                    m_dist = modVector(dist); // Find distance module

                    mh = (particles->particle[i].h + particles->particle[j].h) / 2.0; // Average of smooth distance

                    if(m_dist < scale_k*mh) { // If is within the SPH radius
                        if(int_pairs->quant+1 > int_pairs->alocated) // If there is no memory alocated for the pairs
                            alocMorePairs(int_pairs); // Alocate more memory

                        int_pairs->int_pair[int_pairs->quant].i = i; // First particle of the interaction pair
                        int_pairs->int_pair[int_pairs->quant].j = j; // Second particle of the interaction pair

                        kernel(m_dist, dist, &(int_pairs->int_pair[int_pairs->quant]), mh); // Calculates weight functoin and derivative

                        int_pairs->quant += 1; // Increase number of pairs
                    }
                }
            }
        }
    }
}