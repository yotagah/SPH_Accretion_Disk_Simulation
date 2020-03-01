#ifndef __FORCE__
#define __FORCE__

// Zerate the acceleration of all particles
void zerateForces(particles_t *parts);

// Calculates the acceleration resulting from internal material forces
void internalForce(particles_t *parts, int_pairs_t *pairs);

// Calculates acceleration due to external forces
void externalForce(particles_t *parts, stars_t *stars);

#endif