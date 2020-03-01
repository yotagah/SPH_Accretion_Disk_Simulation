#ifndef __INTEGRATION__
#define __INTEGRATION__

// Calculate the values for a time step
void singleStep(particles_t *particles, int_pairs_t *int_pairs, stars_t *stars);

// Find the particles that will make the interaction (that are within the SPH radius) and calculate the weight function for each pair
void directFind(particles_t *particles, int_pairs_t *int_pairs);

#endif