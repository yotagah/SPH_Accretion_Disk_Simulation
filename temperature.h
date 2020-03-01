#ifndef __TEMPERATURE__
#define __TEMPERATURE__

// Calculates energy dissipated by artificial heat (Eq. 28 Libersky at. Al 1993, High Strain Lagrangian Hydrodynamics)
void artificialHeat(particles_t *particles, int_pairs_t *int_pairs);

// Dissipation of thermal energy by Beta-Cooling radiation, Rice at all 2012
void betaCooling(particles_t *parts, int_pairs_t *pairs, stars_t *stars, double beta);

#endif