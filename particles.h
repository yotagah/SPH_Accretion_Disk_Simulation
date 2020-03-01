#ifndef __PARTICLES__
#define __PARTICLES__

#define ALOC_STEP 100 // Allocates memory to each quantity of particles / pairs

// Presets for particles
particles_t newParticles();

// Add new particle
int addParticle(particles_t *parts, queue_t *queue, double x, double y, double z, double vx, double vy, double vz, double m, double e_0, double h_size);

// Allocates more memory for particles
void alocMoreParticles(particles_t *parts);

// Presets for interaction pairs
int_pairs_t newIntPairs();

// Allocates more memory for interaction pairs
void alocMorePairs(int_pairs_t *pairs);

// Presets for particles in the queue
queue_t newQueue();

// Allocates more memory for the particles in the queue
void alocMoreQueue(queue_t *queue);

// Add an inactive particle to the inactive queue
void pushInactiveToQueue(queue_t *queue, int part);

// Remove an inactive particle from the inactive queue
int popInactiveFromQueue(queue_t *queue);

// Tests whether any particles fell on a star or left the system
void testInactive(particle_t *part, int id, queue_t *queue, stars_t *stars, double t);

#endif
