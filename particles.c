#include <stdio.h>
#include <stdlib.h>
#include "datatypes.h"
#include "particles.h"
#include "vector.h"

// Presets for particles
particles_t newParticles() {
	particles_t new;

	new.quant = 0;
	new.alocated = ALOC_STEP;
	new.particle = (particle_t*) malloc(sizeof(particle_t)*ALOC_STEP);

	return new;
}

// Add new particle
int addParticle(particles_t *parts, queue_t *queue, double x, double y, double z, double vx, double vy, double vz, double m, double e_0, double h_size) {

	int inactive;

	particle_t new;

	new.r[0] = x;
	new.r[1] = y;
	new.r[2] = z;

	new.v[0] = vx;
	new.v[1] = vy;
	new.v[2] = vz;

	new.m = m;
	new.e = e_0;

	new.h = h_size;

	new.active = 1;

	inactive = popInactiveFromQueue(queue);

	if(inactive < 0) { // If there are no inactive particles to reuse
		inactive = parts->quant;
		parts->quant += 1;

	    if(parts->quant > parts->alocated) { // If there is no memory allocated for that amount
	    	alocMoreParticles(parts); // Allocate more memory
	    }

		parts->particle[parts->quant-1] = new;
	} else {
		inactive = queue->particle[inactive];
		parts->particle[inactive] = new;
	}

	return inactive;
}

// Allocates more memory for particles
void alocMoreParticles(particles_t *parts) {
	parts->alocated += ALOC_STEP;
	parts->particle = (particle_t*) realloc(parts->particle, sizeof(particle_t)*parts->alocated);
}

// Presets for interaction pairs
int_pairs_t newIntPairs() {
	int_pairs_t new;

	new.quant = 0;
	new.alocated = ALOC_STEP;
	new.int_pair = (int_pair_t*) malloc(sizeof(int_pair_t)*ALOC_STEP);

	return new;
}

// Allocates more memory for interaction pairs
void alocMorePairs(int_pairs_t *pairs) {
	pairs->alocated += ALOC_STEP;
	pairs->int_pair = (int_pair_t*) realloc(pairs->int_pair, sizeof(int_pair_t)*pairs->alocated);
}

// Presets for particles in the queue
queue_t newQueue() {
	queue_t new;

	new.end = 0;
	new.first = 0;
	new.alocated = ALOC_STEP;
	new.particle = (int*) malloc(sizeof(int)*ALOC_STEP);

	return new;
}

// Allocates more memory for the particles in the queue
void alocMoreQueue(queue_t *queue) {
	queue->alocated += ALOC_STEP;
	queue->particle = (int*) realloc(queue->particle, sizeof(int)*queue->alocated);
}

// Add an inactive particle to the inactive queue
void pushInactiveToQueue(queue_t *queue, int part) {
	int k, ant;
	int *copy;

	ant = queue->end;
	queue->end = (queue->end+1) % queue->alocated;
    if(queue->end == queue->first) { // If there is no memory allocated for that amount
		copy = (int*) malloc(sizeof(int)*queue->alocated);
		for(k = 0; k < queue->alocated; k++) {
			copy[k] = queue->particle[k];
		}
		for(k = 0; k < queue->alocated; k++) {
			queue->particle[k] = copy[(k+queue->first)%queue->alocated];
		}
		free(copy);

		queue->first = 0;
		queue->end = queue->alocated;
		ant = queue->alocated-1;

		alocMoreQueue(queue); // Allocate more memory
    }
    queue->particle[ant] = part;
}

// Remove an inactive particle from the inactive queue
int popInactiveFromQueue(queue_t *queue) {
	int aux;

	if(queue->end == queue->first) // If the queue is empty
		return -1;
	else { // If not
		aux = queue->first;
		queue->first = (aux+1) % queue->alocated;
		return aux;
	}
}

// Tests whether any particles fell on a star or left the system
void testInactive(particle_t *part, int id, queue_t *queue, stars_t *stars, double t) {

    double dist1[DIM], dist2[DIM];
    
    subVector(part->r, stars->p1, dist1);
   	subVector(part->r, stars->p2, dist2);
    if(modVector(dist1) < stars->rad1) { // Fell on the primary star
    	FILE *accret = fopen("accret1.dat", "a");
	    fprintf(accret, "%g\n", t); // Append the value of current time to the accretion ratio file data of the primary
    	fclose(accret);
	}
	else if(modVector(dist2) < stars->rad2) { // Fell on the secondary star
		FILE *accret = fopen("accret2.dat", "a");
	    fprintf(accret, "%g\n", t); // Append the value of current time to the accretion ratio file data of the secondary
	    fclose(accret);
	}
	else if(abs(part->r[2]) < stars->a/4.0 && abs(part->r[0]) < stars->a && abs(part->r[1]) < stars->a) { // Tests if the particle is "inside" the system
		return; // If it is, finished the function
	}

	// If the particle leave the system, inactivate and move it to the queue
	pushInactiveToQueue(queue, id);
	part->active = 0;
}