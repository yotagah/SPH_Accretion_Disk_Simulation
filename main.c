#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "datatypes.h"
#include "particles.h"
#include "integration.h"
#include "vector.h"

double norm; // Normalization of the weight function
double dt; // Integration interval
double ini_e; // Initial energy of the particles

double determineL1(double p, double m1, double m2, double dps, double dcm);
void initParticle(particles_t *parts, queue_t *queue, int_pairs_t *pairs, stars_t *stars, double mpart, double h_size, double t);

int main() {

	int i, alfa;

	double t, end;
	double h_size; // Size of the interaction area between particles
	double mpart; // Particles mass
	double pem; // Next mass ejection
	double dm2; // Secondary mass loss rate

	particles_t particles; // Particles
	queue_t queue; // Inactive particles
	int_pairs_t int_pairs; // Interaction pairs
	stars_t stars; // Information about the stars

	int c; // Counter for selection of "frames to print"

	particles = newParticles();
	queue = newQueue();
	int_pairs = newIntPairs();

	stars.Gm1 = 0.84*MASSA_SOL*G; // Primary mass
	stars.Gm2 = 0.125*MASSA_SOL*G; // Secondary mass

	stars.a = 5.12289676e8; // Distance between primary and secondary (Calculated with Kepler's 3rd Law)

	stars.rad1 = 0.078e8; // Primary radius
	stars.rad2 = (0.38 + 0.2*log10(stars.Gm2/stars.Gm1)) * stars.a; // Secondary radius

	stars.mia = stars.Gm2 * stars.a / (stars.Gm1 + stars.Gm2); // Center of mass

	stars.w_orb = sqrt((stars.Gm1 + stars.Gm2) / pow(stars.a, 3)); // Orbital angular speed

	stars.l1 = determineL1(2*M_PI / stars.w_orb, stars.Gm1/G, stars.Gm2/G, stars.a, stars.mia) - stars.mia; // Lagrangian point L1

	// Primary position
	stars.p1[0] = -stars.mia;
	stars.p1[1] = 0.0;
	stars.p1[2] = 0.0;

	// Secondary position
	stars.p2[0] = stars.a-stars.mia;
	stars.p2[1] = 0.0;
	stars.p2[2] = 0.0;

	// Random seeds
	srand48(time(NULL));
	srand(time(NULL));

	// Remove old simulation data files
	system("rm -f accret1.dat");
	system("rm -f accret2.dat");

	h_size = 0.02*stars.a;

	// The normalization constant depends on the number of dimensions of the problem and the weight function choosed
	if(DIM == 2)
		norm = 15.0/(7.0*M_PI*h_size*h_size);
	else if(DIM == 3)
		norm = 3.0/(2.0*M_PI*h_size*h_size*h_size);

	dt = 3.0; // Define the integration step
	t = 0.0; // Initial time

	dm2 = 2.0e12; // Secondary mass loss rate
	mpart = dm2 * (2*M_PI/stars.w_orb) / 250; // Calcule of the particles mass (250 is the particles per orbit wanted)
	pem = 0.0; // Time of the next mass ejection
	ini_e = 8.31 * 2000 / (1.7e-3 * 0.01);  // Calcule of initial energy of the particles (2000 K)

	end = dt/2.0 + (2*M_PI/stars.w_orb) * 10; // Time of the end of simulation (100 orbits)

	c = 0; // Counter for print
	for(t = 0.0; t <= end; t += dt) {

		// LeapFrog integration, Liu 2003 (pg 209)

		for(i = 0; i < particles.quant; i++) {
			if(particles.particle[i].active) {
				particles.particle[i].e_prev = particles.particle[i].e;
				particles.particle[i].e += dt*particles.particle[i].dedt/2.0;
				if(particles.particle[i].e < 0.0) particles.particle[i].e = 0;
				for(alfa = 0; alfa < DIM; alfa++) {
					particles.particle[i].v_prev[alfa] = particles.particle[i].v[alfa];
					particles.particle[i].v[alfa] += dt*particles.particle[i].a[alfa]/2.0;
				}
			}
		}

		if(t*dm2 >= pem) { // Is time for new particle?

			initParticle(&particles, &queue, &int_pairs, &stars, mpart, h_size, t); // Add new particle and proceed with one step (slightly different for the new particle)
            pem += mpart; // Adjust for next ejection

        } else {

        	// Calculate step
			singleStep(&particles, &int_pairs, &stars);

			for(i = 0; i < particles.quant; i++) {
				if(particles.particle[i].active) {
					particles.particle[i].e = particles.particle[i].e_prev + dt*particles.particle[i].dedt;
					if(particles.particle[i].e < 0.0) particles.particle[i].e = 0;
					for(alfa = 0; alfa < DIM; alfa++) {
						particles.particle[i].v[alfa] = particles.particle[i].v_prev[alfa] + dt*particles.particle[i].a[alfa];
						particles.particle[i].r[alfa] += dt*particles.particle[i].v[alfa];
					}
					testInactive(&particles.particle[i], i, &queue, &stars, t); // Check if has accreted or lost particles
				}
			}
		}

		if(c%33 == 0) {
			fprintf(stderr, "%0.2f %% %d\n", (t/end)*100, particles.quant); // Progress and particle count to Standard Error
			for(i = 0; i < particles.quant; i++) {
				if(particles.particle[i].active) { // Print active particles position to Standard Output
					printf("%g %g %g ", particles.particle[i].r[0], particles.particle[i].r[1], particles.particle[i].r[2]);
				}
			}
			printf("\n");
		}
		c++;
	}

	return 0;
}

// Determine the position of L1 (point where the acceleration is zero)
double determineL1(double p, double m1, double m2, double dps, double dcm) {
    // Use Newton's method to find the position of L1
    // p - orbital period
    // m - masses
    // dps - distance between primary and secondary (secondary position)
    // dcm - position of center of mass

    double r;
    double r_p;
    double der;
    double func;

    r = dps / 2.0;

    while(abs(r-r_p) > 1e-200) { // Lower the number, higher the precision
        r_p = r;
        func = G * (m1/pow(r,2) - m2/pow(dps-r,2)) - pow(2*M_PI/p,2) * (r - dcm);
        der = G * (m1*(-2.0)/pow(r,3) - m2*(2.0)/pow(dps-r,3)) - pow(2*M_PI/p,2);
        r = r_p - func/der;
    }
    
    return r;
}

// Initiate a particle in a random position inside a cube with h_size size close to L1 point
void initParticle(particles_t *parts, queue_t *queue, int_pairs_t *pairs, stars_t *stars, double mpart, double h_size, double t) {

	int last, i, alfa;
	double rand_size, pa[DIM];

	// Set a random position inside a cube with size 1
    rand_size = (double) rand() / RAND_MAX;
    pa[0] = h_size * rand_size;
    rand_size = (double) rand() / RAND_MAX - 0.5;
    pa[1] = h_size * rand_size;
	rand_size = (double) rand() / RAND_MAX - 0.5;
	pa[2] = h_size * rand_size;

    last = addParticle(parts, queue, stars->l1*0.99 - pa[0], pa[1], pa[2], 0.0, 0.0, 0.0, mpart, ini_e, h_size);

	singleStep(parts, pairs, stars); // Proceed with a step

	parts->particle[last].e += dt*parts->particle[last].dedt/2.0;
	for(alfa = 0; alfa < DIM; alfa++) {
		parts->particle[last].v[alfa] += dt*parts->particle[last].a[alfa]/2.0;
		parts->particle[last].r[alfa] += dt*parts->particle[last].v[alfa];
	}

	for(i = 0; i < parts->quant; i++) {
		if(i != last && parts->particle[i].active) { // Active particles only and also exclude the last added
			parts->particle[i].e = parts->particle[i].e_prev + dt*parts->particle[i].dedt;
			if(parts->particle[i].e < 0.0) parts->particle[i].e = 0;
			for(alfa = 0; alfa < DIM; alfa++) {
				parts->particle[i].v[alfa] = parts->particle[i].v_prev[alfa] + dt*parts->particle[i].a[alfa];
				parts->particle[i].r[alfa] += dt*parts->particle[i].v[alfa];
			}
			testInactive(&parts->particle[i], i, queue, stars, t); // Check if has accreted or lost particles
		}
	}
}