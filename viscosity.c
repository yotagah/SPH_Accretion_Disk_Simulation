#include <stdio.h>

#include "datatypes.h"

// Calculates acceleration due to artificial viscosity Eq. 23 Simpson 1995
void artificialViscosity(particles_t * parts, int_pairs_t * pairs)
{

	int k, d;
	int i, j;			// Interaction pair

	double vintr;		// Internal product of v and r
	double modr2;		// Module of r squared
	double dv[DIM];		// Speed difference
	double dr;			// Position difference

	double phi;			// Eq. 4.67 Liu 2003
	double cmed;		// Eq. 4.68 Liu 2003
	double rhomed;		// Eq. 4.69 Liu 2003
	double pi;			// Eq. 23 Simpson 1995
	double aux;
	double h_size;

	double alfa = 0.4;

	for (k = 0; k < pairs->quant; k++) {
		i = pairs->int_pair[k].i;
		j = pairs->int_pair[k].j;

		// Calculates v internal r and module of r squared
		vintr = 0.0;
		modr2 = 0.0;
		for (d = 0; d < DIM; d++) {
			dv[d] = parts->particle[i].v[d] - parts->particle[j].v[d];
			dr = parts->particle[i].r[d] - parts->particle[j].r[d];
			vintr += dv[d] * dr;
			modr2 += dr * dr;
		}

		if (vintr <= 0.0) { // Eq. 23 Simpson 1995
			h_size = (parts->particle[i].h + parts->particle[j].h) / 2.0;
			phi = h_size * vintr / (modr2 + 0.01 * h_size * h_size); // Eq. 24 Simpson 1995 -> 0.1*h_size to prevent division by zero
			cmed = (parts->particle[i].c + parts->particle[j].c) / 2.0;
			rhomed = (parts->particle[i].rho + parts->particle[j].rho) / 2.0;

			pi = (-alfa * cmed * phi) / rhomed; // Eq. 23 Simpson 1995

			// Sum acceleration due to viscosity
			for (d = 0; d < DIM; d++) {
				aux =  -pi * pairs->int_pair[k].dwdx[d];
				parts->particle[i].a[d] += aux * parts->particle[j].m;
				parts->particle[j].a[d] -= aux * parts->particle[i].m;

				parts->particle[i].dedt -= 0.5 * aux * dv[d] * parts->particle[j].m; // Eq. 8.28 Liu 2003
   				parts->particle[j].dedt -= 0.5 * aux * dv[d] * parts->particle[i].m;; // Eq. 8.28 Liu 2003
			}
		}
	}
}
