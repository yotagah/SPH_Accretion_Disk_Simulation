#include "datatypes.h"
#include "force.h"
#include "vector.h"

// Zerate the acceleration of all particles
void zerateForces(particles_t *parts) {

	int k, alfa;

	for(k = 0; k < parts->quant; k++) {
        if(parts->particle[k].active) {
    		for(alfa = 0; alfa < DIM; alfa++) {
    			parts->particle[k].a[alfa] = 0;
    		}
    		parts->particle[k].dedt = 0;
        }
	}
}

// Calculates the acceleration resulting from internal material forces
void internalForce(particles_t *parts, int_pairs_t *pairs) {
    int k, alfa;
    int i, j; // Particles

    double aux_a, aux_e;

    double dv; // Speed difference
    double rhoirhoj;

	for(k = 0; k < pairs->quant; k++) { // All interaction pairs

        i = pairs->int_pair[k].i;
        j = pairs->int_pair[k].j;

        rhoirhoj = parts->particle[i].rho * parts->particle[j].rho;

        aux_e = 0.0;
    	for(alfa = 0; alfa < DIM; alfa++) {

   			aux_a = (parts->particle[i].p + parts->particle[j].p) * pairs->int_pair[k].dwdx[alfa] / rhoirhoj; // Eq. 4.60 Liu 2003

            dv = parts->particle[i].v[alfa] - parts->particle[j].v[alfa];
            aux_e += aux_a * dv; // Eq. 4.62 Liu 2003

    		parts->particle[i].a[alfa] -= parts->particle[j].m * aux_a; // Eq. 4.60 Liu 2003
    		parts->particle[j].a[alfa] += parts->particle[i].m * aux_a; // Eq. 4.60 Liu 2003
        }

   		parts->particle[i].dedt += 0.5 * parts->particle[j].m * aux_e; // Eq. 4.62 Liu 2003
   		parts->particle[j].dedt += 0.5 * parts->particle[i].m * aux_e; // Eq. 4.62 Liu 2003
    }
}

// Calculates acceleration due to external forces
void externalForce(particles_t *parts, stars_t *stars) {
	int k, alfa;
    double dist1[DIM], dist2[DIM];
    double m_dist1_cub, m_dist2_cub;

	for(k = 0; k < parts->quant; k++) { // All particles
        if(parts->particle[k].active) { // Active particles only

            subVector(parts->particle[k].r, stars->p1, dist1); // Distance to the primary (result to dist1)
            subVector(parts->particle[k].r, stars->p2, dist2); // Distance to the secondary (result to dist2)

            m_dist1_cub = modVector(dist1); // Module of the distance
            m_dist1_cub = m_dist1_cub*m_dist1_cub*m_dist1_cub; // Module^3

            m_dist2_cub = modVector(dist2); // Module of the distance
            m_dist2_cub = m_dist2_cub*m_dist2_cub*m_dist2_cub; // Module^3

            for(alfa = 0; alfa < DIM; alfa++) {

                parts->particle[k].a[alfa] -= stars->Gm1 * dist1[alfa] / m_dist1_cub; // Gravitational acceleration to the primary star
                parts->particle[k].a[alfa] -= stars->Gm2 * dist2[alfa] / m_dist2_cub; // Gravitational acceleration to the secondary star

                switch(alfa) { // Calculates the centripetal and coriolis acceleration (only 2 dimensions)
                    case 0:
                        parts->particle[k].a[0] += stars->w_orb * (stars->w_orb * parts->particle[k].r[0] + 2 * parts->particle[k].v[1]); // Centripetal and coriolis acceleration
                        break;
                    case 1:
                        parts->particle[k].a[1] += stars->w_orb * (stars->w_orb * parts->particle[k].r[1] - 2 * parts->particle[k].v[0]); // Centripetal and coriolis acceleration
                        break;
                }
            }
        }
	}
}