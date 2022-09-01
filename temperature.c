#include <stdlib.h>
#include "datatypes.h"

// Calculates energy dissipated by artificial heat (Eq. 28 Libersky at. Al 1993, High Strain Lagrangian Hydrodynamics)
void artificialHeat(particles_t *parts, int_pairs_t *pairs) {

    int i, j, k, alfa;
    double dr, dv, aux, mrho, aux_divv, mui, muj, muij, rdwdx, g1, g2, modr2;
    double *divv;
    double h_size;

    // Heat conduction parameters
    g1 = 0.5;
    g2 = 1.0;

    divv = (double*) calloc(parts->quant, sizeof(double));

    for(k = 0; k < pairs->quant; k++) {
        i = pairs->int_pair[k].i;
        j = pairs->int_pair[k].j;

        aux_divv = 0;
        for(alfa = 0; alfa < DIM; alfa++) {
            dv = parts->particle[j].v[alfa] - parts->particle[i].v[alfa];
            aux_divv += dv * pairs->int_pair[k].dwdx[alfa];
        }

        divv[i] += parts->particle[j].m * aux_divv/parts->particle[j].rho;
        divv[j] += parts->particle[i].m * aux_divv/parts->particle[i].rho;
    }

    for(k = 0; k < pairs->quant; k++) {
        i = pairs->int_pair[k].i;
        j = pairs->int_pair[k].j;

        mrho = (parts->particle[i].rho + parts->particle[j].rho) / 2.0;
        h_size = (parts->particle[i].h + parts->particle[j].h) / 2.0;

        modr2 = 0.0;
        rdwdx = 0.0;
        for(alfa = 0; alfa < DIM; alfa++) {
            dr = parts->particle[i].r[alfa] - parts->particle[j].r[alfa];
            modr2 += dr * dr;
            rdwdx += dr * pairs->int_pair[k].dwdx[alfa];
        }

        mui = g1 * h_size * parts->particle[i].c + g2 * h_size * h_size * (abs(divv[i]) - divv[i]);
        muj = g1 * h_size * parts->particle[j].c + g2 * h_size * h_size * (abs(divv[j]) - divv[j]);
        muij = (mui + muj);

        aux = muij * rdwdx / (mrho * (modr2 + 0.01*h_size*h_size));

        parts->particle[i].dedt += parts->particle[j].m * aux * (parts->particle[i].e - parts->particle[j].e);
        parts->particle[j].dedt += parts->particle[i].m * aux * (parts->particle[j].e - parts->particle[i].e);
    }
    
    free(divv);
}

// Dissipation of thermal energy by Beta-Cooling radiation, Rice at all 2012
void betaCooling(particles_t *parts, int_pairs_t *pairs, stars_t *stars, double beta) {
    int k, i, j;
    double inv_tau_c;

    inv_tau_c = (stars->w_orb / beta);

    for(k = 0; k < pairs->quant; k++) {
        i = pairs->int_pair[k].i;
        j = pairs->int_pair[k].j;

        parts->particle[i].dedt -= inv_tau_c * parts->particle[i].m * parts->particle[i].e * pairs->int_pair[k].w / parts->particle[i].rho; // Eq. 24 Rice 2012
        parts->particle[j].dedt -= inv_tau_c * parts->particle[j].m * parts->particle[j].e * pairs->int_pair[k].w / parts->particle[j].rho; // Eq. 24 Rice 2012
    }
}
