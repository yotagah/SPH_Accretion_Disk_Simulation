#include <math.h>
#include "datatypes.h"
#include "kernel.h"

// Assign value of the weight function and gradient to an interaction pair
void kernel(double m_dist, double dist[DIM], int_pair_t *int_pair, double h_size) {

    double q; // "Percentage" of the distance in relation to the SPH radius

    q = m_dist/h_size;

    int_pair->w = weightAndGrad(q, dist, m_dist, int_pair->dwdx, h_size); // Weight function and gradient
}

// Calculates the weight function and gradient
double weightAndGrad(double q, double dist[DIM], double m_dist, double grad[DIM], double h_size) {

    int d;

    // Function: cubic spline kernel by W4 - Spline (Monaghan 1985)

    if (q >= 0.0 && q <= 1.0) {
        for(d = 0; d < DIM; d++) { // Gradient
            grad[d] = norm * (-2.0 + 3.0/2.0*q) / (h_size*h_size) * dist[d];
        }
        return norm * (2.0/3.0 - q*q + (q*q*q)/2.0);
    } else if (q > 1.0 && q <= 2.0) {
        for(d = 0; d < DIM; d++) { // Gradient
            grad[d] = -norm * 1.0/6.0 * 3.0*pow(2.0-q,2)/h_size * (dist[d]/m_dist);
        }
        return norm * 1.0/6.0 * pow(2.0-q,3);
    } else {
        for(d=0; d < DIM; d++) { // Gradient
            grad[d] = 0.0;
        }
        return 0.0;
    }
}

// Calculates only the weight function
double weight(double q) {

    // Function: cubic spline kernel by W4 - Spline (Monaghan 1985)
    
    if (q >= 0.0 && q <= 1.0) {
        return norm * (2.0/3.0 - q*q + (q*q*q)/2.0);
    } else if (q > 1.0 && q <= 2.0) {
        return norm * 1.0/6.0 * pow(2.0-q,3);
    } else {
        return 0.0;
    }
}