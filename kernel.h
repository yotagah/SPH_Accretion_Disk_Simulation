#ifndef __KERNEL__
#define __KERNEL__

// Assign value of the weight function and gradient to an interaction pair
void kernel(double m_dist, double dist[DIM], int_pair_t *int_pair, double h_size);

// Calculates the weight function and gradient
double weightAndGrad(double q, double dist[DIM], double m_dist, double grad[DIM], double h_size);

// Calculates only the weight function
double weight(double q);

#endif