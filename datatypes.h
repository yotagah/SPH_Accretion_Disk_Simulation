#ifndef __DATATYPES__
#define __DATATYPES__

#define DIM 3					// Number of dimensions of the problem
#define G 6.67428e-11			// [m^3/kg*s] universal constant of gravitation
#define MASSA_SOL 1.98892e30	// [kg] mass of the sun

extern double norm; // Weight function normalization
extern double dt; // Integration interval

// Particle
typedef struct particle_t {
	double a[DIM];				// Acceleration		(m/s^2)
	double v[DIM];				// Speed			(m/s)
	double r[DIM];				// Position			(m)

	double v_prev[DIM];			// Previous speed (for LeapFrog integration)

	double m; 					// Mass				(kg)
	double p;					// Pressure			(N/m^3)
	double rho;					// Density			(kg/m^3)
	double e;					// Energy			(J)
	double dedt;				// Energy var.		(J/s)
	double c;					// Sound speed 		(m/s)

	double e_prev;				// Previous energy

	double h;					// Smoothing size

	short active;				// Active particle?
} particle_t;

// Particles
typedef struct particles_t {
	int quant; 				// Particles quantity
	particle_t *particle; 	// Particles pointer
	int alocated;			// Alocated memory quantity
} particles_t;

// Inactive particles
typedef struct queue_t {
	int end; 				// Last of queue
	int first;				// First of queue
	int *particle; 			// Particles pointer
	int alocated;			// Alocated memory quantity
} queue_t;

// Interaction pair
typedef struct int_pair_t {
	double w;				// Weight function (0.0 to 1.0)
	double dwdx[DIM];		// Derivative weight function (-1.0 to 1.0)
	int i;					// First particle of the pair
	int j;					// Second particle of the pair
} int_pair_t;

// Interaction pairs
typedef struct int_pairs_t {
	int quant; 				// Pairs quantity
	int_pair_t *int_pair;	// Pairs pointer
	int alocated;			// Alocated memory quantity
} int_pairs_t;

// Stars
typedef struct stars_t {
	double Gm1, Gm2;			// Mass times constant of gravitation
	double rad1, rad2;			// Primary and secondary radius (m)
	double p1[DIM], p2[DIM];	// Positions (m)
	double w_orb;				// Orbital angular speed (rad/s) 2PI/p_orb
	double mia;					// Center of mass (m) (mi * a)
	double l1;					// L1 point position (m)
	double a;					// Distance between primary and secondary
} stars_t;

#endif
