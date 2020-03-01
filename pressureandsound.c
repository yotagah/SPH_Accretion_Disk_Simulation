#include <math.h>
#include "datatypes.h"
#include "pressureandsound.h"

// Calculates the pressure and speed of sound for the particles (Equation 21 and 22, Simpson 1995)
void pressureAndSoundSpeed(particles_t *parts, double gamma) {

	// gamma: adiabatic index

	int k;

	for(k = 0; k < parts->quant; k++) { // All particles
		if(parts->particle[k].active) {

			parts->particle[k].p = (gamma - 1.0) * parts->particle[k].rho * parts->particle[k].e; // Eq. 21 Simpson 1995

			parts->particle[k].c = sqrt(gamma * (gamma - 1.0) * parts->particle[k].e); // Eq. 22 Simpson 1995
		}
	}

}