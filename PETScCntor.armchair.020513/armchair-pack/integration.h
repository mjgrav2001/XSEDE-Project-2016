#ifndef INTEGRATION_H_
#define INTEGRATION_H_
#include <cmath>
using namespace std;

double integrate(Parameters p, double IVvalues[], double vsd)
{
	//Initializations and declarations
	double currtofe[p.numEnergySteps],
		  c1 = 14./45, c2 = 64./45, c3 = 24./45, c4 = c2, c5 = c1,
		  corrsd = 0, therm = .03,
		  energy, fermi1, fermi2, emu1 = -vsd/2, emu2 = vsd/2;
	
	for(int i=0; i<p.numEnergySteps; i++) {
		energy = p.energyMin + p.energyStep*i;
		fermi1 = 1 / (exp((energy-emu1) / therm) + 1);
		fermi2 = 1 / (exp((energy-emu2) / therm) + 1);
		if(abs(energy) < p.energyStep) currtofe[i] = 0;
		else currtofe[i] = 2 * IVvalues[i] * (fermi1-fermi2);
	}

	for(int i=0; i<p.numEnergySteps-4; i+=4) {
		corrsd += p.energyStep * (c1*currtofe[i] + c2*currtofe[i+1] + c3*currtofe[i+2] + c4*currtofe[i+3] + c5*currtofe[i+4]);
	}

	return corrsd;
} //end integrate

#endif /* INTEGRATION_H_ */

