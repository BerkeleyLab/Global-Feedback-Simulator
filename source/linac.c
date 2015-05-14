#include "linac.h"

Linac_dp Linac_Allocate_Array(int n)
{
  Linac_dp linac_net = calloc(n, sizeof(Linac *));
  return linac_net;
}

void Linac_Append(Linac** linac_arr, Linac* linac, int index)
{
  // XXX Add some check!!
  linac_arr[index] = linac;
}

void Linac_Allocate_In(Linac *linac, Cryomodule_dp cryo_net, int n_cryos,
	double dE, double R56, double T566, double phi, double lam, double s0, double a, double L)
{
	linac -> cryo_net = cryo_net;
	linac -> n_cryos = n_cryos;
	linac -> dE = dE;
	linac -> R56 = R56;
	linac -> T566 = T566;
	linac -> phi = phi;
	linac -> lam = lam;
	linac -> s0 = s0;
	linac -> a = a;
	linac -> L = L;

}

Linac *Linac_Allocate_New(Cryomodule_dp cryo_net, int n_cryos,
	double dE, double R56, double T566, double phi, double lam, double s0, double a, double L)
{
	Linac *linac = calloc(1,sizeof(Linac));
	Linac_Allocate_In(linac, cryo_net, n_cryos,
		dE, R56, T566, phi, lam, s0, a, L);
	return linac;
}

void Linac_Deallocate(Linac *linac)
{
	for(int i=0;i<linac->n_cryos;i++){
		Cryomodule_Deallocate(linac->cryo_net[i]);
	}
	free(linac->cryo_net);

	linac->dE = 0.0;
	linac->R56 = 0.0;
	linac->T566 = 0.0;
	linac->phi = 0.0;
	linac->lam = 0.0;
	linac->s0 = 0.0;
	linac->a = 0.0;
	linac->L = 0.0;
}

Cryomodule *Get_Cryomodule(Linac *linac, int index)
{
	if(linac->cryo_net[index]) return linac->cryo_net[index];
	else return NULL;
}

Cryomodule_State *Get_Cryo_State(Linac_State *linac_state, int index)
{
	if(linac_state->cryo_state_net[index]) return linac_state->cryo_state_net[index];
	else return NULL;
}

void Linac_State_Allocate(Linac_State *linac_state, Linac *linac)
{

	// Allocate memory for the array of states
	linac_state->cryo_state_net = calloc(linac->n_cryos, sizeof(Cryomodule_State *));

	// Then allocate memory for the States themselves
	int i;
 	for(i=0;i<linac->n_cryos;i++) {
 		linac_state->cryo_state_net[i] = (Cryomodule_State*)calloc(1,sizeof(Cryomodule_State));
 		Cryomodule_State_Allocate(linac_state->cryo_state_net[i], linac->cryo_net[i]);
 	}

 	linac_state->linac_V = 0.0;
 	linac_state->linac_Kg = 0.0;
}

void Linac_State_Deallocate(Linac_State *linac_state, Linac *linac)
{
	for(int i=0;i<linac->n_cryos;i++){
		Cryomodule_State_Deallocate(linac_state->cryo_state_net[i], linac->cryo_net[i]);
	}
	free(linac_state->cryo_state_net);
}

double complex Linac_Step(Linac *linac, Linac_State *linac_state, double delta_tz, double beam_charge,
	double *amp_error, double *phase_error)
{
	// Overall Linac accelerating voltage (and its projection on the beam)
	double complex linac_V=0.0, linac_V_beam=0.0;
	double complex linac_Kg=0.0;

	// Run state-space simulation step
	// Iterate over the Cryomodules
	for(int i=0;i<linac->n_cryos;i++){
		linac_V += Cryomodule_Step(linac->cryo_net[i], linac_state->cryo_state_net[i], delta_tz, beam_charge);
		linac_Kg += linac_state->cryo_state_net[i]->cryo_Kg;
	}

	// Project the Linac voltage over the beam (phi radians relative RF phase to the beam),
	linac_V_beam = linac_V*cos(linac->phi);
	// Any deviations into the imaginary component of this vector is a phase error
	*phase_error = carg(linac_V);

	// Now calculate relative amplitude error
	// - in relative units with respect to the Linac Energy increase,
	// - eV considered equivalent to V here since we are dealing with Electrons.

	if(linac->dE!=0.0){ // Avoid division by 0
		*amp_error = (cabs(linac_V_beam)-linac->dE)/linac->dE;
	} else{
		*amp_error = cabs(linac_V_beam)-linac->dE;
	}

	// Store total Linac drive signal (vector sum of all RF Station drive signals)
	linac_state->linac_Kg = linac_Kg;

	// Total Linac Accelerating (store and return)
	linac_state->linac_V = linac_V;
	return linac_V;
}

void Gun_Allocate_In(Gun *gun, double E, double sz0, double sd0, double Q)
{
	gun->E = E;
	gun->sz0 = sz0;
	gun->sd0 = sd0;
	gun->Q = Q;
}

Gun *Gun_Allocate_New(double E, double sz0, double sd0, double Q)
{
	Gun *gun = calloc(1,sizeof(Gun));
	Gun_Allocate_In(gun, E, sz0, sd0, Q);
	return gun;
}
