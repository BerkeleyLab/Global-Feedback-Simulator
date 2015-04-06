#include "linac.h"

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
}

void Linac_State_Deallocate(Linac_State *linac_state, Linac *linac)
{
	for(int i=0;i<linac->n_cryos;i++){
		Cryomodule_State_Deallocate(linac_state->cryo_state_net[i], linac->cryo_net[i]);
	}
	free(linac_state->cryo_state_net);
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

void Gun_Deallocate(Gun *gun)
{
	gun->E = 0.0;
	gun->sz0 = 0.0;
	gun->sd0 = 0.0;
	gun->Q = 0.0;

}
