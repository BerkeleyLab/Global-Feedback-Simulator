#ifndef LINAC_H
#define LINAC_H


/*
 *  linac.h
 *
 * Definitions for a Linac section model
 *
 */

#include <complex.h>
#include "cryomodule.h"

/*
 * Data struture storing the parameters for a single Linac
 * section comprised of multiple Cryomodules.
 */
typedef struct str_Linac {
	// Paramters used by double compress
	double dE, R56, T566, phi, lam, s0, a, L;

	// Array of Cryomodules in a Linac
	int n_cryos;
	Cryomodule **cryo_net;

} Linac;


typedef struct str_Linac_State {
	Cryomodule_State **cryo_state_net;
} Linac_State;


// Rename pointer for SWIG work around
typedef Linac* Linac_p;
typedef Linac ** Linac_dp;

Cryomodule *Get_Cryomodule(Linac *linac, int index);
void Linac_Allocate_In(Linac *linac, Cryomodule_dp cryo_net, int n_cryos,
	double dE, double R56, double T566, double phi, double lam, double s0, double a, double L);
Linac * Linac_Allocate_New(Cryomodule_dp cryo_net, int n_cryos,
	double dE, double R56, double T566, double phi, double lam, double s0, double a, double L);
void Linac_Deallocate(Linac *linac);
void Linac_State_Allocate(Linac_State *linac_state, Linac *linac);
// void Linac_Step(Linac *linac, Linac_State * linac_state,
	// double delta_tz, double complex beam_charge, double complex probe_ns, double complex rev_ns,
	// int openloop);
void Linac_State_Deallocate(Linac_State * linac_state, Linac *linac);

/*
 * Data struture storing the beam parameters on 
 * exit from gun for Double Compress
 */

typedef struct str_Gun {
 	// Paramters used by double compress
	// (and for calculating dE_E)
	double E, sz0, sd0;
	double Q;
} Gun;

void Gun_Allocate_In(Gun *gun, double E, double sz0, double sd0, double Q);
Gun *Gun_Allocate_New(double E, double sz0, double sd0, double Q);
void Gun_Deallocate(Gun *gun);
// Model for Gun is currently Phase-space only and therefore no state support is present

#endif
