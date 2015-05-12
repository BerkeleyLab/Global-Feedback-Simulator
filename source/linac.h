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
 * Data structure storing the parameters for a single Linac
 * section comprised of multiple Cryomodules.
 */
typedef struct str_Linac {
	// Parameters used by double compress
	double dE, R56, T566, phi, lam, s0, a, L;

	// Array of Cryomodules in a Linac
	int n_cryos;
	Cryomodule **cryo_net;

} Linac;


typedef struct str_Linac_State {
	Cryomodule_State **cryo_state_net;
	double amp_error, phase_error;
} Linac_State;

// Rename pointer for SWIG work around
typedef Linac* Linac_p;
typedef Linac ** Linac_dp;

Cryomodule *Get_Cryomodule(Linac *linac, int index);
Cryomodule_State *Get_Cryo_State(Linac_State *linac, int index);
void Linac_Allocate_In(Linac *linac, Cryomodule_dp cryo_net, int n_cryos,
	double dE, double R56, double T566, double phi, double lam, double s0, double a, double L);
Linac * Linac_Allocate_New(Cryomodule_dp cryo_net, int n_cryos,
	double dE, double R56, double T566, double phi, double lam, double s0, double a, double L);
void Linac_Deallocate(Linac *linac);
Linac_dp Linac_Allocate_Array(int n);
void Linac_Append(Linac** linac_arr, Linac* linac, int index);
void Linac_State_Allocate(Linac_State *linac_state, Linac *linac);
void Linac_State_Deallocate(Linac_State * linac_state, Linac *linac);

double complex Linac_Step(Linac *linac, Linac_State *linac_state, double delta_tz, double beam_charge,\
	double *amp_error, double *phase_error);

/*
 * Data structure storing the beam parameters on
 * exit from gun for Double Compress
 */

typedef struct str_Gun {
 	// Parameters used by double compress
	// (and for calculating dE_E)
	double E, sz0, sd0;
	double Q;
} Gun;

void Gun_Allocate_In(Gun *gun, double E, double sz0, double sd0, double Q);
Gun *Gun_Allocate_New(double E, double sz0, double sd0, double Q);
// Model for Gun is currently Phase-space only and therefore no state support is present
// Also no Deallocate routine since it does not contain any pointers

#endif
