/**
 * @file linac.c
 * @brief Linac Model: Allocation, configuration and stepping functions for the Linac model, which iterates over an arbitrary number of Cryomodules.
 * @author Carlos Serrano (Cserrano@lbl.gov)
 */

#include "linac.h"

/** Allocates memory for an array of Linac Sections.
  * Linacs themselves need to be allocated and filled individually,
  * and appended to this array using Linac_Append.*/
Linac_dp Linac_Allocate_Array(int n 	///< Size of the array of Linacs
	)
{
  Linac_dp linac_net = calloc(n, sizeof(Linac *));
  return linac_net;
}

/** Append a Linac previously allocated and filled to the array given as an argument. */
void Linac_Append(
	Linac** linac_arr,	///< Pointer to array of Linacs
	Linac* linac,				///< Pointer to Linac to be appended
	int index						///< Position of Linac in the array
	)
{
  // XXX Add some check!!
  linac_arr[index] = linac;
}

/** Takes a pointer to a Linac struct (representing a Linac Section)
  * which has been previously allocated and fills it in with the values passed as arguments. */
void Linac_Allocate_In(
	Linac *linac,							///< Pointer to Linac
	Cryomodule_dp cryo_net,		///< Array of Cryomodules in a Linac
	int n_cryos,							///< Number of Cryomodules in Linac
	double dE,								///< Energy increase in Linac Section in eV
	double R56,								///< Longitudinal dispersion (if any)
	double T566,							///< 2nd-order longitudinal dispersion (if any)
	double phi,								///< Nominal Linac RF phase in radians
	double lam,								///< RF wavelength in meters
	double s0,								///< Wakefield characteristic length (Sband=0.105m, Xband=0.02625m) in meters
	double a,									///< Mean iris radius (Sband=11.654mm,Xband=4.72mm) in meters
	double L 									///< Total Linac electrical length in meters
	)
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

/** Allocates memory for a Linac struct (representing a Linac Section)
  * and fills it in with the values passed as arguments. Returns a pointer to the newly allocated struct. */
Linac *Linac_Allocate_New(
	Cryomodule_dp cryo_net,		///< Array of Cryomodules in a Linac
	int n_cryos,							///< Number of Cryomodules in Linac
	double dE,								///< Energy increase in Linac Section in eV
	double R56,								///< Longitudinal dispersion (if any)
	double T566,							///< 2nd-order longitudinal dispersion (if any)
	double phi,								///< Nominal Linac RF phase in radians
	double lam,								///< RF wavelength in meters
	double s0,								///< Wakefield characteristic length (Sband=0.105m, Xband=0.02625m) in meters
	double a,									///< Mean iris radius (Sband=11.654mm,Xband=4.72mm) in meters
	double L 									///< Total Linac electrical length in meters
	)
{
	Linac *linac = calloc(1,sizeof(Linac));
	Linac_Allocate_In(linac, cryo_net, n_cryos,
		dE, R56, T566, phi, lam, s0, a, L);
	return linac;
}

/** Frees memory of a Linac struct (representing a Linac Section). */
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

/** Helper routine to get a reference to a given Cryomodule given the Linac struct. */
Cryomodule *Get_Cryomodule(
	Linac *linac,		///< Pointer to Linac struct
	int index				///< Position of Cryomodule in Linac
	)
{
	if(linac->cryo_net[index]) return linac->cryo_net[index];
	else return NULL;
}
/** Helper routine to get a reference to a given Cryomodule State given the Linac State struct. */
Cryomodule_State *Get_Cryo_State(
	Linac_State *linac_state,	///< Pointer to Linac State struct
	int index									///< Position of Cryomodule in Linac
	)
{
	if(linac_state->cryo_state_net[index]) return linac_state->cryo_state_net[index];
	else return NULL;
}

/** Takes a previously configured Linac and allocates its State struct accordingly.
  * It allocates states for the Linac's Cryomodules recursively. */
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

/** Frees memory of Linac State struct. */
void Linac_State_Deallocate(Linac_State *linac_state, Linac *linac)
{
	for(int i=0;i<linac->n_cryos;i++){
		Cryomodule_State_Deallocate(linac_state->cryo_state_net[i], linac->cryo_net[i]);
	}
	free(linac_state->cryo_state_net);
}

/** Step function for Linac:
  * Calculates the state for the next simulation step.
  * Returns vector sum of all Cryomodule accelerating voltages
  * (as seen by the beam) and stores current state in State struct. */
double complex Linac_Step(
	Linac *linac,								///< Pointer to Linac struct
	Linac_State *linac_state,		///< Pointer to Linac State
	// Inputs
	double delta_tz,						///< Timing jitter in seconds (RF reference noise)
	double beam_charge,					///< Beam charge in Coulombs
	// Outputs
	double *amp_error,					///< Pointer to RF amplitude error (output, relative to overall Linac nominal voltage)
	double *phase_error					///< Pointer to RF phase error (output, in radians)
	)
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

	// Total Linac Accelerating voltage (store and return)
	linac_state->linac_V = linac_V;
	return linac_V;
}

/** Takes a pointer to a Gun struct which has been previously allocated
  * and fills it in with the values passed as arguments. */
void Gun_Allocate_In(
	Gun *gun,			///< Pointer to Gun struct
	double E,			///< Nominal gun exit energy in eV
	double sz0,		///< Nominal initial RMS bunch length in meters
	double sd0,		///< Nominal initial incoh. energy spread at nominal gun exit energy (fraction)
	double Q			///< Nominal beam charge in Coulombs
	)
{
	gun->E = E;
	gun->sz0 = sz0;
	gun->sd0 = sd0;
	gun->Q = Q;
}

/** Allocates memory for a Gun struct and fills it in with the values passed as arguments.
  * Returns a pointer to the newly allocated struct. */
Gun *Gun_Allocate_New(
	double E,			///< Nominal gun exit energy in eV
	double sz0,		///< Nominal initial RMS bunch length in meters
	double sd0,		///< Nominal initial incoh. energy spread at nominal gun exit energy (fraction)
	double Q			///< Nominal beam charge in Coulombsdouble E, double sz0, double sd0, double Q)
	)
{
	Gun *gun = calloc(1,sizeof(Gun));
	Gun_Allocate_In(gun, E, sz0, sd0, Q);
	return gun;
}
