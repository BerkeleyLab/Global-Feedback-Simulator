/**
 * @file cryomodule.c
 * @brief Cryomodule Model: Allocation, configuration and stepping functions for the cryomodule model,
 * which iterates over an arbitrary number of cavities, mechanical modes and electro-mechanical interactions.
 * @author Carlos Serrano (Cserrano@lbl.gov)
 */

#include "cryomodule.h"

#include <math.h>

/** Allocates memory for an array of Mechanical Eigenmodes.
  * Mechanical modes themselves need to be allocated and filled individually,
  * and appended to this array using MechMode_Append.*/
MechMode_dp MechMode_Allocate_Array(int n 	///< Size of the array of Mechanical Eigenmodes
	)
{
  MechMode_dp mechMode_net = calloc(n, sizeof(MechMode *));

  return mechMode_net;
}

/** Append a Mechanical Eigenmode previously allocated and filled to the array given as an argument. */
void MechMode_Append(MechMode** mechMode_arr,	///< Pointer to array of Mechanical Eigenmodes
	MechMode* mechMode,													///< Pointer to Mechanical Eigenmode to be appended
	int index 																	///< Position of Mechanical Eigenmode in the array
	)
{
  // XXX Add some check!!
  mechMode_arr[index] = mechMode;
}

/** Helper routine to get a reference to a given RF Station State given the Cryomodule State struct. */
RF_State *Get_RF_State(Cryomodule_State *cryo_state,	///< Pointer to Cryomodule State
	int index 																					///< Position of RF Station in the Cryomodule
	)
{
	if(cryo_state->rf_state_net[index]) return cryo_state->rf_state_net[index];
	else return NULL;
}

/** Helper routine to get a reference to a given Mechanical mode State given the Cryomodule State struct. */
MechMode_State *Get_MechMode_State(Cryomodule_State *cryo_state,	///< Pointer to Cryomodule State
	int index 																											///< Position of Mechanical mode in the Cryomodule
	)
{
	if(cryo_state->mechMode_state_net[index]) return cryo_state->mechMode_state_net[index];
	else return NULL;
}

/** Takes a pointer to a MechMode struct (representing a mechanical eigenmode)
  * which has been previously allocated and fills it in with the values passed as arguments. */
void MechMode_Allocate_In(
	MechMode *mechMode,		///< Pointer to MechMode struct
	double f_nu,					///< Mechanical mode's resonant frequency in Hz
	double Q_nu,					///< Mechanical mode's Quality factor (unitless)
	double k_nu,					///< Mechanical mode's spring constant in N/m
	double Tstep					///< Simulation time-step in seconds
	)
{
	// Pre-calculate angular frequency
	double omega_nu = 2.0*M_PI*f_nu;

	// Calculate Filter coefficients
	mechMode->a_nu = -omega_nu/(2.0*Q_nu);
	mechMode->b_nu = omega_nu*sqrt(1-1/(4.0*pow(Q_nu,2.0)));
	// c_nu is in the physics equations in order to normalize the filter gain at DC
	// This normalization already takes place in the Filter routine (see filter.h),
	// and therefore c_nu is included here to match the equations but is set to unity.
	mechMode->c_nu = 1.0/k_nu;
	mechMode->Tstep = Tstep;

	// Append modes to Mechanical Mode Filter
	// 2nd-order low-pass filter with two conjugate poles at a_nu +/- j*b_nu
	double complex mode_p1 = mechMode->a_nu + _Complex_I*mechMode->b_nu;
	double complex mode_p2 = mechMode->a_nu - _Complex_I*mechMode->b_nu;

	Filter_Append_Modes(&mechMode->fil, &mode_p1, 1, Tstep);
	Filter_Append_Modes(&mechMode->fil, &mode_p2, 1, Tstep);
}

/** Allocates memory for a MechMode struct (representing a mechanical eigenmode)
  * and fills it in with the values passed as arguments. Returns a pointer to the newly allocated struct. */
MechMode *MechMode_Allocate_New(
	double f_nu,					///< Mechanical mode's resonant frequency in Hz
	double Q_nu,					///< Mechanical mode's Quality factor (unitless)
	double k_nu,					///< Mechanical mode's spring constant in N/m
	double Tstep					///< Simulation time-step in seconds
	)
{
	// Allocate MechMode
	MechMode * mechMode = calloc(1,sizeof(MechMode));

	// Allocate Filter
	Filter_Allocate_In(&mechMode->fil,2,2);

	MechMode_Allocate_In(mechMode, f_nu, Q_nu, k_nu, Tstep);

	// Return allocated MechMode
	return mechMode;
}

/** Frees memory of a MechMode mode struct (representing an mechanical eigenmode). */
void MechMode_Deallocate(MechMode *mechMode)
{
	Filter_Deallocate(&mechMode->fil);
	mechMode->f0 = 0.0;
	mechMode->Q = 0.0;
	mechMode->a_nu = 0.0;
	mechMode->b_nu = 0.0;
	mechMode->c_nu = 0.0;
	mechMode->Tstep = 0.0;
}

/** Takes a previously configured Mechanical mode and allocates its State struct accordingly. */
void MechMode_State_Allocate(MechMode_State *mechMode_State, MechMode *mechMode)
{
	// Initialize attributes
	mechMode_State->x_nu = 0.0;

	// Allocate Filter State
	Filter_State_Allocate(&mechMode_State->fil_state, &mechMode->fil);

}

/** Frees memory of Mechanical Eigenmode State struct. */
void MechMode_State_Deallocate(MechMode_State *mechMode_State)
{
	Filter_State_Deallocate(&mechMode_State->fil_state);
	mechMode_State->x_nu = 0.0;
}

/** Step function for Mechanical Eigenmode:
  * Calculates the state for the next simulation step, which comes down to stepping the 2nd-order LPF
  * It stores the current state in State struct (displacement x_nu and filter state). */
void MechMode_Step(MechMode *mechMode,	///< Pointer to MechMode struct
	MechMode_State *mechMode_State,				///< Pointer to MechMode State
	double complex F_nu 									///< Lorentz force in Newtons
	)
{
	// Run the 2nd-order low-pass filter and store the result in x_nu (displacement)
 	mechMode_State->x_nu = Filter_Step(&(mechMode->fil), mechMode->c_nu*F_nu, &(mechMode_State->fil_state));
}

/** Allocates memory for an array of Cryomodules.
  * Cryomodules themselves need to be allocated and filled individually,
  * and appended to this array using Cryomodule_Append.*/
Cryomodule_dp Cryomodule_Allocate_Array(int n 	///< Size of the array of Cryomodules
	)
{
  Cryomodule_dp cryomodule_net = calloc(n, sizeof(Cryomodule *));
  return cryomodule_net;
}

/** Append a Cryomodule previously allocated and filled to the array given as an argument. */
void Cryomodule_Append(Cryomodule** cryo_arr, Cryomodule* cryo, int index)
{
  // XXX Add some check!!
  cryo_arr[index] = cryo;
}

/** Takes a pointer to a Cryomodule struct which has been previously allocated and fills it in with the values passed as arguments. */
void Cryomodule_Allocate_In(Cryomodule *cryo,		///< Pointer to Cryomodule struct
	RF_Station_dp rf_station_net,									///< Array of RF Stations in a Cryomodule
	int n_rf_stations,														///< Number of RF Stations
	MechMode_dp mechMode_net,											///< Array of Mechanical modes in a Cryomodule
	int n_mechModes																///< Number of Mechanical modes
	)
{
	cryo -> rf_station_net = rf_station_net;
	cryo -> mechMode_net = mechMode_net;
	cryo -> n_rf_stations = n_rf_stations;
	cryo -> n_mechModes = n_mechModes;
}

/** Allocates memory for a Cryomodule struct and fills it in with the values passed as arguments.
	* Returns a pointer to the newly allocated struct. */
Cryomodule * Cryomodule_Allocate_New(
	RF_Station_dp rf_station_net,									///< Array of RF Stations in a Cryomodule
	int n_rf_stations,														///< Number of RF Stations
	MechMode_dp mechMode_net,											///< Array of Mechanical modes in a Cryomodule
	int n_mechModes																///< Number of Mechanical modesRF_Station **rf_station_net, int n_rf_stations, MechMode **mechMode_net, int n_mechModes)
	)
{
	Cryomodule * cryo = calloc(1,sizeof(Cryomodule));
	Cryomodule_Allocate_In(cryo, rf_station_net, n_rf_stations, mechMode_net, n_mechModes);
	return cryo;
}

/** Frees memory of a Cryomodule struct.*/
void Cryomodule_Deallocate(Cryomodule* cryo)
{
	for(int i=0;i<cryo->n_rf_stations;i++){
		RF_Station_Deallocate(cryo->rf_station_net[i]);
	}
	for (int i=0;i<cryo->n_mechModes;i++){
		MechMode_Deallocate(cryo->mechMode_net[i]);
	}

	free(cryo->rf_station_net);
	free(cryo->mechMode_net);

	cryo->n_rf_stations = 0.0;
	cryo->n_mechModes = 0.0;
}

/** Helper routine to get a reference to a given RF Station given the Cryomodule struct. */
RF_Station *Get_RF_Station(Cryomodule *cryo, int index)
{
	if(cryo->rf_station_net[index]) return cryo->rf_station_net[index];
	else return NULL;
}

/** Takes a previously configured Cryomodule and allocates its State struct accordingly.
  * It allocates states for the Cryomodule's RF Stations and Mechanical modes recursively. */
void Cryomodule_State_Allocate(Cryomodule_State *cryo_state, Cryomodule *cryo)
{

	// Allocate memory for the array of states (RF Stations and MechModes)
	cryo_state->rf_state_net = calloc(cryo->n_rf_stations, sizeof(RF_State *));
	cryo_state->mechMode_state_net = calloc(cryo->n_mechModes, sizeof(MechMode_State *));

	// Then allocate memory for the States themselves
	// Allocate RF Station States
	int i;
 	for(i=0;i<cryo->n_rf_stations;i++) {
 		cryo_state->rf_state_net[i] = (RF_State*)calloc(1,sizeof(RF_State));
 		RF_State_Allocate(cryo_state->rf_state_net[i], cryo->rf_station_net[i]);
 	}

	// Allocate Mechanical Mode States
 	for(i=0;i<cryo->n_mechModes;i++) {
 		cryo_state->mechMode_state_net[i] = (MechMode_State*)calloc(1,sizeof(MechMode_State));
 		MechMode_State_Allocate(cryo_state->mechMode_state_net[i], cryo->mechMode_net[i]);
 	}

 	// Allocate vectors to store Lorentz forces
 	cryo_state->F_nu = calloc(cryo->n_mechModes,sizeof(double));

}

/** Frees memory of Cryomodule State struct. */
void Cryomodule_State_Deallocate(Cryomodule_State *cryo_state, Cryomodule *cryo)
{
	for(int i=0;i<cryo->n_rf_stations;i++){
		RF_State_Deallocate(cryo_state->rf_state_net[i], cryo->rf_station_net[i]);
	}
	for (int i=0;i<cryo->n_mechModes;i++){
		MechMode_State_Deallocate(cryo_state->mechMode_state_net[i]);
	}

	free(cryo_state->rf_state_net);
	free(cryo_state->mechMode_state_net);
	free(cryo_state->F_nu);
}

/** Step function for Cryomodule:
  * Calculates the state for the next simulation step.
  * Returns vector sum of all cavity accelerating voltages
  * (as seen by the beam) and stores current state in State struct. */
double complex Cryomodule_Step(
	Cryomodule *cryo,								///< Pointer to Cryomodule
	Cryomodule_State * cryo_state,	///< Pointer to Cryomodule State
	double delta_tz,								///< Timing jitter in seconds (RF reference noise)
	double beam_charge							///< Beam charge in Coulombs
	)
{
	// Electrical state-space model
	// ----------------------------------------------------

	int i;
	// Total Cryomodule accelerating voltage
	double complex cryo_V=0.0;
	// Total Cryomodule drive signal
	double complex cryo_Kg=0.0;

	// Iterate over RF Stations in a Cryomodule
	for(i=0;i<cryo->n_rf_stations;i++) {
		 cryo_V += RF_Station_Step(cryo->rf_station_net[i], delta_tz, beam_charge, 0.0, cryo_state->rf_state_net[i]);
		 cryo_Kg += cryo_state->rf_state_net[i]->cav_state.Kg;
	} // End iterate over i

	// Electro-mechanical state-space model
	// ----------------------------------------------------

	// Indexes (as used in equations in the documentation)
		// nu: Mechanical Eigenmodes index,
		// mu: Electrical Eigenmode index,
	int nu, mu;
	int n_Emodes;
	// Temporary variables for summation components
	double V_2_mu=0.0, A_nu_mu=0.0, x_nu=0.0, C_mu_nu=0.0, delta_omega_now=0.0, F_nu_now=0.0;

	// Calculate Electrical to Mechanical couplings
	// and calculate displacements for each Mechanical Eigenmode
	// Iterate over nu
	for(nu=0;nu<cryo->n_mechModes;nu++) {
		// Iterate over RF Stations in a Cryomodule
		for(i=0;i<cryo->n_rf_stations;i++) {
			// Iterate over Electrical Eigenmodes (mu index) in a RF Station's Cavity
			n_Emodes = cryo->rf_station_net[i]->cav->n_modes;
			for(mu=0;mu<n_Emodes;mu++) {
				// Grab square of accelerating voltage of Electrical Mode (mu)
				V_2_mu = cryo_state->rf_state_net[i]->cav_state.elecMode_state_net[mu]->V_2;
				// Grab its coupling to Mechanical Mode (nu,mu)
				A_nu_mu = cryo->rf_station_net[i]->cav->elecMode_net[mu]->A[nu];
				// Calculate contribution of Electrical Eigenmode to Lorentz force
				F_nu_now += A_nu_mu*V_2_mu;
			} // End iterate over mu
		} // End iterate over i
		cryo_state->F_nu[nu] = F_nu_now;

		// Once F_nu is known for a Mechanical Eigenmode,
		// apply state-space simulation step to Mechanical Mode (2nd-order LPF)
		MechMode_Step(cryo->mechMode_net[nu], cryo_state->mechMode_state_net[nu], cryo_state->F_nu[nu]);

	} // End iterate over nu

	// Displacements are now available for each Mechanical Mode
	// Apply matrix to translate into Lorentz-force detuning

	// Iterate over RF Stations in a Cryomodule
	for(i=0;i<cryo->n_rf_stations;i++) {
		// Iterate over Electrical Eigenmodes (mu index) in a RF Station's Cavity
		n_Emodes = cryo->rf_station_net[i]->cav->n_modes;
		for(mu=0;mu<n_Emodes;mu++) {
			// Iterate over nu
			for(nu=0;nu<cryo->n_mechModes;nu++) {
				// Grab displacement for Mechanical Eigenmode (nu)
				x_nu = cryo_state->mechMode_state_net[nu]->x_nu;
				// Grab Mechanical to Electrical coupling coefficient (mu,nu)
				C_mu_nu = cryo->rf_station_net[i]->cav->elecMode_net[mu]->C[nu];
				// Calculate contribution of mechanical displacement to Lorentz-force detune frequency
				delta_omega_now += C_mu_nu*x_nu;
			} // End iterate over nu
			cryo_state->rf_state_net[i]->cav_state.elecMode_state_net[mu]->delta_omega = delta_omega_now;
		} // End iterate over mu
	} // End iterate over i

	// Store total Cryomodule drive signal (vector sum of all RF Station drive signals)
	 cryo_state->cryo_Kg = cryo_Kg;

	// Return vector sum of all cavity accelerating voltages
	// (as seen by the beam)
	return cryo_V;

}
