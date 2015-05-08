#include "cryomodule.h"

#include <math.h>

MechMode_State_dp MechMode_State_Allocate_Array(int n)
{
  MechMode_State_dp mechMode_State_net = calloc(n, sizeof(MechMode_State *));
  
  return mechMode_State_net;
}

MechMode_dp MechMode_Allocate_Array(int n)
{
  MechMode_dp mechMode_net = calloc(n, sizeof(MechMode *));
  
  return mechMode_net;
}

void MechMode_Append(MechMode** mechMode_arr, MechMode* mechMode, int index)
{
  // XXX Add some check!!
  mechMode_arr[index] = mechMode;
}

RF_State *Get_RF_State(Cryomodule_State *cryo_state, int index)
{
	if(cryo_state->rf_state_net[index]) return cryo_state->rf_state_net[index];
	else return NULL;
}

MechMode_State *Get_MechMode_State(Cryomodule_State *cryo_state, int index)
{
	if(cryo_state->mechMode_state_net[index]) return cryo_state->mechMode_state_net[index];
	else return NULL;
}

void MechMode_Allocate_In(MechMode *mechMode, double f_nu, double Q_nu, double k_nu, double Tstep) 
{
	// Pre-calculate angular frequency
	double omega_nu = 2.0*M_PI*f_nu;

	// Calculate Filter coefficients
	mechMode->a_nu = -omega_nu/(2.0*Q_nu);
	mechMode->b_nu = omega_nu*sqrt(1-1/(4.0*pow(Q_nu,2.0)));
	mechMode->c_nu = omega_nu/k_nu/mechMode->b_nu;
	mechMode->Tstep = Tstep;

	// Append modes to Mechanical Mode Filter
	// 2nd-order low-pass filter with two conjugate poles at a_nu +/- j*b_nu
	double complex mode_p1 = mechMode->a_nu + _Complex_I*mechMode->b_nu;
	double complex mode_p2 = mechMode->a_nu - _Complex_I*mechMode->b_nu;

	Filter_Append_Modes(&mechMode->fil, &mode_p1, 1, Tstep);
	Filter_Append_Modes(&mechMode->fil, &mode_p2, 1, Tstep);
}

MechMode *MechMode_Allocate_New(double f_nu, double Q_nu, double k_nu, double Tstep)
{
	// Allocate MechMode
	MechMode * mechMode = calloc(1,sizeof(MechMode));

	// Allocate Filter
	Filter_Allocate_In(&mechMode->fil,2,2);

	MechMode_Allocate_In(mechMode, f_nu, Q_nu, k_nu, Tstep);

	// Return allocated MechMode
	return mechMode;
}

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

void MechMode_State_Allocate(MechMode_State *mechMode_State, MechMode *mechMode)
{
	// Initialize attributes
	mechMode_State->x_nu = 0.0;

	// Allocate Filter State
	Filter_State_Allocate(&mechMode_State->fil_state, &mechMode->fil);

}

void MechMode_State_Deallocate(MechMode_State *mechMode_State)
{
	Filter_State_Deallocate(&mechMode_State->fil_state);
	mechMode_State->x_nu = 0.0;
}

void MechMode_Step(MechMode *mechMode, MechMode_State *mechMode_State, double complex F_nu)
{
	// Run the 2nd-order low-pass filter and store the result in x_nu (displacement)
 	mechMode_State->x_nu = Filter_Step(&(mechMode->fil), mechMode->c_nu*F_nu, &(mechMode_State->fil_state));
}

Cryomodule_dp Cryomodule_Allocate_Array(int n)
{
  Cryomodule_dp cryomodule_net = calloc(n, sizeof(Cryomodule *));
  return cryomodule_net;
}

void Cryomodule_Append(Cryomodule** cryo_arr, Cryomodule* cryo, int index)
{
  // XXX Add some check!!
  cryo_arr[index] = cryo;
}

void Cryomodule_Allocate_In(Cryomodule *cryo, RF_Station_dp rf_station_net, int n_rf_stations, MechMode_dp mechMode_net, int n_mechModes)
{
	cryo -> rf_station_net = rf_station_net;
	cryo -> mechMode_net = mechMode_net;
	cryo -> n_rf_stations = n_rf_stations;
	cryo -> n_mechModes = n_mechModes;
}

Cryomodule * Cryomodule_Allocate_New(RF_Station **rf_station_net, int n_rf_stations, MechMode **mechMode_net, int n_mechModes)
{
	Cryomodule * cryo = calloc(1,sizeof(Cryomodule));
	Cryomodule_Allocate_In(cryo, rf_station_net, n_rf_stations, mechMode_net, n_mechModes);
	return cryo;
}
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

RF_Station *Get_RF_Station(Cryomodule *cryo, int index)
{
	if(cryo->rf_station_net[index]) return cryo->rf_station_net[index];
	else return NULL;
}

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

void Cryomodule_Step(Cryomodule *cryo, Cryomodule_State * cryo_state,
	double delta_tz, double complex beam_charge, double complex probe_ns, double complex rev_ns,
	int openloop)
{
	// Electrical state-space model
	// ----------------------------------------------------

	int i;
	double complex RF_out;
	// Iterate over RF Stations in a Cryomodule
	for(i=0;i<cryo->n_rf_stations;i++) {
		 RF_out = RF_Station_Step(cryo->rf_station_net[i], delta_tz, beam_charge, probe_ns, rev_ns, openloop, cryo_state->rf_state_net[i]);
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

}