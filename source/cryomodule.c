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
	mechMode->c_nu = -omega_nu/k_nu/mechMode->b_nu;

	// Append modes to Filters for x_nu and y_nu
	// Mode's single-pole, low-pass filter allocation
  	// Calculate poles: Same pole for both Filters (a_nu)
  	double complex mode_p = mechMode->a_nu;
	Filter_Append_Modes(&mechMode->fil_x, &mode_p, 1, Tstep);
	Filter_Append_Modes(&mechMode->fil_y, &mode_p, 1, Tstep);
}

MechMode *MechMode_Allocate_New(double f_nu, double Q_nu, double k_nu, double Tstep)
{
	// Allocate MechMode
	MechMode * mechMode = calloc(1,sizeof(MechMode));

	// Allocate Filters for x_nu and y_nu
	Filter_Allocate_In(&mechMode->fil_x,1,1);
	Filter_Allocate_In(&mechMode->fil_y,1,1);

	MechMode_Allocate_In(mechMode, f_nu, Q_nu, k_nu, Tstep);

	// Return allocated MechMode
	return mechMode;
}

void MechMode_State_Allocate(MechMode_State *mechMode_State, MechMode *mechMode)
{
	// Initialize attributes
	mechMode_State->x_nu = 0.0;
	mechMode_State->y_nu = 0.0;

	// Allocate Filter State for x_nu
	Filter_State_Allocate(&mechMode_State->fil_state_x, &mechMode->fil_x);
	
	// Allocate Filter State for y_nu
	Filter_State_Allocate(&mechMode_State->fil_state_y, &mechMode->fil_y);

}

void MechMode_Step(MechMode *mechMode, MechMode_State *mechMode_State, double complex F_nu)
{
	// Calculate drive term for x_nu filter
	double complex x_in = -mechMode->b_nu*mechMode_State->y_nu;
	// Step x_nu filter (displacement coordinate)
 	mechMode_State->x_nu = Filter_Step(&(mechMode->fil_x), x_in, &(mechMode_State->fil_state_x));

	// Calculate drive term for y_nu filter
 	double complex y_in = mechMode->b_nu*mechMode_State->x_nu + mechMode->c_nu*F_nu;
	// Step y_nu filter (velocity coordinate)
 	mechMode_State->y_nu = Filter_Step(&(mechMode->fil_y), y_in, &(mechMode_State->fil_state_y));

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

	// Indeces (as used in equations in the documentation)
		// nu: Mechanical Eigenmodes index,
		// mu: Electrical Eigenmode index,
	int nu, mu;
	int n_Emodes;
	// Temporary variables for summation components
	double V_2_mu=0.0, A_nu_mu=0.0, x_nu=0.0, C_mu_nu=0.0;

	// Calculate Electrical to Mecahnical couplings
	// and calculate displacements for each Mechanical Eigenmode
	// Iterate over nu
	for(nu=0;nu<cryo->n_mechModes;nu++) {
		// Iterate over RF Stations in a Cryomodule
		for(i=0;i<cryo->n_rf_stations;i++) {
			// Iterate over Electrical Eignenmodes (mu index) in a RF Station's Cavity
			n_Emodes = cryo->rf_station_net[i]->cav->n_modes; 
			for(mu=0;mu<n_Emodes;mu++) {
				// Grab square of accelerating voltage of Electrical Mode (mu)
				V_2_mu = cryo_state->rf_state_net[i]->cav_state.elecMode_state_net[mu]->V_2;
				// Grab its coupling to Mechanical Mode (nu,mu)
				A_nu_mu = cryo->rf_station_net[i]->cav->elecMode_net[mu]->A[nu];
				
				// Calculate contribution of Electrical Eigenmode to Lorentz force
				cryo_state->F_nu[nu] += A_nu_mu*V_2_mu;
			} // End iterate over mu
		} // End iterate over i

		// Once F_nu is known for a Mechanical Eigenmode,
		// apply state-space simulation step to Mechanical Mode (2nd-order LPF)
		MechMode_Step(cryo->mechMode_net[nu], cryo_state->mechMode_state_net[nu], cryo_state->F_nu[nu]);

	} // End iterate over nu

	// Displacements are now available for each Mechanical Mode
	// Apply matrix to transalate into Lorentz-force detuning

	// Iterate over RF Stations in a Cryomodule
	for(i=0;i<cryo->n_rf_stations;i++) {
		// Iterate over Electrical Eignenmodes (mu index) in a RF Station's Cavity
		n_Emodes = cryo->rf_station_net[i]->cav->n_modes; 
		for(mu=0;mu<n_Emodes;mu++) {
			// Iterate over nu
			for(nu=0;nu<cryo->n_mechModes;nu++) {
				// Grab displacement for Mechanical Eigenmode (mu)
				x_nu = cryo_state->mechMode_state_net[nu]->x_nu;
				// Grab Mechanical to Electrical coupling coefficient (mu,nu)
				C_mu_nu = cryo->rf_station_net[i]->cav->elecMode_net[mu]->C[nu];
				// Calculate contribution of mechanical displacement to Lorentz-force detuning 
				cryo_state->rf_state_net[i]->cav_state.elecMode_state_net[mu]->delta_omega += C_mu_nu*V_2_mu;
			} // End iterate over nu
		} // End iterate over mu
	} // End iterate over i

}