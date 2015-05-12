#ifndef SIMULATION_TOP_H
#define SIMULATION_TOP_H

/*
 * simulation_top.h/c
 *
 * The Top Level program for the beam and accelerator simulation.
 *
 */

#include <complex.h>

#include "linac.h"
#include "doublecompress.h"
// #include "dynamic_noise.h"
// #include "beam_based_feedback.h"

/*
 * Data struture storing the parameters for a full Accelerator Simulation
 * comprised of a Gun and multiple Linacs.
 */
typedef struct str_Simulation {

	// Simulation parameters
	double Tstep;	// Simulation time-step size
	int time_steps;	// Total number of Simulation steps
	
	// Electron Gun
	Gun *gun;

	// Array of Linacs in an Accelerator
	int n_linacs;
	Linac **linac_net;

	// Beam-based feedback parameters
	// BBF *bbf;

	// Noise parameters
	Noise_Srcs *noise_src;

} Simulation;

typedef struct str_Simulation_State {
	// Array of Linac States
	Linac_State **linac_state_net;

	// Longitudinal beam dynamics noise sources
	Noise_Srcs *noise_srcs;

	// Doublecompress State
	Doublecompress_State *dc_state;

	// Array of Linac accelerating voltage errors (amplitude and phase)
	// (amplitude normalized by Linac increase in Energy in eV)
	double *amp_error_net, *phase_error_net; 

} Simulation_State;

void Sim_Allocate_In(Simulation *sim, double Tstep, int time_steps,
	Gun *gun, Linac ** linac_net, int n_linacs);
Simulation *Sim_Allocate_New(double Tstep, double time_steps,
	Gun *gun, Linac ** linac_net, int n_linacs);
void Sim_Deallocate(Simulation *sim);

void Sim_State_Allocate(Simulation_State *sim_state, Simulation *sim);
void Sim_State_Deallocate(Simulation_State *sim_state, Simulation *sim);

/*
 * Performs sim.time_steps simulation time-steps (top-level of the entire Simulation Engine)
 */
void Simulation_Run(Simulation *sim, Simulation_State *sim_state,char * fname, int OUTPUTFREQ);

#endif
