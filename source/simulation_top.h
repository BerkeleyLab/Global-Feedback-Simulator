/**
  * @file simulation_top.h
  * @brief Header file for simulation_top.c
  * The Top Level program for the beam and accelerator simulation.
  * @author Carlos Serrano (CSerrano@lbl.gov)
*/

#ifndef SIMULATION_TOP_H
#define SIMULATION_TOP_H

#include <complex.h>
#include <stdio.h>

#include "linac.h"
#include "doublecompress.h"
#include "noise.h"
// #include "beam_based_feedback.h"



/**
 * Data structure storing the parameters for a full Accelerator Simulation
 * comprised of a Gun and multiple Linacs.
 */
typedef struct str_Simulation {

	// Simulation parameters
	double Tstep;	///< Simulation time-step size
	int time_steps;	///< Total number of Simulation steps

	// Electron Gun
	Gun *gun;

	// Array of Linacs in an Accelerator
	int n_linacs;
	Linac **linac_net;

	// Beam-based feedback parameters
	// BBF *bbf;


} Simulation;

typedef struct str_Simulation_State {

	Linac_State **linac_state_net;	///< Array of Linac States
	Noise_Srcs *noise_srcs;///< Longitudinal beam dynamics noise sources

	Doublecompress_State *dc_state; ///< Doublecompress State

	// Array of Linac accelerating voltage errors (amplitude and phase)
	// (amplitude normalized by Linac increase in Energy in eV)
	double *amp_error_net, *phase_error_net;

} Simulation_State;

void Sim_Allocate_In(Simulation *sim, double Tstep, int time_steps,
	Gun *gun, Linac ** linac_net, int n_linacs);
Simulation *Sim_Allocate_New(double Tstep, int time_steps, Gun *gun, Linac ** linac_net, int n_linacs);
void Sim_Deallocate(Simulation *sim);

void Sim_State_Allocate(Simulation_State *sim_state, Simulation *sim, Noise_Srcs* noise_srcs);
void Sim_State_Deallocate(Simulation_State *sim_state, Simulation *sim);

void Apply_Correlated_Noise(int t_now, double Tstep, Noise_Srcs * noise_srcs);
void Write_Sim_Step( FILE * fp, double time, Simulation *sim, Simulation_State *sim_state);

/**
 * Performs sim.time_steps simulation time-steps (top-level of the entire Simulation Engine)
 */
void Simulation_Run(Simulation *sim, Simulation_State *sim_state,char * fname, int OUTPUTFREQ);

#endif
