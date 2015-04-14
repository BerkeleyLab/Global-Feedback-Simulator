#include "simulation_top.h"

#include <stdio.h>
#include <stdlib.h>

void Sim_Allocate_In(Simulation *sim, double Tstep, int time_steps,
	Gun *gun, Linac ** linac_net, int n_linacs)
{
	sim->Tstep = Tstep;
	sim->time_steps = time_steps;
	sim->gun = gun;
	sim->linac_net = linac_net;
	sim->n_linacs = n_linacs;
}

Simulation *Sim_Allocate_New(double Tstep, double time_steps,
	Gun *gun, Linac ** linac_net, int n_linacs)
{
	Simulation *sim = calloc(1,sizeof(Simulation));
	Sim_Allocate_In(sim, Tstep, time_steps, gun, linac_net, n_linacs);

	return sim;
}
void Sim_Deallocate(Simulation *sim)
{
	for(int i=0;i<sim->n_linacs;i++){
		Linac_Deallocate(sim->linac_net[i]);
	}
	free(sim->linac_net);
	free(sim->gun);
	
	sim->Tstep = 0.0;
	sim->time_steps = 0;
	sim->time_steps = 0;
	sim->n_linacs = 0;
}

void Sim_State_Allocate(Simulation_State *sim_state, Simulation *sim)
{
	// Allocate memory for the array of states
	sim_state->linac_state_net = calloc(sim->n_linacs, sizeof(Linac_State *));

	// Then allocate memory for the States themselves
 	for(int i=0;i<sim->n_linacs;i++) {
 		sim_state->linac_state_net[i] = (Linac_State*)calloc(1,sizeof(Linac_State));
 		Linac_State_Allocate(sim_state->linac_state_net[i], sim->linac_net[i]);
 	}
}

void Sim_State_Dellocate(Simulation_State *sim_state, Simulation *sim)
{
 	for(int i=0;i<sim->n_linacs;i++) {
 		Linac_State_Deallocate(sim_state->linac_state_net[i], sim->linac_net[i]);
	}
	free(sim_state->linac_state_net);
}

/*
 * Performs sim.time_steps simulation time-steps (top-level of the entire Simulation Engine)
 */
// void Simulation_Run(Simulation *sim, Simulation_State *sim_state,
	// char * fname, int OUTPUTFREQ);