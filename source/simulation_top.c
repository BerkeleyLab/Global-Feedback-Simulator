#include "simulation_top.h"

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

void Sim_State_Allocate(Simulation_State *sim_state, Simulation *sim, Noise_Srcs* noise_srcs)
{
	// Allocate memory for the array of Linac states
	sim_state->linac_state_net = calloc(sim->n_linacs, sizeof(Linac_State *));

	// Then allocate memory for the States themselves
 	for(int i=0;i<sim->n_linacs;i++) {
 		sim_state->linac_state_net[i] = (Linac_State*)calloc(1,sizeof(Linac_State));
 		Linac_State_Allocate(sim_state->linac_state_net[i], sim->linac_net[i]);
 	}
	// Allocate memory for the array of Linac amplitude and phase errors
	sim_state->amp_error_net = (double*)calloc(sim->n_linacs,sizeof(double));
	sim_state->phase_error_net = (double*)calloc(sim->n_linacs,sizeof(double));

	// Allocate memory for Longitudinal beam dynamics noise sources
	sim_state->noise_srcs = noise_srcs;

	// Allocate Doublecompress State
	sim_state->dc_state =(Doublecompress_State*)calloc(1,sizeof(Doublecompress_State));
	Doublecompress_State_Allocate(sim_state->dc_state, sim->n_linacs);
}

void Sim_State_Deallocate(Simulation_State *sim_state, Simulation *sim)
{
 	for(int i=0;i<sim->n_linacs;i++) {
 		Linac_State_Deallocate(sim_state->linac_state_net[i], sim->linac_net[i]);
	}
	free(sim_state->linac_state_net);
	free(sim_state->amp_error_net);
	free(sim_state->phase_error_net);
	free(sim_state->noise_srcs);
	Doublecompress_State_Deallocate(sim_state->dc_state);
	free(sim_state->dc_state);
}

void Apply_Correlated_Noise(int t_now, double Tstep, Noise_Srcs * noise_srcs)
{
	Noise_Step(t_now, Tstep, noise_srcs->type[0], noise_srcs->settings+N_NOISE_SETTINGS*0, &noise_srcs->dQ_Q);
	Noise_Step(t_now, Tstep, noise_srcs->type[1], noise_srcs->settings+N_NOISE_SETTINGS*1, &noise_srcs->dtg);
	Noise_Step(t_now, Tstep, noise_srcs->type[2], noise_srcs->settings+N_NOISE_SETTINGS*2, &noise_srcs->dE_ing);
	Noise_Step(t_now, Tstep, noise_srcs->type[3], noise_srcs->settings+N_NOISE_SETTINGS*3, &noise_srcs->dsig_z);
	Noise_Step(t_now, Tstep, noise_srcs->type[4], noise_srcs->settings+N_NOISE_SETTINGS*4, &noise_srcs->dsig_E);
	Noise_Step(t_now, Tstep, noise_srcs->type[5], noise_srcs->settings+N_NOISE_SETTINGS*5, &noise_srcs->dchirp);
}

#define CPRINT(c) {if(cimag(c)<0) fprintf(fp,"%10.16e%10.16ej ",creal(c),cimag(c)); else fprintf(fp,"%10.16e+%10.16ej ",creal(c),cimag(c)); }

void Write_Sim_Step(FILE * fp, double time, Simulation *sim, Simulation_State *sim_state)
{
 	fprintf(fp, "%10.16e   ",time);					// Current simulation time [s]
 	fprintf(fp, "%10.16e %10.16e %10.16e %10.16e %10.16e %10.16e  ",
		sim_state->noise_srcs->dQ_Q,				// Beam charge jitter [relative to nominal beam charge]
		sim_state->noise_srcs->dtg,					// Timing error of gun wrt RF (<0 is an early bunch) [s]
		sim_state->noise_srcs->dE_ing,			// Energy deviation at end of injector [eV]
		sim_state->noise_srcs->dsig_z,			// Deviation of bunch length from nominal length [m]
		sim_state->noise_srcs->dsig_E,			// Deviation of energy spread from nominal energy spread [fraction of nominal energy]
		sim_state->noise_srcs->dchirp);			// <Ez> correlation [m]

	for(int l=0;l<sim->n_linacs;l++) {
		fprintf(fp,"%10.16e %10.16e %10.16e %10.16e %10.16e %10.16e %10.16e %10.16e %10.16e %10.16e %10.16e  ",
			sim_state->amp_error_net[l],
			sim_state->phase_error_net[l],
			sim_state->dc_state->dE_E[l],
			sim_state->dc_state->dt[l],
			sim_state->dc_state->sz[l],
			sim_state->dc_state->dE_Ei[l],
			sim_state->dc_state->dE_Ei2[l],

			cabs(sim_state->linac_state_net[l]->linac_V),
			carg(sim_state->linac_state_net[l]->linac_V)*180.0/M_PI,

			cabs(sim_state->linac_state_net[l]->linac_Kg),
			carg(sim_state->linac_state_net[l]->linac_Kg)*180.0/M_PI
		);

	}
	// End simulation step entry with a new line
	fprintf(fp, "\n");
}

#undef CPRINT

/*
 * Performs sim.time_steps simulation time-steps (top-level of the entire Simulation Engine)
 */
void Simulation_Run(Simulation *sim, Simulation_State *sim_state, char * fname, int OUTPUTFREQ)
{
	// Total Linac accelerating voltage in Volts
	double complex linac_V=0.0;
	// Current simulation time [s]
	double time=0.0;
	// Temporary signals
	double delta_tz=0.0, beam_charge=0.0;

	// Open the file to store simulation results in
	FILE * fp = NULL;
	if(fname != NULL) fp = fopen(fname,"wb");

	// Iterate over time steps
	for(int t=0;t<sim->time_steps;t++){

		// Calculate current simulation time
		time = (t+1)*sim->Tstep;

		Apply_Correlated_Noise(t, sim->Tstep,sim_state->noise_srcs);

		// Iterate over Linacs
		for(int l=0;l<sim->n_linacs;l++){
			// Timing jitter
			delta_tz = *sim_state->dc_state->dt;
			// Add charge jitter to nominal beam charge to obtain instantaneous value
			beam_charge = sim->gun->Q*(1.0+sim_state->noise_srcs->dQ_Q);

			// Run Linac simulation step
			linac_V = Linac_Step(sim->linac_net[l], sim_state->linac_state_net[l],
			// Input errors
			delta_tz, beam_charge,
			// Output errors
			&sim_state->amp_error_net[l], &sim_state->phase_error_net[l]);

		} // End of iterate over Linacs

	// Calculate Longitudinal Beam Dynamics parameters
	Doublecompress(sim->gun, sim->linac_net, sim->n_linacs,
		// Input noise sources
		sim_state->noise_srcs,
		// Linac amplitude and phase errors
		sim_state->phase_error_net, sim_state->amp_error_net,
		// Doublecompress State
		sim_state->dc_state);

	if(t%OUTPUTFREQ == 0){
		if(fp!=NULL) Write_Sim_Step(fp, time, sim, sim_state);
    }

	// Apply Beam-based feedback
	// BBF_Step(sim->bbf, sim_state->dc_state, sim->linac_net, sim->n_linacs);

	} // End iteration over time-steps


	// Close the file to store simulation results in
	if(fp!=NULL) fclose(fp);
}

