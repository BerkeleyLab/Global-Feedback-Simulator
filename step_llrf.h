#ifndef STEP_LLRF_H
#define STEP_LLRF_H


/*
 *  step_llrf.h
 *
 * Routines for timestepping Low-Level RF model
 *
 * Alejandro F Queiruga
 * Daniel "Scott" Driver
 * 2013 LBL
 *
 */

/*
 * TODO: Fill in license
 */
#include <complex.h>
#include "linac_param.h"


typedef struct str_fpga_state {
  double complex drive, state, err;
} FPGA_State;

typedef struct str_cavity_state {
  double complex beam;
  double complex voltage;
} Cavity_State;

typedef struct str_linac_state {

  Filter_State RXF, TRF1, TRF2, Cav_Fil;

  FPGA_State fpga;
  Cavity_State cav;

  double complex RXF_out;
} Linac_State;

typedef Linac_State * LINS;


Linac_State** make_state_arrays(int n);
void cycle_buffer(Linac_State** linss, int n);

/*
 * Routine for allocating the buffers inside the Linac_State.
 * Mostly for mallocing the filter states, but it also zeroes 
 * other things out for good practice.
 */
void Linac_State_Allocate(Linac_State * lins, Linac_Param * linp);
void Linac_State_Deallocate(Linac_State * lins);

/*
 * step_llrf: Integrate the linac section forward in time
 *
 * linp is the paramter structure for this linac
 *
 * linss is a buffer of all the states,
 *  as a pointer to an array of pointers to Linac_States:
 *  linss[0] is the current timestep to be computed t=t;
 *  linss[1] is the previous timestep, t=t-dt used to integrate;
 *  linss[2...n] is further back in time, t=t-n*dt, used 
 *  for delays. The amount of history data depends on how long
 *  of delays are present, and should be kept at a minimum to
 *  preserve resources. The pointers are to be pingponged in a circle.
 */

double complex step_llrf(Linac_Param *linp,
			 double simt, double delta_tz,
			 double complex beam_charge, double complex adc_noise,
			 int openloop,
			 Linac_State ** linss);




/*
 * Helper routines
 */

/*
 * Exactly what it sounds like, apply a phase shift to a complex signal.
 * TODO: Can be optimized.
 */
double complex phase_shift(double complex in, double theta);

/*
 * Apply an FPGA's control forward in time. Input is the measured cavity voltage.
 * Output state to stnow.
 * returns sig_error which does not get propograted through the block diagram.
 */
double complex step_fpga(FPGA_Param * fpga, double complex cavity_vol,
			 FPGA_State * stnow, FPGA_State * stpast, 
			 int openloop);

double complex step_PI_fpga(FPGA_Param * fpga,
			    double dt, double complex cavity_vol,
			 FPGA_State * stnow, FPGA_State * stpast);

double complex saturate(double complex in, double harshness);

/*
 * Step a linac's triode configuration in time, 
 * drive_in -> [[ TRF1 -> saturate_c -> TRF2 ]] -> TRF2_OUTPUT_D
 *                  `------.  triode ,------'
 */
double complex step_triode(Linac_Param *linp, double complex drive_in,
		 Linac_State * linnow, Linac_State * linpast);


double complex step_cavity(Linac_Param *linp, double delta_tz,
			   double complex drive_in, double complex beam_charge,
			   Linac_State * linnow, Linac_State * linpast);

#endif
