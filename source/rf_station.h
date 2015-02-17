#ifndef RF_STATION_H
#define RF_STATION_H

/*
 *  rf_station.h
 *
 * Routines for timestepping RF Station model
 *
 */

#include <complex.h>
#include "filter.h"
#include "cavity.h"

/*
 * FPGA controller
 */

typedef struct str_FPGA {
  double kp;
  double ki;
  double complex set_point;
  double Tstep;
} FPGA;

typedef struct str_fpga_state {
  double complex drive, state, err;
} FPGA_State;

void FPGA_Allocate_In(FPGA * fpga, double kp, double ki, double complex set_point, double Tstep);
// double complex FPGA_Step(FPGA * fpga, double complex cavity_vol, FPGA_State * stnow, int openloop);
void FPGA_Clear(FPGA_State * stnow);
double complex step_PI_fpga(FPGA *fpga, double complex cavity_vol, FPGA_State *stnow, int openloop);

/*
 * Delay
 */
 
typedef struct str_Delay {
	int size;
} Delay;

typedef struct str_Delay_State {
	complex double *buffer;
	int index;
} Delay_State;

void Delay_State_Allocate(Delay *delay, Delay_State *delay_state);
void Delay_State_Deallocate(Delay *delay, Delay_State *delay_state);
double complex Delay_Step(double complex in, Delay *delay, Delay_State *delay_state);
void Delay_Clear(Delay *delay, Delay_State *delay_state);

/*
 * RF Station
 */

typedef struct str_RF_Station {

  double nom_grad;	// Nominal Cavity gradient
  double saturate_c;	// Saturation parameter

  Filter RXF;	// Anti-alias filter
  Filter TRF1, TRF2;	// Triode Filters
  FPGA fpga;	// FPGA Controller
  Cavity *cav;	// Cavity

  Delay loop_delay; // Loop Delay

} RF_Station;

typedef struct str_RF_Station_State {

  Filter_State RXF, TRF1, TRF2;
  Cavity_State cav_state;
  FPGA_State fpga_state;
  Delay_State loop_delay_state;

} RF_State;


typedef RF_State * RF_State_p;

RF_State** make_rf_state_array(int n);

void RF_Station_Allocate_In(RF_Station * rf_station,
   double Tstep, double saturate_c,  double kly_max_v,  
   double complex * p_TRF1, double complex * p_TRF2,  double complex * p_RXF,  
   Cavity *cav, double stable_gbw,  int loop_delay_size);

void RF_State_Allocate(RF_State *rf_state, RF_Station *rf_station);
void RF_State_Deallocate(RF_State *rf_state, RF_Station *rf_station);

double complex RF_Station_Step(
  RF_Station *rf_station,
  double delta_tz, double complex beam_charge,
  double complex probe_ns, double complex rev_ns,
  int openloop,
  RF_State *rf_state);

void RF_Station_Clear(RF_Station *rf_station, RF_State *rf_state);

/*
 * Helper routines
 */

/*
 * Exactly what it sounds like, apply a phase shift to a complex signal.
 * TODO: Can be optimized.
 */
double complex Phase_Shift(double complex in, double theta);

double complex Saturate(double complex in, double harshness);

/*
 * Step a linac's triode configuration in time, 
 * drive_in -> [[ TRF1 -> saturate_c -> TRF2 ]] -> TRF2_OUTPUT_D
 *                  `------.  triode ,------'
 */
double complex Triode_Step(RF_Station *rf_station, double complex drive_in,
		 RF_State *rf_state);

void Triode_Clear(RF_Station *rf_station, RF_State *rf_state);

#endif
