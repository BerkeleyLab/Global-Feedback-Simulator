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
  double out_sat;
  double state_sat;
  double Tstep;
} FPGA;

typedef struct str_fpga_state {
  double complex drive, state, err;
  // FPGA control boolean to open and close the feedback loop
  int openloop;
} FPGA_State;

void FPGA_Allocate_In(FPGA * fpga, double kp, double ki, double complex set_point, double out_sat, double Tstep);
void FPGA_Deallocate(FPGA *fpga);
void FPGA_Clear(FPGA_State * stnow);
double complex FPGA_Step(FPGA *fpga, double complex cavity_vol, FPGA_State *stnow);

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

void Delay_State_Deallocate(Delay_State *delay_state);
void Delay_State_Allocate(Delay *delay, Delay_State *delay_state);
void Delay_State_Deallocate(Delay_State *delay_state);
double complex Delay_Step(double complex in, Delay *delay, Delay_State *delay_state);
void Delay_Clear(Delay *delay, Delay_State *delay_state);
void Delay_Deallocate(Delay *delay);

/*
 * RF Station
 */

typedef struct str_RF_Station {

  double nom_grad;	// Nominal Cavity gradient
  double Clip;  // Saturation parameter
  double PAscale;	// Amplifier scaling (from unitless to sqrt(W))
  double PAmax;

  Filter RXF;	// Anti-alias filter
  Filter TRF1, TRF2;	// SSA Filters
  FPGA fpga;	// FPGA Controller
  Cavity *cav;	// Cavity

  Delay loop_delay; // Loop Delay

} RF_Station;

typedef RF_Station* RF_Station_p;
typedef RF_Station** RF_Station_dp;

RF_Station_dp RF_Station_Allocate_Array(int n);
void RF_Station_Append(RF_Station** rf_station_arr, RF_Station* rf_station, int index);

typedef struct str_RF_Station_State {

  Filter_State RXF, TRF1, TRF2;
  Cavity_State cav_state;
  FPGA_State fpga_state;
  Delay_State loop_delay_state;

} RF_State;

typedef RF_State* RF_State_p;
typedef RF_State** RF_State_dp;

RF_State_dp RF_State_Allocate_Array(int n);

void RF_Station_Allocate_In(RF_Station * rf_station,
  double Tstep,
  double Clip,
  double PAmax,
  double PAscale,
  double complex * p_TRF1, double complex * p_TRF2,
  double complex * p_RXF,
  Cavity *cav,
  double stable_gbw,
  double FPGA_out_sat,
  int loop_delay_size);

RF_Station * RF_Station_Allocate_New(
  double Tstep,
  double Clip,
  double PAmax,
  double PAscale,
  double complex * p_TRF1, double complex * p_TRF2,
  double complex * p_RXF,
  Cavity *cav,
  double stable_gbw,
  double FPGA_out_sat,
  int loop_delay_size);

void RF_Station_Deallocate(RF_Station *rf_station);

void RF_State_Allocate(RF_State *rf_state, RF_Station *rf_station);
void RF_State_Deallocate(RF_State *rf_state, RF_Station *rf_station);

double complex RF_Station_Step(
  RF_Station *rf_station,
  double delta_tz, double complex beam_charge,
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
 * Step a linac's SSA configuration in time, 
 * drive_in -> [[ TRF1 -> saturate_c -> TRF2 ]] -> TRF2_OUTPUT_D
 *                  `------.  triode ,------'
 */
double complex SSA_Step(RF_Station *rf_station, double complex drive_in,
		 RF_State *rf_state);

void SSA_Clear(RF_Station *rf_station, RF_State *rf_state);

#endif
