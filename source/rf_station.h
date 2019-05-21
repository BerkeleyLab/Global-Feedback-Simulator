/**
  * @file rf_station.h
  * @brief Header file for rf_station.c
  * @author Carlos Serrano (CSerrano@lbl.gov)
*/

#ifndef RF_STATION_H
#define RF_STATION_H

#include <complex.h>
#include "filter.h"
#include "cavity.h"
#include "noise.h"

/*
 * FPGA
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
  int drive_off; ///< FPGA control boolean to open and close the feedback loop
} FPGA_State;

void FPGA_Allocate_In(FPGA * fpga, double kp, double ki, double complex set_point, double out_sat, double Tstep);
void FPGA_Deallocate(FPGA *fpga);
void FPGA_Clear(FPGA_State * stnow);
double complex FPGA_Step(FPGA *fpga, double complex cavity_vol, double complex feed_forward, FPGA_State *stnow);

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

  double nom_grad;  ///< Nominal Cavity gradient
  double Clip;  ///< Saturation parameter
  double PAscale; ///< Amplifier scaling (from unitless to sqrt(W))
  double PAmax;

  Filter noise_shape_fil; ///< Noise-shaping low-pass filter
  Filter SSA_fil; ///< SSA Filters
  FPGA fpga;  ///< FPGA Controller
  Cavity *cav;  ///< Cavity

  Delay loop_delay; ///< Loop Delay

  // RMS Noise setting for each ADC
  double probe_ns_rms, rev_ns_rms, fwd_ns_rms;

} RF_Station;

typedef RF_Station* RF_Station_p;
typedef RF_Station** RF_Station_dp;

RF_Station_dp RF_Station_Allocate_Array(int n);
void RF_Station_Append(RF_Station** rf_station_arr, RF_Station* rf_station, int index);

typedef struct str_RF_Station_State {

  Filter_State noise_shape_fil, SSA_fil;
  Cavity_State cav_state;
  FPGA_State fpga_state;
  Delay_State loop_delay_state;

  // Noise signal for each sampled signal
  double complex probe_ns, rev_ns, fwd_ns;

} RF_State;

typedef RF_State* RF_State_p;
typedef RF_State** RF_State_dp;

void RF_Station_Allocate_In(RF_Station * rf_station,
  double Tstep,
  double Clip,
  double PAmax,
  double PAscale,
  double PAbw,
  double noise_shape_bw,
  Cavity *cav,
  double stable_gbw,
  double control_zero,
  double FPGA_out_sat,
  int loop_delay_size,
  double probe_ns_rms,
  double rev_ns_rms,
  double fwd_ns_rms);

RF_Station * RF_Station_Allocate_New(
  double Tstep,
  double Clip,
  double PAmax,
  double PAscale,
  double PAbw,
  double noise_shape_bw,
  Cavity *cav,
  double stable_gbw,
  double control_zero,
  double FPGA_out_sat,
  int loop_delay_size,
  double probe_ns_rms,
  double rev_ns_rms,
  double fwd_ns_rms);

void RF_Station_Deallocate(RF_Station *rf_station);

void RF_State_Allocate(RF_State *rf_state, RF_Station *rf_station);
void RF_State_Deallocate(RF_State *rf_state, RF_Station *rf_station);

double complex RF_Station_Step(
  RF_Station *rf_station,
  double delta_tz, double complex beam_current, double complex feed_forward,
  RF_State *rf_state);

/**
 * Step a Solid-State Amplifier (SSA) in time
 */
double complex SSA_Step(RF_Station *rf_station, double complex drive_in,
  RF_State *rf_state);

void RF_Station_Clear(RF_Station *rf_station, RF_State *rf_state);

/**
 * Exactly what it sounds like, apply a phase shift to a complex signal.
 */
double complex Phase_Shift(double complex in, double theta);
double complex Saturate(double complex in, double harshness);
void Apply_LLRF_Noise(RF_Station *rf_station, RF_State *rf_state);

void SSA_Clear(RF_Station *rf_station, RF_State *rf_state);

#endif
