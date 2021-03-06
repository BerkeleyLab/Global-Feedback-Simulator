/**
  * @file cavity.h
  * @brief Header file for cavity.c
  * @author Carlos Serrano (CSerrano@lbl.gov)
*/

#ifndef CAVITY_H
#define CAVITY_H

#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include "filter.h"

typedef struct str_elecmode{
  double LO_w0;
  double omega_f, omega_d_0;
	double complex k_drive, k_beam, k_probe, k_em;
	double Tstep;
	Filter fil;
  double *A;  ///< Matrix coefficients: Convert square of accelerating voltage to Force
	double *C;  ///< Matrix coefficients: Convert displacement to frequency shift
} ElecMode;

typedef ElecMode * ElecMode_p;
typedef ElecMode ** ElecMode_dp;

typedef struct str_elecmode_state{
	Filter_State fil_state;
	double delta_omega;
	double d_phase;
	double V_2;
} ElecMode_State;

typedef struct str_cavity{
	double rf_phase, design_voltage;
	int fund_index;
	double L;
	double nom_grad;
	int n_modes;
	ElecMode_dp elecMode_net;
} Cavity;

typedef struct str_cavity_state{
	double complex E_probe, E_reverse, E_fwd, V;
  double complex Kg;
	ElecMode_State **elecMode_state_net;
} Cavity_State;

void Cavity_Allocate_In(Cavity *cav,
  ElecMode_dp elec_mode_net, int n_modes,
  double L, double nom_grad,
  double rf_phase, double design_voltage,
  int fund_index);

Cavity * Cavity_Allocate_New(ElecMode_dp elec_mode_net, int n_modes,
  double L, double nom_grad,
  double rf_phase, double design_voltage,
  int fund_index);

void Cavity_Deallocate(Cavity *cav);

double complex Cavity_Step(Cavity *cav, double delta_tz, double complex drive_in, double complex beam_current, Cavity_State *cav_state);
void Cavity_Clear(Cavity *cav, Cavity_State *cav_state);

void Cavity_State_Allocate(Cavity_State *cav_state, Cavity *cav);
void Cavity_State_Deallocate(Cavity_State *cav_state, Cavity *cav);
ElecMode_State *ElecMode_State_Get(Cavity_State *cav_state, int idx);

void ElecMode_Allocate_In(ElecMode *elecMode,
  double RoverQ, double foffset, double LO_w0,
  double Q_0, double Q_drive, double Q_probe,
  double rf_phase,  double phase_rev, double phase_probe,
  double Tstep,
  double *mech_couplings, int n_mech
  );

ElecMode *ElecMode_Allocate_New(
  double RoverQ, double foffset, double LO_w0,
  double Q_0, double Q_drive, double Q_probe,
  double rf_phase,  double phase_rev, double phase_probe,
  double Tstep,
  double *mech_couplings, int n_mech
  );
void ElecMode_Deallocate(ElecMode * elecMode);

ElecMode_dp ElecMode_Allocate_Array(int n);
void ElecMode_Append(ElecMode** elecMode_arr, ElecMode* elecMode, int index);

void ElecMode_State_Allocate(ElecMode_State *elecMode_state, ElecMode *elecMode);
void ElecMode_State_Deallocate(ElecMode_State *elecMode_state);

double complex ElecMode_Step(ElecMode *elecMode, double complex Kg_fwd,
	double beam_current, double delta_tz,
  ElecMode_State *elecMode_state,
  double complex *v_probe, double complex *v_em);

#endif
