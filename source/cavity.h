#ifndef CAVITY_H
#define CAVITY_H

#include <complex.h>

#include "filter.h"


// typedef struct str_elecmode{
// 	double omega_0, omega_f;
// 	double k_drive, k_beam, k_probe;
// 	double d_phase;
// 	double *mech_couplings;
// } ElecMode;

// ElecMode * ElecMode_p;
// ElecMode ** ElecMode_dp;

// typedef struct str_cavity{
// 	double L;
// 	double nom_grad;
// 	ElecMode_dp elec_mode_net;
// } Cavity;

typedef struct str_cavity {
  // Seemed like fixed inputs?
  double psd_llrf, w0, bunch_rep, Q_L, R_Q;
  double beta_in, beta_out, beta_beam;

  // Calculated Quantities by cavity init
  double bandw, noise_rms, bw_ol, k;

  // Added in by the LLRF
  double nom_beam_phase, rf_phase, design_voltage, unity_voltage;

  Filter cav_fil;

} Cavity;

typedef struct str_cavity_state{
	double complex voltage;
	Filter_State fil_state;
} Cavity_State;

void Cavity_Config(Cavity * cav,
		   double dt,int n_cav,
		   double psd_llrf, double w0, double bunch_rep,
		   double Q_L, double R_Q,
		   double beta_in, double beta_out, double beta_beam,
		   double complex cav_pole);

double complex Cavity_Step(Cavity *cav, double delta_tz, double complex drive_in, double complex beam_charge, Cavity_State *cav_state);
void Cavity_Clear(Cavity *cav, Cavity_State *cav_state);

#endif