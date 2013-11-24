#ifndef LINAC_H
#define LINAC_H


/*
 *  linac_param.h
 *
 * Definitions for a linear accelerator model
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
#include "filter.h"

/*
 * Data structure storing the paramters of the FPGA controller
 */
typedef struct str_FPGA {
  double gain;
  double int_gain;
  double complex set_point;
} FPGA_Param;

/*
 * Data structure storing the paramters of the Cavity array in a linac
 */
typedef struct str_cavity {
  // Seemed like fixed inputs?
  double psd_llrf, w0, bunch_rep, Q_L, R_Q;
  double beta_in, beta_out, beta_beam;

  // Calculated Quantities by cavity init
  double bandw, noise_rms, bw_ol, k;

  // Added in by the LLRF
  double nom_beam_phase, rf_phase, design_voltage,
    unity_voltage;
} Cavity;


/*
 * Data struture storing the parameters for a single linac
 * section comprised of multiple cavities.
 */
typedef struct str_linac_param {
  /* Paramters that are used by double compress */
  double dE, R56, T566, phi, lam, s0, a, L;

  /* Parameters used to set up my cavity */
  int n_cav;
  double nom_grad;

  double saturate_c;

  /* Sub objects */
  Filter RXF, TRF1, TRF2, Cav_Fil;
  FPGA_Param fpga;
  
  Cavity cav;

  double drift[3];
} Linac_Param;


/*Rename pointer for SWIG work around*/
typedef Linac_Param* LINP;
typedef Linac_Param ** LINdp;

/*
 * Data struture storing the beam parameters on 
 * exit from gun for Double Compress
 */

typedef struct str_gun_param {
  /* Paramters that are used by double compress */
  /* and for calculating dE_E */
  double E, sz0, sd0;
  double Q;
} Gun_Param;



/*
 * This routine takes all of the values read in from a file (TBD in python)
 * and fills out a Linac_Param structure with the calculated values.
 * Does not allocate.
 */
void Linac_Config(Linac_Param * linp,
		  // General (nonphysical) simulation paramters
		  double dt,
		  // Linac paramters
		  double dE, double R56, double T566, double phi,
		  double lam, double s0, double a, double L,
		  // Double compress parameters
		  // double Eg,
		  // Properties of the triode system
		  double saturate_c,
		  double complex * p_TRF1, double complex * p_TRF2,
		  // Properties of the RX Filter
		  double complex * p_RXF,
		  //Properties of the cavity 
		  int n_cav, double nom_grad,
		  double psd_llrf, double w0, double bunch_rep,
		  double Q_L, double R_Q,
		  double beta_in, double beta_out, double beta_beam,
		  double kly_max_v, double stable_gbw
		  );

void Linac_Deallocate(Linac_Param * linp);


void Cavity_Config(Cavity * cav,
		   double dt,int n_cav,
		   double psd_llrf, double w0, double bunch_rep,
		   double Q_L, double R_Q,
		   double beta_in, double beta_out, double beta_beam);


void FPGA_Config(FPGA_Param * fpga,
		 double kp, double ki, double complex set_point);

void Gun_Config(Gun_Param * gun, 
		 double E, double sz0, double sd0, double Q);

#endif
