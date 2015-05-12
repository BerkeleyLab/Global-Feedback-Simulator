#ifndef DOUBLECOMPRESS_H
#define DOUBLECOMPRESS_H

#include "linac.h"

typedef struct str_doublecompress_state {
  // Variables of length Nlinac
  double *Ipk, *sz, *dE_E, *sd, *dt, *sdsgn, *k,
    *Eloss, *dE_Ei, *dE_Ei2;
  // Of length Nlinac
  double *cor;
} Doublecompress_State;

/*
 * A data structure to store dynamically set simulation variables,
 * i.e. they take on a new value each time-step
 */

typedef struct str_noise_srcs {
  
  double dQ_Q;
  double dtg;
  double dE_ing;
  double dsig_z;
  double dsig_E;
  double dchirp; //go into gun

} Noise_Srcs;

void Doublecompress_State_Allocate(Doublecompress_State * dcs, int Nlinac);
void Doublecompress_State_Deallocate(Doublecompress_State * dcs);
void Doublecompress_State_Attach(Doublecompress_State * dcs, int Nlinac, double * payload);

void Doublecompress(Gun * gun, Linac ** linac_array, int Nlinac, 
	//Inputs which change with time potentially
	Noise_Srcs *noise_srcs, double * dphivr, double * dV_Vvr,
	//double_compress output states
	Doublecompress_State * dcs
);

void Doublecompress_Octave_Benchmark(Gun * gun, Linac ** linac_array, int Nlinac, 
	//Inputs which change with time potentially
	Noise_Srcs *noise_srcs, double * dphivr, double * dV_Vvr,
	//double_compress output states
	Doublecompress_State * dcs
);

#endif
