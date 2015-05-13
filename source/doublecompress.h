#ifndef DOUBLECOMPRESS_H
#define DOUBLECOMPRESS_H

#include "linac.h"
#include "noise.h"

/*
 * A data structure to store dynamically set simulation variables,
 * i.e. they take on a new value each time-step
 */

typedef struct str_doublecompress_state {
  // Variables of length Nlinac
  double *Ipk, *sz, *dE_E, *sd, *dt, *sdsgn, *k,
    *Eloss, *dE_Ei, *dE_Ei2;
  // Of length Nlinac
  double *cor;
} Doublecompress_State;

/*
 * Data structure for storing the noise configurations
 * for the entries in Noise_Srcs
 *
 * The indices into the arrays map to entries in Noise_Srcs like so:
 * 0 -> dQ_Q
 * 1 -> dtg
 * 2 -> dE_ing
 * 3 -> dsig_z
 * 4 -> dsig_E
 * 5 -> dchirp
 *
 * The types are keys for what type of white noise to apply:
 * 0 -> None  = 0
 * 1 -> White = s0 * rand[-1/2,1/2]
 * 2 -> Sine  = s0 * sin( s1 * t );
 * 3 -> Chirp = s0 * sin( s1 * t^2 );
 * 4 -> Step  = 0.0 if time<s0 else s1
 * default -> None
 */

#define N_NOISE_SRCS 6
#define N_NOISE_ATTRIBUTES 2

typedef struct str_noise_srcs {

  double dQ_Q;
  double dtg;
  double dE_ing;
  double dsig_z;
  double dsig_E;
  double dchirp; //go into gun

  int type[N_NOISE_SRCS];
  double  attributes[N_NOISE_SRCS*N_NOISE_ATTRIBUTES];

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
