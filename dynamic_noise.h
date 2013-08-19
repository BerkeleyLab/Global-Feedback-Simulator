#ifndef DYNAMIC_NOISE_H
#define DYNAMIC_NOISE_H


#include "doublecompress.h"


#define N_NOISE_PARAM 8
#define N_NOISE_SET 1
/*
 * Data structure for storing the noise configurations
 * for the entries in Dynamic_Param.
 *
 * The indices into the arrays map to entries in Dynamic_Param like so:
 * 0 -> dQ_Q
 * 1 -> dtg
 * 2 -> dE_ing
 * 3 -> dsig_z
 * 4 -> dsig_E
 * 5 -> dchirp
 * 6,7 -> real,imag of adc_noise
 *
 * The types are keys for what type of whitenoise to apply:
 * 0 -> None
 * 1 -> White
 * default -> None
 */

typedef struct str_Noise_Source {
  int type[N_NOISE_PARAM];
  double  settings[N_NOISE_PARAM*N_NOISE_SET];
} Noise_Source;

void Apply_Noise(int ti,double simdt, Noise_Source * ns,
		 Dynamic_Param * dynp);
#endif
