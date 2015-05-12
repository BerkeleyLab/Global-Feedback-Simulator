#ifndef NOISE_H
#define NOISE_H

#define N_NOISE_PARAM 8
#define N_NOISE_SET 2
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
 * 0 -> None  = 0
 * 1 -> White = s0 * rand[-1/2,1/2] 
 * 2 -> Sine  = s0 * sin( s1 * t );
 * 3 -> Chirp = s0 * sin( s1 * t^2 );
 * 4 -> Step  = 0.0 if time<s0 else s1
 * default -> None
 */

typedef struct str_Noise_Source {
  int type[N_NOISE_PARAM];
  double  settings[N_NOISE_PARAM*N_NOISE_SET];
} Noise_Source;

void Apply_Noise(int ti,double simdt, Noise_Source * ns,
		 Dynamic_Param * dynp);
#endif
