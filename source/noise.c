/**
 * @file noise.c
 * @brief Noise Generator Model: Stepping functions for configurable noise sources.
 * @author Carlos Serrano (Cserrano@lbl.gov)
 */

#include "noise.h"

#include "stdlib.h"
#include "math.h"

/** Step function for Noise:
  * Calculates the Noise value for next simulation step.
  * It can be configured to generate several types of noise. */
void Noise_Step(
  int t_now,          ///< Current simulation time in seconds
  double Tstep,       ///< Simulation time step in seconds
  int type,           ///< 1 (White Noise), 2 (Sine Wave), 3 (Chirp), 4 (Step), 0 (do nothing)
  double *settings,   ///< Pointer to noise settings (different meaning according to type of noise)
  double *val         ///< Pointer to current generated noise value (output)
  )
{
  switch(type) {
  case 1:
    /* White Noise */
    *val = randn(0.0, settings[0]);
    break;
  case 2:
    /* Sine Wave */
    *val = settings[0]*( sin(settings[1] * t_now*Tstep ) );
    break;
  case 3:
    /* Chirp */
    *val = settings[0]*( sin(settings[1] * (t_now*Tstep)*(t_now*Tstep)) );
    break;
  case 4:
    /* Step */
    *val = settings[0] > (t_now+1)*Tstep ? 0.0 : settings[1];
    break;
  case 0:
  default:
    /* Do Nothing */
    break;
  }
}

/** Gaussian distribution pseudo-random number generator.
  * Returns noise value provided Mean and Standard Deviation.
  * Credit: http://phoxis.org/2013/05/04/generating-random-numbers-from-normal-distribution-in-c/ */
double randn(
  double mu,    ///< Mean
  double sigma  ///< Standard deviation
  )
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1){
    call = !call;
    return (mu + sigma * (double) X2);
  }

  do{
    U1 = -1 + ((double) rand () / RAND_MAX) * 2;
    U2 = -1 + ((double) rand () / RAND_MAX) * 2;
    W = pow (U1, 2) + pow (U2, 2);
  }

  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}
