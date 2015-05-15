#include "noise.h"

#include "stdlib.h"
#include "math.h"

void Noise_Step(int t_now, double Tstep, int type, double *settings, double *val)
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

double randn(double mu, double sigma)
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
