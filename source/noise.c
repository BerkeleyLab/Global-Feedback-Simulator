#include "noise.h"

#include "stdlib.h"
#include "math.h"

void Noise_Step(int t_now, double Tstep, int type, double *settings, double *val)
{
  switch(type) {
  case 1:
    /* White Noise */
    *val = settings[0]*(0.5- rand() / (double)RAND_MAX);
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
