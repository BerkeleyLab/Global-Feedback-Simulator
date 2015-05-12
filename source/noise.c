#include "dynamic_noise.h"

#include "stdlib.h"
#include "math.h"

void Noise_Step(int ti, double simdt, int type, double *settings, double *val) {
  switch(type) {
  case 1:
    /* White Noise */
    *val = settings[0]*(0.5- rand() / (double)RAND_MAX);
    break;
  case 2:
    /* Sine Wave */
    *val = settings[0]*( sin(settings[1] * ti*simdt ) );
    break;
  case 3:
    /* Chirp */
    *val = settings[0]*( sin(settings[1] * (ti*simdt)*(ti*simdt)) );
    break;
  case 4:
    /* Step */
    *val = settings[0] > (ti+1)*simdt ? 0.0 : settings[1];
    break;
  case 0:
  default:
    /* Do Nothing */
    break;
  }
}

void Apply_Noise(int ti,double simdt, Noise_Source * ns,
		 Dynamic_Param * dynp)
{
  Noise_Step(ti,simdt, ns->type[0],ns->settings+N_NOISE_SET*0, &dynp->dQ_Q);
  Noise_Step(ti,simdt, ns->type[1],ns->settings+N_NOISE_SET*1, &dynp->dtg);
  Noise_Step(ti,simdt, ns->type[2],ns->settings+N_NOISE_SET*2, &dynp->dE_ing);
  Noise_Step(ti,simdt, ns->type[3],ns->settings+N_NOISE_SET*3, &dynp->dsig_z);
  Noise_Step(ti,simdt, ns->type[4],ns->settings+N_NOISE_SET*4, &dynp->dsig_E);
  Noise_Step(ti,simdt, ns->type[5],ns->settings+N_NOISE_SET*5, &dynp->dchirp);
  // Do error on the real component of ADC and then the imag component
  double adc_real = creal(dynp->adc_noise);
  double adc_imag = cimag(dynp->adc_noise);
  Noise_Step(ti,simdt, ns->type[6],ns->settings+N_NOISE_SET*6, 
	  (double *)&(adc_real) );
  Noise_Step(ti,simdt, ns->type[7],ns->settings+N_NOISE_SET*7, 
	  (double *)&(adc_imag) );
  dynp->adc_noise = adc_real + I*adc_imag;

}
