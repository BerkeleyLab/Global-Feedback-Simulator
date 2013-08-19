#include "dynamic_noise.h"

#include "stdlib.h"

void do_item(int ti,double simdt, int type,double *settings, double *val) {
  switch(type) {
  case 1:
    /* White Noise */
    *val = settings[0]*(0.5- rand() / (double)RAND_MAX);
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
  do_item(ti,simdt, ns->type[0],ns->settings+N_NOISE_SET*0, &dynp->dQ_Q);
  do_item(ti,simdt, ns->type[1],ns->settings+N_NOISE_SET*1, &dynp->dtg);
  do_item(ti,simdt, ns->type[2],ns->settings+N_NOISE_SET*2, &dynp->dE_ing);
  do_item(ti,simdt, ns->type[3],ns->settings+N_NOISE_SET*3, &dynp->dsig_z);
  do_item(ti,simdt, ns->type[4],ns->settings+N_NOISE_SET*4, &dynp->dsig_E);
  do_item(ti,simdt, ns->type[5],ns->settings+N_NOISE_SET*5, &dynp->dchirp);
  // ummm
  //do_item(ti,simdt, ns->type[6],ns->settings+N_NOISE_SET*6, (double *)&dynp->adc_noise);
  //do_item(ti,simdt, ns->type[7],ns->settings+N_NOISE_SET*7, (double *)&dynp->adc_noise);

}
