/*
 * Unit tests for linac_params
 *
 * Alejandro F Queiruga
 * Daniel "Scott" Driver
 * 2013 LBL
 *
 */
#include <complex.h>

#include "linac_param.h"


int main(int argc, char ** argv) {
  Linac_Param linp;

  complex double p_TRF1[2] = { 1.0+1.0*I, 2.0-1.0*I };
  complex double p_TRF2[1] = { 3.0-3.0*I };
  Linac_Config(&linp,
	       1.0,2.0,3.0,4.0,
	       5.0,6.0,7.0,8.0,
	       9,10.0,11.0,
	       p_TRF1,p_TRF2);

  return 0;
}
