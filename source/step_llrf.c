#include "step_llrf.h"

#include <math.h>

#include <stdlib.h>

Linac_State** make_state_arrays(int n) {
  Linac_State** lins = calloc(n, sizeof(Linac_State *));
  int i;
  for(i=0;i<n;i++) {
    lins[i] = calloc(1,sizeof(Linac_State));
  }
  return lins;
}
void cycle_buffer(Linac_State** linss, int n) {
  Linac_State * tmp;
  int i;
  tmp = linss[n-1];
  for(i=n-1;i>0;i--) {
    linss[i]=linss[i-1];
  }
  linss[0]=tmp;
}

void Linac_State_Allocate(Linac_State * lins, Linac_Param * linp) {
  Filter_State_Allocate(&lins->RXF,    &linp->RXF);
  Filter_State_Allocate(&lins->TRF1,   &linp->TRF1);
  Filter_State_Allocate(&lins->TRF2,   &linp->TRF2);
  Filter_State_Allocate(&lins->cav_state.fil_state, &linp->cav.cav_fil);

  lins->cav_state.voltage = 0.0 + 0.0*_Complex_I;
  
  lins->fpga.drive = 0.0 + 0.0*_Complex_I;
  lins->fpga.state = 0.0 + 0.0*_Complex_I;

  lins->RXF_out = 0.0 + 0.0*_Complex_I;
}

void Linac_State_Deallocate(Linac_State * lins) {
  Filter_State_Deallocate(&lins->RXF);
  Filter_State_Deallocate(&lins->TRF1);
  Filter_State_Deallocate(&lins->TRF2);
  Filter_State_Deallocate(&lins->cav_state.fil_state);
}

double complex phase_shift(double complex in, double theta) {
  return (creal(in)*cos(theta) - cimag(in)*sin(theta))
    + I*(creal(in)*sin(theta) + cimag(in)*cos(theta));
}

double complex step_fpga(FPGA_Param * fpga, double complex cavity_vol,
			 FPGA_State * stnow, int openloop)
{
  double complex sig_error = cavity_vol - fpga->set_point;
  if( openloop ) {
    stnow->state = fpga->set_point;
    stnow->drive = fpga->set_point;
  } else {
    double complex correction = sig_error*fpga->gain;
    
    stnow->state = stnow->state + fpga->int_gain*correction;
    if( cabs(stnow->state) > 1.0) stnow->state /= cabs(stnow->state);
    stnow->drive = stnow->state + correction;
    if( cabs(stnow->drive) > 1.0) stnow->drive /= cabs(stnow->drive);
  }
  stnow->err = sig_error;
  return sig_error;
}

void clear_fpga(FPGA_State * stnow)
{
  stnow-> drive = 0.0+0.0*_Complex_I;
  stnow-> state = 0.0+0.0*_Complex_I;
  stnow-> err = 0.0+0.0*_Complex_I;

}

double complex step_PI_fpga(FPGA_Param * fpga, double dt,
			    double complex cavity_vol,
			    FPGA_State * stnow, FPGA_State * stpast)
{
  double complex err = fpga->set_point - cavity_vol;
  if( 0 ) {
    stnow->state = fpga->set_point;
    stnow->drive = fpga->set_point;
  } else {
    stnow->state = stpast->state + dt*err;
    stnow->drive = stnow->state*fpga->gain*fpga->int_gain + fpga->gain*err;
  }
  return err;
}

double complex saturate(double complex in, double harshness) {
  return in*cpow( 1.0+cpow(cabs(in),harshness) , -1.0/harshness);
}

double complex step_triode(Linac_Param *linp, double complex drive_in, Linac_State * linnow)
{
  double complex trf1out, satout, trf2out;
  trf1out = Filter_Step(&linp->TRF1, drive_in, &linnow->TRF1);
  satout = saturate(trf1out,linp->saturate_c);
  trf2out = Filter_Step(&linp->TRF2, satout, &linnow->TRF2);
  return trf2out;
}

void clear_triode(Linac_Param *linp, Linac_State * linnow)
{
  Filter_State_Clear(&linp->TRF1, &linnow->TRF1);
  Filter_State_Clear(&linp->TRF2, &linnow->TRF2);
}

#define CPRINTs(c) {if(cimag(c)<0) printf("%10.16e%10.16ej ",creal(c),cimag(c)); else printf("%10.16e+%10.16ej ",creal(c),cimag(c)); }

double complex step_llrf(Linac_Param *linp,
			 double dt, double delta_tz,
			 double complex beam_charge,
			 double complex adc_noise,
			 int openloop,
			 Linac_State ** linss)
{
  Linac_State * linnow = linss[0];
  double complex sig_error,cavity_meas,fpga_drive_d,
    triode_out, triode_out_d, cav_out, cav_out_d,
    rxfvoltage;
  int i;
  int RXFDELAY = 1; //this is a delay number for the feedback...

  /*for(i=0;i<5;i++) {
    printf(" state is %lf,%lf\n",
	   creal(linss[i]->fpga.state),cimag(linss[i]->fpga.state));
	   }*/
 

 /*
   * Apply FPGA control
   */
  /*
  printf("linss[RXFDELAY]->RXF_out");
  CPRINTs(linss[RXFDELAY]->RXF_out);
  printf("\n");

  printf(" adc_noise");
  CPRINTs( adc_noise);
  printf("\n");
  */

  cavity_meas = linss[RXFDELAY]->RXF_out + adc_noise;
  /*
  printf("cavity_meas ");
  CPRINTs(cavity_meas);
  printf("\n");
  */
  sig_error = step_fpga(&linp->fpga, cavity_meas,
  			&linnow->fpga, openloop);
  /*
  printf("sig_error ");
  CPRINTs(sig_error);
  printf("\n");

  printf("linnow->fpga.driver ");
  CPRINTs(linnow->fpga.drive);
  printf("\n");
  */
  // drift on fpga->drive
  fpga_drive_d = phase_shift( linnow->fpga.drive, linp->drift[0] );

  /*
  printf("fpga_drive_d ");
  CPRINTs(fpga_drive_d);
  printf("\n");
  //Call Triode
   */
  triode_out = step_triode(linp, fpga_drive_d,
	      linnow);

  /*
  printf("triode_out ");
  CPRINTs(triode_out);
  printf("\n");
  */
  // drift on triode_output
  triode_out_d = phase_shift( triode_out, linp->drift[1]);
  
  /*
  printf("triode_out_d ");
  CPRINTs(triode_out_d);
  printf("\n");
  */

  /*
   * Drive the Cavity
   */
  cav_out = Cavity_Step(&linp->cav, delta_tz, triode_out_d, beam_charge, &linnow->cav_state);
  
  /*
  printf("cav_out ");
  CPRINTs(cav_out);
  printf("\n");
  */
  linnow->cav_state.voltage = cav_out;
  // drift
  cav_out_d = phase_shift( cav_out, linp->drift[2]);

  /*
  printf("cav_out_d ");
  CPRINTs(cav_out_d);
  printf("\n");
  */
  /*
   * RXF
   */
  linnow->RXF_out = Filter_Step(&linp->RXF, cav_out_d,
  		&linnow->RXF);
  
  /*
  printf("linnow->RXF_out ");
  CPRINTs(linnow->RXF_out);
  printf("\n");
  */
  
  return linnow->RXF_out;
}

void clear_linac(Linac_Param *linp, Linac_State * linnow)
{
  clear_fpga(&linnow->fpga);
  Cavity_Clear(&linp->cav, &linnow->cav_state);
  Filter_State_Clear(&linp->RXF, &linnow->RXF);
  clear_triode(linp, linnow);
  linnow->RXF_out = 0.0+0.0*_Complex_I;
}

#undef CPRINTs
