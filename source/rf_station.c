#include "rf_station.h"

#include <math.h>
#include <stdlib.h>

RF_State** make_rf_state_array(int n) {
  RF_State** rf_state_net = calloc(n, sizeof(RF_State *));
  int i;
  for(i=0;i<n;i++) {
    rf_state_net[i] = calloc(1,sizeof(RF_State));
  }
  return rf_state_net;
}

/*
 * Delay
 */
void Delay_State_Allocate(Delay *delay, Delay_State *delay_state)
{
  delay_state -> buffer = (double complex*) calloc(delay->size,sizeof(double complex));
}
void Delay_State_Deallocate(Delay *delay, Delay_State *delay_state)
{
  free(delay_state->buffer);
}

double complex Delay_Step(double complex in, Delay *delay, Delay_State *delay_state)
{
  double complex out = delay_state -> buffer[delay_state->index];
  delay_state -> buffer[delay_state->index] = in;
  delay_state->index += 1 % (delay->size);
  return out;
}

void Delay_Clear(Delay *delay, Delay_State *delay_state){
  for(int i=0 ; i<delay->size; i++){
    delay_state -> buffer[delay_state->index] = (double complex) 0.0;
  }
  delay_state->index = 0;
}

/*
 * RF Station
 */

void RF_Station_Allocate_In(RF_Station * rf_station,
  // General (nonphysical) simulation parameters
  double Tstep,
  // Properties of the triode system
  double saturate_c,
  double kly_max_v,
  double complex * p_TRF1, double complex * p_TRF2,
  // Properties of the RX Filter
  double complex * p_RXF,
  
  // Cavity
  Cavity *cav,
  // FPGA controller
  double stable_gbw, // Gain-Bandwidth product
  
  // Loop Delay
  int loop_delay_size
  )
{

  rf_station->saturate_c = saturate_c;
  
  /*
  * Configure the Filters using their poles
  */
  
  Filter_Allocate_In(&rf_station->TRF1,2,2);
  for(int i=0;i<2;i++) {
    Filter_Append_Modes(&rf_station->TRF1, p_TRF1+i, 1,Tstep);
  }
  
  Filter_Allocate_In(&rf_station->TRF2,1,1);
  Filter_Append_Modes(&rf_station->TRF2, p_TRF2, 1, Tstep);
  
  Filter_Allocate_In(&rf_station->RXF,3,3);
  for(int i=0;i<3;i++) {
    Filter_Append_Modes(&rf_station->RXF,p_RXF+i,1,Tstep);
  }
  
  /*
  * Assign previously configured Cavity
  */
  rf_station->cav = cav;

  /*
  * Configure FPGA controller
  */
  double cav_open_loop_bw = cav -> open_loop_bw;
  // Calculate proportional gain (kp) based on cavity's open-loop bandwith and controllers Gain-Bandwidth product
  double kp =  -stable_gbw/cav_open_loop_bw;

  // Calculate (normalized) controller's set-point
  double complex set_point = cav->design_voltage*cexp(I*cav->rf_phase)
                         / cav->unity_voltage;

  // Configure FPGA: note that ki (integral gain) is defined as kp/10.
  FPGA_Allocate_In(&rf_station->fpga,kp,kp*0.1, set_point, Tstep);

  // Coarse Loop Delay
  rf_station -> loop_delay.size = loop_delay_size;

}

void RF_State_Allocate(RF_State *rf_state, RF_Station *rf_station){

  Filter_State_Allocate(&rf_state->RXF,    &rf_station->RXF);
  Filter_State_Allocate(&rf_state->TRF1,   &rf_station->TRF1);
  Filter_State_Allocate(&rf_state->TRF2,   &rf_station->TRF2);
  for(int i=0; i<rf_station->cav->n_modes ; i++){
    Filter_State_Allocate(&rf_state->cav_state.elecMode_state_net[i]->fil_state, &rf_station->cav->elecMode_net[i]->fil);
  }

  rf_state->cav_state.E_probe = (double complex) 0.0;
  rf_state->cav_state.E_reverse = (double complex) 0.0;
  
  rf_state->fpga_state.drive = (double complex) 0.0;
  rf_state->fpga_state.state = (double complex) 0.0;

  Delay_State_Allocate(&rf_station->loop_delay, &rf_state->loop_delay_state);
}

void RF_State_Deallocate(RF_State *rf_state, RF_Station *rf_station) {
  Filter_State_Deallocate(&rf_state->RXF);
  Filter_State_Deallocate(&rf_state->TRF1);
  Filter_State_Deallocate(&rf_state->TRF2);
  Cavity_State_Deallocate(&rf_state->cav_state, rf_station->cav);
}

// double complex Phase_Shift(double complex in, double theta) {
//   return (creal(in)*cos(theta) - cimag(in)*sin(theta))
//     + I*(creal(in)*sin(theta) + cimag(in)*cos(theta));
// }

double complex Phase_Shift(double complex in, double theta) {
  return in*cexp(I*theta);
}

void FPGA_Allocate_In(FPGA * fpga,
     double kp, double ki, double complex set_point, double Tstep)
{

  fpga->kp = kp;
  fpga->ki = ki;

  fpga->set_point = set_point;
  fpga->Tstep = Tstep;
}

// double complex FPGA_Step(FPGA *fpga, double complex cavity_vol,
// 			 FPGA_State * stnow, int openloop)
// {
//   double complex sig_error = cavity_vol - fpga->set_point;
//   if( openloop ) {
//     stnow->state = fpga->set_point;
//     stnow->drive = fpga->set_point;
//   } else {
//     double complex correction = sig_error*fpga->kp;
    
//     stnow->state = stnow->state + fpga->int_gain*correction;
//     if( cabs(stnow->state) > 1.0) stnow->state /= cabs(stnow->state);
//     stnow->drive = stnow->state + correction;
//     if( cabs(stnow->drive) > 1.0) stnow->drive /= cabs(stnow->drive);
//   }
//   stnow->err = sig_error;
//   return sig_error;
// }

void FPGA_Clear(FPGA_State * stnow)
{
  stnow-> drive = 0.0+0.0*_Complex_I;
  stnow-> state = 0.0+0.0*_Complex_I;
  stnow-> err = 0.0+0.0*_Complex_I;
}

double complex step_PI_fpga(FPGA *fpga,
			    double complex cavity_vol,
			    FPGA_State *stnow, int openloop)
{
  double complex err = fpga->set_point - cavity_vol;
  if( openloop ) {
    stnow->state = fpga->set_point;
    stnow->drive = fpga->set_point;
  } else {
    stnow->state = stnow->state + fpga->Tstep*err;
    stnow->drive = stnow->state*fpga->ki + fpga->kp*err;
  }
  return err;
}

double complex Saturate(double complex in, double harshness) {
  return in*cpow( 1.0+cpow(cabs(in),harshness) , -1.0/harshness);
}

double complex Triode_Step(RF_Station *rf_station, double complex drive_in, RF_State *rf_state)
{
  double complex trf1out, satout, trf2out;
  trf1out = Filter_Step(&rf_station->TRF1, drive_in, &rf_state->TRF1);
  satout = Saturate(trf1out,rf_station->saturate_c);
  trf2out = Filter_Step(&rf_station->TRF2, satout, &rf_state->TRF2);
  return trf2out;
}

void Triode_Clear(RF_Station *rf_station, RF_State * rf_state)
{
  Filter_State_Clear(&rf_station->TRF1, &rf_state->TRF1);
  Filter_State_Clear(&rf_station->TRF2, &rf_state->TRF2);
}

#define CPRINTs(c) {if(cimag(c)<0) printf("%10.16e%10.16ej ",creal(c),cimag(c)); else printf("%10.16e+%10.16ej ",creal(c),cimag(c)); }


double complex RF_Station_Step(
  // Configuration parameters
  RF_Station *rf_station,
  // Beam amplitude and phase noise (beam loading)
  double delta_tz, double complex beam_charge,
  // LLRF noise (probe and reverse signals)
  double complex probe_ns, double complex rev_ns,
  // FPGA control boolean to open and close the feedback loop
  int openloop,
  // State vectors
	RF_State *rf_state)
{

  // Intermediate signals
  double complex E_probe_delayed, E_probe_lp, sig_error, Kg, V_acc;

  // Apply LLRF Noise
  rf_state->cav_state.E_probe =+ probe_ns;
  rf_state->cav_state.E_reverse =+ rev_ns;

  // Apply (overall) loop delay to cavity probe signal
  E_probe_delayed = Delay_Step(rf_state->cav_state.E_probe, &rf_station->loop_delay, &rf_state->loop_delay_state);

  // Apply Low-Pass filter to probe signal (Rx Filter)
  E_probe_lp = Filter_Step(&rf_station->RXF, E_probe_delayed, &rf_state->RXF);
  
  // FPGA
  sig_error = step_PI_fpga(&rf_station->fpga, E_probe_lp, &rf_state->fpga_state, openloop);

  // Triode (includes Tx filters and saturation)
  Kg = Triode_Step(rf_station, rf_state->fpga_state.drive, rf_state);

  // Cavity
  V_acc = Cavity_Step(rf_station->cav, delta_tz, Kg, beam_charge, &rf_state->cav_state);
  
  // Overall cavity accelerating voltage (as seen by the beam)
  return V_acc;
}

void RF_Station_Clear(RF_Station *rf_station, RF_State * rf_state)
{
  FPGA_Clear(&rf_state->fpga_state);
  Cavity_Clear(rf_station->cav, &rf_state->cav_state);
  Filter_State_Clear(&rf_station->RXF, &rf_state->RXF);
  Triode_Clear(rf_station, rf_state);
  Delay_Clear(&rf_station->loop_delay, &rf_state->loop_delay_state);
}

#undef CPRINTs
