/**
 * @file rf_station.c
 * @brief RF Station Model: Allocation, configuration and stepping functions for the RF Station model,
 * including a cavity (with an arbitrary number of Eigenmodes), FPGA controller,
 * Solid-State Amplifier, loop delay and LLRF noise.
 * @author Carlos Serrano (Cserrano@lbl.gov)
 */

#include "rf_station.h"

#include <math.h>

/** Allocates memory for an array of RF Stations.
  * RF Stations themselves need to be allocated and filled individually,
  * and appended to this array using RF_Station_Append.*/
RF_Station_dp RF_Station_Allocate_Array(int n)
{
  RF_Station_dp RF_Station_net = calloc(n, sizeof(RF_Station *));
  return RF_Station_net;
}

/** Append an RF Station previously allocated and filled to the array given as an argument. */
void RF_Station_Append(RF_Station** rf_station_arr, RF_Station* rf_station, int index)
{
  // XXX Add some check!!
  rf_station_arr[index] = rf_station;
}

/** Takes a previously configured Delay and allocates its State struct accordingly. */
void Delay_State_Allocate(Delay *delay, Delay_State *delay_state)
{
  // Allocate memory buffer
  delay_state -> buffer = (double complex*) calloc(delay->size,sizeof(double complex));
  // Initialize buffer
  for(int i=0 ; i<delay->size; i++){
    delay_state -> buffer[delay_state->index] = (double complex) 0.0;
  }
  // Initialize index
  delay_state -> index = 0;
}

/** Frees memory of Delay State struct. */
void Delay_State_Deallocate(Delay_State *delay_state)
{
  free(delay_state->buffer);
}

/** Frees memory of Delay struct. */
void Delay_Deallocate(Delay *delay)
{
  delay->size = 0;
}

/** Step function for Delay: Takes a complex signal in, iterates a circular buffer and returns delayed value. */
double complex Delay_Step(
  double complex in,        ///< Input (complex) vector
  Delay *delay,             ///< Pointer to Delay struct
  Delay_State *delay_state  ///< Pointer to Delay State
  )
{
  double complex out;

  if(delay->size==0)
    out = in;
  else{
    out = delay_state -> buffer[delay_state->index];
    delay_state -> buffer[delay_state->index] = in;
    delay_state->index = (delay_state->index+1) % (delay->size);
  }

  return out;
}

/** Helper routine to zero out Delay buffer. Useful to restore initial state in unit tests. */
void Delay_Clear(Delay *delay, Delay_State *delay_state){
  for(int i=0 ; i<delay->size; i++){
    delay_state -> buffer[delay_state->index] = (double complex) 0.0;
  }
  delay_state->index = 0;
}

/** Takes a pointer to a RF Station struct which has been previously allocated and fills it in with the values passed as arguments. */
 void RF_Station_Allocate_In(
  RF_Station * rf_station,    ///< Pointer to RF Station
  // General (nonphysical) simulation parameters
  double Tstep,               ///< Simulation time-step in seconds
  // Properties of the SSA
  double Clip,                ///< Harshness parameter of SSA clipping function
  double PAmax,               ///< Maximum SSA output power sqrt(W)
  double PAscale,             ///< SSA scaling (from unitless to sqrt(W))
  double PAbw,                ///< SSA bandwidth (Hz)

  // Properties of the Noise Shaping filter
  double noise_shape_bw,     ///< Noise Shaping low-pass filter bandwidth in Hz

  // Cavity
  Cavity *cav,                ///< RF cavity
  // FPGA controller
  double stable_gbw,          ///< FPGA controller Gain-Bandwidth product in Hz
  double control_zero,        ///< FPGA controller zero location in Hz
  double FPGA_out_sat,        ///< FPGA output saturation limit

  // Loop Delay
  int loop_delay_size,        ///< RF feedback loop delay in simulation time steps (multiply by Tstep to get seconds)

  // LLRF Noise Sources
  double probe_ns_rms,      ///< LLRF noise: Cavity field probe port
  double rev_ns_rms,        ///< LLRF noise: Reverse port
  double fwd_ns_rms         ///< LLRF noise: Forward port
  )
{
  rf_station->Clip = Clip;
  rf_station->PAscale = PAscale;
  rf_station->PAmax = PAmax;

  /*
  * Configure the Filters using their poles
  */

  // SSA
  double complex p_SSA_fil = -2.0*M_PI*PAbw;
  Filter_Allocate_In(&rf_station->SSA_fil,1,1);
  Filter_Append_Modes(&rf_station->SSA_fil, &p_SSA_fil,1,Tstep);

  // Noise-shaping filter
  double complex p_noise_shape_fil = -2.0*M_PI*noise_shape_bw;
  Filter_Allocate_In(&rf_station->noise_shape_fil,1,1);
  Filter_Append_Modes(&rf_station->noise_shape_fil,&p_noise_shape_fil,1,Tstep);

  /*
  * Assign previously configured Cavity
  */
  rf_station->cav = cav;

  /*
  * Configure FPGA controller
  */
  // Grab cavity's open loop bandwidth and convert into Hz
  double cav_open_loop_bw = cav -> elecMode_net[cav->fund_index]->omega_f/(2.0*M_PI);
  // Calculate proportional gain (kp) based on cavity's open-loop bandwidth and controllers Gain-Bandwidth product
  double kp =  -stable_gbw/cav_open_loop_bw;
  double ki = kp*2.0*M_PI*control_zero;

  // Find fundamental mode couplings (could emulate some sort of calibration procedure here)
  double complex fund_k_probe = cav -> elecMode_net[cav->fund_index]-> k_probe;
  double complex fund_k_drive = cav -> elecMode_net[cav->fund_index]-> k_drive;

  // Calculate the FPGA controller set-point
  // (scale using probe and drive couplings to the fundamental mode to convert to FPGA units)
  double complex set_point = cav->design_voltage*cexp(I*cav->rf_phase*M_PI/180)*fund_k_probe;

  // Configure FPGA: note that ki (integral gain) is defined as kp/10.
  FPGA_Allocate_In(&rf_station->fpga, kp, ki, set_point, FPGA_out_sat, Tstep);

  // Coarse Loop Delay
  rf_station->loop_delay.size = loop_delay_size;

  rf_station->probe_ns_rms = probe_ns_rms;
  rf_station->rev_ns_rms = rev_ns_rms;
  rf_station->fwd_ns_rms = fwd_ns_rms;

}

/** Allocates memory for an RF Station and fills it in with the values passed as arguments. Returns a pointer to the newly allocated struct. */
RF_Station * RF_Station_Allocate_New(
  // General (nonphysical) simulation parameters
  double Tstep,               ///< Simulation time-step in seconds
  // Properties of the SSA
  double Clip,                ///< Harshness parameter of SSA clipping function
  double PAmax,               ///< Maximum SSA output power sqrt(W)
  double PAscale,             ///< SSA scaling (from unitless to sqrt(W))
  double PAbw,                ///< SSA bandwidth (Hz)

  // Properties of the RX Noise Shaping filter
  double noise_shape_bw,     ///< Noise Shaping low-pass filter bandwidth in Hz

  // Cavity
  Cavity *cav,                ///< RF cavity
  // FPGA controller
  double stable_gbw,          ///< FPGA controller Gain-Bandwidth product
  double control_zero,        ///< FPGA controller zero location
  double FPGA_out_sat,        ///< FPGA output saturation limit

  // Loop Delay
  int loop_delay_size,        ///< RF feedback loop delay in simulation time steps (multiply by Tstep to get seconds)

  // LLRF Noise Sources
  double probe_ns_rms,      ///< LLRF noise: Cavity field probe port
  double rev_ns_rms,        ///< LLRF noise: Reverse port
  double fwd_ns_rms         ///< LLRF noise: Forward port
  )
{

  RF_Station * rf_station;
  rf_station = calloc(1,sizeof(RF_Station));

  RF_Station_Allocate_In(rf_station, Tstep, Clip, PAmax, PAscale, PAbw,
    noise_shape_bw, cav, stable_gbw, control_zero, FPGA_out_sat, loop_delay_size,
    probe_ns_rms, rev_ns_rms, fwd_ns_rms);

  return rf_station;
}

/** Frees memory of an RF Station struct. */
void RF_Station_Deallocate(RF_Station *rf_station)
{

  Filter_Deallocate(&rf_station->noise_shape_fil);
  Filter_Deallocate(&rf_station->SSA_fil);

  FPGA_Deallocate(&rf_station->fpga);
  Cavity_Deallocate(rf_station->cav);
  Delay_Deallocate(&rf_station->loop_delay);

  rf_station->nom_grad = 0.0;
  rf_station->Clip = 0.0;
  rf_station->PAscale = 0.0;
  rf_station->PAmax = 0.0;

}

/** Takes a previously configured RF Station and allocates its State struct accordingly. */
void RF_State_Allocate(RF_State *rf_state, RF_Station *rf_station){

  Filter_State_Allocate(&rf_state->noise_shape_fil,    &rf_station->noise_shape_fil);
  Filter_State_Allocate(&rf_state->SSA_fil,   &rf_station->SSA_fil);
  Cavity_State_Allocate(&rf_state->cav_state, rf_station->cav);

  rf_state->cav_state.E_probe = (double complex) 0.0;
  rf_state->cav_state.E_reverse = (double complex) 0.0;

  rf_state->fpga_state.drive = (double complex) 0.0;
  rf_state->fpga_state.state = (double complex) 0.0;
  rf_state->fpga_state.openloop = (int) 0;

  Delay_State_Allocate(&rf_station->loop_delay, &rf_state->loop_delay_state);
}

/** Frees memory of RF Station State struct. */
void RF_State_Deallocate(RF_State *rf_state, RF_Station *rf_station) {
  Filter_State_Deallocate(&rf_state->noise_shape_fil);
  Filter_State_Deallocate(&rf_state->SSA_fil);
  Cavity_State_Deallocate(&rf_state->cav_state, rf_station->cav);
}

/** Apply a phase shift of theta radians to complex signal in*/
double complex Phase_Shift(double complex in, double theta) {
  return in*cexp(I*theta);
}

/** Takes a pointer to an FPGA struct which has been previously allocated and fills it in with the values passed as arguments. */
void FPGA_Allocate_In(
  FPGA * fpga,                ///< Pointer to FPGA
  double kp,                  ///< Proportional Gain
  double ki,                  ///< Integral Gain
  double complex set_point,   ///< Set-point
  double out_sat,             ///< Output saturation limit
  double Tstep                ///< Simulation time-step in seconds
  )
{

  fpga->kp = kp;
  fpga->ki = ki;

  fpga->set_point = set_point;

  // Saturation limits
  fpga->out_sat = out_sat;
  fpga->state_sat = out_sat;

  fpga->Tstep = Tstep;
}

/** Frees memory of an FPGA struct. */
void FPGA_Deallocate(FPGA *fpga)
{
  fpga->kp = 0.0;
  fpga->ki = 0.0;
  fpga->set_point = 0.0;
  fpga->out_sat = 0.0;
  fpga->state_sat = 0.0;
  fpga->Tstep = 0.0;
}

/** Helper routine to zero out FPGA State. Useful to restore initial state in unit tests. */
void FPGA_Clear(FPGA_State * stnow)
{
  stnow-> drive = 0.0+0.0*_Complex_I;
  stnow-> state = 0.0+0.0*_Complex_I;
  stnow-> err = 0.0+0.0*_Complex_I;
  stnow-> openloop = 0;
}

/** Step function for FPGA:
  * Calculates the state for the next simulation step.
  * Returns the error signal and stores current state in State struct. */
double complex FPGA_Step(
  FPGA *fpga,                   ///< Pointer to FPGA
  double complex cavity_vol,    ///< Measured cavity field
  double complex feed_forward,  ///< Measured cavity field
  FPGA_State *stnow             ///< Pointer to FPGA State
  )
{
  double complex state, drive;
  double scale;

  // Calculate error signal
  double complex err = cavity_vol - fpga->set_point;

  if(stnow->openloop == 1) { // Open loop
    stnow->drive = fpga->set_point;
    stnow->state = fpga->set_point;
  } else {  //Closed loop
    // Integrator state
    state = stnow->state + feed_forward + fpga->Tstep*err*fpga->ki;
    // state = stnow->state + fpga->Tstep*fpga->ki*0.5*(err+stnow->err);
    // Compare state magnitude and saturation limit
    scale = cabs(state)/fpga->state_sat;
    // Clip integrator state if above limit
    stnow->state = (scale > 1.0) ? state/scale : state;

    // Drive signal
    drive = stnow->state + fpga->kp*err;
    // Compare output drive magnitude with saturation limit
    scale = cabs(drive)/fpga->out_sat;
    // Clip drive output if above limit
    stnow->drive = (scale > 1.0) ? drive/scale : drive;

    // Store error state
    stnow->err = err;
  }
  // Return FPGA error signal
  return err;
}

/** Implement Saturation of input signal for SSA. */
double complex Saturate(
  double complex in,    ///< Input (complex) signal
  double harshness      ///< Saturation harshness parameter
  )
{
  return in*cpow(1.0+cpow(cabs(in),harshness), -1.0/harshness);
}

/** Step function for Solid-State Amplifier: Saturation and band limit.
  * Calculates the state for the next simulation step.
  * Returns the error signal and stores current state in State struct. */
double complex SSA_Step(
  RF_Station *rf_station,   ///< Pointer to RF Station (contains SSA struct)
  double complex drive_in,  ///< Input (complex) signal
  RF_State *rf_state        ///< Pointer to RF Station State (contains SSA State)
  )
{
  double complex fil_out, satout;

  // Scale input signal (sqrt(W) -> Normalized units)
  drive_in = drive_in/rf_station->PAscale;
  // Apply low-pass filter to limit bandwidth (SSA_fil)
  fil_out = Filter_Step(&rf_station->SSA_fil, drive_in, &rf_state->SSA_fil);
  // Clip
  satout = Saturate(fil_out,rf_station->Clip);
  // Scale output signal (sqrt(W) -> Normalized units)
  satout = satout*rf_station->PAscale;

  // Return SSA output
  return satout;
}

/** Helper routine to zero out SSA State. Useful to restore initial state in unit tests. */
void SSA_Clear(RF_Station *rf_station, RF_State * rf_state)
{
  Filter_State_Clear(&rf_station->SSA_fil, &rf_state->SSA_fil);
}

/** Update Noise signals for each LLRF port. LLRF noise is defined by its Power Spectral Density
  * and considered Gaussian, where 1/f noise is ignored (see Physics documentation) */
void Apply_LLRF_Noise(RF_Station *rf_station, RF_State *rf_state)
{
  rf_state->probe_ns = randn(0.0, rf_station->probe_ns_rms) + _Complex_I*randn(0.0, rf_station->probe_ns_rms);
  rf_state->rev_ns = randn(0.0, rf_station->rev_ns_rms) + _Complex_I*randn(0.0, rf_station->rev_ns_rms);
  rf_state->fwd_ns = randn(0.0, rf_station->fwd_ns_rms) + _Complex_I*randn(0.0, rf_station->fwd_ns_rms);
}

/** Step function for RF Station:
  * Calculates the state for the next simulation step.
  * Returns the cavity accelerating voltage (as seen by the beam)
  * and stores current state in State struct. */
double complex RF_Station_Step(
  // Configuration parameters
  RF_Station *rf_station,       ///< Pointer to RF Station
  // Beam amplitude and phase noise (beam loading)
  double delta_tz,              ///< Timing jitter in seconds (RF reference noise)
  double complex beam_current,   ///< Beam current in Amps
  double complex feed_forward,   ///< Feed-forward signal
  // State vectors
	RF_State *rf_state            ///< Pointer to RF State
  )
{

  // Intermediate signals
  double complex E_probe_delayed, E_probe_lp, sig_error, Kg, V_acc;

  // LLRF noise (probe, forward and reverse signals)
  // Generate Gaussian Noise
  Apply_LLRF_Noise(rf_station, rf_state);

  // Apply LLRF Noise
  rf_state->cav_state.E_probe += rf_state->probe_ns;
  rf_state->cav_state.E_reverse += rf_state->rev_ns;

  // // Apply (overall) loop delay to cavity probe signal
  E_probe_delayed = Delay_Step(rf_state->cav_state.E_probe, &rf_station->loop_delay, &rf_state->loop_delay_state);

  // Apply noise-shaping low-pass filter
  E_probe_lp = Filter_Step(&rf_station->noise_shape_fil, E_probe_delayed, &rf_state->noise_shape_fil);

  // FPGA
  sig_error = FPGA_Step(&rf_station->fpga, E_probe_lp, feed_forward, &rf_state->fpga_state);

  // SSA (includes band-limiting low-pass filter and saturation)
  Kg = SSA_Step(rf_station, rf_state->fpga_state.drive, rf_state);

  // Cavity
  V_acc = Cavity_Step(rf_station->cav, delta_tz, Kg, beam_current, &rf_state->cav_state);

  // Store digitized forward signal
  rf_state->cav_state.E_fwd = Kg + rf_state->fwd_ns;

  //  Overall cavity accelerating voltage (as seen by the beam)
  return V_acc;
}

/** Helper routine to zero out RF Station State. Useful to restore initial state in unit tests. */
void RF_Station_Clear(RF_Station *rf_station, RF_State * rf_state)
{
  FPGA_Clear(&rf_state->fpga_state);
  Cavity_Clear(rf_station->cav, &rf_state->cav_state);
  Filter_State_Clear(&rf_station->noise_shape_fil, &rf_state->noise_shape_fil);
  SSA_Clear(rf_station, rf_state);
  Delay_Clear(&rf_station->loop_delay, &rf_state->loop_delay_state);
}
