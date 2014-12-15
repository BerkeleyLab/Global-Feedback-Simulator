
#include "cavity.h"

#include <math.h>
#include <stdlib.h>


// void ElecMode_Allocate_In(ElecMode *elecMode,
//   double RoverQ, double omega0, double foffset,
//   double Q_0, double Q_drive, double Q_probe,
//   double beam_phase,  double phase_refl, double phase_probe,
//   double Tstep, double d_phase,
//   double *mech_couplings)
// {

//   // Loaded Q
//   Q_L = 1/(1/Q_0 + 1/Q_drive + 1/Q_probe);

//   // Beam impedance/derivator (conversion factor between charge and beam induced voltage)
//   // Includes shift to take beam phase relative to the RF into account
//   elecMode -> k_beam = Tstep*RoverQ*Q_L*cexp(I*beam_phase);
  
//   // Drive port imperdance
//   elecMode -> k_drive = 2*sqrt(Q_drive*RoverQ)

//   // Probe port impedance (includes phase shift between cavity cell and probe ADC)
//   elecMode -> k_probe = cexp(-I*phase_probe)/sqrt(Q_probe*RoverQ)

//   // Reflected port impedance (includes phase shift between cavity cell and reflected ADC)
//   elecMode -> k_refl = cexp(-I*phase_refl)/sqrt(Q_probe*RoverQ)

//   // Mode's angular frequency (accelerator angular frequency + offset)
//   elecMode -> omega_mode = omega0 + 2*PI*foffset;

//   // Mode's open-loop bandwidth 
//   elecMode -> omega_f = omega0/2*Q_L;

//   // Mode's single-pole, low-pass filter allocation
//   // Calculate pole (mode's bandwidth)
//   double complex mode_p = -elecMode -> omega_f;

//   // Allocate single-pole filter
//   Filter_Allocate_In(&elecMode->fil,1,1);
//   Filter_Append_Modes(&elecMode->fil, &mode_p, 1, Tstep);

//   // Mechanical couplings
//   elecMode -> mech_couplings = mech_couplings;
// }

// ElecMode *ElecMode_Allocate_New(
//   double RoverQ, double omega0, double foffset,
//   double Q_0, double Q_drive, double Q_probe,
//   double beam_phase,  double phase_refl, double phase_probe,
//   double Tstep, double *mech_couplings)
// {
//   ElecMode * elecMode;
//   elecMode = calloc(1,sizeof(ElecMode));

//   ElecMode_Allocate_In(elecMode, RoverQ, omega0, foffset, Q_0, Q_drive, Q_probe, beam_phase,  phase_refl, phase_probe, Tstep, *mech_couplings);

//   return elecMode;
// }

// void step_ElecMode(ElecMode *elecMode,
//   // Inputs
//   double complex K1, double beam_charge, double delta_omega, double delta_tz, double Tstep,
//   // States
//   ElecMode_State *elecMode_state,
//   // Outputs
//   double complex *v_2, double complex *v_probe, double complex *v_refl)
// {
  
//   // Intermediate signals
//   double complex v_beam, v_drive, v_in, v_out;

//   // Instantaneous mode's angular frequency (including time-varying detuning)
//   omega_now = elecMode -> omega_mode + delta_omega;

//   // Beam-induced voltage (convert charge to voltage and add timing noise)
//   // k_beam = Tstep * (R/Q) * Q_L
//   v_beam = beam_charge * elecMode -> k_beam * cexp(I*omega_now*delta_tz);

//   // Drive term (k_drive = 2*sqrt(Q_drive*(R/Q))
//   v_drive = K1 * elecMode -> k_drive;

//   // Integrate detuning angular frequency to obtain phase
//   phase_now = elecMode_state-> d_phase + delta_omega * Tstep;
//   // Store phase state
//   elecMode_state-> d_phase = phase_now;

//   // Calculate mode's driving term (drive + beam)
//   v_in = (v_drive - v_beam)*cexp(I*phase_now)*elecMode->omega_f;

//   // Apply low-pass filter
//   v_out = Filter_Step(&(elecMode->fil), v_in, &(elecMode_state->fil_state));

//   // Calculate outputs based on v_vout
//   // Voltage squared
//   *v_2 = pow(cabs(v_out), 2);

//   // Probe
//   *v_probe = v_out * elecMode -> k_probe;

//   // Reflected
//   *v_refl = v_out * elecMode -> k_refl;
// }

// void Cavity_Allocate_In()
// {

// }

// Cavity *Cavity_Allocate_New()
// {

// }

// void Cavity_State_Allocate()
// {

// }

void Cavity_Config(Cavity * cav,
       double dt,int n_cav,
       double psd_llrf, double w0, double bunch_rep,
       double Q_L, double R_Q,
       double beta_in, double beta_out, double beta_beam,
       double complex cav_pole)
{
  /*
   * Input paramters
   */
  cav->psd_llrf = psd_llrf;
  cav->w0 = w0;
  cav->bunch_rep = bunch_rep;
  cav->Q_L = Q_L;
  cav->R_Q = R_Q;
  cav->beta_in = beta_in;
  cav->beta_out = beta_out;
  cav->beta_beam = beta_beam;
  /*
   * Calculated quantities
   */
  cav->bandw = 0.5/dt;
  cav->noise_rms = 1.5*sqrt( 0.5*psd_llrf*cav->bandw / (double)n_cav );
  cav->bw_ol = w0 / Q_L;

  cav->k = bunch_rep * R_Q * Q_L;

  Filter_Allocate_In(&cav->cav_fil,1,1);
  Filter_Append_Modes(&cav->cav_fil,&cav_pole,1,dt);

}

double complex Cavity_Step(Cavity *cav, double delta_tz,
     double complex drive_in, double complex beam_charge,
     Cavity_State *cav_state)
{
  double complex beam_charge_in, beam_induced_vol, voltage_in, cav_out;
  /*
   * FROM OCTAVE CODE, Carlos Serrano Dec2011:
   % Calculates beam charge. Note that amplitude and phase errors
   % for beam induced voltage are added in this
   % line. (cavity_params.nom_beam_phase is not an error, but a
   % nominal parameter)
  */
  beam_charge_in = beam_charge * 
    cexp(I*(cav->nom_beam_phase + delta_tz*cav->w0));

  /* printf("%lf %lf   %lf %lf   %lf %lf   %lf %lf\n", */
  /*   creal(beam_charge_in),cimag(beam_charge_in), */
  /*   creal(linp->cav->nom_beam_phase),cimag(linp->cav->nom_beam_phase), */
  /*   creal(linp->cav->w0),cimag(linp->cav->w0), */
  /*   creal(beam_charge),cimag(beam_charge)); */
  /*
   % Calculate beam induced voltage
   */
  beam_induced_vol = cav->k*beam_charge_in;
  /*
   % Add drive signal and beam induced voltage
   % Negative sign expressing energy absorption...
   */
  voltage_in = drive_in*cav->beta_in 
    + beam_induced_vol*cav->beta_beam;
  /*
   % Passes signal through cavity filter, output is smoothed out
  */

  cav_out = Filter_Step(&cav->cav_fil, voltage_in, &cav_state->fil_state);

  return cav_out;
}

void Cavity_Clear(Cavity *cav, Cavity_State *cav_state)
{
  // cav_state->voltage = 0.0+0.0*_Complex_I;
  Filter_State_Clear(&cav-> cav_fil, &cav_state->fil_state);
}