/**
 * @file cavity.c
 * @brief Cavity Model: Allocation, configuration and stepping functions for the cavity model, which iterates over an arbitrary number of electrical eigenmodes.
 * @author Carlos Serrano (Cserrano@lbl.gov)
 */

#include "cavity.h"
#include <stdio.h>
#include <string.h>

/** Allocates memory for an array of Electrical Eigenmodes.
  * Electrical modes themselves need to be allocated and filled individually,
  * and appended to this array using ElecMode_Append.*/
ElecMode_dp ElecMode_Allocate_Array(int n   ///< Size of the array of Electrical Eigenmodes
  )
{
  ElecMode_dp elecMode_net = calloc(n, sizeof(ElecMode *));

  return elecMode_net;
}

/** Append an Electrical Eigenmode previously allocated and filled to the array given as an argument. */
void ElecMode_Append(
  ElecMode** elecMode_arr,  ///< Pointer to array of Electrical Eigenmodes
  ElecMode* elecMode,       ///< Pointer to Electrical Eigenmode to be appended
  int idx                   ///< Position of Electrical Eigenmode in the array
  )
{
  // XXX Add some check!!
  elecMode_arr[idx] = elecMode;
}

/** Takes a pointer to a ElecMode struct (representing an electrical eigenmode)
  * which has been previously allocated and fills it in with the values passed as arguments. */
void ElecMode_Allocate_In(
  ElecMode *elecMode,     ///< Pointer to ElecMode struct
  double RoverQ,          ///< Mode's (R/Q) in Ohms
  double foffset,         ///< Mode's frequency offset (with respect to the RF reference frequency) in Hz
  double LO_w0,           ///< Local Oscillator frequency in rad/s
  double Q_0,             ///< Represents losses in the cavity walls
  double Q_drive,         ///< Represents coupling to the input coupler
  double Q_probe,         ///< Represents coupling to the field probe
  double rf_phase,        ///< Beam phase relative to the RF in radians
  double phase_rev,       ///< Phase shift between Cavity cells and reverse ADC
  double phase_probe,     ///< Phase shift between Cavity cells and probe ADC
  double Tstep,           ///< Simulation time-step in seconds
  double *mech_couplings, ///< Pointer to array of mechanical couplings in (rad/s)/V^2
  int n_mech              ///< Number of mechanical eigenmodes
  )
{

  // Store time step size in order to integrate detune frequency and obtain phase
  elecMode -> Tstep = Tstep;

  // Loaded Q
  double Q_L = 1/(1/Q_0 + 1/Q_drive + 1/Q_probe);

  // Beam impedance/derivator (conversion factor between charge and beam induced voltage)
  // Includes shift to take beam phase relative to the RF into account
  elecMode -> k_beam = RoverQ*Q_L*cexp(-I*rf_phase)/Tstep*1e-12;

  // Drive port impedance
  elecMode -> k_drive = 2*sqrt(Q_drive*RoverQ);

  // Probe port impedance (includes phase shift between cavity cell and probe ADC)
  elecMode -> k_probe = cexp(I*phase_probe)/sqrt(Q_probe*RoverQ);

  // Emitted port impedance (includes phase shift between cavity cell and reverse ADC)
  elecMode -> k_em = cexp(I*phase_rev)/sqrt(Q_drive*RoverQ);

  // Linac nominal angular frequency
  elecMode -> LO_w0 = LO_w0;

  // Mode's resonance frequency (RF frequency + offset)
  double omega_0_mode = LO_w0 + 2*M_PI*foffset;
  // Mode's open-loop bandwidth
  elecMode -> omega_f = omega_0_mode/(2*Q_L);

  // Mode's baseline frequency offset
  elecMode -> omega_d_0 = 2*M_PI*foffset;

  // Mode's single-pole, low-pass filter allocation
  // Calculate pole (mode's bandwidth)
  double complex mode_p = -elecMode -> omega_f;

  // Append mode to cavity filter
  Filter_Append_Modes(&elecMode->fil, &mode_p, 1, Tstep);

  // Electro-Mechanical couplings
  // First allocate memory for A and C matrices
  elecMode -> C = calloc(n_mech, sizeof(double));
  elecMode -> A = calloc(n_mech, sizeof(double));

  // Fill in coefficients given mechanical couplings
  int i;
  for(i=0; i<n_mech;i++){
    // mech_couplings ((rad/s)/V^2) are always negative but sometimes referred to as positive quantities
    // Take the absolute value to let user configuration the freedom of expressing it either way
    elecMode -> A[i] = sqrt(fabs(mech_couplings[i])/RoverQ)/omega_0_mode;
    elecMode -> C[i] = -(omega_0_mode)*sqrt(fabs(mech_couplings[i]*RoverQ));
  }
}

/** Allocates memory for a ElecMode struct (representing an electrical eigenmode)
  * and fills it in with the values passed as arguments. Returns a pointer to the newly allocated struct. */
ElecMode * ElecMode_Allocate_New(
  double RoverQ,          ///< Mode's (R/Q) in Ohms
  double foffset,         ///< Mode's frequency offset (with respect to the RF reference frequency) in Hz
  double LO_w0,           ///< Local Oscillator frequency in rad/s
  double Q_0,             ///< Represents losses in the cavity walls
  double Q_drive,         ///< Represents coupling to the input coupler
  double Q_probe,         ///< Represents coupling to the field probe
  double rf_phase,        ///< Beam phase relative to the RF in radians
  double phase_rev,       ///< Phase shift between Cavity cells and reverse ADC
  double phase_probe,     ///< Phase shift between Cavity cells and probe ADC
  double Tstep,           ///< Simulation time-step in seconds
  double *mech_couplings, ///< Pointer to array of mechanical couplings in (rad/s)/V^2
  int n_mech              ///< Number of mechanical eigenmodes
  )
{
  ElecMode * elecMode = calloc(1,sizeof(ElecMode));

  // Allocate single-pole, n_mode Filter
  Filter_Allocate_In(&elecMode->fil,1,1);
  ElecMode_Allocate_In(elecMode, RoverQ, foffset, LO_w0, Q_0, Q_drive, Q_probe, rf_phase,  phase_rev, phase_probe, Tstep ,mech_couplings, n_mech);

  return elecMode;
}
/** Frees memory of a Electrical mode struct (representing an electrical eigenmode).*/
void ElecMode_Deallocate(ElecMode * elecMode)
{
  free(elecMode->A);
  free(elecMode->C);
  Filter_Deallocate(&elecMode->fil);
  free(elecMode);
}

/** Takes a previously configured Electrical mode and allocates its State struct accordingly. */
void ElecMode_State_Allocate(ElecMode_State *elecMode_state, ElecMode *elecMode)
{
  elecMode_state -> delta_omega = 0.0;
  Filter_State_Allocate(&elecMode_state->fil_state, &elecMode->fil);
}

/** Frees memory of Electrical Eigenmode state struct. */
void ElecMode_State_Deallocate(ElecMode_State *elecMode_state)
{
  Filter_State_Deallocate(&elecMode_state->fil_state);
  free(elecMode_state);
}

/** Helper routine to get a reference to a given Electrical mode given the Cavity struct. */
ElecMode_State *ElecMode_State_Get(Cavity_State *cav_state, int idx)
{
  return cav_state->elecMode_state_net[idx];
}

/** Takes a previously configured Cavity and allocates its State struct accordingly.
  * It allocates states for the Cavity's Electrical modes recursively. */
void Cavity_State_Allocate(Cavity_State *cav_state, Cavity *cav)
{

  cav_state -> E_probe = (double complex) 0.0;
  cav_state -> E_reverse = (double complex) 0.0;
  cav_state -> Kg = (double complex) 0.0;
  cav_state -> elecMode_state_net = (ElecMode_State**)calloc(cav->n_modes,sizeof(ElecMode_State*));

  int i;
  for(i=0;i<cav->n_modes;i++) {
      cav_state -> elecMode_state_net[i] = (ElecMode_State*)calloc(1,sizeof(ElecMode_State));
      ElecMode_State_Allocate(cav_state -> elecMode_state_net[i], cav->elecMode_net[i]);
  }
}

/** Frees memory of Cavity state struct. */
void Cavity_State_Deallocate(Cavity_State *cav_state, Cavity *cav)
{
  for(int i=0;i<cav->n_modes;i++) {
      ElecMode_State_Deallocate(cav_state -> elecMode_state_net[i]);
  }
  free(cav_state -> elecMode_state_net);
  free(cav_state);
}

/** Step function for Electrical Eigenmode:
  * Calculates the state for the next simulation step.
  * Returns mode's accelerating voltage and stores current state in State struct. */
double complex ElecMode_Step(
  ElecMode *elecMode,             ///< Pointer to ElecMode struct

  // Inputs
  double complex Kg_fwd,          ///< Drive input in sqrt(W)
  double beam_charge,             ///< Beam charge in Coulombs
  double delta_tz,                ///< Timing jitter in seconds (RF reference noise)

  // States
  ElecMode_State *elecMode_state, ///< Pointer to ElecMode State

  // Outputs (V^2 stored in elecMode_state)
  double complex *v_probe,        ///< Field probe signal in Volts (including LLRF noise on that port)
  double complex *v_em            ///< Emitted signal in Volts (including LLRF noise on that port)
  )
{
  // Intermediate signals
  double complex v_beam=0.0, v_drive=0.0, v_in=0.0, v_out=0.0;
  double omega_now=0.0, d_phase_now=0.0;

  // Beam-induced voltage (convert charge to voltage and add timing noise)
  v_beam = beam_charge * elecMode -> k_beam * cexp(-I*elecMode->LO_w0*delta_tz);  // k_beam = Tstep * (R/Q) * Q_L

  // RF drive term
  v_drive = Kg_fwd * elecMode -> k_drive; // Drive term (k_drive = 2*sqrt(Q_drive*(R/Q))

  // Integrate mode's offset frequency to obtain phase
  // Add baseline frequency offset to perturbation (delta_omega) to obtain total frequency offset
  omega_now = elecMode->omega_d_0 + elecMode_state->delta_omega;
  d_phase_now = elecMode_state-> d_phase + omega_now * elecMode->Tstep;
  elecMode_state-> d_phase = d_phase_now; // Store phase state

  // Calculate mode's driving term (drive + beam)
  // Note the absence of omega_f on this term with respect to the governing equations
  // That term implies the unity gain at DC: normalization takes place in Filter_Step.
  v_in = (v_drive + v_beam)*cexp(-I*d_phase_now);

  // Apply first-order low-pass filter
  v_out = Filter_Step(&(elecMode->fil), v_in, &(elecMode_state->fil_state))*cexp(I*d_phase_now);

  // Calculate outputs based on v_vout
  elecMode_state->V_2 = pow(cabs(v_out), 2.0);  // Voltage squared

  *v_probe = v_out * elecMode -> k_probe; // Probe (k_probe = exp(j*phase_probe)/sqrt(Q_probe*(R/Q)) )

  *v_em = v_out * elecMode -> k_em;  // Emitted (k_em = exp(j*phase_emm)/sqrt(Q_drive*(R/Q)) )

  // Return mode's accelerating voltage
  return v_out;
}

/** Takes a pointer to a Cavity struct which has been previously allocated
  * and fills it in with the values passed as arguments. */
void Cavity_Allocate_In(
  Cavity *cav,                ///< Pointer to ElecMode struct
  ElecMode_dp elecMode_net,   ///< Pointer to an array of Electrical modes (already allocated)
  int n_modes,                ///< Number of electrical eigenmodes
  double L,                   ///< Cavity electrical length in meters
  double nom_grad,            ///< Cavity nominal gradient in V/m
  double rf_phase,            ///< Cavity nominal phase with respect to the beam
  double design_voltage,      ///< Nominal voltage operating point in Volts
  int fund_index              ///< Index of the fundamental mode in the array of Electrical modes
  )
{
  cav -> rf_phase = rf_phase;
  cav -> design_voltage = design_voltage;
  cav -> fund_index = fund_index;
  cav -> L = L;
  cav -> nom_grad = nom_grad;
  cav -> n_modes = n_modes;
  cav-> elecMode_net = elecMode_net;
}
/** Allocates memory for a Cavity struct and fills it in with the values passed as arguments.
  * Returns a pointer to the newly allocated struct. It assumes that the Electrical Eigenmodes have been previously allocated. */
Cavity * Cavity_Allocate_New(
  ElecMode_dp elecMode_net,   ///< Pointer to an array of Electrical modes (already allocated)
  int n_modes,                ///< Number of electrical eigenmodes
  double L,                   ///< Cavity electrical length in meters
  double nom_grad,            ///< Cavity nominal gradient in V/m
  double rf_phase,            ///< Cavity nominal phase with respect to the beam
  double design_voltage,      ///< Nominal voltage operating point in Volts
  int fund_index              ///< Index of the fundamental mode in the array of Electrical modes
  )
{
  Cavity *cav;
  cav = calloc(1,sizeof(Cavity));

  Cavity_Allocate_In(cav, elecMode_net, n_modes, L,
    nom_grad, rf_phase, design_voltage,
    fund_index);

  return cav;
}

/** Frees memory of a Cavity struct and recursively frees its Electrical Eigenmodes.*/
void Cavity_Deallocate(Cavity *cav)
{
  for(int i=0;i<cav->n_modes;i++) {
      ElecMode_Deallocate(cav->elecMode_net[i]);
  }
  free(cav);
}

/** Step function for Cavity model:
  * Calculates the state for the next simulation step.
  * Returns the overall cavity accelerating voltage (as seen by the beam) and stores current state in State struct. */
double complex Cavity_Step(
  Cavity *cav,                ///< Pointer to Cavity struct
  double delta_tz,            ///< Timing jitter in seconds (RF reference noise)
  double complex Kg,          ///< Drive input in sqrt(W)
  double complex beam_charge, ///< Beam charge in Coulombs (complex since it includes RF phase information)
  Cavity_State *cav_state     ///< Pointer to the Cavity State
  )
{

  // Intermediate signals
  double complex v_out=0.0;
  double complex v_probe_now=0.0, v_probe_sum=0.0;
  double complex v_em_now=0.0, v_em_sum=0.0;

  // Propagate high-power drive signal to cavity coupler through waveguide
  double complex Kg_fwd = Kg; // Instantaneous propagation through perfect waveguide for now
  cav_state->Kg = Kg; // Instantaneous propagation through perfect waveguide for now

  int i;
  // Iterate over Electrical Modes and add up
  // mode contributions to probe and reflected signals
  for(i=0;i<cav->n_modes;i++){

    // Sum of mode's accelerating voltages (Seen by the beam, no port couplings)
    v_out += ElecMode_Step(cav->elecMode_net[i], Kg_fwd, beam_charge, delta_tz, cav_state->elecMode_state_net[i], &v_probe_now, &v_em_now);

    // Sum of cavity probe signals (including probe coupling and phase shift between cavity port and probe ADC)
    v_probe_sum += v_probe_now;

    // Sum of emitted voltages (including emitted coupling and phase shift between cavity port and reverse ADC)
    v_em_sum += v_em_now;

  } // End of Electrical mode iteration

  // Re-apply propagation through waveguide between cavity port and directional coupler on the reverse path
  double complex Kg_rfl = Kg; // Instantaneous propagation through perfect waveguide for now

  // Propagate new values into Cavity state
  cav_state -> E_probe = v_probe_sum;
  cav_state -> E_reverse = v_em_sum-Kg_rfl;
  cav_state -> V = v_out;

  // Return overall accelerating voltage (as seen by the beam)
  return v_out;
}

/** Helper routine to zero out Cavity state. Useful to restore initial state in unit tests. */
void Cavity_Clear(Cavity *cav, Cavity_State *cav_state)
{
  // Zero out signals
  cav_state -> E_probe = 0.0;
  cav_state -> E_reverse = 0.0;
  cav_state -> V = 0.0;

  // Clear filter states for each electrical mode
  for(int i=0;i<cav->n_modes;i++){
    Filter_State_Clear(&cav-> elecMode_net[i]->fil, &cav_state->elecMode_state_net[i]->fil_state);
  }
}
