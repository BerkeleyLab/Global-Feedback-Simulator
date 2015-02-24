
#include "cavity.h"
#include <stdio.h>
#include <string.h>

#define CPRINTs(c) {if(cimag(c)<0) printf("%10.16e%10.16ej ",creal(c),cimag(c)); else printf("%10.16e+%10.16ej ",creal(c),cimag(c)); }

ElecMode_dp ElecMode_Allocate_Array(int n)
{
  ElecMode_dp elecMode_net = calloc(n, sizeof(ElecMode *));
  
  return elecMode_net;
}

void ElecMode_Append(ElecMode** elecMode_arr, ElecMode* elecMode, int index)
{
  // XXX Add some check!!
  elecMode_arr[index] = elecMode;
}


void ElecMode_Allocate_In(ElecMode *elecMode,
  double RoverQ, double foffset, double omega_0_mode,
  double Q_0, double Q_drive, double Q_probe,
  double beam_phase,  double phase_rev, double phase_probe,
  double Tstep
  // ,double *mech_couplings
  )
{

  // Store time step size in order to integrate detune frequency and obtain phase
  elecMode -> Tstep = Tstep;

  // Loaded Q
  double Q_L = 1/(1/Q_0 + 1/Q_drive + 1/Q_probe);

  // Beam impedance/derivator (conversion factor between charge and beam induced voltage)
  // Includes shift to take beam phase relative to the RF into account
  elecMode -> k_beam = RoverQ*Q_L*cexp(-I*beam_phase)/Tstep*1e-12;
  
  // Drive port imperdance
  elecMode -> k_drive = 2*sqrt(Q_drive*RoverQ);

  // Probe port impedance (includes phase shift between cavity cell and probe ADC)
  elecMode -> k_probe = cexp(I*phase_probe)/sqrt(Q_probe*RoverQ);

  // Emitted port impedance (includes phase shift between cavity cell and reverse ADC)
  elecMode -> k_em = cexp(I*phase_rev)/sqrt(Q_drive*RoverQ);

  // Mode's angular frequency (accelerator angular frequency + offset)
  elecMode -> omega_0_mode = omega_0_mode;

  // Mode's open-loop bandwidth 
  elecMode -> omega_f = omega_0_mode/(2*Q_L);
  
  // Mode's single-pole, low-pass filter allocation
  // Calculate pole (mode's bandwidth)
  double complex mode_p = -elecMode -> omega_f - I*2*M_PI*foffset;

  // Append mode to cavity filter
  Filter_Append_Modes(&elecMode->fil, &mode_p, 1, Tstep);

  // Mechanical couplings
  // elecMode -> mech_couplings = mech_couplings;
}

ElecMode * ElecMode_Allocate_New(
  double RoverQ, double foffset, double omega_0_mode,
  double Q_0, double Q_drive, double Q_probe,
  double beam_phase,  double phase_rev, double phase_probe,
  double Tstep
  // ,double *mech_couplings
  )
{
  ElecMode * elecMode = calloc(1,sizeof(ElecMode));

  // Allocate single-pole, n_mode Filter
  Filter_Allocate_In(&elecMode->fil,1,1);
  ElecMode_Allocate_In(elecMode, RoverQ, foffset, omega_0_mode, Q_0, Q_drive, Q_probe, beam_phase,  phase_rev, phase_probe, Tstep
    // ,mech_couplings
    );

  return elecMode;
}

void ElecMode_Deallocate(ElecMode * elecMode)
{
  free(elecMode->mech_couplings);
  Filter_Deallocate(&elecMode->fil);
  free(elecMode);
}


void ElecMode_State_Allocate(ElecMode_State *elecMode_state, ElecMode *elecMode)
{
  elecMode_state -> delta_omega = 0.0;
  Filter_State_Allocate(&elecMode_state->fil_state, &elecMode->fil);
}

void ElecMode_State_Deallocate(ElecMode_State *elecMode_state, ElecMode * elecMode)
{
  Filter_State_Deallocate(&elecMode_state->fil_state);
  free(elecMode_state);
}

void Cavity_State_Allocate(Cavity_State *cav_state, Cavity *cav)
{

  cav_state -> E_probe = (double complex) 0.0;
  cav_state -> E_reverse = (double complex) 0.0;
  cav_state -> elecMode_state_net = (ElecMode_State**)calloc(cav->n_modes,sizeof(ElecMode_State*));

  int i;
  for(i=0;i<cav->n_modes;i++) {
      cav_state -> elecMode_state_net[i] = (ElecMode_State*)calloc(1,sizeof(ElecMode_State));
      ElecMode_State_Allocate(cav_state -> elecMode_state_net[i], cav->elecMode_net[i]);
  }
}

void Cavity_State_Deallocate(Cavity_State *cav_state, Cavity *cav)
{
  for(int i=0;i<cav->n_modes;i++) {
      ElecMode_State_Deallocate(cav_state -> elecMode_state_net[i], cav->elecMode_net[i]);
  }
  free(cav_state -> elecMode_state_net);
  free(cav_state);
}


double complex ElecMode_Step(ElecMode *elecMode,
  // Inputs
  double complex Kg_fwd, double beam_charge, double delta_tz,
  // States
  ElecMode_State *elecMode_state,
  // Outputs (V^2 stored in elecMode_state)
  double complex *v_probe, double complex *v_em)
{
  // Intermediate signals
  double complex v_beam, v_drive, v_in, v_out;
  double omega_now, d_phase_now;

  // Instantaneous mode's angular frequency (including time-varying detuning)
  omega_now = elecMode -> omega_0_mode + elecMode_state->delta_omega;

  // Beam-induced voltage (convert charge to voltage and add timing noise)
  v_beam = beam_charge * elecMode -> k_beam * cexp(-I*omega_now*delta_tz);  // k_beam = Tstep * (R/Q) * Q_L

  // RF drive term
  v_drive = Kg_fwd * elecMode -> k_drive; // Drive term (k_drive = 2*sqrt(Q_drive*(R/Q))

  // Integrate detuning angular frequency to obtain phase
  d_phase_now = elecMode_state-> d_phase + elecMode_state->delta_omega * elecMode->Tstep;
  elecMode_state-> d_phase = d_phase_now; // Store phase state

  // Calculate mode's driving term (drive + beam)
  // Note the absence of omega_f on this term with respect to the governing equations
  // That term implies the unity gain at DC: normalization takes place in Filter_Step.
  v_in = (v_drive + v_beam)*cexp(-I*d_phase_now);

  // Apply first-order low-pass filter
  v_out = Filter_Step(&(elecMode->fil), v_in, &(elecMode_state->fil_state))*cexp(I*d_phase_now);

  // Calculate outputs based on v_vout
  elecMode_state->v_2 = pow(cabs(v_out), 2);  // Voltage squared

  *v_probe = v_out * elecMode -> k_probe; // Probe (k_probe = exp(j*phase_probe)/sqrt(Q_probe*(R/Q)) )

  *v_em = v_out * elecMode -> k_em;  // Emitted (k_em = exp(j*phase_emm)/sqrt(Q_drive*(R/Q)) )

  // Return mode's accelerating voltage
  return v_out;
}

void Cavity_Allocate_In(Cavity *cav, 
  ElecMode_dp elecMode_net, int n_modes,
  double L, double nom_grad,
  // XXX
  double nom_beam_phase, double rf_phase, double design_voltage, double unity_voltage,
  double open_loop_bw)
  // XXX
{
  // XXX Check if used in the future: Inherited properties
  cav -> nom_beam_phase = nom_beam_phase;
  cav -> rf_phase = rf_phase;
  cav -> design_voltage = design_voltage;
  cav -> unity_voltage = unity_voltage;
  cav -> open_loop_bw = open_loop_bw;
  /// XXX

  cav -> L = L;
  cav -> nom_grad = nom_grad;
  cav -> n_modes = n_modes;
  cav-> elecMode_net = elecMode_net;
}

Cavity * Cavity_Allocate_New(ElecMode_dp elecMode_net, int n_modes, 
  double L, double nom_grad,
   // XXX
  double nom_beam_phase, double rf_phase, double design_voltage, double unity_voltage,
  double open_loop_bw)
  // XXX
{
  Cavity *cav;
  cav = calloc(1,sizeof(Cavity));

  Cavity_Allocate_In(cav, elecMode_net, n_modes, L, nom_grad, 
    nom_beam_phase, rf_phase, design_voltage, unity_voltage, 
    open_loop_bw);

  return cav;
}

void Cavity_Deallocate(Cavity *cav)
{
  for(int i=0;i<cav->n_modes;i++) {
      ElecMode_Deallocate(cav->elecMode_net[i]);
  }
  free(cav);
}

double complex Cavity_Step(Cavity *cav, double delta_tz,
     double complex Kg, double complex beam_charge,
     Cavity_State *cav_state)
{

  // Intermediate signals
  double complex v_out=0.0;
  double complex v_probe_now=0.0, v_probe_sum=0.0;
  double complex v_em_now=0.0, v_em_sum=0.0;

  // Propagate high-power drive signal to cavity coupler through waveguide
  double complex Kg_fwd = Kg; // Instantaneous propagation through perfect waveguide for now

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

  // Return overall accelerating voltage (as seen by the beam)
  return v_out;
}

double complex ElecMode_Get(Cavity *cav, int idx, const char *k_type)
{
  // Check index bound
   if(idx > cav->n_modes-1){
      fprintf(stderr, "ElecMode_Get: Requested ElecMode index is too large\n");
      exit(-1);
  }
  
  if (!strcmp(k_type, "k_drive")){ 
    printf("ElecMode_Get: k_drive = %f\n", cabs(cav->elecMode_net[idx]->k_drive));
    return cav->elecMode_net[idx]->k_drive;
  }
  else if(!strcmp(k_type, "k_beam")){
    printf("ElecMode_Get: k_beam = %f\n", cabs(cav->elecMode_net[idx]->k_beam));
    return cav->elecMode_net[idx]->k_beam;
  }
  else if(!strcmp(k_type, "k_probe")) return cav->elecMode_net[idx]->k_probe;
  else if(!strcmp(k_type, "k_em")) return cav->elecMode_net[idx]->k_em;
  else if(!strcmp(k_type, "omega_0_mode")) return cav->elecMode_net[idx]->omega_0_mode;
  else if(!strcmp(k_type, "omega_f")) return cav->elecMode_net[idx]->omega_f;
  else{
    printf("ElecMode_Get: Invalid k_type"); 
    return -1;
  }

}

void ElecMode_Set(Cavity *cav, int idx, const char *k_type, double complex new_value)
{
  // Check index bound
   if(idx > cav->n_modes-1){
      fprintf(stderr, "ElecMode_Set: Requested ElecMode index is too large\n");
      exit(-1);
  }
  
  if (!strcmp(k_type, "k_drive")) cav->elecMode_net[idx]->k_drive = new_value;
  else if(!strcmp(k_type, "k_beam")) cav->elecMode_net[idx]->k_beam = new_value;
  else{
    printf("ElecMode_Set: Invalid k_type"); 
  }

}

void Cavity_Clear(Cavity *cav, Cavity_State *cav_state)
{
  // Zero out signals
  cav_state -> E_probe = 0.0;
  cav_state -> E_reverse = 0.0;

  // Clear filter states for each electrical mode
  for(int i=0;i<cav->n_modes;i++){
    Filter_State_Clear(&cav-> elecMode_net[i]->fil, &cav_state->elecMode_state_net[i]->fil_state);
  }
}