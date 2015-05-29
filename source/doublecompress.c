/**
 * @file doublecompress.c
 * @brief Longitudinal Phase-Space Beam Dynamics Model (doublecompress).
 * Allocation, configuration and beam dynamics calculcation functions, which iterates over an arbitrary number of Linac Sections.
 * Logic initially wrote by Paul Emma (SLAC) in Matlab, and translated into C by Daniel Driver (UC Berkeley).
 * @author Carlos Serrano (Cserrano@lbl.gov)
 */

#include "doublecompress.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/** Allocates Doublecompress State struct accordingly for Linac of length Nlinac. */
void Doublecompress_State_Allocate(
  Doublecompress_State * dcs,  ///< Pointer to State struct
  int Nlinac                   ///< Number of Linac Sections
  )
{
  dcs->Ipk =    (double*)calloc(Nlinac,sizeof(double));
  dcs->sz =     (double*)calloc(Nlinac,sizeof(double));
  dcs->dE_E =   (double*)calloc(Nlinac,sizeof(double));
  dcs->sd =     (double*)calloc(Nlinac,sizeof(double));
  dcs->dt =     (double*)calloc(Nlinac,sizeof(double));
  dcs->sdsgn =  (double*)calloc(Nlinac,sizeof(double));
  dcs->k =      (double*)calloc(Nlinac,sizeof(double));
  dcs->Eloss =  (double*)calloc(Nlinac,sizeof(double));
  dcs->dE_Ei =  (double*)calloc(Nlinac,sizeof(double));
  dcs->dE_Ei2 = (double*)calloc(Nlinac,sizeof(double));
  dcs->cor =    (double*)calloc(Nlinac,sizeof(double));
}

/** Frees memory of Doublecompress State struct. */
void Doublecompress_State_Deallocate(Doublecompress_State * dcs)
{
  free(dcs->Ipk);
  free(dcs->sz);
  free(dcs->dE_E);
  free(dcs->sd);
  free(dcs->dt);
  free(dcs->sdsgn);
  free(dcs->k);
  free(dcs->Eloss);
  free(dcs->dE_Ei);
  free(dcs->dE_Ei2);
  free(dcs->cor);
}

void Doublecompress_State_Attach(Doublecompress_State * dcs, int Nlinac, double * payload)
{
  dcs->Ipk = payload;
  dcs->sz = payload + Nlinac;
  dcs->dE_E = payload + 2*Nlinac;
  dcs->sd = payload+ 3*Nlinac;
  dcs->dt = payload+ 4*Nlinac;
  dcs->sdsgn = payload + 5*Nlinac;
  dcs->Eloss = payload + 6*Nlinac;
  dcs->dE_Ei = payload + 7*Nlinac;
  dcs->dE_Ei2 = payload + 8*Nlinac;
  dcs->cor = payload + 9*Nlinac;

  // Payload better have been of size 10*Nlinac!
}

#ifndef c
#define c 299792458.0       ///< Speed of light in m/s
#endif
#ifndef mu0
#define mu0 4*M_PI*1e-7     ///< Vacuum permeability in N/A^2
#endif
#ifndef Z0
#define Z0 120.0*M_PI       ///< Impedance of free space in Ohms
#endif
#ifndef e
#define e 1.60217656535E-19 ///< Charge of the electron in Coulombs
#endif
/** Square of x. */
#define SQ(x) ((x)*(x))

/** Function to calculate bunch length and energy spread after a
  * compressor system (just like LCLS) where the E-z correlations are generated
  * by Linacs at off crest RF phases (not zero crossing).  The wakefield is
  * included in linear form assuming rectangular z-distributions.  The RMS
  * bunch lengths and energy spreads calculated here are linear in that they
  * do not directly include the T566, or RF curvature non-linearities along the
  * bunch.  The calculation does, however, include the E-z correlation dependence
  * on incoming timing jitter (due to RF curvature) and charge jitter, and the
  * T566 effect on R56 for mean-off-energy beams.  The bunch head is at z<0
  * (same as LiTrack), so a chicane compressor has both R56 < 0 and phi < 0.
  * A reference for most equations used in this functions can be found in the
  * Longitudinal beam dynamics section of the Physics document, in the doc/ directory.
  * All units are in SI units except for Energy, which is expressed in eV.
  -# Inputs:
    -# gun.sz0: Nominal initial RMS bunch length [m]
    -# gun.sd0:  Nominal initial incoh. energy spread at Eg [fraction]
    -# gun.E:  Nominal gun exit energy [eV]
    -# gun.N:        Bunch population [e.g. 6.25E9]

    -# linac_array[j].dE: Energy difference of bunch from start to finish of Linac j [eV]
    -# linac_array[j].R56:  Nominal R56 values for Linac j (chicane-R56 < 0) [m]
    -# linac_array[j].T566: Nominal T566 values for Linac j (always > 0) [m]
    -# linac_array[j].phi:  Nominal Linac RF phase for Linac j (-30 deg accelerates and puts head energy lower than tail) [radians]
    -# linac_array[j].L:  Linac length (scales wake) [m]
    -# linac_array[j].lam:  RF wavelength for each Linac (Sband=0.105m, Xband=0.02625m) [m]
    -# linac_array[j].s0: Wakefield characteristic length (Sband=1.322mm), Xband=0.77mm) [m]
    -# linac_array[j].a:  Mean iris radius (Sband=11.654mm,Xband=4.72mm) [m]

    -# dN_Nf:  Relative bunch population error (fraction) (e.g. -2 => 2-low in bunch charge)
    -# dtg:    Timing error of gun wrt RF (<0 is an early bunch) [sec]
    -# dEg:    Energy deviation at end of injector [eV]
    -# dsig_z: Deviation of bunch length from nominal length [m]
    -# dsig_E: Deviation of energy spread from nominal energy
    -# spread [fraction of nominal energy]
    -# chirp:  <Ez> correlation [m]
    -# dphivr: Vector of 5 Linac RF phase errors (<0 is early bunch arrival) [rad] {Nlinacs}
    -# dV_Vvr: Vector of 5 Linac RF relative voltage errors [fraction] {Nlinacs}

  -# Outputs:

    -# Ipk:  Peak current at end of j-th Linac [A] {Nlinacs}
    -# sz: rms bunch length after j-th Linac [m] {Nlinacs}
    -# dE_E: Relative energy offset after j-th Linac [fraction] {Nlinacs}
    -# sd: rms rel. energy spread after j-th Linac [fraction] {Nlinacs}
    -# dt: Timing error at end of j-th Linac [sec] {Nlinacs}
    -# sdsgn:  Signed, correlated E-sprd per Linac (dE_E-z slope * sigz) [fraction] {Nlinacs}
    -# k:  <Ez> correlation const. of j-th Linac [1/m] {Nlinacs}
    -# Eloss:  Energy loss per Linac due to wake [eV] {Nlinacs}
    -# dE_Ei:  Relative energy offset of JUST j-th Linac [fraction] {Nlinacs}
    -# dE_Ei2:Energy offset error relative to final energy of JUST j-th Linac [fraction] {Nlinacs}
    -# cor:  Signed-correlated energy spread (slope*sigz) [fraction] {Nlinacs+1}
  *
*/
void Doublecompress(
  Gun * gun,                   ///< Electron Gun
  Linac ** linac_array,        ///< Array of Linac Sections
  int Nlinac,                  ///< Number of Linac Sections
	// Inputs which change with time potentially
	Noise_Srcs * noise_srcs,     ///< Noise Sources
  double * dphivr,             ///< Array of RF Amplitude errors (fraction relative to Linac nominal voltage)
  double * dV_Vvr,             ///< Array of RF phase errors in radians
	// double_compress Outputs
	Doublecompress_State * dcs   ///< Pointer to Doublecompress State struct
  )
{

  // Rename for some consistency with the original double compress
  double dN_Nf = noise_srcs->dQ_Q;

  // Values that get updated with the Linac, here initial conditions are set
  double Eprev = gun->E;                          // Energy after leaving gun [ev]
  double szprev = gun->sz0 + noise_srcs->dsig_z;  // RMS bunch length [m]
  double sdprev = gun->sd0 + noise_srcs->dsig_E;  // incoh. energy spread [fraction]
  double dE_Eprev = noise_srcs->dE_ing/gun->E;    // Relative energy error at start
  double dtprev = noise_srcs->dtg;                // Timing error of previous Linac
  double corprev = noise_srcs->dchirp;

  // Declarations just to reduce calculations and typing
  Linac * lin;
  double sqrt12=sqrt(12.0);

  double E;     // Energy at the end of the current Linac
  double Er;    // Energy ratio
  double ds;    // FW bunch length for uniform dist. [m]
  double dphi;  // Total local phase error (gun, prev-R56 + local-RF) [rad]
  double s0;    // Wakefield characteristic length [m]
  double lambar;// lam/2*pi
  double R56;   // R56 value adjust for energy error
  double C;     // Temp for cos(phir)
  double kR561; double sd2; double sz2;
  double Nec;   // Intermediate in calc of kw
  double kw;    // wake's effect on linear correlation factor (<0) [1/m]
  double kn;    // RF phase induced linear correlation factor [1/m]
  double k;     // k for the loop to reduce typing stored in dcs->k[j]
  for(int j=0;j<Nlinac;j++){

    lin = linac_array[j]; // Pull out current Linac from pointer array
    E = (Eprev+lin->dE);  // Calculate the Energy after Linac
    Er = Eprev/E;         // Calculate Energy ratio
    ds = szprev*sqrt(12.0); // Convert from RMS to FWHM

    // Pre-calculate a couple of factors for the calculation of kw
    s0 = lin->s0;
    Nec = 2.0*fabs(gun->Q)*s0*Z0*c/M_PI/SQ(lin->a);

    // Wake's effect on linear correlation factor (<0) [1/m]
    kw = -(1.0+dN_Nf)*(Nec*lin->L/(SQ(ds)*E))*
      (1.0-(1.0+sqrt(ds/s0))*exp(-sqrt(ds/s0)));

    // Pre-calculate a couple of factors for the calculation of kn
    lambar = lin->lam/2.0/M_PI;
    C = cos(lin->phi);
    dphi = dtprev*c/lambar + dphivr[j];

    // RF phase induced linear correlation factor [1/m]
    kn = (Er-1.0)*sin(lin->phi + dphi)/(lambar*C);

    k = kw + kn;
    dcs->k[j] = k;

    // Relative Energy error, but only individual
    // Linac contribution (not original,
    // changed by Stefan Paret and introduced dV_Vvr)
    dcs->dE_Ei[j] = (1.0-Er)*((1.0+dV_Vvr[j])*cos(lin->phi + dphi)/C - 1.0);

    // Relative energy error due to dphase and dN error [ ]
    // (changed to use dcs->dE_Ei by Daniel Driver)
    dcs->dE_E[j] = dE_Eprev*Er + dcs->dE_Ei[j] + kw*dN_Nf/(1.0+dN_Nf)*ds/2.0;

    // Relative energy error, but only individual Linac, relative to
    // final E (not original) (Note see below for division by final energy)
    // (changed to use dcs->dE_Ei by Daniel Driver)
    dcs->dE_Ei2[j] = dcs->dE_Ei[j]*E;

    // R56 value changed by T566*dE/E [m];
    // (changed by S. Paret, multiplication by 1 instead of 2)
    R56 = lin->R56 + 1.0*dcs->dE_E[j]*lin->T566;

    // Approximate energy loss due to wake (>0) [GeV]
    dcs->Eloss[j] = -E*kw*ds/2.0;
    kR561 = 1.0+k*R56;  // save computation time
    sd2 = SQ(sdprev);   // save computation time
    sz2 = SQ(szprev);   // save computation time

    // RMS bunch length after Linac and R56 #(j-1) [m]
    dcs->sz[j] = sqrt(SQ(kR561)*sz2 + SQ(R56*Er*sdprev) + 2.0*Er*R56*kR561*corprev);

    // RMS Energy spread after Linac and R56 #(j-1)
    dcs->sd[j] = sqrt(SQ(k)*sz2 + SQ(Er)*sd2 + 2.0*Er*k*corprev);

    // save new E-z correlation [m]
    dcs->cor[j] = k*kR561*sz2 + SQ(Er)*R56*sd2 + Er*(1.0+2.0*k*R56)*corprev;

    // Signed-correlated energy spread (slope*sigz) [ ]
    dcs->sdsgn[j] = dcs->cor[j]/dcs->sz[j];

    // Calculate peak current
    dcs->Ipk[j] = (1.0+dN_Nf)*fabs(gun->Q)*c/sqrt12/dcs->sz[j];

    // Timing error
    dcs->dt[j] = dtprev + dcs->dE_E[j]*R56/c;  // timing error after k-th Linac [s]

    // Move current to previous for the next Linac step
    Eprev=E;
    dtprev=dcs->dt[j];
    szprev=dcs->sz[j];
    sdprev=dcs->sd[j];
    dE_Eprev=dcs->dE_E[j];
    corprev=dcs->cor[j];
  }

  /*
  Finish calculating dcs->dE_Ei2, which requires normalization by the final Energy.
  I would like calculation to be all in one line above but we stored dE and added to energy at each
  step to get the final energy is so it is unknown until the end of the loop.
  The E here is different from above. Before it changed with each loop and was the Energy
  of the beam after Linac j. Now it is static and is the final energy after all Linacs
  */

  for(int j=0;j<Nlinac;j++){
    dcs->dE_Ei2[j]=dcs->dE_Ei2[j]/E;

  }

}

void Doublecompress_Octave_Benchmark(Gun * gun, Linac ** linac_array, int Nlinac,
			// Time-varying inputs
			Noise_Srcs * noise_srcs, double * dphivr, double * dV_Vvr,
			// double_compress output states
			Doublecompress_State * dcs
			)
{

  // Rename for some consistency with the original double_compress
  double dN_Nf = noise_srcs->dQ_Q;
   //e used in octave code
  double e_oct = 1.602177E-19;

  // Variables which get updated with the Linac, here initial conditions are set
  double Eprev=gun->E;                        // Energy after leaving gun [ev]
  double szprev=gun->sz0+noise_srcs->dsig_z;  // RMS bunch length [m]
  double sdprev=gun->sd0+noise_srcs->dsig_E;  // incoh. energy spread [fraction]
  double dE_Eprev=noise_srcs->dE_ing/gun->E;  // Relative energy error at start
  double dtprev=noise_srcs->dtg;              // Timing error of previous Linac NOT THE TIME STEP!!!!
  double corprev=noise_srcs->dchirp;

  // Variable definition just to reduce calculations and typing
  Linac * lin;
  double sqrt12=sqrt(12.0);

  double N = fabs(gun->Q)*6.241E18;   // Number of particles from gun
  double E;                           // Energy at the end of the current Linac
  double Er;                          // Energy ratio
  double ds;                          // FW bunch length for uniform dist. [m]
  double dphi;                        // total local phase error (gun, prev-R56 + local-RF) [rad]
  double s0;                          // Wakefield characteristic length [m]
  double lambar;                      // lam/2*pi
  double R56;                         // R56 value adjust for Energy error
  double C;                           // Temp for cos(phir)
  double kR561; double sd2; double sz2;
  double Nec;                         // Intermediate in calc of kw
  double kw;                          // Wake's effect on linear correlation factor (<0) [1/m]
  double kn;                          // RF phase induced linear correlation factor [1/m]
  double k;                           // k for the loop to reduce typing stored in dcs->k[j]
  for(int j=0;j<Nlinac;j++){

    lin=linac_array[j];               // Pull out current Linac from pointer array
    E=(Eprev+lin->dE);                // calculate the energy after Linac
    Er=Eprev/E;                       // Calculate energy ratio
    ds=szprev*sqrt(12.0);
    lambar=lin->lam/2.0/M_PI;
    C = cos(lin->phi);

    dphi= dtprev*c/lambar + dphivr[j];
    s0=lin->s0;

    Nec=2.0*N*e_oct*s0*Z0*c/M_PI/SQ(lin->a);

    // More Physics calculations
    // Wake's effect on linear correlation factor (<0) [1/m]
    kw= -(1.0+dN_Nf)*(Nec*lin->L/(SQ(ds)*E))*
      (1.0-(1.0+sqrt(ds/s0))*exp(-sqrt(ds/s0)));

    // RF phase induced linear correlation factor [1/m]
    kn = (Er-1.0)*sin(lin->phi + dphi)/(lambar*C);
    k = kw + kn;
    dcs->k[j]=k;

    // Relative Energy error, but only individual
    // Linac contribution (not original,
    // changed by Stefan Paret and introduced dV_Vvr)
    dcs->dE_Ei[j]    = (1.0-Er)*((1.0+dV_Vvr[j])*cos(lin->phi + dphi)/C - 1.0);

    // Relative energy error due to dphase and dN error [ ]
    // (changed to use dcs->dE_Ei by Daniel Driver)
    dcs->dE_E[j]= dE_Eprev*Er + dcs->dE_Ei[j]+kw*dN_Nf/(1.0+dN_Nf)*ds/2.0;

    //relative energy error, but only individual Linac, relative to
    //final E (not original) (Note see below for division by final energy)
    //changed to use dcs->dE_Ei by Daniel Driver
    dcs->dE_Ei2[j]   = dcs->dE_Ei[j]*E;

    // R56 value changed by T566*dE/E [m];
    // (changed by S. Paret, multiplication by 1 instead of 2)
    R56 = lin->R56 + 1.0*dcs->dE_E[j]*lin->T566;

    // Approximate energy loss due to wake (>0) [GeV]
    dcs->Eloss[j]    = -E*kw*ds/2.0;
    kR561       = 1.0+k*R56;    // save computation time
    sd2         = SQ(sdprev);   // save computation time
    sz2         = SQ(szprev);   // save computation time

    // RMS bunch length after Linac and R56 #(j-1) [m]
    dcs->sz[j] = sqrt(SQ(kR561)*sz2 + SQ(R56*Er*sdprev) +
		 2.0*Er*R56*kR561*corprev);

    // RMS energy spread after Linac and R56 #(j-1) [
    dcs->sd[j]       = sqrt(SQ(k)*sz2 + SQ(Er)*sd2 + 2.0*Er*k*corprev);

    // Save new E-z correlation [m]
    dcs->cor[j]   = k*kR561*sz2 + SQ(Er)*R56*sd2 +
      Er*(1.0+2.0*k*R56)*corprev;

    // Signed-correlated energy spread (slope*sigz) [ ]
    dcs->sdsgn[j]    = dcs->cor[j]/dcs->sz[j];

    // Calculate peak current
    dcs->Ipk[j]=(1.0+dN_Nf)*N*e_oct*c/sqrt12/dcs->sz[j];

    // Timing error
    dcs->dt[j]=dtprev + dcs->dE_E[j]*R56/c;  // timing error after Linac k [s]

    // Move current to previous for the next Linac step
    Eprev=E;
    dtprev=dcs->dt[j];
    szprev=dcs->sz[j];
    sdprev=dcs->sd[j];
    dE_Eprev=dcs->dE_E[j];
    corprev=dcs->cor[j];
  }

 /*
  Finish calculating dcs->dE_Ei2, which requires normalization by the final Energy.
  I would like calculation to be all in one line above but we stored dE and added to energy at each
  step to get the final energy is so it is unknown until the end of the loop.
  The E here is different from above. Before it changed with each loop and was the Energy
  of the beam after Linac j. Now it is static and is the final energy after all Linacs
  */

  for(int j=0;j<Nlinac;j++){
    dcs->dE_Ei2[j]=dcs->dE_Ei2[j]/E;

  }

}
