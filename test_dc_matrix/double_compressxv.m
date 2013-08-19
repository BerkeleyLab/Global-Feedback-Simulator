function [Ipk,sz,dE_E,sd,dt,sdsgn,k,Eloss,dE_Ei,dE_Ei2,cor1] = double_compressxv(params,dN_N,dtg,dEg,dsig_z,...
                                                  dsig_E,chirp,dphiv,dV_Vv,Q)
%	
%	Function to calculate bunch length and energy spread after a five-stage
%	compressor system (just like LCLS) where the E-z correlations are generated
%	by five linacs at off crest rf phases (not zero crossing).  The wakefield is
%	included in linear form assuming rectangular z-distributions.  The rms
%	bunch lengths and energy spreads calculated here are linear in that they
%	do not directly include the T566, or rf curvature non-linearities along the
%	bunch.  The calculation does, however, include the E-z correlation dependence
%	on incoming timing jitter (due to rf curvature) and charge jitter, and the
%	T566 effect on R56 for mean-off-energy beams.  The bunch head is at z<0
%	(same as LiTrack), so a chicane compressor has both R56 < 0 and phi < 0.
%       Written by P. Emma
%
%       Modifications:
%       A reference for most equations used below is 'doc/presentations/PE_slides/COMP_OPT1.pdf'
%       Note that there are a few known differences from the above mentioned reference. These
%       differences are alright and have been marked as 'changed by author'.
%       Deviations of the mean energy, bunch length, energy spread and chirp of the injected bunch
%       were added to the input parameters and are taken into account in the computations.
%       The sign of the charge is ignored.
%       Nov. 2012, S. Paret
%
%       Inputs: sz0:	Nominal initial rms bunch length [mm]
%               sd0:	Nominal initial incoh. energy spread at Eg [%]
%               Eg:	Nominal gun exit energy [GeV]
%               Ev:	Vector of 5 nominal post-linac energies [GeV]
%               R56v:	Vector of 5 nominal R56 values (chicane-R56 < 0) [m]
%               T566v:	Vector of 5 nominal T566 values (always > 0) [m]
%               phiv:	Vector of 5 nominal linac RF phases (-30 deg accelerates 
%		          and puts head energy lower than tail) [degrees]
%               Lv:	Vector of 5 linac lengths (scales wake) [m]
%               N:	Bunch population [e.g. 6.25E9]
%               lamv:	RF wavelength for each linac (Sband=0.105m, Xband=0.02625m) [m]
%               s0v:	Wakefield characteristic length (Sband=1.322mm), Xband=0.77mm) [mm]
%               av:	Mean iris radius (Sband=11.654mm,Xband=4.72mm) [mm]
%               dN_N:	Relative bunch population error [%] (e.g. -2 => 2%-low in bunch charge)
%               dtg:    Timing error of gun wrt RF (<0 is an early bunch) [psec]
%               dEg:    Energy deviation at end of injector [GeV]
%               dsig_z: Deviation of bunch length from nominal length [mm]
%               dsig_E: Deviation of energy spread from nominal energy spread [% of nominal energy]
%               chirp:  <Ez> correlation [m]
%               dphiv:	Vector of 5 linac RF phase errors (<0 is early bunch arrival) [deg]
%               dV_Vv: 	Vector of 5 linac RF relative voltage errors [%]
%		Q:	absolute value of nominal bunch charge [Coulombs]

%
%	Outputs:
%               Ipk:	Peak current at end of i-th linac [A]
%		sz:	rms bunch length after i-th linac [mm]
%               dE_E:	Relative energy offset after i-th linac [%]
%               sd:	rms rel. energy spread after i-th linac [%]
%               dt:	Timing error at end of i-th linac [psec]
%               sdsgn:	Signed, correlated E-sprd per linac (dE_E-z slope * sigz) [%]
%               k:	<Ez> correlation const. of i-th linac [1/m]
%               Eloss:	Energy loss per linac due to wake [GeV]



%===============================================================================

% (sz0,sd0,Eg,Ev,R56v,T566v,phiv,Lv,lamv,s0v,av) are provided thru params, the rest are input individually, except N
if (length(dphiv) ~= length(params.phiv))
  error('wrong size to dphiv input');
end
format long
%N      = abs(Q)*6.241E18                              % bunch population [particles/bunch]

%lam	=params.lam;
sz0	=params.sz0;
sd0	=params.sd0;
Eg	=params.Eg;
Ev	=params.Ev;
R56v	=params.R56v;
T566v	=params.T566v;
phiv	=params.phiv;
lamv	=params.lamv;
s0v	=params.s0v;
av	=params.av;
Lv	=params.Lv;

temp=0;
if(temp==1)
  dN_N
  dtg
  dEg
  dsig_z
  dsig_E
  chirp
  dphiv
  dV_Vv
end

Nv        = length(Ev);
%e         =1.60217656535E-19;                         % constants and conversions...
c         = 2.99792458E8;
Z0        = 120*pi;
dN_Nf     = dN_N*1e-2;		                  % percent -> unitless
s0m       = s0v*1E-3;
am        = av*1E-3;
deg_2_rad = pi/180;
lambar    = lamv/2/pi;
sz0m      = (sz0+dsig_z)*1E-3;
sd0f      = (sd0+dsig_E)*1E-2;
sqrt12    = sqrt(12);
Nec       = 2*Q*s0m*Z0*c/pi./am.^2/1E9;         %equation in doc looks like e^2, missing e is in the conversion from ev->joules

%Nec       = 2*abs(Q)*s0m*Z0*c/pi./am.^2/1E9;         %equation in doc looks like e^2, missing e is in the conversion from ev->joules
phivr     = phiv*deg_2_rad;
dphivr    = dphiv*deg_2_rad;
dV_Vvr    = dV_Vv/100;                            % voltage errors in % converted to dimensionless
C         = cos(phivr);
sz1       = [sz0m zeros(1,Nv)];                   % bunch length 6-vector (includes gun value)
sd1       = [sd0f zeros(1,Nv)];                   % energy spread 6-vector (includes gun value)
cor1      = [chirp zeros(1,Nv)];                  % E-z correlation 6-vector
dE_E1     = [dEg/Eg zeros(1,Nv)];                 % 6-vector for energy deviation
dt1       = [dtg*1E-12 zeros(1,Nv)];              % 6-vector for timing errors [sec]
Ev1       = [Eg Ev(:)'];                          % 6-vector nominal of energies (so all user inputs are 5-vectors)

sz        = zeros(1,Nv);
sd        = zeros(1,Nv);
sdsgn     = zeros(1,Nv);
k         = zeros(1,Nv);
Ipk       = zeros(1,Nv);
dE_E      = zeros(1,Nv);
dE_Ei     = zeros(1,Nv);
dE_Ei2    = zeros(1,Nv);
dt        = zeros(1,Nv);
Eloss     = zeros(1,Nv);
Er        = zeros(1,Nv);

for j = 1:Nv                                      % loop over all linac-compressor segents
  Er(j)       = Ev1(j)/Ev1(j+1);                 % energy ratio [ ]
  ds          = sz1(j)*sqrt12;                    % FW bunch length for uniform dist. [m]
  dphi        = dt1(j)*c/lambar(j) + dphivr(j);   % total local phase error (gun, prev-R56 + local-rf) [rad]
  % wake's effect on linear correlation factor (<0) [1/m]
  kw          = -(1+dN_Nf)*(Nec(j)*Lv(j)/(ds^2*Ev(j)))*(1-(1+sqrt(ds/s0m(j)))*exp(-sqrt(ds/s0m(j))));
  % rf phase induced linear correlation factor [1/m]
  kn          = (Er(j)-1)*sin(phivr(j) + dphi)/(lambar(j)*C(j));
  k(j)        = kw + kn;                          % total linear correlation factor [1/m]
  % relative energy error due to dphase and dN error [ ]
  % changed by author (introduced dV_Vvr)
  dE_E(j)     = dE_E1(j)*Er(j) + (1-Er(j))*((1+dV_Vvr(j))*cos(phivr(j) + dphi)/C(j) - 1) +...
	 kw*dN_Nf/(1+dN_Nf)*ds/2;
  % relative energy error, but only individual linac contribution (not original)
  dE_Ei(j)    = (1-Er(j))*((1+dV_Vvr(j))*cos(phivr(j) + dphi)/C(j) - 1);
  % relative energy error, but only individual linac, relative to L3 final E (not original)
  dE_Ei2(j)   = (Ev1(j+1)-Ev1(j))*((1+dV_Vvr(j))*cos(phivr(j) + dphi)/C(j) - 1)/Ev(Nv);

  % R56 value changed by T566*dE/E [m]; changed by author (multiplication by 1 instead of 2)
  R56         = R56v(j) + 1*dE_E(j)*T566v(j);
  Eloss(j)    = -Ev(j)*kw*ds/2;                   % approximate energy loss due to wake (>0) [GeV]
  kR561       = 1+k(j)*R56;                       % save computation time
  sd2         = sd1(j)^2;                         % save computation time
  sz2         = sz1(j)^2;                         % save computation time

  % rms bunch length after linac and R56 #(j-1) [m]
  sz(j)       = sqrt(kR561^2*sz2 + (R56*Er(j)*sd1(j))^2 + 2*Er(j)*R56*kR561*cor1(j));
  % rms energy spread after linac and R56 #(j-1) [ ]
  sd(j)       = sqrt(k(j)^2*sz2 + Er(j)^2*sd2 + 2*Er(j)*k(j)*cor1(j));
  % save new E-z correlation [m]
  cor1(j+1)   = k(j)*kR561*sz2 + Er(j)^2*R56*sd2 + Er(j)*(1+2*k(j)*R56)*cor1(j);
  sdsgn(j)    = cor1(j+1)/sz(j);                  % signed-correlated energy spread (slope*sigz) [ ]
  dE_E1(j+1)  = dE_E(j);                          % save dE_E [ ]
  sz1(j+1)    = sz(j);                            % save new rms bunch length [m]
  sd1(j+1)    = sd(j);                            % save new total rms energy spread [ ]
  dt1(j+1)    = dt1(j) + dE_E(j)*R56/c;           % timing error after R56 [sec]
  dt(j)       = 1E12*dt1(j+1);                    % save timing error to this point [psec]
end

Ipk   = (1+dN_Nf)*Q*c/sqrt12./sz;               % calculate peak current over uniform bunch dist. [A]
%Ipk   = (1+dN_Nf)*abs(Q)*c/sqrt12./sz;               % calculate peak current over uniform bunch dist. [A]
sz    = sz*1E3;                                   % m -> mm
sd    = sd*1E2;                                   % [ ] -> %
sdsgn = sdsgn*1E2;                                % [ ] -> % (dE/E-z slope times sigz)
dE_E  = dE_E*1E2;                                 % [ ] -> %
dE_Ei = dE_Ei*1E2;                                % [ ] -> % (not original)
dE_Ei2 = dE_Ei2*1E2;                              % [ ] -> % (not original)