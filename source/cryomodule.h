#ifndef CRYOMODULE_H
#define CRYOMODULE_H

/*
 *  cryomudule.h
 *
 * Routines for time-stepping Cryomodule model
 *
 */

#include <complex.h>
#include "rf_station.h"

 typedef struct str_MechMode {
	double f0, Q;
	Filter fil;
	double a_nu, b_nu, c_nu;
	double Tstep;

} MechMode;

typedef MechMode* MechMode_p;
typedef MechMode** MechMode_dp;

typedef struct str_MechMode_State {
	Filter_State fil_state;
	double x_nu;

} MechMode_State;

typedef MechMode_State* MechMode_State_p;
typedef MechMode_State** MechMode_State_dp;

typedef struct str_Cryomodule {
	int n_rf_stations, n_mechModes;
	RF_Station **rf_station_net;
	MechMode **mechMode_net;

} Cryomodule;

typedef Cryomodule* Cryomodule_p;
typedef Cryomodule** Cryomodule_dp;

typedef struct str_Cryomodule_State {
	RF_State **rf_state_net;
	MechMode_State **mechMode_state_net;
	double *F_nu;

} Cryomodule_State;

MechMode_State_dp MechMode_State_Allocate_Array(int n);
MechMode_dp MechMode_Allocate_Array(int n);
void MechMode_Append(MechMode** mechMode_arr, MechMode* mechMode, int index);

void MechMode_Allocate_In(MechMode *mechMode, double f_nu, double Q_nu, double k_nu, double Tstep);
MechMode *MechMode_Allocate_New(double f_nu, double Q_nu, double k_nu, double Tstep);
void MechMode_Deallocate(MechMode *mechMode);
void MechMode_State_Allocate(MechMode_State *mechMode_State, MechMode *mechMode);
void MechMode_State_Deallocate(MechMode_State *mechMode_State);
void MechMode_Step(MechMode *mechMode, MechMode_State *mechMode_State, double complex F_nu);

RF_Station *Get_RF_Station(Cryomodule *cryo, int index);
RF_State *Get_RF_State(Cryomodule_State *cryo_state, int index);
MechMode_State *Get_MechMode_State(Cryomodule_State *cryo_state, int index);

Cryomodule_dp Cryomodule_Allocate_Array(int n);
void Cryomodule_Append(Cryomodule** cryo_arr, Cryomodule* cryo, int index);
void Cryomodule_Allocate_In(Cryomodule *cryo, RF_Station_dp rf_station_net, int n_rf_stations, MechMode_dp mechMode_net, int n_mechModes);
Cryomodule * Cryomodule_Allocate_New(RF_Station **rf_station_net, int n_rf_stations, MechMode **mechMode_net, int n_mechModes);
void Cryomodule_Deallocate(Cryomodule* cryo);
void Cryomodule_State_Allocate(Cryomodule_State *cryo_state, Cryomodule *cryo);
void Cryomodule_State_Deallocate(Cryomodule_State *cryo_state, Cryomodule *cryo);
void Cryomodule_Step(Cryomodule *cryo, Cryomodule_State * cryo_state,
	double delta_tz, double complex beam_charge, double complex probe_ns, double complex rev_ns,
	int openloop);

#endif