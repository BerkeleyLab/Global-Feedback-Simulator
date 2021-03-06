/**
  * @file filter.h
  * @brief Header file for filter.c. Generic filter.
  * Can handle a multi-order filter, with an arbitrary number
  * of modes in parallel. Stored in a simulator format to
  * Compressed-Row-Sparse.
  * @author Carlos Serrano (CSerrano@lbl.gov)
*/

#ifndef FILTER_H
#define FILTER_H
#include "complex.h"

typedef struct str_filter {
  int alloc_order;
  int alloc_coeffs;
  int n_coeffs;
  int order;
  int *modes;
  int *coeff_start;
  double complex * coeffs;
  double complex * poles;
} Filter;

typedef struct str_filter_state {
  double complex * state; // of length n_coeffs
  double complex * input; // of length order
} Filter_State;

Filter * Filter_Allocate_New(int alloc_order, int alloc_coeffs);
void Filter_Allocate_In(Filter * fil, int alloc_order, int alloc_coeffs);
void Filter_Deallocate(Filter * fil);

void Filter_Append_Modes(Filter * fil, double complex * poles,int ord,double dt);

void Filter_State_Allocate(Filter_State * sf, Filter * fil);
void Filter_State_Deallocate(Filter_State * sf);
void Filter_State_Clear(Filter * fil, Filter_State * sf);
double complex Filter_Step(Filter * fil, double complex innow, Filter_State * fil_state);
void Filter_Set_State(Filter * fil, Filter_State * fil_state, double complex state);

#endif
