#ifndef FILTER_H
#define FILTER_H
#include "complex.h"

/*
 * Generic Filter Data structure
 *
 * Can handle a multi-order filter, with an arbitrary number 
 * of modes in parallel. Stored in a similator format to 
 * Compressed-Row-Sparse. 
 *
 *
 *
 */

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
double complex Filter_Step(Filter * fil, double complex innow,
			   Filter_State * filnow, Filter_State * filpast);

#endif
