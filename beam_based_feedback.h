#ifndef BEAM_BASED_FEEDBACK_H
#define BEAM_BASED_FEEDBACK_H

/*
 * beam_based_feedback.h
 *
 * Routines for Beam Based Feedback.
 *
 * Stores and applies the pseudoinverse of the Jacobian
 * to the results from doublecompress. The Jacobian and SVD
 * are computed in python, and then the result is packed into
 * the C structure densely.
 * We could use LAPACK here... but whatev's.
 *
 * AFQ LBL Summer 2013
 */

#include "linac_param.h"
#include "doublecompress.h"


/*
 * A structure to store the pseudoinverse and the selectors for
 * BBF.
 */
typedef struct str_bbf_param {
  /*
   * Index selectors for the control and measured sides
   * of length 2*U_control|measured.
   * These arrays index into the in/out variables of Doub.Comp. 
   * as pairs (L,A), where A is a switch on which array (there
   * are various for the ins/outs), and L is an index into 
   * that array. The vectors that multiply and return from a dense
   * Mpinv are of length U_control|measured and are assembled by
   *    V[k] = array[A[k]][L[k]], 
   * where
   *    A[k] = idx[2*k+1], L[k] = idx[2*k+0]
   * for k = [0,U-1]. The selection array[A] is performed via a switch.
   * 
   */
  int U_control;
  int * idx_control;

  int  U_measured;
  int * idx_measured;

  /*
   * The preallocated input and output Vectors, 
   * of length U_control|measured.
   */

  double * V_control;
  double * V_measured;

  /*
   * The pseudoinverse, stored dense of size Um * Uc. 
   * U_measured is the inside dimension, so 
   *    [M]_i,j = Mpinv[Um*i + j]
   * Copied in by scipy.
   */
  double * Mpinv;
} BBF_Param;


/*
 * Do the mallocing/freeing for the data structures.
 * Configuration is done in python.
 */
void BBF_Param_Alloc(BBF_Param * bbf, int Uc, int Um);
void BBF_Param_Free(BBF_Param * bbf);

void BBF_Step(BBF_Param * bbf, Doublecompress_State * dcs,
	      Linac_Param ** linp, int Nlinac);
#endif
