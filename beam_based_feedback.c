#include "beam_based_feedback.h"

#include <stdlib.h>
#include <stdio.h>

void BBF_Param_Alloc(BBF_Param * bbf, int Uc, int Um) {
  bbf->U_control = Uc;
  bbf->U_measured = Um;

  bbf->idx_control = (int*)calloc(2*Uc, sizeof(int));
  bbf->idx_measured = (int*)calloc(2*Um, sizeof(int));

  bbf->V_control = (double*)calloc(Uc,sizeof(double));
  bbf->V_measured = (double*)calloc(Um,sizeof(double));

  bbf->Mpinv = (double*)calloc(Uc*Um,sizeof(double));
}

void BBF_Param_Free(BBF_Param * bbf) {
  free(bbf->idx_control);
  free(bbf->idx_measured);
  free(bbf->V_control);
  free(bbf->V_measured);
  free(bbf->Mpinv);
  // Just for good measure.
  bbf->U_control = 0;
  bbf->U_measured = 0;
}

void BBF_pack_meas(BBF_Param * bbf, Doublecompress_State * dcs) {
  int k, A, L;
  for(k=0;k<bbf->U_measured;k++) {
    L = bbf->idx_measured[2*k+0];
    A = bbf->idx_measured[2*k+1];
    switch(A)
      {
      case 0:
	bbf->V_measured[k] = dcs->dE_E[L];
	break;
      case 1:
	bbf->V_measured[k] = dcs->sz[L];
	break;
      case 2:
	bbf->V_measured[k] = dcs->dt[L];
	break;
      case 3:
	bbf->V_measured[k] = dcs->sd[L];
	break;
      default:
	printf("Invalid BBF key %d...\n",A);
	bbf->V_measured[k] = 0.0;
      }
  }
}

void BBF_Step(BBF_Param * bbf, Doublecompress_State * dcs,
	      Linac_Param ** linp_arr, int Nlinac)
{
  int i,j,Um,Uc;
  Um = bbf->U_measured;
  Uc = bbf->U_control;
  
  /*
   * Pack the new state in the vector
   */
  BBF_pack_meas(bbf,dcs);

  /*
   * Multiply by the pseudoinverse
   */
  for(i=0;i<Uc;i++) {
    //bbf->V_control[i]=0.0;
    for(j=0;j<Um;j++) {
      bbf->V_control[i]+=bbf->Mpinv[i*Um+j]*bbf->V_measured[j];
    }
  }

  /*
   * Now adjust the setpoints on the Linacs.
   * One pass assuming that the list is sorted by linacs.
   */
  double complex sp_modulator = 1.0;
  int L;
  int Llast = bbf->idx_control[0];
  int A = bbf->idx_control[1];
  if(A==0) { //V
    sp_modulator += bbf->V_control[0];
  } else { //Phi
    sp_modulator += cexp(I*bbf->V_control[0]) - 1.0;
  }

  for(i=1; i<bbf->U_control;i++) {
    L = bbf->idx_control[2*i+0];
    A = bbf->idx_control[2*i+1];
    if(L!=Llast) {
      linp_arr[Llast]->fpga.set_point *= sp_modulator;
      sp_modulator = 1.0;
      Llast = L;
    }
    if(A==0) { //V
      sp_modulator += bbf->V_control[i];
    } else { //Phi
      sp_modulator += cexp(I*bbf->V_control[i]) - 1.0;
    }
  }

  
}
