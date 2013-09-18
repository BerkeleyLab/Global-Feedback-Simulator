#include "state_space_top.h"

#include <stdio.h>
#include <stdlib.h>



Linac_State *** allocate_states(Linac_Param ** linp_array, 
				int Nlinac, int Nhist)
{
  Linac_State *** linss_array;
  int l, i;
  linss_array = (Linac_State***)calloc(Nlinac,sizeof(Linac_State**));
  for(l=0;l<Nlinac;l++) {
    linss_array[l] = (Linac_State**)calloc(Nhist,sizeof(Linac_State*));
    for(i=0;i<Nhist;i++) {
      linss_array[l][i] = (Linac_State*)calloc(1,sizeof(Linac_State));
      Linac_State_Allocate(linss_array[l][i],linp_array[l]);
    }
  }
  return linss_array;
}
void deallocate_states(Linac_State *** linss_array, int Nlinac, int Nhist)
{
  int l, i;
  for(l=0;l<Nlinac;l++) {
    for(i=0;i<Nhist;i++) {
      free(linss_array[l][i]);
    }
    free(linss_array[l]);
  }
  free(linss_array);
}

#define CPRINT(c) {if(cimag(c)<0) fprintf(fp,"%10.16e%10.16ej ",creal(c),cimag(c)); else fprintf(fp,"%10.16e+%10.16ej ",creal(c),cimag(c)); }

void filewrite(FILE * fp, double time, Gun_Param * gun, 
	       Dynamic_Param *dynp, Doublecompress_State * dcs,
	       Linac_Param ** linp_array, Linac_State *** linss_array, int Nlinac,
	       double * error_vol_a, double * error_vol_p) 
{
  int l;
  fprintf(fp,"%10.16e   ",time);
  fprintf(fp,"%10.16e %10.16e %10.16e   ",gun->Q, dynp->dQ_Q, 
	  dynp->dtg,dcs->dE_E[Nlinac-1],dcs->dt[Nlinac-1]);

  CPRINT(dynp->adc_noise);

  for(l=0;l<Nlinac;l++) {
    fprintf(fp,"%10.16e %10.16e %10.16e %10.16e %10.16e %10.16e %10.16e   ",error_vol_a[l],error_vol_p[l],dcs->dE_E[l],dcs->dt[l], dcs->sz[l], dcs->dE_Ei[l], dcs->dE_Ei2[l]);
    
    CPRINT(linss_array[l][0]->cav.voltage);
    CPRINT(linss_array[l][0]->fpga.err);
    CPRINT(linp_array[l]->fpga.set_point);
    CPRINT(linss_array[l][0]->fpga.drive);
    CPRINT(linss_array[l][0]->fpga.state);
    
  }
  
  fprintf(fp,"\n");
  
}
#undef CPRINT

void state_space_top(Gun_Param * gun, Linac_Param ** linp_array, int Nlinac,
		     Linac_State *** linss_array, int Nhist,
		     BBF_Param * bbf, Noise_Source * nsrc, 
		     double simdt, int NT, int openloop,
		     char * fname, int OUTPUTFREQ
		     )
{
  double time;
  int ti,l;
  Dynamic_Param dynp;
  // Outputs of step_llrf that go into double compress
  double * error_vol_p, *error_vol_a;
  // Outputs of double compress
  Doublecompress_State dcs; //TODO: THESE NEED AN INITIAL VALUE

  dynp.dQ_Q = 0.0;
  dynp.dtg = 0.0;
  dynp.dE_ing = 0.0;
  dynp.dsig_z = 0.0;
  dynp.dsig_E = 0.0;
  dynp.dchirp = 0.0;
  dynp.adc_noise = 0.0;

  /*
   * Allocate temporary memory
   */
  error_vol_p = (double*)calloc(Nlinac,sizeof(double));
  error_vol_a = (double*)calloc(Nlinac,sizeof(double));
  Doublecompress_State_Alloc(&dcs, Nlinac);
  
  /*
   * Open the file
   */
  FILE * fp = NULL;
  if(fname != NULL) {
    fp = fopen(fname,"w");
  }
  /*******************************************
   * Start the mainloop
   *******************************************/
  for(ti=0; ti<NT; ti++) {
    time = (ti+1)*simdt;

    // printf("%i\n",ti);
    /***
     * Call apply noise, where noise and stuff is introduced
     ***/
    if(nsrc!=NULL) {
      Apply_Noise(ti,simdt, nsrc, &dynp);
    }


    /***
     * Loop over each linac and call step_llrf
     ***/
    for(l=0; l<Nlinac; l++) {
      // Cycle the buffer of histories
      cycle_buffer(linss_array[l], Nhist);

      // Call step_llrf: passing in simdt is unneccessary it seems
      step_llrf(linp_array[l], simdt, dcs.dt[l], gun->Q*(1.0+dynp.dQ_Q),
		dynp.adc_noise,
	        openloop,
		linss_array[l]);

      // calculate those errors
      error_vol_p[l] = carg( linss_array[l][0]->cav.voltage );
      if(linp_array[l]->fpga.set_point!=0.0) {
	error_vol_a[l] =
	  (cabs(linss_array[l][0]->cav.voltage )
	   - linp_array[l]->fpga.set_point)
	  / (linp_array[l]->fpga.set_point);
      } else {
	error_vol_a[l] =
	  (cabs(linss_array[l][0]->cav.voltage )
	   - linp_array[l]->fpga.set_point);
      }
    } // End loop over linacs
    

    /***
     * Call double_compress
     ***/
    doublecompress_octave_ourtypes(gun,linp_array,Nlinac,
		       &dynp,error_vol_p,error_vol_a,
		       &dcs);


    /***
     * Do beam based feedback
     ***/
    BBF_Step(bbf,&dcs, linp_array,Nlinac);


    /***
     * Output transient data
     ***/
    if(ti%OUTPUTFREQ==0) {
      if(fp!=NULL) {
	filewrite(fp, ti*simdt, gun, &dynp, &dcs,
		  linp_array, linss_array, Nlinac,
		  error_vol_a, error_vol_p);
      }
    }


  }
  /**********************************************
   * End mainloop over ti
   **********************************************/


  /*
   * Free allocated memory
   */
  free(error_vol_p);
  free(error_vol_a);
  Doublecompress_State_Dealloc(&dcs);
  if(fp!=NULL) fclose(fp);
}
