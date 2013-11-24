#ifndef STATE_SPACE_TOP_H
#define STATE_SPACE_TOP_H

/*
 * state_space_top.h/c
 *
 * The mainloop for the beam and linac simulation.
 *
 * Alejandro F Queiruga
 * Daniel S Driver
 * LBL 2013
 *
 */

/*
 * TODO: Fill in license
 */

#include <complex.h>
#include "filter.h"
#include "linac_param.h"
#include "step_llrf.h"

#include "doublecompress.h"
#include "dynamic_noise.h"

#include "beam_based_feedback.h"

Linac_State *** allocate_states(Linac_Param ** linp_array, 
				int Nlinac, int Nhist);
void deallocate_states(Linac_State *** linss_array, int Nlinac, int Nhist);

/*
 * Performs NT timesteps
 */
void state_space_top(Gun_Param * gun, Linac_Param ** linp_array, int Nlinac,
		     Linac_State *** linss_array, int Nhist,
		     BBF_Param * bbf, Noise_Source * nrsc,
		     double simdt, int NT, int openloop,
		     char * fname, int OUTPUTFREQ
		     );

#endif
