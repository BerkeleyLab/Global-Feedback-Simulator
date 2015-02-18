#!/usr/bin/python

import readjson.parse_accelerator as parseAcc

def Compute_ElecModes(cav, Tstep):
    import numpy as np
    
    # Initialize an empty list to return
    modes_out = []

    beam_phase = cav.nom_beam_phase['value']
    for idx, mode in enumerate(cav.elec_modes):
        mode_name = mode.mode_name
        w0 = mode.omega_0_mode['value']
        foffset = mode.foffset['value']
        RoverQ = mode.RoverQ['value']

        k_probe = np.exp(1j*mode.phase_probe['value'])/np.sqrt(mode.Q_probe['value']*RoverQ);
        k_em = np.exp(1j*mode.phase_rev['value'])/np.sqrt(mode.Q_drive['value']*RoverQ);

        Q_L = 1/(1/mode.Q_0['value'] + 1/mode.Q_drive['value'] + 1/mode.Q_probe['value'])
        bw = w0/(2.0*Q_L);
        k_beam = RoverQ*Q_L*np.exp(-1j*beam_phase)/Tstep;
        k_drive = 2*np.sqrt(mode.Q_drive['value']*RoverQ);
        mode_dict = {"mode_name": mode_name,"w0": w0, "beam_phase": beam_phase, "RoverQ": RoverQ, "foffset": foffset, "Q_L": Q_L, "bw": bw, "k_beam": k_beam, "k_drive": k_drive, "k_probe": k_probe, "k_em": k_em}
        modes_out.append(mode_dict)

    return modes_out

def Get_SWIG_Cavity(cavity_test_file, Verbose=True):

    if Verbose: print "\nLoading JSON configuration files ..."
    #Get the simulation and accelerator objects from the JSON parser
    file_list =  [
        "source/configfiles/unit_tests/default_accelerator.json",
        "source/configfiles/unit_tests/LCLS-II_accelerator.json",
        cavity_test_file]

    simulation, accelerator = parseAcc.ParseAccelerator(file_list, Verbose)

    Tstep = simulation.Tstep['value']
    cavity_object = accelerator.linac_list[0].module_list[0].station_list[0].cavity
    cavity, cavity_state = cavity_object.Configure_C()
        
    modes_config = Compute_ElecModes(cavity_object, Tstep)

    return cavity, cavity_state, Tstep, modes_config

if __name__=="__main__":

    # Convert user-defined configuration into SWIG-wrapped C handlers
    Get_SWIG_Cavity()