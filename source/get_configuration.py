#!/usr/bin/python

import readjson.parse_accelerator as parseAcc



def Get_SWIG_Cavity(cavity_test_file, Verbose=True):

    if Verbose: print "\nLoading JSON configuration files ..."
    #Get the simulation and accelerator objects from the JSON parser
    file_list =  [
        "source/configfiles/unit_tests/default_accelerator.json",
        "source/configfiles/unit_tests/LCLS-II_accelerator.json",
        cavity_test_file]

    simulation, accelerator = parseAcc.ParseAccelerator(file_list, Verbose)

    Tstep = simulation.Tstep['value']
    cavity_object = accelerator.linac_list[0].cryomodule_list[0].station_list[0].cavity
    cavity = cavity_object.Get_C_Pointer()
    cavity_state = cavity_object.Get_State_Pointer(cavity) 
    
    modes_config = []
    for idx, mode in enumerate(cavity_object.elec_modes):
        mode_dict = mode.Compute_ElecMode(Tstep, cavity_object.nom_beam_phase['value'])
        modes_config.append(mode_dict)

    return cavity, cavity_state, Tstep, modes_config

def Get_SWIG_RF_Station(rf_station_test_file, Verbose=True):

    if Verbose: print "\nLoading JSON configuration files ..."
    #Get the simulation and accelerator objects from the JSON parser
    file_list =  [
        "source/configfiles/unit_tests/default_accelerator.json",
        "source/configfiles/unit_tests/LCLS-II_accelerator.json"]

    if rf_station_test_file:
        file_list.append(rf_station_test_file)

    simulation, accelerator = parseAcc.ParseAccelerator(file_list, Verbose)

    Tstep = simulation.Tstep['value']
    rf_station_object = accelerator.linac_list[0].cryomodule_list[0].station_list[0]

    fund_index = rf_station_object.cavity.fund_index['value']
    nom_beam_phase = rf_station_object.cavity.nom_beam_phase['value']

    fund_mode_dict = rf_station_object.cavity.elec_modes[fund_index].Compute_ElecMode(Tstep, nom_beam_phase)

    rf_station = rf_station_object.Get_C_Pointer()
    rf_state = rf_station_object.Get_State_Pointer(rf_station)

    return rf_station, rf_state, Tstep, fund_mode_dict

if __name__=="__main__":

    # Convert user-defined configuration into SWIG-wrapped C handlers
    Get_SWIG_Cavity()