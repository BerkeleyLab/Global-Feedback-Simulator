#!/usr/bin/python

import readjson.parse_simulation as parseSim


def Get_SWIG_Cavity(cavity_test_file, Verbose=True):

    # Get the Simulation object from the JSON parser
    if Verbose: print "\nLoading JSON configuration files ..."

    file_list =  [
        "source/configfiles/unit_tests/default_accelerator.json",
        "source/configfiles/unit_tests/LCLS-II_accelerator.json",
        cavity_test_file]

    simulation = parseSim.ParseSimulation(file_list, Verbose)

    Tstep = simulation.Tstep['value']
    cavity_object = simulation.linac_list[0].cryomodule_list[0].station_list[0].cavity
    cavity = cavity_object.Get_C_Pointer()
    cavity_state = cavity_object.Get_State_Pointer() 
    
    modes_config = []
    for idx, mode in enumerate(cavity_object.elec_modes):
        mode_dict = mode.Compute_ElecMode(Tstep, cavity_object.rf_phase['value'])
        modes_config.append(mode_dict)

    return cavity_object, Tstep, modes_config

def Get_SWIG_RF_Station(rf_station_test_file, Verbose=True):

    # Get the Simulation object from the JSON parser
    if Verbose: print "\nLoading JSON configuration files ..."

    file_list =  [
        "source/configfiles/unit_tests/default_accelerator.json",
        "source/configfiles/unit_tests/LCLS-II_accelerator.json"]

    if rf_station_test_file:
        file_list.append(rf_station_test_file)

    simulation = parseSim.ParseSimulation(file_list, Verbose)

    Tstep = simulation.Tstep['value']
    rf_station_object = simulation.linac_list[0].cryomodule_list[0].station_list[0]

    fund_index = rf_station_object.cavity.fund_index['value']
    rf_phase = rf_station_object.cavity.rf_phase['value']

    fund_mode_dict = rf_station_object.cavity.elec_modes[fund_index].Compute_ElecMode(Tstep, rf_phase)

    # Overwrite pre-computed settings for Set-point from Linac configuration
    nom_grad = rf_station_object.cavity.nom_grad['value']
    L = rf_station_object.cavity.L['value']
    rf_station_object.cavity.design_voltage['value'] = nom_grad*L

    rf_station = rf_station_object.Get_C_Pointer()
    rf_state = rf_station_object.Get_State_Pointer()

    return rf_station_object, Tstep, fund_mode_dict

def Get_SWIG_Cryomodule(cryo_test_file, Verbose=True):

    # Get the Simulation object from the JSON parser
    if Verbose: print "\nLoading JSON configuration files ..."

    file_list =  [
        "source/configfiles/unit_tests/default_accelerator.json",
        "source/configfiles/unit_tests/LCLS-II_accelerator.json"]

    if cryo_test_file:
        file_list.append(cryo_test_file)

    simulation = parseSim.ParseSimulation(file_list, Verbose)

    Tstep = simulation.Tstep['value']
    cryo_object = simulation.linac_list[0].cryomodule_list[0]

    cryo = cryo_object.Get_C_Pointer()
    cryo_state = cryo_object.Get_State_Pointer()

    fund_mode_dicts = []
    station_list = simulation.linac_list[0].cryomodule_list[0].station_list
    for rf_station in station_list:
        rf_phase = rf_station.cavity.rf_phase['value']
        fund_index = rf_station.cavity.fund_index['value']
        fund_mode_dict = rf_station.cavity.elec_modes[fund_index].Compute_ElecMode(Tstep, rf_phase)
        fund_mode_dicts.append(fund_mode_dict)

    return cryo_object, Tstep, fund_mode_dicts

def Get_SWIG_Simulation(sim_test_file=None, Verbose=True):

    # Get the Simulation object from the JSON parser
    if Verbose: print "\nLoading JSON configuration files ..."

    file_list =  [
        "source/configfiles/unit_tests/default_accelerator.json"]

    if sim_test_file:
        file_list.append(sim_test_file)

    sim_object = parseSim.ParseSimulation(file_list, Verbose)

    # Extract simulation time-step
    Tstep = sim_object.Tstep['value']

    sim_C_Pointer = sim_object.Get_C_Pointer()
    sim_State_Pointer = sim_object.Get_State_Pointer()

    return sim_object

if __name__=="__main__":

    # Convert user-defined configuration into SWIG-wrapped C handlers
    Get_SWIG_Cavity()