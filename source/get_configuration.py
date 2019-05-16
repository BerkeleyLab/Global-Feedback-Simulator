"""
Wrapper functions to get references to SWIG-wrapped C structures (used by unit tests).
Functions are provided for different levels of the model hierarchy to run
unit tests on individual components, where the configuration is parsed from
JSON files.
"""

import readjson.parse_simulation as parseSim

def Get_SWIG_Cavity(cavity_test_file, Verbose=True):

    # Get the Simulation object from the JSON parser
    if Verbose:
        print "\nLoading JSON configuration files ..."

    file_list = [
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

def Get_SWIG_RF_Station(rf_station_test_file, Verbose=True, set_point=None):

    # Get the Simulation object from the JSON parser
    if Verbose:
        print "\nLoading JSON configuration files ..."

    file_list = [
        "source/configfiles/unit_tests/default_accelerator.json"]

    if isinstance(rf_station_test_file, list):
        for file in rf_station_test_file:
            file_list.append(file)
    else:
        file_list.append(rf_station_test_file)

    simulation = parseSim.ParseSimulation(file_list, Verbose)

    Tstep = simulation.Tstep['value']
    rf_station_object = simulation.linac_list[0].cryomodule_list[0].station_list[0]

    fund_index = rf_station_object.cavity.fund_index['value']
    rf_phase = rf_station_object.cavity.rf_phase['value']

    fund_mode_dict = rf_station_object.cavity.elec_modes[fund_index].Compute_ElecMode(Tstep, rf_phase)

    # Overwrite pre-computed settings for Set-point from Linac configuration
    if set_point == None:
        nom_grad = rf_station_object.cavity.nom_grad['value']
        L = rf_station_object.cavity.L['value']
        rf_station_object.cavity.design_voltage['value'] = nom_grad*L
    else:
        rf_station_object.cavity.design_voltage['value'] = set_point

    rf_station = rf_station_object.Get_C_Pointer()
    rf_state = rf_station_object.Get_State_Pointer()

    return rf_station_object, Tstep, fund_mode_dict

def Get_SWIG_Cryomodule(cryo_test_file, Verbose=True):

    # Get the Simulation object from the JSON parser
    if Verbose:
        print "\nLoading JSON configuration files ..."

    file_list = [
        "source/configfiles/unit_tests/default_accelerator.json",
        "source/configfiles/unit_tests/LCLS-II_accelerator.json"]

    if cryo_test_file:
        file_list.append(cryo_test_file)

    simulation = parseSim.ParseSimulation(file_list, Verbose)

    Tstep = simulation.Tstep['value']
    cryo_object = simulation.linac_list[0].cryomodule_list[0]

    fund_mode_dicts = []
    station_list = simulation.linac_list[0].cryomodule_list[0].station_list
    for rf_station in station_list:
        rf_phase = rf_station.cavity.rf_phase['value']
        fund_index = rf_station.cavity.fund_index['value']
        fund_mode_dict = rf_station.cavity.elec_modes[fund_index].Compute_ElecMode(Tstep, rf_phase)
        fund_mode_dicts.append(fund_mode_dict)

        # Overwrite set-point values to maximum accelerating voltage for test
        nom_grad = rf_station.cavity.nom_grad['value']
        L = rf_station.cavity.L['value']
        rf_station.cavity.design_voltage['value'] = nom_grad*L

    cryo = cryo_object.Get_C_Pointer()
    cryo_state = cryo_object.Get_State_Pointer()

    return cryo_object, Tstep, fund_mode_dicts

def Get_SWIG_Simulation(sim_test_files=None, Verbose=True):

    # Get the Simulation object from the JSON parser
    if Verbose:
        print "\nLoading JSON configuration files ..."

    file_list = [
        "source/configfiles/unit_tests/default_accelerator.json"]

    if sim_test_files:
        for sim_test_file in sim_test_files:
            file_list.append(sim_test_file)

    sim_object = parseSim.ParseSimulation(file_list, Verbose)

    # Extract simulation time-step
    Tstep = sim_object.Tstep['value']

    sim_C_Pointer = sim_object.Get_C_Pointer()
    sim_State_Pointer = sim_object.Get_State_Pointer()

    return sim_object

if __name__ == "__main__":

    # Convert user-defined configuration into SWIG-wrapped C handlers
    Get_SWIG_Cavity()
