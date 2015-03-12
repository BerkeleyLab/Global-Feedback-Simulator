#!/usr/bin/python

#
# Accelerator-specific configuration file:
#    Defines all the Python classes involved in the simulation.
#    Parses configuration information from a dictionary and instantiates Python objects with the configuration values.
#
# Classes defined:
#   Simulation
#   Accelerator
#   Cavity
#   ElecMode
#   MechMode
#   Piezo
#
# Call the readConfiguration function in order to get a full collection of instances for a simulation run.
#
# readjson_accelerator.py
#

from readjson import readentry
from math import pi

# Define Simulation time step as global
Tstep_global = 0.0

class Simulation:
    """ Simulation class: contains parameters specific to a simulation run"""

    def __init__(self, name, comp_type, Tstep, time_steps, nyquist_sign, synthesis):
        self.name = name
        self.type = comp_type
        self.Tstep = Tstep
        self.time_steps = time_steps
        self.nyquist_sign = nyquist_sign
        self.synthesis = synthesis

        # Assign Tstep to the global variable
        global Tstep_global
        Tstep_global = Tstep['value']

    def __str__ (self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--Simulation Object--\n"
         + "name: " + self.name  + "\n"
         + "type: " + self.type  + "\n"
         + "Tstep: " + str(self.Tstep) + "\n"
         + "time_steps: " + str(self.time_steps) + "\n"
         + "nyquist_sign: " + str(self.nyquist_sign) + "\n"
         + "synthesis: " + str(self.synthesis) + "\n")

def readSimulation(confDict):
    """ readSimulation: Takes the global configuration dictionary and
    returns a Simulation object with the configuration values filled in"""

    # Read name and component type (strings, no need to go through readentry)
    name = confDict["Simulation"]["name"]
    comp_type = confDict["Simulation"]["type"]

    # Read rest of configuration parameters
    Tstep = readentry(confDict, confDict["Simulation"]["Tstep"])
    time_steps = readentry(confDict, confDict["Simulation"]["time_steps"])
    nyquist_sign = readentry(confDict,confDict["Simulation"]["nyquist_sign"])

    # Check if simulation dictionary has a Synthesis entry, and if so parse it
    if confDict["Simulation"].has_key("Synthesis"):
        synthesis = readSynthesis(confDict["Simulation"])
    else:
        synthesis = None

    # Instantiate Simulation object and return
    simulation = Simulation(name, comp_type, Tstep, time_steps, nyquist_sign, synthesis)

    return simulation

class Synthesis:
    """ Synthesis class: contains parameters specific to a Synthesis run.
    The parameters in this class are not run-time configurable. They therefore
    need to mirror the synthesizable Verilog and are used to compute FPGA
    register settings.
    """

    def __init__(self, name, comp_type, n_mech_modes, df_scale):
        self.name = name
        self.type = comp_type
        self.n_mech_modes = n_mech_modes
        self.df_scale = df_scale

    def __str__ (self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--Synthesis Object--\n"
         + "name: " + self.name  + "\n"
         + "type: " + self.type  + "\n"
         + "n_mech_modes: " + str(self.n_mech_modes) + "\n"
         + "df_scale: " + str(self.df_scale) + "\n")

def readSynthesis(simDict):
    """ readSynthesis: Takes the global configuration dictionary and
    returns a Synthesis object with the configuration values filled in"""

    # Read name and component type (strings, no need to go through readentry)
    name = simDict["Synthesis"]["name"]
    comp_type = simDict["Synthesis"]["type"]

    # Read rest of configuration parameters
    n_mech_modes = readentry(simDict, simDict["Synthesis"]["n_mech_modes"])
    df_scale = readentry(simDict, simDict["Synthesis"]["df_scale"])

    # Instantiate Synthesis object and return
    synthesis = Synthesis(name, comp_type, n_mech_modes, df_scale)

    return synthesis

class Cavity:
    """ Cavity class: contains parameters specific to a cavity,
            including a nested list of electrical modes"""
    def __init__(self, name, comp_type, param_dic, elec_modes):
        """ Cavity constructor:
            Inputs:
                name: Cavity instance name,
                comp_type: component type (Cavity),
                param_dic: dictionary containing all the Cavity parameters.
                elec_modes: list of ElecMode objects (one per electrical mode).
            Output: Cavity object."""

        # Name and component type (strings, no need to go through readentry)
        self.name = name
        self.type = comp_type

        # Rest of configuration parameters
        self.L = param_dic["L"]
        self.nom_grad = param_dic["nom_grad"]
        self.nom_beam_phase = param_dic["nom_beam_phase"]
        self.rf_phase = param_dic["rf_phase"]
        self.design_voltage = param_dic["design_voltage"]
        self.unity_voltage = param_dic["unity_voltage"]

        # List of electrical mode (ElecMode) instances (objects)
        self.elec_modes = elec_modes

        # Find the fundamental mode based on coupling to the beam
        ## Criterium here is the fundamental mode being defined as that with the highest shunt impedance (R/Q)
        RoverQs = map(lambda x: x.RoverQ['value'],self.elec_modes)

        fund_index = RoverQs.index(max(RoverQs))

        # Store the index of the fundamental mode
        self.fund_index = {"value" : fund_index, "units" : "N/A", "description" : "Index of the fundamental mode in array"}

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--Cavity Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "L: " + str(self.L) + "\n"
        + "nom_grad: " + str(self.nom_grad) + "\n"
        + "nom_beam_phase: " + str(self.nom_beam_phase) + "\n"
        + "rf_phase: " + str(self.rf_phase) + "\n"
        + "design_voltage: " + str(self.design_voltage) + "\n"
        + "unity_voltage: " + str(self.unity_voltage) + "\n"
        + "electrical modes: " + '\n'.join(str(x) for x in self.elec_modes))

    def Get_C_Pointer(self):
        import accelerator as acc
        # First count number of Electrical Modes and Allocate Array
        n_modes = len(self.elec_modes)
        elecMode_net = acc.ElecMode_Allocate_Array(n_modes)

        # Allocate each Electrical Mode and append it to the elecMode_net
        for idx, mode in enumerate(self.elec_modes):
            n_mech = len(mode.mech_couplings_list)
            mech_couplings = acc.double_Array(n_mech)
            for m in xrange(n_mech):
                mech_couplings[m] = mode.mech_couplings_list[m]

            elecMode = acc.ElecMode_Allocate_New(mode.RoverQ['value'], \
                mode.foffset['value'], mode.omega_0_mode['value'], \
                mode.Q_0['value'], mode.Q_drive['value'], mode.Q_probe['value'], \
                self.nom_beam_phase['value'],  mode.phase_rev['value'], mode.phase_probe['value'], \
                Tstep_global, mech_couplings, n_mech)

            acc.ElecMode_Append(elecMode_net, elecMode, idx)

        L = self.L['value']
        nom_grad = self.nom_grad['value']
        nom_beam_phase = self.nom_beam_phase['value']
        rf_phase = self.rf_phase['value']
        design_voltage = self.design_voltage['value']
        unity_voltage = self.unity_voltage['value']
        fund_index = self.fund_index['value']

        # Get a C-pointer to a Cavity structure
        cavity = acc.Cavity_Allocate_New(elecMode_net, n_modes, L, nom_grad, \
            nom_beam_phase, rf_phase, design_voltage, unity_voltage, \
            fund_index)

        return cavity

    @staticmethod
    def Get_State_Pointer(cavity_C_pointer):
        import accelerator as acc

        cavity_state = acc.Cavity_State()
        acc.Cavity_State_Allocate(cavity_state, cavity_C_pointer)

        return cavity_state


class ElecMode:
    def __init__(self, name, comp_type, mode_name, param_dic, mech_couplings_list):
        """ ElecMode class: contains parameters specific to an electrical mode,
            including a dictionary specifying the mechanical couplings.
            Note the absence of a readElecMode method, the process for parsing
            the global configuration dictionary and creating ElecMode objects
            is done recursively in Cavity.readCavity(...)"""

        self.name = name
        self.type = comp_type
        self.mode_name = mode_name

        self.RoverQ = param_dic["RoverQ"]
        self.foffset = param_dic["foffset"]
        self.peakV = param_dic["peakV"]
        self.Q_0 = param_dic["Q_0"]
        self.Q_drive = param_dic["Q_drive"]
        self.Q_probe = param_dic["Q_probe"]
        self.phase_rev = param_dic["phase_rev"]
        self.phase_probe = param_dic["phase_probe"]
        ## Add (replicate) a parameter that will be filled after object instance
        self.omega_0_mode = {"value" : 0.0, "units" : "rad/s", "description" : "Linac's Nominal resonance angular frequency"}

        self.mech_couplings_list = mech_couplings_list

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--ElecMode Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "mode_name: " + str(self.mode_name) + "\n"
        + "RoverQ: " + str(self.RoverQ) + "\n"
        + "foffset: " + str(self.foffset) + "\n"
        + "peakV: " + str(self.peakV) + "\n"
        + "Q_0: " + str(self.Q_0) + "\n"
        + "Q_drive: " + str(self.Q_drive) + "\n"
        + "Q_probe: " + str(self.Q_probe) + "\n"
        + "phase_rev: " + str(self.phase_rev) + "\n"
        + "phase_probe: " + str(self.phase_probe) + "\n"
        + "mech_couplings_list: " + str(self.mech_couplings_list))

    def Compute_ElecMode(self, Tstep, nom_beam_phase):
        import numpy as np

        # Initialize an empty list to return
        modes_out = []

        beam_phase = nom_beam_phase

        mode_name = self.mode_name
        w0 = self.omega_0_mode['value']
        foffset = self.foffset['value']
        RoverQ = self.RoverQ['value']

        k_probe = np.exp(1j*self.phase_probe['value'])/np.sqrt(self.Q_probe['value']*RoverQ);
        k_em = np.exp(1j*self.phase_rev['value'])/np.sqrt(self.Q_drive['value']*RoverQ);

        Q_L = 1/(1/self.Q_0['value'] + 1/self.Q_drive['value'] + 1/self.Q_probe['value'])
        bw = w0/(2.0*Q_L);
        k_beam = RoverQ*Q_L*np.exp(-1j*beam_phase)/Tstep;
        k_drive = 2*np.sqrt(self.Q_drive['value']*RoverQ);
        mode_dict = {"mode_name": mode_name,"w0": w0, "beam_phase": beam_phase, "RoverQ": RoverQ, "foffset": foffset, "Q_L": Q_L, "bw": bw, "k_beam": k_beam, "k_drive": k_drive, "k_probe": k_probe, "k_em": k_em}

        return mode_dict

class MechMode:
    """ MechMode class: contains parameters specific to a mechanical mode.
        Information concerning couplings with electrical modes and Piezos is
        contained in ElecMode and Piezo objects respectively."""

    def __init__(self, name, comp_type, param_dic):
        self.name = name
        self.type = comp_type

        self.f0 = param_dic["f0"]
        self.Q = param_dic["Q"]
        self.full_scale= param_dic["full_scale"]

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--MechMode Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "f0: " + str(self.f0) + "\n"
        + "Q: " + str(self.Q) + "\n"
        + "full_scale: " + str(self.full_scale) + "\n")

class Piezo:
    """ Piezo class: contains couplings between the Piezo and each
    one of the mechanical modes (MechMode instances)."""

    def __init__(self, name, comp_type, mech_couplings_list, VPmax):
        self.name = name
        self.type = comp_type

        self.VPmax = VPmax
        self.mech_couplings_list = mech_couplings_list

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--Piezo Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "VPmax: " + str(self.VPmax) + "\n"
        + "mech_couplings_list: " + str(self.mech_couplings_list))

def readMechMode(confDict, mech_mode_entry):
    """ readMechMode: Takes the global configuration dictionary and
    returns a MechMode object with the configuration values filled in
    Inputs:
        confDict: Global configuration dictionary,
        mech_mode_entry: Name of the mechanical mode to be read (string)
    Output:
        mech_mode: MechMode object"""

    # Read name and component type
    name = confDict[mech_mode_entry]['name']
    comp_type = confDict[mech_mode_entry]['type']

    # Read the rest of the configuration parameters and store in a dictionary
    param_dic = {}

    param_dic["f0"] = readentry(confDict,confDict[mech_mode_entry]["f0"])
    param_dic["Q"] = readentry(confDict,confDict[mech_mode_entry]["Q"])
    param_dic["full_scale"] = readentry(confDict,confDict[mech_mode_entry]["full_scale"])

    # Create and return a mechanical mode (MechMode) instance
    mech_mode = MechMode(name, comp_type, param_dic)

    return mech_mode

def readCouplings(confDict, mech_couplings, module_entry):
    """ readCouplings: Takes the global configuration dictionary and a dictionary containing non-zero
    values for the mechanical couplings (i.e. coupling between electrical modes or piezos and mechanical modes).
    The module_entry input is necessary in order to access the proper module's mechanical mode list.
    The length of the coupling vector used in the simulation code must be equal to the number of mechanical modes (M).
    The mech_couplings input is supposed to contain non-zero values from the configuration file,
    and readCouplings always returns a dictionary of M elements, where the couplings not specified in
    the input file are filled with 0s.
    Inputs:
        confDict: Global configuration dictionary,
        mech_couplings: dictionary containing couplings defined in the configuration file (0 to M elements).
        module_entry: Module entry in global dictionary in order to access the proper modules mechanical mode list.
    Output:
        mech_couplings_list: ordered list containing mechanical couplings for an electrical mode or piezo (length M).
            Order corresponds to the order of appearance of the mechanical mode in mech_net."""

    # Grab the full list of mechanical modes in the Module
    mech_net = confDict[module_entry]["mechanical_mode_connect"]

    # Make an ordered list of size M (where M is the total number of mechanical modes, length of mech_net)
    # Fill with 0s if no coupling is specified in mech_couplings by the user
    mech_couplings_list = [mech_couplings[m] if m in mech_couplings else 0.0 for m in mech_net]

    return mech_couplings_list

def readCavity(confDict, cav_entry, module_entry):
    """ readCavity: Takes the global configuration dictionary and returns a Cavity object
    with the configuration values filled in. The process includes a recursive read of
    electrical modes in each cavity, where ElecMode objects are created for each electrical mode
    and contained as a list of ElecMode objects in the Cavity object.
    Inputs:
        confDict: Global configuration dictionary,
        cav_entry: Name of the cavity to be read (string).
        module_entry: Module entry in global dictionary in order to access the proper module's
            mechanical mode list, wich is used as a consistency check to generate mechanical
            coupling vectors for each electrical mode.
    Output:
        cavity: Cavity object."""

    # Read name and component type
    name = confDict[cav_entry]['name']
    cav_comp_type = confDict[cav_entry]['type']

    # Read and store the rest of the parameters in a dictionary
    cav_param_dic = {}

    cav_param_dic["L"] = readentry(confDict,confDict[cav_entry]["L"])
    cav_param_dic["nom_grad"] = readentry(confDict,confDict[cav_entry]["nom_grad"])

    cav_param_dic["nom_beam_phase"] = readentry(confDict,confDict[cav_entry]["nom_beam_phase"])
    cav_param_dic["rf_phase"] = readentry(confDict,confDict[cav_entry]["rf_phase"])
    cav_param_dic["design_voltage"] = readentry(confDict,confDict[cav_entry]["design_voltage"])
    cav_param_dic["unity_voltage"] = readentry(confDict,confDict[cav_entry]["unity_voltage"])


    # Grab the list of electrical modes
    elec_mode_connect = confDict[cav_entry]["elec_mode_connect"]
    n_elec_modes = len(elec_mode_connect) # Number of electrical modes

    elec_mode_list = []

    ## Start of loop through electrical modes
    # Cycle through electrical modes, read parameters from global dictionary and append to list of modes.
    for m in range(n_elec_modes):
        # Take mth element of mode list
        elec_mode_now = elec_mode_connect[m]

        # Read component name and type
        elec_mode_name = confDict[elec_mode_now]['name']
        elec_comp_type = confDict[elec_mode_now]['type']

        # Identifier for mode
        mode_name = confDict[elec_mode_now]['mode_name']

        # Create dictionary to store electrical mode parameters
        elec_param_dic = {}

        # Read rest of parameters and store in dictionary
        elec_param_dic["RoverQ"] = readentry(confDict,confDict[elec_mode_now]["RoverQ"])
        elec_param_dic["foffset"] = readentry(confDict,confDict[elec_mode_now]["foffset"])
        elec_param_dic["peakV"] = readentry(confDict,confDict[elec_mode_now]["peakV"])
        elec_param_dic["Q_0"] = readentry(confDict,confDict[elec_mode_now]["Q_0"])
        elec_param_dic["Q_drive"] = readentry(confDict,confDict[elec_mode_now]["Q_drive"])
        elec_param_dic["Q_probe"] = readentry(confDict,confDict[elec_mode_now]["Q_probe"])
        elec_param_dic["phase_rev"] = readentry(confDict,confDict[elec_mode_now]["phase_rev"])
        elec_param_dic["phase_probe"] = readentry(confDict,confDict[elec_mode_now]["phase_probe"])

        # Read dictionary of couplings from global configuration dictionary
        mech_couplings = readentry(confDict,confDict[elec_mode_now]["mech_couplings"]["value"])

        # Get a coupling list of length M (number of mechanical modes),
        # filled with 0s if no coupling is specified by user
        mech_couplings_list = readCouplings(confDict, mech_couplings, module_entry)

        # Instantiate ElecMode object
        elec_mode = ElecMode(elec_mode_name, elec_comp_type, mode_name, elec_param_dic, mech_couplings_list)

        # Append to list of electrical modes
        elec_mode_list.append(elec_mode)
    ## End of loop through electrical modes

    # Instantiate Cavity object and return
    cavity = Cavity(name, cav_comp_type, cav_param_dic, elec_mode_list)

    return cavity

def readPiezo(confDict, piezo_entry, module_entry):
    """ readPiezo: Takes the global configuration dictionary and returns a Piezo object
    with the configuration parameters and the couplings with each one of the the mechanical
    modes values filled in.
    Inputs:
        confDict: Global configuration dictionary,
        piezo_entry: Name of the Piezo to be read (string).
        module_entry: Module entry in global dictionary in order to access the proper module's
            mechanical mode list, wich is used as a consistency check to generate mechanical
            coupling vectors for each Piezo.
    Output:
        piezo: Piezo object."""

    # Read name and component type
    name = confDict[piezo_entry]['name']
    piezo_comp_type = confDict[piezo_entry]['type']

    # Read dictionary of couplings from global configuration dictionary
    mech_couplings = readentry(confDict,confDict[piezo_entry]["mech_couplings"]["value"])

    # Check consistency of coupling entries with the list of mechanical modes,
    # and get a coupling dictionary of length M (number of mechanical modes)
    mech_couplings_list = readCouplings(confDict, mech_couplings, module_entry)

    # Read rest of parameters
    VPmax = readentry(confDict,confDict[piezo_entry]["VPmax"])

    # Create a Piezo instance and return
    piezo = Piezo(name, piezo_comp_type, mech_couplings_list, VPmax)

    return piezo

def readList(confDict, list_in, readFunction, module_entry=None):
    """ readList: Generic function to read list of componentns.
    Takes the global configuration dictionary, cycles through the list of components
    (list of names, list_in), uses the names to identify the configuration entries in
    the global dictionary, calls the proper read function for each component (readFunction),
    and returns a list of instances (objects).
    Inputs:
        confDict: Global configuration dictionary,
        list_in: list of components to cycle through (list of strings),
        readFunction: name of the read function for the component,
        module_entry: necessary in some cases in order to pass along module entry information,
            needed by readStation and readPiezo in order to find the mechanical modes in
            their corresponding Module.
    Output:
        list_out: List of component objects"""

    # Create empty list for component instances
    list_out = []

    # Cycle through list of component names
    for k in range(len(list_in)):
        # Read component configuration and create component instance
        if module_entry == None:
            component = readFunction(confDict, list_in[k])
        else:
            component = readFunction(confDict, list_in[k], module_entry)
        # Append object to the component list
        list_out.append(component)

    # Return list
    return list_out

class Controller:
    """ Controller class: contains parameters specific to a Controller configuration"""

    def __init__(self, name, comp_type, stable_gbw):
        self.name = name
        self.type = comp_type

        self.stable_gbw = stable_gbw

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--Controller Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "stable_gbw: " + str(self.stable_gbw) + "\n")

def readController(confDict, controller_entry):

     # Read name and component type
    name = confDict[controller_entry]['name']
    controller_comp_type = confDict[controller_entry]['type']

    # Read the rest of parameters
    stable_gbw = readentry(confDict,confDict[controller_entry]["stable_gbw"])

    # Create a Controller instance and return
    controller = Controller(name, controller_comp_type, stable_gbw)

    return controller

class ZFilter:
    """ ZFilter class: contains parameters specific to a Filter configuration"""

    def __init__(self, name, comp_type, order, nmodes, poles):
        self.name = name
        self.type = comp_type

        self.order = order
        self.nmodes = nmodes
        self.poles = poles

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--ZFilter Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "order: " + str(self.order) + "\n"
        + "nmodes: " + str(self.nmodes) + "\n"
        + "poles: " + str(self.poles) + "\n")

def readZFilter(confDict, zfilter_entry):

    # Read name and component type
    name = confDict[zfilter_entry]['name']
    zfilter_comp_type = confDict[zfilter_entry]['type']

    # Read the rest of parameters
    order = readentry(confDict,confDict[zfilter_entry]["order"])
    nmodes = readentry(confDict,confDict[zfilter_entry]["nmodes"])
    poles = readentry(confDict,confDict[zfilter_entry]["poles"])

    # Create a ZFilter instance and return
    zfilter =  ZFilter(name, zfilter_comp_type, order, nmodes, poles)

    return zfilter

class ADC:
    """ ADC class: contains parameters specific to a ADC configuration"""

    def __init__(self, name, comp_type, adc_max, adc_off, psd_llrf):
        self.name = name
        self.type = comp_type

        self.adc_max = adc_max
        self.adc_off = adc_off
        self.psd_llrf = psd_llrf

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--ADC Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "adc_max: " + str(self.adc_max) + "\n"
        + "adc_off: " + str(self.adc_off) + "\n"
        + "psd_llrf: " + str(self.psd_llrf) + "\n")

def readADC(confDict, adc_entry):

    # Read name and component type
    name = confDict[adc_entry]['name']
    adc_comp_type = confDict[adc_entry]['type']

    # Read the rest of parameters
    adc_max = readentry(confDict,confDict[adc_entry]["adc_max"])
    adc_off = readentry(confDict,confDict[adc_entry]["adc_off"])
    psd_llrf = readentry(confDict,confDict[adc_entry]["psd_llrf"])

    # Create an ADC instance and return
    adc = ADC(name, adc_comp_type, adc_max, adc_off, psd_llrf)

    return adc


class Amplifier:
    """ Amplifier class: contains parameters specific to a Amplifier configuration"""

    def __init__(self, name, comp_type, PAmax, PAbw, Clip, top_drive):
        self.name = name
        self.type = comp_type

        self.PAmax = PAmax
        self.PAbw = PAbw
        self.Clip = Clip
        self.top_drive = top_drive

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--Amplifier Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "PAmax: " + str(self.PAmax) + "\n"
        + "PAbw: " + str(self.PAbw) + "\n"
        + "Clip: " + str(self.Clip) + "\n"
        + "top_drive: " + str(self.top_drive) + "\n")

    def Get_Saturation_Limit(self):
        """Get_Saturation_Limit: Calculate (measure) the output drive limit from the FPGA controller
        based on the clipping parameter of the saturation function and the maximum output power
        percentage level setting. Returns maximum input value to the saturation function in order
        to reach the percentage of the maximum amplifier output power indicated by top_drive."""

        import numpy as np
        import accelerator as acc
        # Input vector
        inp = np.arange(0.0,10.0,0.1,dtype=np.complex)
        # Output vector
        oup = np.zeros(inp.shape,dtype=np.complex)

        # Boolean indicating finding percentile reach
        found = False
        V_sat = max(inp)

        c = self.Clip['value']
        top_drive = float(self.top_drive['value'])/100

        # Sweep input
        for i in xrange(len(inp)):
            oup[i] = acc.Saturate(inp[i],c)
            if (found == False) and (oup[i].real >= top_drive):
                V_sat = inp[i]
                found = True

        return V_sat.real

def readAmplifier(confDict, amplifier_entry):

    # Read name and component type
    name = confDict[amplifier_entry]['name']
    amplifier_comp_type = confDict[amplifier_entry]['type']

    # Read the rest of parameters
    PAmax = readentry(confDict,confDict[amplifier_entry]["PAmax"])
    PAbw = readentry(confDict,confDict[amplifier_entry]["PAbw"])
    Clip = readentry(confDict,confDict[amplifier_entry]["Clip"])
    top_drive = readentry(confDict,confDict[amplifier_entry]["top_drive"])

    # Create an Amplifier instance and return
    amplifier = Amplifier(name, amplifier_comp_type, PAmax, PAbw, Clip, top_drive)

    return amplifier

class Station:
    """ Station class: contains parameters specific to a Station configuration"""

    def __init__(self, name, comp_type, amplifier, cavity, rx_filter, tx_filter1, tx_filter2, controller, loop_delay_size, cav_adc, fwd_adc, rfl_adc, piezo_list):
        self.name = name
        self.type = comp_type

        self.amplifier = amplifier
        self.cavity = cavity
        self.rx_filter = rx_filter
        self.tx_filter1 = tx_filter1
        self.tx_filter2 = tx_filter2
        self.controller = controller
        self.loop_delay_size = loop_delay_size
        self.cav_adc = cav_adc
        self.fwd_adc = fwd_adc
        self.rfl_adc = rfl_adc
        self.piezo_list = piezo_list

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--Station Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"

        + "amplifier: " + str(self.amplifier) + "\n"
        + "cavity: " + str(self.cavity) + "\n"
        + "rx_filter: " + str(self.rx_filter) + "\n"
        + "tx_filter1: " + str(self.tx_filter1) + "\n"
        + "tx_filter2: " + str(self.tx_filter2) + "\n"
        + "controller: " + str(self.controller) + "\n"
        + "loop_delay_size: " + str(self.loop_delay_size) + "\n"
        + "cav_adc: " + str(self.cav_adc) + "\n"
        + "fwd_adc: " + str(self.fwd_adc) + "\n"
        + "rfl_adc: " + str(self.rfl_adc) + "\n"
        + "piezo_list: " + '\n'.join(str(x) for x in self.piezo_list))

    def Get_C_Pointer(self):

        import accelerator as acc
        import numpy as np

        p_RXF = acc.complexdouble_Array(3)
        p_RXF[0] = complex(self.rx_filter.poles['value'][0][0])*1e6
        p_RXF[1] = complex(self.rx_filter.poles['value'][1][0])*1e6
        p_RXF[2] = complex(self.rx_filter.poles['value'][2][0])*1e6

        p_TRF1 = acc.complexdouble_Array(2)
        p_TRF1[0] = complex(self.tx_filter1.poles['value'][0][0])*1e6
        p_TRF1[1] = complex(self.tx_filter1.poles['value'][1][0])*1e6

        p_TRF2 = acc.complexdouble_Array(1)
        p_TRF2[0] = complex(self.tx_filter2.poles['value'][0][0])*1e6

        Clip = self.amplifier.Clip['value']
        PAmax = self.amplifier.PAmax['value']
        PAscale = np.sqrt(PAmax)*(float(self.amplifier.top_drive['value'])/100) # [sqrt(W)]

        FPGA_out_sat = self.amplifier.Get_Saturation_Limit()*PAscale

        stable_gbw = self.controller.stable_gbw['value']
        loop_delay_size = self.loop_delay_size['value']

        cavity_pointer = self.cavity.Get_C_Pointer()

        rf_station = acc.RF_Station_Allocate_New(Tstep_global, Clip, PAmax, PAscale, p_TRF1, p_TRF2, p_RXF, cavity_pointer, stable_gbw, FPGA_out_sat, loop_delay_size)

        return rf_station

    def Get_State_Pointer(self, RF_Station_C_Pointer):
        import accelerator as acc

        rf_state = acc.RF_State()

        acc.RF_State_Allocate(rf_state, RF_Station_C_Pointer)

        return rf_state

def readStation(confDict, station_entry, module_entry):

    # Read name and component type
    name = confDict[station_entry]['name']
    station_comp_type = confDict[station_entry]['type']

    # Read all the station components
    amplifier_entry = confDict[station_entry]['Amplifier']
    amplifier = readAmplifier(confDict, amplifier_entry)

    cavity_entry = confDict[station_entry]['Cavity']
    cavity = readCavity(confDict, cavity_entry, module_entry)

    rx_filter_entry = confDict[station_entry]['Rx_filter']
    rx_filter = readZFilter(confDict, rx_filter_entry)

    tx_filter1_entry = confDict[station_entry]['Tx_filter1']
    tx_filter1 = readZFilter(confDict, tx_filter1_entry)

    tx_filter2_entry = confDict[station_entry]['Tx_filter2']
    tx_filter2 = readZFilter(confDict, tx_filter2_entry)

    controller_entry = confDict[station_entry]['Controller']
    controller = readController(confDict, controller_entry)

    loop_delay_size = readentry(confDict, confDict[station_entry]['loop_delay_size'])

    cav_adc_entry = confDict[station_entry]['cav_adc']
    cav_adc = readADC(confDict, cav_adc_entry)

    rfl_adc_entry = confDict[station_entry]['rfl_adc']
    rfl_adc = readADC(confDict, rfl_adc_entry)

    fwd_adc_entry = confDict[station_entry]['fwd_adc']
    fwd_adc = readADC(confDict, fwd_adc_entry)

    piezo_connect = confDict[station_entry]['piezo_connect']
    piezo_list = readList(confDict, piezo_connect, readPiezo, module_entry)

    # Create a Station instance and return
    station = Station(name, station_comp_type, amplifier, cavity, rx_filter, tx_filter1, tx_filter2, controller, loop_delay_size, cav_adc, fwd_adc, rfl_adc, piezo_list)

    return station

class Chicane:
    """ Chicane class: contains parameters specific to a Chicane configuration"""

    def __init__(self, comp_type, name, R56, T566):
        self.name = name
        self.type = comp_type

        self.R56 = R56
        self.T566 = T566

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--Module Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"

        + "R56: " + str(self.R56) + "\n"
        + "T566: " + str(self.T566) + "\n")

def readChicane(confDict, chicane_entry):

    # Read name and component type
    name = confDict[chicane_entry]['name']
    chicane_comp_type = confDict[chicane_entry]['type']

    # Read parameters
    R56 = readentry(confDict,confDict[chicane_entry]["R56"])
    T566 = readentry(confDict,confDict[chicane_entry]["T566"])

    # Create a Chicane instance and return
    chicane = Chicane(name, chicane_comp_type, R56, T566)

    return chicane

class Module:
    """ Module class: contains parameters specific to a Module configuration"""

    def __init__(self, name, comp_type, station_list, mechanical_mode_list, lp_shift):
        self.name = name
        self.type = comp_type

        self.station_list = station_list
        self.mechanical_mode_list = mechanical_mode_list
        self.lp_shift = lp_shift

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--Module Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "station_list: " + '\n'.join(str(x) for x in self.station_list)
        + "mechanical_mode_list: " + '\n'.join(str(x) for x in self.mechanical_mode_list)
        + "lp_shift: " + str(self.lp_shift) + "\n")

def readModule(confDict, module_entry):

    # Read name and component type
    name = confDict[module_entry]['name']
    module_comp_type = confDict[module_entry]['type']

    # Read the station and mechanical mode connectivity
    station_connect = confDict[module_entry]['station_connect']
    mechanical_mode_connect = confDict[module_entry]['mechanical_mode_connect']

    # Read list of stations and mechanical modes recursively
    station_list = readList(confDict, station_connect, readStation, module_entry)
    mechanical_mode_list = readList(confDict, mechanical_mode_connect, readMechMode)

    # Read lp_shift
    lp_shift = readentry(confDict,confDict[module_entry]["lp_shift"])

    # Create a Module instance and return
    module = Module(name, module_comp_type, station_list, mechanical_mode_list, lp_shift)

    return module

class Linac:
    """ Linac class: contains parameters specific to a Linac configuration"""

    def __init__(self, name, comp_type, param_dic, module_list, chicane):
        self.name = name
        self.type = comp_type

        self.f0 = param_dic["f0"]
        self.lam = param_dic["lam"]
        self.s0 = param_dic["s0"]
        self.dds_numerator = param_dic["dds_numerator"]
        self.dds_denominator = param_dic["dds_denominator"]

        self.module_list = module_list
        self.chicane = chicane

        # Need to manually propagate the value of f0 down to the Electrical Mode level
        for module in module_list:
            for station in module.station_list:
                for mode in station.cavity.elec_modes:
                    mode.omega_0_mode["value"] = 2*pi*(self.f0["value"] + mode.foffset["value"])

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--Linac Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"

        + "f0: " + str(self.f0) + "\n"
        + "lam: " + str(self.lam) + "\n"
        + "s0: " + str(self.s0) + "\n"
        + "dds_numerator: " + str(self.dds_numerator) + "\n"
        + "dds_denominator: " + str(self.dds_denominator) + "\n"
        + "chicane: " + str(self.chicane) + "\n"

        + "module_list: " + '\n'.join(str(x) for x in self.module_list))

def readLinac(confDict, linac_entry):

    # Read name and component type
    name = confDict[linac_entry]['name']
    linac_comp_type = confDict[linac_entry]['type']

    # Read and store the rest of the parameters in a dictionary
    linac_param_dic = {}

    linac_param_dic["f0"] = readentry(confDict,confDict[linac_entry]["f0"]) # Resonance frequency [Hz]
    linac_param_dic["lam"] = readentry(confDict,confDict[linac_entry]["lam"])
    linac_param_dic["lam"] = readentry(confDict,confDict[linac_entry]["lam"])
    linac_param_dic["s0"] = readentry(confDict,confDict[linac_entry]["s0"])
    linac_param_dic["iris_rad"] = readentry(confDict,confDict[linac_entry]["iris_rad"])
    linac_param_dic["dds_numerator"] = readentry(confDict,confDict[linac_entry]["dds_numerator"])
    linac_param_dic["dds_denominator"] = readentry(confDict,confDict[linac_entry]["dds_denominator"])

    # Read the module connectivity
    module_connect = confDict[linac_entry]['module_connect']

    # Read list of modules recursively
    module_list = readList(confDict, module_connect, readModule)

    # Read the chicane
    chicane_name = confDict[linac_entry]["Chicane"]
    chicane = readChicane(confDict, chicane_name)

    # Instantiate Linac object and return
    linac = Linac(name, linac_comp_type, linac_param_dic, module_list, chicane)

    return linac

class Accelerator:
    """ Accelerator class: contains parameters specific to an accelerator configuration"""

    def __init__(self, name, comp_type, bunch_rate, gun, linac_list):
        self.name = name
        self.type = comp_type
        self.bunch_rate = bunch_rate
        self.gun = gun
        self.linac_list = linac_list

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--Accelerator Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "bunch_rate: " + str(self.bunch_rate) + "\n"
        + "gun: " + str(self.gun) + "\n"
        + "linac_list: " + '\n'.join(str(x) for x in self.linac_list))

class Gun:
    """ Gun class: contains parameters specific to an Gun configuration"""

    def __init__(self, name, comp_type, Q, sz0, sd0, E):
        self.name = name
        self.type = comp_type

        self.Q = Q
        self.sz0 = sz0
        self.sd0 = sd0
        self.E = E

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--Gun Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "Q: " + str(self.Q) + "\n"
        + "sz0: " + str(self.sz0) + "\n"
        + "sd0: " + str(self.sd0) + "\n"
        + "E: " + str(self.E) + "\n")

def readGun(confDict):

    # Read name and component type
    gun_entry = confDict["Accelerator"]['gun']

    name = confDict[gun_entry]['name']
    gun_comp_type = confDict[gun_entry]['type']

    Q = readentry(confDict,confDict[gun_entry]["Q"])
    sz0 = readentry(confDict,confDict[gun_entry]["sz0"])
    sd0 = readentry(confDict,confDict[gun_entry]["sd0"])
    E = readentry(confDict,confDict[gun_entry]["E"])

    # Instantiate Linac object and return
    gun = Gun(name, gun_comp_type, Q, sz0, sd0, E)

    return gun


def readAccelerator(confDict):
    """ readAccelerator : Takes the global configuration dictionary and
    returns an Accelerator object with the configuration values filled in"""

    # Read name and component type
    name = confDict["Accelerator"]["name"]
    comp_type = confDict["Accelerator"]["type"]

    # Read other parameters
    bunch_rate = readentry(confDict,confDict["Accelerator"]["bunch_rate"])

    # Read Accelerator components (Gun + series of linacs)
    # # Read gun
    gun = readGun(confDict)

    # # Read connectivity of linacs
    linac_connect = confDict["Accelerator"]["linac_connect"]
    # # Read linacs recursively
    linac_list = readList(confDict, linac_connect, readLinac)

    # Get Accelerator instance and return
    accelerator = Accelerator(name, comp_type, bunch_rate, gun, linac_list)

    return accelerator

def readConfiguration(confDict):
    """ readConfiguration : Takes the global configuration dictionary and
    returns Simulation and Accelerator objects. This routine is to be called
    from upper level programs in order to obtain the necessary instances
    to run a full simulation """

    # Read Simulation parameters and create an instance
    simulation = readSimulation(confDict)

    # Read Accelerator parameters (and recursively all of its components) and create an instance
    accelerator = readAccelerator(confDict)

    # Return Simulation and accelerator instances
    return simulation, accelerator