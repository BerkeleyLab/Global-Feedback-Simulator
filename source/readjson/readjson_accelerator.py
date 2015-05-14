#!/usr/bin/python

#
# Accelerator-specific configuration file:
#    Defines all the Python classes involved in the simulation.
#    Parses configuration information from a dictionary and instantiates
#    Python objects with the configuration values.
#
# Classes defined:
#   Simulation
#   Noise
#   Gun
#   Linac
#   Cryomodule
#   Station
#   Controller
#   Cavity
#   ElecMode
#   MechMode
#   Piezo
#   ADC
#   ZFilter
#
# Get a Simulation Class instance in order to get a full collection of instances for a simulation run.
#
# readjson_accelerator.py
#

from readjson import readentry
from math import pi

# Define Simulation time step as global
Tstep_global = 0.0

class Synthesis:
    """ Synthesis class: contains parameters specific to a Synthesis run.
    The parameters in this class are not run-time configurable. They therefore
    need to mirror the synthesizable Verilog and are used to compute FPGA
    register settings.
    """

    def __init__(self, confDict):
        """ Synthesis Constructor:
            Input:
                confDict: Global configuration dictionary."""

         # Read name and component type (strings, no need to go through readentry)
        self.name = confDict["Synthesis"]["name"]
        self.comp_type = confDict["Synthesis"]["type"]

        # Read rest of configuration parameters
        self.n_mech_modes = readentry(confDict, confDict["Synthesis"]["n_mech_modes"])
        self.df_scale = readentry(confDict, confDict["Synthesis"]["df_scale"])

    def __str__ (self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--Synthesis Object--\n"
         + "name: " + self.name  + "\n"
         + "type: " + self.type  + "\n"
         + "n_mech_modes: " + str(self.n_mech_modes) + "\n"
         + "df_scale: " + str(self.df_scale) + "\n")

class Cavity:
    """ Cavity class: contains parameters specific to a cavity,
            including a nested list of electrical modes"""

    def __init__(self, confDict, cav_entry, cryomodule_entry):
        """ Cavity constructor: includes a recursive read of
            electrical modes in each cavity, where ElecMode objects are created for each electrical mode
            and contained as a list of ElecMode objects in the Cavity object.

            Inputs:
                confDict: Global configuration dictionary,
                cav_entry: Name of the cavity to be read (string).
                cryomodule_entry: Cryomodule entry in global dictionary in order to access the proper cryomodules

            Attributes:
                name: Cavity instance name,
                comp_type: component type (Cavity),
                elec_modes: list of ElecMode objects (one per electrical mode).

            (The mechanical mode list, wich is used as a consistency check to generate mechanical
            coupling vectors for each electrical mode)"""

        # Read name and component type
        self.name = confDict[cav_entry]['name']
        self.type = confDict[cav_entry]['type']

        # Read and store the rest of the parameters in a dictionary
        cav_param_dic = {}

        self.L = readentry(confDict,confDict[cav_entry]["L"])
        self.nom_grad = readentry(confDict,confDict[cav_entry]["nom_grad"])

        # Grab the list of electrical modes
        elec_mode_connect = confDict[cav_entry]["elec_mode_connect"]
        n_elec_modes = len(elec_mode_connect) # Number of electrical modes

        elec_mode_list = []

        ## Start of loop through electrical modes
        # Cycle through electrical modes, read parameters from global dictionary and append to list of modes.
        for m in range(n_elec_modes):
            # Take mth element of mode list
            elecMode_entry = elec_mode_connect[m]

            # Instantiate ElecMode object
            elec_mode = ElecMode(confDict, elecMode_entry, cryomodule_entry)

            # Append to list of electrical modes
            elec_mode_list.append(elec_mode)
        ## End of loop through electrical modes

        # Make the List of Electrical Modes an attribute of the Cavity object
        self.elec_modes = elec_mode_list

        # Find the fundamental mode based on coupling to the beam
        ## Criterium here is the fundamental mode being defined as that with the highest shunt impedance (R/Q)
        RoverQs = map(lambda x: x.RoverQ['value'],self.elec_modes)

        fund_index = RoverQs.index(max(RoverQs))

        # Store the index of the fundamental mode
        self.fund_index = {"value" : fund_index, "units" : "N/A", "description" : "Index of the fundamental mode in array"}

        ## Add (replicate) parameters that will be filled after object instance
        # rf_phase corresponds to Linac's phi parameter
        self.rf_phase = {"value" : 0.0, "units" : "deg", "description" : "Nominal Linac RF phase (-30 deg accelerates and puts head energy lower than tail)"}
        # design_voltage is related to the Cavity set-point (Default at max)

        self.design_voltage = {"value" : self.nom_grad["value"]*self.L["value"], "units" : "V", "description" : "Design operating Cavity voltage"}

    def __str__(self):
        """str: Convenient concatenated string output for printout"""

        return ("\n--Cavity Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "L: " + str(self.L) + "\n"
        + "nom_grad: " + str(self.nom_grad) + "\n"
        + "rf_phase: " + str(self.rf_phase) + "\n"
        + "design_voltage: " + str(self.design_voltage) + "\n"
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
                mode.foffset['value'], mode.LO_w0['value'], \
                mode.Q_0['value'], mode.Q_drive['value'], mode.Q_probe['value'], \
                self.rf_phase['value'],  mode.phase_rev['value'], mode.phase_probe['value'], \
                Tstep_global, mech_couplings, n_mech)

            acc.ElecMode_Append(elecMode_net, elecMode, idx)
            mode.C_Pointer = elecMode

        L = self.L['value']
        nom_grad = self.nom_grad['value']
        rf_phase = self.rf_phase['value']
        design_voltage = self.design_voltage['value']
        fund_index = self.fund_index['value']

        # Get a C-pointer to a Cavity structure
        cavity = acc.Cavity_Allocate_New(elecMode_net, n_modes, L, nom_grad, \
            rf_phase, design_voltage, \
            fund_index)

        self.C_Pointer = cavity

        return cavity

    def Get_State_Pointer(self):
        import accelerator as acc

        cavity_state = acc.Cavity_State()
        acc.Cavity_State_Allocate(cavity_state, self.C_Pointer)

        self.State = cavity_state

        return cavity_state

class ElecMode:
    def __init__(self, confDict, elecMode_entry, cryomodule_entry):
        """ ElecMode class: contains parameters specific to an electrical mode,
            including a dictionary specifying the mechanical couplings.
            Note the absence of a readElecMode method, the process for parsing
            the global configuration dictionary and creating ElecMode objects
            is done recursively in Cavity.readCavity(...)"""

        # Read component name and type
        self.name = confDict[elecMode_entry]['name']
        self.type = confDict[elecMode_entry]['type']

        # Identifier for mode (e.g pi, 8pi/9, etc.)
        self.mode_name = confDict[elecMode_entry]['mode_name']

        # Read rest of parameters and store in dictionary
        self.RoverQ = readentry(confDict,confDict[elecMode_entry]["RoverQ"])
        self.foffset = readentry(confDict,confDict[elecMode_entry]["foffset"])
        self.peakV = readentry(confDict,confDict[elecMode_entry]["peakV"])
        self.Q_0 = readentry(confDict,confDict[elecMode_entry]["Q_0"])
        self.Q_drive = readentry(confDict,confDict[elecMode_entry]["Q_drive"])
        self.Q_probe = readentry(confDict,confDict[elecMode_entry]["Q_probe"])
        self.phase_rev = readentry(confDict,confDict[elecMode_entry]["phase_rev"])
        self.phase_probe = readentry(confDict,confDict[elecMode_entry]["phase_probe"])

        # Read dictionary of couplings from global configuration dictionary
        mech_couplings = readentry(confDict,confDict[elecMode_entry]["mech_couplings"]["value"])
        # Get a coupling list of length M (number of mechanical modes),
        # filled with 0s if no coupling is specified by user
        self.mech_couplings_list = readCouplings(confDict, mech_couplings, cryomodule_entry)

        ## Add (replicate) a parameter that will be filled after object instance
        self.LO_w0 = {"value" : 0.0, "units" : "rad/s", "description" : "Linac's Nominal resonance angular frequency"}

    def __str__(self):
        """str: Convenient concatenated string output for printout"""

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

    def Compute_ElecMode(self, Tstep, rf_phase):
        import numpy as np

        # Initialize an empty list to return
        modes_out = []

        beam_phase = rf_phase

        mode_name = self.mode_name
        LO_w0 = self.LO_w0['value']
        foffset = self.foffset['value']
        w0 = LO_w0 + 2.0*pi*foffset
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
    """ MechMode Class: contains parameters specific to a mechanical mode.
        Information concerning couplings with electrical modes and Piezos is
        contained in ElecMode and Piezo objects respectively."""

    def __init__(self, confDict, mechMode_entry):
        """ MechMode Constructor:
                Inputs:
                    confDict: Global configuration dictionary,
                    mech_mode_entry: Name of the mechanical mode to be read (string)"""

        # Read name and component type
        self.name = confDict[mechMode_entry]['name']
        self.type = confDict[mechMode_entry]['type']

        # Read the rest of the configuration parameters and store in a dictionary
        self.f0 = readentry(confDict,confDict[mechMode_entry]["f0"])
        self.Q = readentry(confDict,confDict[mechMode_entry]["Q"])
        self.full_scale = readentry(confDict,confDict[mechMode_entry]["full_scale"])

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--MechMode Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "f0: " + str(self.f0) + "\n"
        + "Q: " + str(self.Q) + "\n"
        + "full_scale: " + str(self.full_scale) + "\n")

    def Get_C_Pointer(self):
        import accelerator as acc

        # Grab attributes from object
        f0 = self.f0['value']
        Q = self.Q['value']
        k = 1.0

        # Allocate Memory for C struct
        mechMode = acc.MechMode_Allocate_New(f0, Q, k, Tstep_global);

        self.C_Pointer = mechMode

        # Return C Pointer
        return mechMode

class Piezo:
    """ Piezo Class: contains couplings between the Piezo and each
    one of the mechanical modes (MechMode instances)."""

    def __init__(self, confDict, piezo_entry, cryomodule_entry):
        """ Piezo Constructor:
                 Inputs:
                    confDict: Global configuration dictionary,
                    piezo_entry: Name of the Piezo to be read (string).
                    cryomodule_entry: Cryomodule entry in global dictionary in order to access
                        the proper cryomodule's mechanical mode list, which is used as a consistency
                        check to generate mechanical coupling vectors for each Piezo."""

        # Read name and component type
        self.name = confDict[piezo_entry]['name']
        self.type = confDict[piezo_entry]['type']

        # Read rest of parameters
        self.VPmax = readentry(confDict,confDict[piezo_entry]["VPmax"])

        # Read dictionary of couplings from global configuration dictionary
        mech_couplings = readentry(confDict,confDict[piezo_entry]["mech_couplings"]["value"])

        # Check consistency of coupling entries with the list of mechanical modes,
        # and get a coupling dictionary of length M (number of mechanical modes)
        self.mech_couplings_list = readCouplings(confDict, mech_couplings, cryomodule_entry)

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--Piezo Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "VPmax: " + str(self.VPmax) + "\n"
        + "mech_couplings_list: " + str(self.mech_couplings_list))


def readCouplings(confDict, mech_couplings, cryomodule_entry):
    """ readCouplings: Takes the global configuration dictionary and a dictionary containing non-zero
    values for the mechanical couplings (i.e. coupling between electrical modes or piezos and mechanical modes).
    The cryomodule_entry input is necessary in order to access the proper cryomodule's mechanical mode list.
    The length of the coupling vector used in the simulation code must be equal to the number of mechanical modes (M).
    The mech_couplings input is supposed to contain non-zero values from the configuration file,
    and readCouplings always returns a dictionary of M elements, where the couplings not specified in
    the input file are filled with 0s.
    Inputs:
        confDict: Global configuration dictionary,
        mech_couplings: dictionary containing couplings defined in the configuration file (0 to M elements).
        cryomodule_entry: Cryomodule entry in global dictionary in order to access the proper cryomodule's mechanical mode list.
    Output:
        mech_couplings_list: ordered list containing mechanical couplings for an electrical mode or piezo (length M).
            Order corresponds to the order of appearance of the mechanical mode in mech_net."""

    # Grab the full list of mechanical modes in the Cryomodule
    mech_net = confDict[cryomodule_entry]["mechanical_mode_connect"]

    # Make an ordered list of size M (where M is the total number of mechanical modes, length of mech_net)
    # Fill with 0s if no coupling is specified in mech_couplings by the user
    mech_couplings_list = [mech_couplings[m] if m in mech_couplings else 0.0 for m in mech_net]

    return mech_couplings_list


def readList(confDict, list_in, constructor, cryomodule_entry=None):
    """ readList: Generic function to read list of components.
    Takes the global configuration dictionary, cycles through the list of components
    (list of names, list_in), uses the names to identify the configuration entries in
    the global dictionary, calls the proper Constructor for each component (constructor),
    and returns a list of instances (objects).
    Inputs:
        confDict: Global configuration dictionary,
        list_in: list of components to cycle through (list of strings),
        constructor: name of Constructor for the component,
        cryomodule_entry: necessary in some cases in order to pass along cryomodule entry information,
            needed by readStation and Piezo in order to find the mechanical modes in
            their corresponding Cryomodule.
    Output:
        list_out: List of component objects"""

    # Create empty list for component instances
    list_out = []

    # Cycle through list of component names
    for k in range(len(list_in)):
        # Read component configuration and create component instance
        if cryomodule_entry == None:
            component = constructor(confDict, list_in[k])
        else:
            component = constructor(confDict, list_in[k], cryomodule_entry)
        # Append object to the component list
        list_out.append(component)

    # Return list
    return list_out

class Controller:
    """ Controller class: contains parameters specific to a Controller configuration"""

    def __init__(self, confDict, controller_entry):
        # Read name and component type
        name = confDict[controller_entry]['name']
        self.type = confDict[controller_entry]['type']

        # Read the rest of parameters
        self.stable_gbw = readentry(confDict,confDict[controller_entry]["stable_gbw"])

    def __str__(self):
        """str: Convenient concatenated string output for printout"""

        return ("\n--Controller Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "stable_gbw: " + str(self.stable_gbw) + "\n")

class ZFilter:
    """ ZFilter class: contains parameters specific to a Filter configuration"""

    def __init__(self, confDict, zfilter_entry):

        # Read name and component type
        self.name = confDict[zfilter_entry]['name']
        self.type = confDict[zfilter_entry]['type']

        # Read the rest of parameters
        self.order = readentry(confDict,confDict[zfilter_entry]["order"])
        self.nmodes = readentry(confDict,confDict[zfilter_entry]["nmodes"])
        self.poles = readentry(confDict,confDict[zfilter_entry]["poles"])

    def __str__(self):
        """str: Convenient concatenated string output for printout"""

        return ("\n--ZFilter Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "order: " + str(self.order) + "\n"
        + "nmodes: " + str(self.nmodes) + "\n"
        + "poles: " + str(self.poles) + "\n")

class ADC:
    """ ADC class: contains parameters specific to a ADC configuration"""

    def __init__(self, confDict, adc_entry):

        # Read name and component type
        self.name = confDict[adc_entry]['name']
        self.type = confDict[adc_entry]['type']

        # Read the rest of parameters
        self.adc_max = readentry(confDict,confDict[adc_entry]["adc_max"])
        self.adc_off = readentry(confDict,confDict[adc_entry]["adc_off"])
        self.noise_psd = readentry(confDict,confDict[adc_entry]["noise_psd"])

    def __str__(self):
        """str: Convenient concatenated string output for printout"""

        return ("\n--ADC Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "adc_max: " + str(self.adc_max) + "\n"
        + "adc_off: " + str(self.adc_off) + "\n"
        + "noise_psd: " + str(self.noise_psd) + "\n")

class Amplifier:
    """ Amplifier class: contains parameters specific to a Amplifier configuration"""

    def __init__(self, confDict, amplifier_entry):

        # Read name and component type
        self.name = confDict[amplifier_entry]['name']
        self.type = confDict[amplifier_entry]['type']

        # Read the rest of parameters
        self.PAmax = readentry(confDict,confDict[amplifier_entry]["PAmax"])
        self.PAbw = readentry(confDict,confDict[amplifier_entry]["PAbw"])
        self.Clip = readentry(confDict,confDict[amplifier_entry]["Clip"])
        self.top_drive = readentry(confDict,confDict[amplifier_entry]["top_drive"])

    def __str__(self):
        """str: Convenient concatenated string output for printout"""

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

class Station:
    """ Station class: contains parameters specific to a Station configuration"""

    def __init__(self, confDict, station_entry, cryomodule_entry):

        # Read name and component type
        self.name = confDict[station_entry]['name']
        self.type = confDict[station_entry]['type']

        # Read all the station components
        amplifier_entry = confDict[station_entry]['Amplifier']
        self.amplifier = Amplifier(confDict, amplifier_entry)

        cavity_entry = confDict[station_entry]['Cavity']
        self.cavity = Cavity(confDict, cavity_entry, cryomodule_entry)

        rx_filter_entry = confDict[station_entry]['Rx_filter']
        self.rx_filter = ZFilter(confDict, rx_filter_entry)

        tx_filter1_entry = confDict[station_entry]['Tx_filter1']
        self.tx_filter1 = ZFilter(confDict, tx_filter1_entry)

        tx_filter2_entry = confDict[station_entry]['Tx_filter2']
        self.tx_filter2 = ZFilter(confDict, tx_filter2_entry)

        controller_entry = confDict[station_entry]['Controller']
        self.controller = Controller(confDict, controller_entry)

        self.loop_delay_size = readentry(confDict, confDict[station_entry]['loop_delay_size'])

        cav_adc_entry = confDict[station_entry]['cav_adc']
        self.cav_adc = ADC(confDict, cav_adc_entry)

        rfl_adc_entry = confDict[station_entry]['rfl_adc']
        self.rfl_adc = ADC(confDict, rfl_adc_entry)

        fwd_adc_entry = confDict[station_entry]['fwd_adc']
        self.fwd_adc = ADC(confDict, fwd_adc_entry)

        piezo_connect = confDict[station_entry]['piezo_connect']
        self.piezo_list = readList(confDict, piezo_connect, Piezo, cryomodule_entry)

        self.N_Stations = confDict[station_entry]['N_Stations']

        ## Add (replicate) parameters that will be filled after object instance
        self.max_voltage = {"value" : 0.0, "units" : "V", "description" : "Maximum accelerating voltage"}

    def __str__(self):
        """str: Convenient concatenated string output for printout"""

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
        + "N_Stations: " + str(self.N_Stations) + "\n"
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

        rf_station = acc.RF_Station()
        acc.RF_Station_Allocate_In(rf_station, Tstep_global, Clip, PAmax, PAscale, p_TRF1, p_TRF2, p_RXF, cavity_pointer, stable_gbw, FPGA_out_sat, loop_delay_size)

        self.C_Pointer = rf_station

        return rf_station

    def Get_State_Pointer(self):
        import accelerator as acc

        rf_state = acc.RF_State()

        acc.RF_State_Allocate(rf_state, self.C_Pointer)
        self.State = rf_state

        return rf_state

class Cryomodule:
    """ Cryomodule class: contains parameters specific to a Cryomodule configuration"""

    def __init__(self, confDict, cryomodule_entry):
        # Read name and component type
        self.name = confDict[cryomodule_entry]['name']
        self.type = confDict[cryomodule_entry]['type']

        # Read the station and mechanical mode connectivity
        station_connect = confDict[cryomodule_entry]['station_connect']
        mechanical_mode_connect = confDict[cryomodule_entry]['mechanical_mode_connect']

        # Read list of stations and mechanical modes recursively
        self.station_list = readList(confDict, station_connect, Station, cryomodule_entry)
        self.mechanical_mode_list = readList(confDict, mechanical_mode_connect, MechMode)

        # Read lp_shift
        self.lp_shift = readentry(confDict,confDict[cryomodule_entry]["lp_shift"])

    def __str__(self):
        """str: Convenient concatenated string output for printout"""

        return ("\n--Cryomodule Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "station_list: " + '\n'.join(str(x) for x in self.station_list)
        + "mechanical_mode_list: " + '\n'.join(str(x) for x in self.mechanical_mode_list)
        + "lp_shift: " + str(self.lp_shift) + "\n")

    def Get_C_Pointer(self):

        import accelerator as acc

        # First count number of Stations and Mechanical Modes and Allocate Arrays
        n_Stations = len(self.station_list)
        n_MechModes = len(self.mechanical_mode_list)

        # Allocate memory for RF Station and Mechanical mode Arrays
        rf_station_net = acc.RF_Station_Allocate_Array(n_Stations)
        mechMode_net = acc.MechMode_Allocate_Array(n_MechModes)

        # Allocate each RF Station and append it to the rf_station_net
        for idx, rf_station in enumerate(self.station_list):
            RF_Station_C_Pointer = rf_station.Get_C_Pointer()
            acc.RF_Station_Append(rf_station_net, RF_Station_C_Pointer, idx)

        # Allocate each MechMode and append it to the mechMode_net
        for idx, mechMode in enumerate(self.mechanical_mode_list):
            mechMode_C_Pointer = mechMode.Get_C_Pointer()
            acc.MechMode_Append(mechMode_net, mechMode_C_Pointer, idx)

        # Instantiate Cryomodule C structure
        cryomodule = acc.Cryomodule()
        # Fill in Cryomodule C data structure
        acc.Cryomodule_Allocate_In(cryomodule, rf_station_net, n_Stations, mechMode_net, n_MechModes)

        self.C_Pointer = cryomodule

        # Return C Pointer for Cryomodule and lists of pointers
        return cryomodule

    def Get_State_Pointer(self, cryo_state=None):
        import accelerator as acc

        # Get C pointer to Cryomodule_State C struct,
        # (if it has not been allocated yet)
        if(cryo_state==None):
            cryo_state = acc.Cryomodule_State()
            # Allocate Memory for cryo_state
            acc.Cryomodule_State_Allocate(cryo_state, self.C_Pointer)

        # Get Pointers to RF States
        for idx, station in enumerate(self.station_list):
            station.State = acc.Get_RF_State(cryo_state, idx)

        # Get Pointers to MechMode States
        for idx, mechMode in enumerate(self.mechanical_mode_list):
            mechMode.State = acc.Get_MechMode_State(cryo_state, idx)

        self.State = cryo_state

        # Return State C pointer
        return cryo_state

class Linac:
    """ Linac class: contains parameters specific to a Linac configuration"""

    def __init__(self, confDict, linac_entry):

        import numpy as np

        # Read name and component type
        name = confDict[linac_entry]['name']
        self.type = confDict[linac_entry]['type']

        self.f0 = readentry(confDict,confDict[linac_entry]["f0"])
        self.E = readentry(confDict,confDict[linac_entry]["E"])
        self.phi = readentry(confDict,confDict[linac_entry]["phi"])
        self.phi['value'] = self.phi['value']*np.pi/180    # Convert degrees to radians
        self.s0 = readentry(confDict,confDict[linac_entry]["s0"])
        self.iris_rad = readentry(confDict,confDict[linac_entry]["iris_rad"])
        self.R56 = readentry(confDict,confDict[linac_entry]["R56"])
        self.dds_numerator = readentry(confDict,confDict[linac_entry]["dds_numerator"])
        self.dds_denominator = readentry(confDict,confDict[linac_entry]["dds_denominator"])

        # Read the cryomodule connectivity
        cryomodule_connect = confDict[linac_entry]['cryomodule_connect']

        # Read list of modules recursively
        self.cryomodule_list = readList(confDict, cryomodule_connect, Cryomodule)

        # Add parameters that will be filled after object instance
        self.dE = {"value" : 0.0, "units" : "eV", "description" : "Energy increase in Linac (final minus initial Energy)"}
        self.max_voltage = {"value" : 0.0, "units" : "V", "description" : "Maximum Accelerating Voltage"}
        self.N_Stations = {"value" : 0.0, "units" : "N/A", "description" : "Total number of RF Stations in Linac"}
        self.L = {"value" : 0.0, "units" : "m", "description" : "Total Linac Length"}

        # Some parameters are deduced from others
        # RF wavelength deduced from f0
        c = 2.99792458e8    # Speed of light [m/s]
        self.lam = {"value" : c/self.f0['value'], "units" : "m", "description" : "RF wavelength for each linac (Sband=0.105m, Xband=0.02625m)"}
        # T566 deduced from R56 (small angle approximation)
        self.T566 = {"value" : -1.5*self.R56['value'], "units" : "m", "description" : "Nominal T566 (always >= 0)"}

        # Need to manually propagate values down to the Cavity and Electrical Mode level
        for cryomodule in self.cryomodule_list:
            for station in cryomodule.station_list:
                # Indicate each cavity the nominal beam phase for the Linac
                station.cavity.rf_phase["value"] = self.phi["value"]

                # Calculate each RF Station's maximum accelerating voltage
                N_Stations = station.N_Stations["value"]
                nom_grad = station.cavity.nom_grad["value"]
                L = station.cavity.L["value"]
                station.max_voltage["value"] = N_Stations*nom_grad*L

                # Add to the Linac's total accelerating voltage, number of RF Stations and Length
                self.max_voltage["value"] += station.max_voltage["value"]
                self.N_Stations["value"] += N_Stations
                self.L["value"] += L*N_Stations

                # Indicate each Electrical Eigenmode the nominal LO frequency for the Linac
                for mode in station.cavity.elec_modes:
                    mode.LO_w0["value"] = 2*pi*self.f0["value"]

    def __str__(self):
        """str: Convenient concatenated string output for printout"""

        return ("\n--Linac Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"

        + "f0: " + str(self.f0) + "\n"
        + "E: " + str(self.E) + "\n"
        + "phi: " + str(self.phi) + "\n"
        + "lam: " + str(self.lam) + "\n"
        + "s0: " + str(self.s0) + "\n"
        + "iris_rad: " + str(self.iris_rad) + "\n"
        + "R56: " + str(self.R56) + "\n"
        + "T566: " + str(self.T566) + "\n"
        + "dE: " + str(self.dE) + "\n"
        + "max_voltage: " + str(self.max_voltage) + "\n"
        + "L: " + str(self.L) + "\n"
        + "dds_numerator: " + str(self.dds_numerator) + "\n"
        + "dds_denominator: " + str(self.dds_denominator) + "\n"

        + "cryomodule_list: " + '\n'.join(str(x) for x in self.cryomodule_list))

    def Get_C_Pointer(self):

        import accelerator as acc

        # First count number of Stations and Mechanical Modes and Allocate Arrays
        n_Cryos = len(self.cryomodule_list)

        # Allocate memory for array of Cryomodules
        cryo_net = acc.Cryomodule_Allocate_Array(n_Cryos)

        # Allocate each Cryomodule and append it to the cryo_net
        for idx, cryo in enumerate(self.cryomodule_list):
            Cryo_C_Pointer = cryo.Get_C_Pointer()
            acc.Cryomodule_Append(cryo_net, Cryo_C_Pointer, idx)

        # Instantiate Linac C structure
        linac = acc.Linac()
        # Fill in Linac C data structure
        acc.Linac_Allocate_In(linac, cryo_net, n_Cryos,
            self.dE["value"], self.R56["value"], self.T566["value"], \
            self.phi["value"], self.lam["value"], self.s0["value"], \
            self.iris_rad["value"], self.L["value"])

        self.C_Pointer = linac

        # Return Linac C Pointer
        return linac

    def Get_State_Pointer(self, linac_state=None):
        import accelerator as acc

        # Get C pointer to Linac_state C struct,
        # (if it has not been allocated yet)
        if(linac_state==None):
            # Get C pointer to Linac_State C struct
            linac_state = acc.Linac_State()
            # Allocate Memory for linac_state
            acc.Linac_State_Allocate(linac_state, self.C_Pointer)

        # Get Pointers to Cryomodule States
        for idx, cryo in enumerate(self.cryomodule_list):
            cryo_state = acc.Get_Cryo_State(linac_state, idx)
            cryo.Get_State_Pointer(cryo_state)

        # Return State C pointer
        return linac_state

class Simulation:
    """ Simulation class: contains parameters specific to a Simulation run,
    as well as all parameters in the Accelerator configuration. This Class is
    to be instantiated from upper level programs in order to obtain all necessary instances
    to run a full simulation"""

    def __init__(self, confDict):
        """Simulation Constructor:
            Inputs:
                confDict: Global configuration dictionary."""

        import numpy as np

        # Read name and component type
        self.name = confDict["Accelerator"]["name"]
        self.type = confDict["Accelerator"]["type"]

        # Read rest of configuration parameters
        self.Tstep = readentry(confDict, confDict["Simulation"]["Tstep"])
        self.time_steps = readentry(confDict, confDict["Simulation"]["time_steps"])
        self.nyquist_sign = readentry(confDict,confDict["Simulation"]["nyquist_sign"])

        # Check if simulation dictionary has a Synthesis entry, and if so parse it
        if confDict["Simulation"].has_key("Synthesis"):
            self.synthesis = Synthesis(confDict["Simulation"])
        else:
            self.synthesis = None

        # Accelerator parameters
        self.bunch_rate = readentry(confDict,confDict["Accelerator"]["bunch_rate"])


        # Read Noise Sources
        self.noise_srcs = Noise(confDict)

        # Read Accelerator components (Gun + series of linacs)
        # # Read gun
        self.gun = Gun(confDict)
        Egun = self.gun.E['value'] # Gun exit Energy

        # # Read connectivity of linacs
        linac_connect = confDict["Accelerator"]["linac_connect"]
        # # Read linacs recursively
        self.linac_list = readList(confDict, linac_connect, Linac)

        # Now that the Array of Linacs has been instantiated and their final Energies are known,
        # fill in the Energy increase parameter for each Linac.
        ## Start with the Energy out of the Gun
        Elast = Egun
        # Iterate over Linacs
        for linac in self.linac_list:
            # Energy increase is the difference between Linac's final and initial Energy
            linac.dE["value"] = linac.E["value"] - Elast
            Elast = linac.E["value"]

            # Check if Energy increase is compatible with Linac configuration
            sp_ratio = linac.dE["value"]/np.cos(linac.phi["value"])/linac.max_voltage["value"]
            if np.abs(sp_ratio) > 1.0:
                error_1 = "Linac "+ linac.name + ": Energy increase higher than tolerated:\n"
                error_2 = "\tEnergy increase requested = %.2f eV"%linac.dE["value"]
                error_3 = "\tAt Beam phase = %.2f deg"%linac.phi["value"]
                error_4 = "\tand Maximum Accelerating Voltage = %.2f V"%linac.max_voltage["value"]
                error_text = error_1 + error_2 + error_3 + error_4
                raise Exception(error_text)
            else:
                # Once Energy increase is known, calculate set-point for each RF Station
                for cryo in linac.cryomodule_list:
                    for station in cryo.station_list:
                        station_max_voltage = station.max_voltage['value']
                        station.cavity.design_voltage['value'] = station_max_voltage*sp_ratio

        # Add a parameter which was not parsed from configuration but deduced
        self.E = {"value" : Elast, "units" : "eV", "description" : "Final Accelerator Energy"}

        # Assign Tstep to the global variable
        global Tstep_global
        Tstep_global = self.Tstep['value']

    def __str__(self):
        """str: Convenient concatenated string output for printout"""

        return ("\n--Simulation Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "Tstep: " + str(self.Tstep) + "\n"
        + "time_steps: " + str(self.time_steps) + "\n"
        + "nyquist_sign: " + str(self.nyquist_sign) + "\n"
        + "synthesis: " + str(self.synthesis) + "\n"
        + "bunch_rate: " + str(self.bunch_rate) + "\n"
        + "noise_srcs: " + str(self.noise_srcs) + "\n"
        + "E: " + str(self.E) + "\n"
        + "gun: " + str(self.gun) + "\n"

        + "linac_list: " + '\n'.join(str(x) for x in self.linac_list))

    def Get_C_Pointer(self):

        import accelerator as acc

        # First count number of Linacs and Allocate Arrays
        n_linacs = len(self.linac_list)

        # Allocate memory for array of Linacs
        linac_net = acc.Linac_Allocate_Array(n_linacs)

        # Allocate memory for Gun
        gun_C_Pointer = self.gun.Get_C_Pointer()

        # Allocate each Linac and append it to the linac_net
        for idx, linac in enumerate(self.linac_list):
            Linac_C_Pointer = linac.Get_C_Pointer()
            acc.Linac_Append(linac_net, Linac_C_Pointer, idx)

        # Instantiate Accelerator C structure
        sim = acc.Simulation()

        # Fill in Linac C data structure
        acc.Sim_Allocate_In(sim, \
            self.Tstep['value'], self.time_steps['value'], \
            gun_C_Pointer, linac_net, n_linacs)

        self.C_Pointer = sim

        # Return C Pointer for Simulation
        return sim

    def Get_State_Pointer(self):
        import accelerator as acc

        # Allocate memory for Noise State
        noise_State_Pointer = self.noise_srcs.Get_State_Pointer()

        sim_state = acc.Simulation_State()
        acc.Sim_State_Allocate(sim_state, self.C_Pointer, noise_State_Pointer)

        self.State = sim_state

        return sim_state

class Gun:
    """ Gun class: contains parameters specific to an Gun configuration"""

    def __init__(self, confDict):

        # Read name and component type
        gun_entry = confDict["Accelerator"]['gun']

        self.name = confDict[gun_entry]['name']
        self.type = confDict[gun_entry]['type']

        self.Q = readentry(confDict,confDict[gun_entry]["Q"])
        self.sz0 = readentry(confDict,confDict[gun_entry]["sz0"])
        self.sd0 = readentry(confDict,confDict[gun_entry]["sd0"])
        self.E = readentry(confDict,confDict[gun_entry]["E"])

    def __str__(self):
        """str: Convinient concatenated string output for printout"""

        return ("\n--Gun Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "Q: " + str(self.Q) + "\n"
        + "sz0: " + str(self.sz0) + "\n"
        + "sd0: " + str(self.sd0) + "\n"
        + "E: " + str(self.E) + "\n")

    def Get_C_Pointer(self):

        import accelerator as acc

        gun = acc.Gun()
        acc.Gun_Allocate_In(gun, self.E['value'], self.sz0['value'], self.sd0['value'], self.Q['value']);

        self.C_Pointer = gun

        return gun

class Noise:
    """ Noise class: contains configuration regarding the correlated noise sources in the Accelerator."""

    def __init__(self, confDict):

        # Read name and component type
        noise_entry = confDict['Noise']

        self.name = confDict[noise_entry]['name']
        self.type = confDict[noise_entry]['type']

        self.dQ_Q = readentry(confDict,confDict[noise_entry]["dQ_Q"])
        self.dtg = readentry(confDict,confDict[noise_entry]["dtg"])
        self.dE_ing = readentry(confDict,confDict[noise_entry]["dE_ing"])
        self.dsig_z = readentry(confDict,confDict[noise_entry]["dsig_z"])
        self.dsig_E = readentry(confDict,confDict[noise_entry]["dsig_E"])
        self.dchirp = readentry(confDict,confDict[noise_entry]["dchirp"])

    def __str__(self):
        """str: Convenient concatenated string output for printout"""

        return ("\n--Noise Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "dQ_Q: " + str(self.dQ_Q) + "\n"
        + "dtg: " + str(self.dtg) + "\n"
        + "dE_ing: " + str(self.dE_ing) + "\n"
        + "dsig_z: " + str(self.dsig_z) + "\n"
        + "dsig_E: " + str(self.dsig_E) + "\n"
        + "dchirp: " + str(self.dchirp) + "\n")

    def Get_State_Pointer(self):

        import accelerator as acc

        noise_srcs = acc.Noise_Srcs()

        type_net = acc.intArray_frompointer(noise_srcs.type)
        setting_net = acc.double_Array_frompointer(noise_srcs.settings)

        # Dictionaries for parameters and indices
        field_dict = ['dQ_Q', 'dtg', 'dE_ing', 'dsig_z', 'dsig_E', 'dchirp']
        type_dict = {'None':0, 'White':1, 'Sine':2, 'Chirp':3, 'Step':4}

        # Loop over all sources of noise in the dictionary,
        # and see if the user specified them
        for i in xrange(len(field_dict)):
            key = field_dict[i]
            # Make sure the key corresponds to a field that the user specified
            try:
                entry = self.key
                # Figure out its type
                type_net[i] = type_dict[entry['Type']]
                # Write the settings for the noise into the C structure
                settings = entry['Settings']
                if not isinstance(settings, list):
                    setting_net[acc.N_NOISE_SETTINGS*i] = float(settings)
                else:
                    for k in xrange(len(settings)):
                        setting_net[acc.N_NOISE_SETTINGS*i+k] = float(settings[k])
            except:
                # Default to no noise and whine to the user
                type_net[i] = 0

        # Return the C Pointer
        return noise_srcs

