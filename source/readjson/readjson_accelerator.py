"""
Accelerator-specific configuration file:

Defines all the Python classes involved in the simulation.
Parses configuration information from a dictionary and instantiates
Python objects with the configuration values.
Get a Simulation Class instance in order to get a full collection of instances for a simulation run.
"""

from readjson import readentry
from math import pi

## Define Simulation time step as global
Tstep_global = 0.0

class Synthesis:
    """ Contains parameters specific to a Synthesis run.
    The parameters in this class are not run-time configurable. They therefore need
    to mirror the synthesizable Verilog and are used to compute FPGA register settings.
    """

    def __init__(self, confDict):

        """ Constructor
        Input:
            - confDict: Global configuration dictionary.
        """

        ## Instance name
        self.name = confDict["Synthesis"]["name"]
        ## Instance type
        self.type = confDict["Synthesis"]["type"]

        # Read rest of configuration parameters
        ## Number of Mechanical modes instantiated in FPGA
        self.n_mech_modes = readentry(confDict, confDict["Synthesis"]["n_mech_modes"])
        ## FPGA scaling factor
        self.df_scale = readentry(confDict, confDict["Synthesis"]["df_scale"])

    def __str__ (self):
        """Convenient concatenated string output for printout."""

        return ("\n--Synthesis Object--\n"
         + "name: " + self.name  + "\n"
         + "type: " + self.type  + "\n"
         + "n_mech_modes: " + str(self.n_mech_modes) + "\n"
         + "df_scale: " + str(self.df_scale) + "\n")

class Cavity:
    """ Contains parameters specific to a cavity, including a nested list of electrical modes"""

    def __init__(self, confDict, cav_entry, cryomodule_entry):
        """
        Cavity constructor: includes a recursive read of electrical modes in each cavity,
        where ElecMode objects are created for each electrical mode and contained as a list of ElecMode objects in the Cavity object.

        Inputs:
            - confDict: Global configuration dictionary,
            - cav_entry: Name of the cavity to be read (string),
            - cryomodule_entry: Cryomodule entry in global dictionary in order to access the proper cryomodules.

        (The mechanical mode list is used as a consistency check to generate mechanical coupling vectors for each electrical mode).
        """

        ## Instance name
        self.name = confDict[cav_entry]['name']
        ## Instance type
        self.type = confDict[cav_entry]['type']

        # Read and store the rest of the parameters in a dictionary
        cav_param_dic = {}

        ## cavity electrical length [m]
        self.L = readentry(confDict,confDict[cav_entry]["L"])

        ## cavity Nominal gradient [V/m]
        self.nom_grad = readentry(confDict,confDict[cav_entry]["nom_grad"])

        # Grab the list of electrical modes
        elec_mode_connect = confDict[cav_entry]["elec_mode_connect"]
        n_elec_modes = len(elec_mode_connect) # Number of electrical modes

        elec_mode_list = []

        # Start of loop through electrical modes
        # Cycle through electrical modes, read parameters from global dictionary and append to list of modes.
        for m in range(n_elec_modes):
            # Take mth element of mode list
            elecMode_entry = elec_mode_connect[m]

            # Instantiate ElecMode object
            elec_mode = ElecMode(confDict, elecMode_entry, cryomodule_entry)

            # Append to list of electrical modes
            elec_mode_list.append(elec_mode)
        # End of loop through electrical modes

        # Make the List of Electrical Modes an attribute of the Cavity object
        ## List of ElecMode objects (one per electrical mode)
        self.elec_modes = elec_mode_list

        # Find the fundamental mode based on coupling to the beam
        # Criterium here is the fundamental mode being defined as that with the highest shunt impedance (R/Q)
        RoverQs = map(lambda x: x.RoverQ['value'],self.elec_modes)

        fund_index = RoverQs.index(max(RoverQs))

        ## Index of the fundamental mode
        self.fund_index = {"value" : fund_index, "units" : "N/A", "description" : "Index of the fundamental mode in array"}

        # Add (replicate) parameters that will be filled after object instance
        ## cavity nominal phase with respect to the beam [deg]
        self.rf_phase = {"value" : 0.0, "units" : "deg", "description" : "Nominal Linac RF phase (-30 deg accelerates and puts head energy lower than tail)"}

        ## Related to the cavity set-point (Default at max) [V]
        self.design_voltage = {"value" : self.nom_grad["value"]*self.L["value"], "units" : "V", "description" : "Design operating Cavity voltage"}

    def __str__(self):
        """Convenient concatenated string output for printout."""

        return ("\n--Cavity Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "L: " + str(self.L) + "\n"
        + "nom_grad: " + str(self.nom_grad) + "\n"
        + "rf_phase: " + str(self.rf_phase) + "\n"
        + "design_voltage: " + str(self.design_voltage) + "\n"
        + "electrical modes: " + '\n'.join(str(x) for x in self.elec_modes))

    def Get_C_Pointer(self):
        """ Return reference to the SWIG-wrapped C structure. """

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

        ## Pointer to the SWIG-wrapped C structure
        self.C_Pointer = cavity

        return cavity

    def Get_State_Pointer(self):
        """ Return reference to the SWIG-wrapped State C structure. """
        import accelerator as acc

        cavity_state = acc.Cavity_State()
        acc.Cavity_State_Allocate(cavity_state, self.C_Pointer)

        ## Pointer to the SWIG-wrapped State C structure
        self.State = cavity_state

        return cavity_state

class ElecMode:
    def __init__(self, confDict, elecMode_entry, cryomodule_entry):
        """
        Contains parameters specific to an electrical mode, including a dictionary specifying the mechanical couplings.
        Note the absence of a readElecMode method, the process for parsing the global configuration dictionary and
        creating ElecMode objects is done recursively in the Cavity constructor.
        """

        ## Instance name
        self.name = confDict[elecMode_entry]['name']
        ## Instance type
        self.type = confDict[elecMode_entry]['type']

        ## Identifier for mode (e.g pi, 8pi/9, etc.)
        self.mode_name = confDict[elecMode_entry]['mode_name']

        # Read rest of parameters and store in dictionary
        ## Mode's (R/Q) [Ohms]
        self.RoverQ = readentry(confDict,confDict[elecMode_entry]["RoverQ"])
        ## Mode's frequency offset (with respect to the RF reference frequency) [Hz]
        self.foffset = readentry(confDict,confDict[elecMode_entry]["foffset"])
        ## Scaling factor for FPGA double precision to fixed point conversion of voltages
        self.peakV = readentry(confDict,confDict[elecMode_entry]["peakV"])
        ## Represents losses in the cavity walls
        self.Q_0 = readentry(confDict,confDict[elecMode_entry]["Q_0"])
        ## Represents coupling to the input coupler
        self.Q_drive = readentry(confDict,confDict[elecMode_entry]["Q_drive"])
        ## Represents coupling to the field probe
        self.Q_probe = readentry(confDict,confDict[elecMode_entry]["Q_probe"])
        ## Phase shift between Cavity cells and reverse ADC
        self.phase_rev = readentry(confDict,confDict[elecMode_entry]["phase_rev"])
        ## Phase shift between Cavity cells and probe ADC
        self.phase_probe = readentry(confDict,confDict[elecMode_entry]["phase_probe"])

        # Read dictionary of couplings from global configuration dictionary
        mech_couplings = readentry(confDict,confDict[elecMode_entry]["mech_couplings"]["value"])
        # Get a coupling list of length M (number of mechanical modes),
        # filled with 0s if no coupling is specified by user
        ## List of mode's mechanical couplings
        self.mech_couplings_list = readCouplings(confDict, mech_couplings, cryomodule_entry)

        # Add (replicate) a parameter that will be filled after object instance
        ## Local Oscillator frequency [rad/s]
        self.LO_w0 = {"value" : 0.0, "units" : "rad/s", "description" : "Linac's Nominal resonance angular frequency"}

    def __str__(self):
        """Convenient concatenated string output for printout."""

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
        """
        Helper function to compute Electrical Mode's parameters (normally computed in C).
        Used in unit tests in order to compared measured properties to physical quantities
        Inputs:
            - Tstep: Simulation time step [s]
            - rf_phase: Beam phase relative to the RF [deg]
        """

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

        Q_L = 1.0/(1.0/self.Q_0['value'] + 1.0/self.Q_drive['value'] + 1.0/self.Q_probe['value'])
        bw = w0/(2.0*Q_L);
        k_beam = RoverQ*Q_L*np.exp(-1j*beam_phase)/Tstep;
        k_drive = 2.0*np.sqrt(self.Q_drive['value']*RoverQ);
        mode_dict = {"mode_name": mode_name,"w0": w0, "beam_phase": beam_phase, "RoverQ": RoverQ, "foffset": foffset, "Q_L": Q_L, "bw": bw, "k_beam": k_beam, "k_drive": k_drive, "k_probe": k_probe, "k_em": k_em}

        return mode_dict

class MechMode:
    """
    Contains parameters specific to a mechanical mode.
    Information concerning couplings with electrical modes and Piezos
    is contained in ElecMode and Piezo objects respectively."""

    def __init__(self, confDict, mechMode_entry):
        """
        MechMode Constructor:
        Inputs:
            - confDict: Global configuration dictionary,
            - mech_mode_entry: Name of the mechanical mode to be read (string).
        """


        ## Instance name
        self.name = confDict[mechMode_entry]['name']
        ## Instance type
        self.type = confDict[mechMode_entry]['type']

        # Read the rest of the configuration parameters and store in a dictionary
        ## Mechanical mode's resonant frequency [Hz]
        self.f0 = readentry(confDict,confDict[mechMode_entry]["f0"])
        ## Mechanical mode's Quality factor [unitless]
        self.Q = readentry(confDict,confDict[mechMode_entry]["Q"])
        ## Scaling factor used by FPGA
        self.full_scale = readentry(confDict,confDict[mechMode_entry]["full_scale"])

    def __str__(self):
        """Convenient concatenated string output for printout."""

        return ("\n--MechMode Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "f0: " + str(self.f0) + "\n"
        + "Q: " + str(self.Q) + "\n"
        + "full_scale: " + str(self.full_scale) + "\n")

    def Get_C_Pointer(self):
        """ Return reference to the SWIG-wrapped C structure. """
        import accelerator as acc

        # Grab attributes from object
        f0 = self.f0['value']
        Q = self.Q['value']
        k = 1.0

        # Allocate Memory for C struct
        mechMode = acc.MechMode_Allocate_New(f0, Q, k, Tstep_global);

        ## Pointer to the SWIG-wrapped C structure
        self.C_Pointer = mechMode

        # Return C Pointer
        return mechMode

class Piezo:
    """ Contains couplings between the Piezo and each
    one of the mechanical modes (MechMode instances)."""

    def __init__(self, confDict, piezo_entry, cryomodule_entry):
        """
        Piezo Constructor:
            Inputs:
                - confDict: Global configuration dictionary,
                - piezo_entry: Name of the Piezo to be read (string).
                - cryomodule_entry: Cryomodule entry in global dictionary in order to access
                    the proper cryomodule's mechanical mode list, which is used as a consistency
                    check to generate mechanical coupling vectors for each Piezo."""


        ## Instance name
        self.name = confDict[piezo_entry]['name']
        ## Instance type
        self.type = confDict[piezo_entry]['type']

        # Read rest of parameters
        ## Scaling factor used by the FPGA
        self.VPmax = readentry(confDict,confDict[piezo_entry]["VPmax"])

        # Read dictionary of couplings from global configuration dictionary
        mech_couplings = readentry(confDict,confDict[piezo_entry]["mech_couplings"]["value"])

        # Check consistency of coupling entries with the list of mechanical modes,
        # and get a coupling dictionary of length M (number of mechanical modes)
        ## Couplings to the Mechanical modes
        self.mech_couplings_list = readCouplings(confDict, mech_couplings, cryomodule_entry)

    def __str__(self):
        """Convenient concatenated string output for printout."""

        return ("\n--Piezo Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "VPmax: " + str(self.VPmax) + "\n"
        + "mech_couplings_list: " + str(self.mech_couplings_list))


def readCouplings(confDict, mech_couplings, cryomodule_entry):
    """
    Takes the global configuration dictionary and a dictionary containing non-zero values for the
    mechanical couplings (i.e. coupling between electrical modes or piezos and mechanical modes).
    The cryomodule_entry input is necessary in order to access the proper cryomodule's mechanical mode list.
    The length of the coupling vector used in the simulation code must be equal to the number of mechanical modes (M).
    The mech_couplings input is supposed to contain non-zero values from the configuration file,
    and readCouplings always returns a dictionary of M elements, where the couplings not specified in
    the input file are filled with 0s.
    Inputs:
        - confDict: Global configuration dictionary,
        - mech_couplings: dictionary containing couplings defined in the configuration file (0 to M elements).
        - cryomodule_entry: Cryomodule entry in global dictionary in order to access the proper cryomodule's mechanical mode list.
    Output:
        - mech_couplings_list: ordered list containing mechanical couplings for an electrical mode or piezo (length M).
            Order corresponds to the order of appearance of the mechanical mode in mech_net.
    """

    # Grab the full list of mechanical modes in the Cryomodule
    mech_net = confDict[cryomodule_entry]["mechanical_mode_connect"]

    # Make an ordered list of size M (where M is the total number of mechanical modes, length of mech_net)
    # Fill with 0s if no coupling is specified in mech_couplings by the user
    mech_couplings_list = [mech_couplings[m] if m in mech_couplings else 0.0 for m in mech_net]

    return mech_couplings_list


def readList(confDict, list_in, constructor, cryomodule_entry=None):
    """
    Generic function to read list of components.
    Takes the global configuration dictionary, cycles through the list of components
    (list of names, list_in), uses the names to identify the configuration entries in
    the global dictionary, calls the proper Constructor for each component (constructor),
    and returns a list of instances (objects).
    Inputs:
        - confDict: Global configuration dictionary,
        - list_in: list of components to cycle through (list of strings),
        - constructor: name of Constructor for the component,
        - cryomodule_entry (optional): necessary in some cases in order to pass along cryomodule entry information,
            needed by readStation and Piezo in order to find the mechanical modes in
            their corresponding Cryomodule.
    Output:
        - list_out: List of component objects.
    """

    # Create empty list for component instances
    list_out = []

    # Cycle through list of Instance names
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
    """ Contains parameters specific to a Controller configuration. """

    def __init__(self, confDict, controller_entry):
        ## Instance name
        self.name = confDict[controller_entry]['name']
        ## Instance type
        self.type = confDict[controller_entry]['type']

        # Read the rest of parameters
        ## FPGA controller Gain-Bandwidth product [Hz]
        self.stable_gbw = readentry(confDict,confDict[controller_entry]["stable_gbw"])

    def __str__(self):
        """Convenient concatenated string output for printout."""

        return ("\n--Controller Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "stable_gbw: " + str(self.stable_gbw) + "\n")

class ZFilter:
    """ Contains parameters specific to a Filter configuration"""

    def __init__(self, confDict, zfilter_entry):

        ## Instance name
        self.name = confDict[zfilter_entry]['name']
        ## Instance type
        self.type = confDict[zfilter_entry]['type']

        # Read the rest of parameters
        ## Filter order
        self.order = readentry(confDict,confDict[zfilter_entry]["order"])
        ## Total number of modes
        self.nmodes = readentry(confDict,confDict[zfilter_entry]["nmodes"])
        ## Filter poles
        self.poles = readentry(confDict,confDict[zfilter_entry]["poles"])

    def __str__(self):
        """Convenient concatenated string output for printout."""

        return ("\n--ZFilter Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "order: " + str(self.order) + "\n"
        + "nmodes: " + str(self.nmodes) + "\n"
        + "poles: " + str(self.poles) + "\n")

class ADC:
    """ Contains parameters specific to a ADC configuration"""

    def __init__(self, confDict, adc_entry):

        ## Instance name
        self.name = confDict[adc_entry]['name']
        ## Instance type
        self.type = confDict[adc_entry]['type']

        # Read the rest of parameters
        ## ADC full scale
        self.adc_max = readentry(confDict,confDict[adc_entry]["adc_max"])
        ## ADC offster
        self.adc_off = readentry(confDict,confDict[adc_entry]["adc_off"])
        ## ADC noise Power-Spectral-Density [dBc/Hz]
        self.noise_psd = readentry(confDict,confDict[adc_entry]["noise_psd"])

    def __str__(self):
        """Convenient concatenated string output for printout."""

        return ("\n--ADC Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "adc_max: " + str(self.adc_max) + "\n"
        + "adc_off: " + str(self.adc_off) + "\n"
        + "noise_psd: " + str(self.noise_psd) + "\n")

class Amplifier:
    """ Contains parameters specific to a Amplifier configuration"""

    def __init__(self, confDict, amplifier_entry):

        ## Instance name
        self.name = confDict[amplifier_entry]['name']
        ## Instance type
        self.type = confDict[amplifier_entry]['type']

        # Read the rest of parameters
        ## Maximum SSA output power [sqrt(W)]
        self.PAmax = readentry(confDict,confDict[amplifier_entry]["PAmax"])
        ## SSA scaling (from unitless to sqrt(W))
        self.PAbw = readentry(confDict,confDict[amplifier_entry]["PAbw"])
        ## Harshness parameter of SSA clipping function
        self.Clip = readentry(confDict,confDict[amplifier_entry]["Clip"])
        ## FPGA drive saturation limit [percentage of PAmax]
        self.top_drive = readentry(confDict,confDict[amplifier_entry]["top_drive"])

    def __str__(self):
        """Convenient concatenated string output for printout."""

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
    """ Contains parameters specific to a Station configuration"""

    def __init__(self, confDict, station_entry, cryomodule_entry):


        ## Instance name
        self.name = confDict[station_entry]['name']
        ## Instance type
        self.type = confDict[station_entry]['type']

        # Read all the station components
        amplifier_entry = confDict[station_entry]['Amplifier']
        ## Amplifier object
        self.amplifier = Amplifier(confDict, amplifier_entry)

        cavity_entry = confDict[station_entry]['Cavity']
        ## cavity object
        self.cavity = Cavity(confDict, cavity_entry, cryomodule_entry)

        rx_filter_entry = confDict[station_entry]['Rx_filter']
        ## filter object: Anti-alias filter
        self.rx_filter = ZFilter(confDict, rx_filter_entry)

        tx_filter1_entry = confDict[station_entry]['Tx_filter1']
        ## filter object: SSA filter 1
        self.tx_filter1 = ZFilter(confDict, tx_filter1_entry)

        tx_filter2_entry = confDict[station_entry]['Tx_filter2']
        ## filter object: SSA filter 2
        self.tx_filter2 = ZFilter(confDict, tx_filter2_entry)

        controller_entry = confDict[station_entry]['Controller']
        ## FPGA controller object
        self.controller = Controller(confDict, controller_entry)

        ## RF feedback loop delay in simulation time steps (multiply by Tstep to get seconds)
        self.loop_delay_size = readentry(confDict, confDict[station_entry]['loop_delay_size'])

        cav_adc_entry = confDict[station_entry]['cav_adc']
        ## Cavity field probe port ADC object
        self.cav_adc = ADC(confDict, cav_adc_entry)

        rev_adc_entry = confDict[station_entry]['rev_adc']
        ## Reverse port ADC object
        self.rev_adc = ADC(confDict, rev_adc_entry)

        fwd_adc_entry = confDict[station_entry]['fwd_adc']
        ## Forward port ADC object
        self.fwd_adc = ADC(confDict, fwd_adc_entry)

        piezo_connect = confDict[station_entry]['piezo_connect']
        ## List of Piezo objects
        self.piezo_list = readList(confDict, piezo_connect, Piezo, cryomodule_entry)

        ## Number of RF Stations
        self.N_Stations = confDict[station_entry]['N_Stations']

        # Add (replicate) parameters that will be filled after object instance
        ## Maximum accelerating voltage [V]
        self.max_voltage = {"value" : 0.0, "units" : "V", "description" : "Maximum accelerating voltage"}

    def __str__(self):
        """Convenient concatenated string output for printout."""

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
        + "rev_adc: " + str(self.rev_adc) + "\n"
        + "N_Stations: " + str(self.N_Stations) + "\n"
        + "piezo_list: " + '\n'.join(str(x) for x in self.piezo_list))

    def Get_C_Pointer(self):
        """ Return reference to the SWIG-wrapped C structure. """

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

        # Translate ADC noise Power Spectral Density (PSD),
        # which is expressed in [dBc/Hz] (where the carrier is the full range of the ADC)
        # into Volts RMS to scale the pseudo-random Gaussian noise
        probe_psd = self.cav_adc.noise_psd['value'] # [dBc/Hz]
        rev_psd = self.rev_adc.noise_psd['value']   # [dBc/Hz]
        fwd_psd = self.fwd_adc.noise_psd['value']   # [dBc/Hz]

        # Pre-calculate some common terms
        cav_V_nominal = self.cavity.nom_grad['value']*self.cavity.L['value']    # Cavity nominal accelerating voltage [V]
        bandwidth = 0.5/Tstep_global                                            # Bandwidth [Hz]

        # Calculate port couplings
        fund_index = self.cavity.fund_index['value']
        fund_Emode = self.cavity.elec_modes[fund_index]
        RoverQ = fund_Emode.RoverQ['value']
        k_probe = 1.0/np.sqrt(fund_Emode.Q_probe['value']*RoverQ);
        k_em = 1.0/np.sqrt(fund_Emode.Q_drive['value']*RoverQ);
        k_drive = 2.0*np.sqrt(fund_Emode.Q_drive['value']*RoverQ);

        # Now calculate ADC noise in Volts RMS (see Physics documentation for details)
        probe_ns_rms = 1.5*np.sqrt(0.5*10.0**(probe_psd/10.0)*bandwidth)*cav_V_nominal*k_probe
        rev_ns_rms = 1.5*np.sqrt(0.5*10.0**(rev_psd/10.0)*bandwidth)*cav_V_nominal*k_em
        fwd_ns_rms = 1.5*np.sqrt(0.5*10.0**(fwd_psd/10.0)*bandwidth)*cav_V_nominal/k_drive

        rf_station = acc.RF_Station()
        acc.RF_Station_Allocate_In(rf_station, Tstep_global, Clip, PAmax, PAscale, p_TRF1, p_TRF2, p_RXF, cavity_pointer, stable_gbw, FPGA_out_sat, loop_delay_size, probe_ns_rms, rev_ns_rms, fwd_ns_rms)

        ## Pointer to the SWIG-wrapped C structure
        self.C_Pointer = rf_station

        return rf_station

    def Get_State_Pointer(self):
        """ Return reference to the SWIG-wrapped State C structure. """
        import accelerator as acc

        rf_state = acc.RF_State()

        acc.RF_State_Allocate(rf_state, self.C_Pointer)
        ## Pointer to the SWIG-wrapped State C structureacc.RF_State_Allocate(rf_state, self.C_Pointer)
        self.State = rf_state

        return rf_state

class Cryomodule:
    """ Contains parameters specific to a Cryomodule configuration"""

    def __init__(self, confDict, cryomodule_entry):

        ## Instance name
        self.name = confDict[cryomodule_entry]['name']
        ## Instance type
        self.type = confDict[cryomodule_entry]['type']

        # Read the station and mechanical mode connectivity
        station_connect = confDict[cryomodule_entry]['station_connect']
        mechanical_mode_connect = confDict[cryomodule_entry]['mechanical_mode_connect']

        # Read list of stations and mechanical modes recursively
        ## List of RF Station objects
        self.station_list = readList(confDict, station_connect, Station, cryomodule_entry)
        ## List of mechanical eigenmodes
        self.mechanical_mode_list = readList(confDict, mechanical_mode_connect, MechMode)

        # Read lp_shift
        ## Scaling factor used in the FPGA
        self.lp_shift = readentry(confDict,confDict[cryomodule_entry]["lp_shift"])

    def __str__(self):
        """Convenient concatenated string output for printout."""

        return ("\n--Cryomodule Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "station_list: " + '\n'.join(str(x) for x in self.station_list)
        + "mechanical_mode_list: " + '\n'.join(str(x) for x in self.mechanical_mode_list)
        + "lp_shift: " + str(self.lp_shift) + "\n")

    def Get_C_Pointer(self):
        """ Return reference to the SWIG-wrapped C structure. """

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

        ## Pointer to the SWIG-wrapped C structure
        self.C_Pointer = cryomodule

        # Return C Pointer for Cryomodule and lists of pointers
        return cryomodule

    def Get_State_Pointer(self, cryo_state=None):
        """
        Return reference to the SWIG-wrapped State C structure.
        Input:
            - cryo_state (optional): if part of a hierarchy and state has already been allocated,
                provide the reference and it will be assigned to the object's State attribute.
        """
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

        ## Pointer to the SWIG-wrapped State C structure
        self.State = cryo_state

        # Return State C pointer
        return cryo_state

class Linac:
    """ Contains parameters specific to a Linac configuration"""

    def __init__(self, confDict, linac_entry):

        import numpy as np

        # Read name and component type
        name = confDict[linac_entry]['name']
        ## Instance type
        self.type = confDict[linac_entry]['type']

        ## Local Oscillator frequency [Hz]
        self.f0 = readentry(confDict,confDict[linac_entry]["f0"])
        ## Energy at the end of Linac Section [eV]
        self.E = readentry(confDict,confDict[linac_entry]["E"])
        ## Nominal Linac RF phase [radians]
        self.phi = readentry(confDict,confDict[linac_entry]["phi"])
        self.phi['value'] = self.phi['value']*np.pi/180    # Convert degrees to radians
        ## Wakefield characteristic length (Sband=0.105m, Xband=0.02625m) [m]
        self.s0 = readentry(confDict,confDict[linac_entry]["s0"])
        ## Mean iris radius (Sband=11.654mm,Xband=4.72mm) [m]
        self.iris_rad = readentry(confDict,confDict[linac_entry]["iris_rad"])
        ## Longitudinal dispersion (if any)
        self.R56 = readentry(confDict,confDict[linac_entry]["R56"])
        ## DDS factor used in FPGA
        self.dds_numerator = readentry(confDict,confDict[linac_entry]["dds_numerator"])
        ## DDS factor used in FPGA
        self.dds_denominator = readentry(confDict,confDict[linac_entry]["dds_denominator"])

        # Read the cryomodule connectivity
        cryomodule_connect = confDict[linac_entry]['cryomodule_connect']

        # Read list of modules recursively
        ## List of cryomodule objects
        self.cryomodule_list = readList(confDict, cryomodule_connect, Cryomodule)

        # Add parameters that will be filled after object instance
        ## Energy increase in Linac Section [eV]
        self.dE = {"value" : 0.0, "units" : "eV", "description" : "Energy increase in Linac (final minus initial Energy)"}
        ## Maximum Accelerating Voltage [V]
        self.max_voltage = {"value" : 0.0, "units" : "V", "description" : "Maximum Accelerating Voltage"}
        ## Total number of RF Stations in Linac
        self.N_Stations = {"value" : 0.0, "units" : "N/A", "description" : "Total number of RF Stations in Linac"}
        ## Total Linac Length [m]
        self.L = {"value" : 0.0, "units" : "m", "description" : "Total Linac Length"}

        # Some parameters are deduced from others
        # RF wavelength deduced from f0
        c = 2.99792458e8    # Speed of light [m/s]

        ## RF wavelength [m]
        self.lam = {"value" : c/self.f0['value'], "units" : "m", "description" : "RF wavelength (Sband=0.105m, Xband=0.02625m)"}
        # T566 deduced from R56 (small angle approximation)
        ## 2nd-order longitudinal dispersion (if any)
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
        """Convenient concatenated string output for printout."""

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
        """ Return reference to the SWIG-wrapped C structure. """

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

        ## Pointer to the SWIG-wrapped C structure
        self.C_Pointer = linac

        # Return Linac C Pointer
        return linac

    def Get_State_Pointer(self, linac_state=None):
        """
        Return reference to the SWIG-wrapped State C structure.
        Input:
            - linac_state (optional): if part of a hierarchy and state has already been allocated,
                provide the reference and it will be assigned to the object's State attribute.
        """
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
    """ Contains parameters specific to a Simulation run,
    as well as all parameters in the Accelerator configuration. This Class is
    to be instantiated from upper level programs in order to obtain all necessary instances
    to run a full simulation"""

    def __init__(self, confDict):
        """Simulation Constructor:
            Inputs:
                confDict: Global configuration dictionary."""

        import numpy as np

        ## Instance name
        self.name = confDict["Accelerator"]["name"]
        ## Instance type
        self.type = confDict["Accelerator"]["type"]

        # Read rest of configuration parameters
        ## Simulation time-step [s]
        self.Tstep = readentry(confDict, confDict["Simulation"]["Tstep"])
        ## Total Simulation time duration in time steps
        self.time_steps = readentry(confDict, confDict["Simulation"]["time_steps"])
        ## Factor used in FPGA
        self.nyquist_sign = readentry(confDict,confDict["Simulation"]["nyquist_sign"])

        # Check if simulation dictionary has a Synthesis entry, and if so parse it
        if confDict["Simulation"].has_key("Synthesis"):
            self.synthesis = Synthesis(confDict["Simulation"])
        else:
            self.synthesis = None

        # Accelerator parameters
        self.bunch_rate = readentry(confDict,confDict["Accelerator"]["bunch_rate"])

        # Read Noise Sources
        ## Simulation correlared noise sources
        self.noise_srcs = Noise(confDict)

        # Read Accelerator components (Gun + series of linacs)
        # Read gun
        ## Gun object
        self.gun = Gun(confDict)
        Egun = self.gun.E['value'] # Gun exit Energy

        # Read connectivity of linacs
        linac_connect = confDict["Accelerator"]["linac_connect"]
        # Read linacs recursively
        ## List of Linac objects
        self.linac_list = readList(confDict, linac_connect, Linac)

        # Now that the Array of Linacs has been instantiated and their final Energies are known,
        # fill in the Energy increase parameter for each Linac.
        # Start with the Energy out of the Gun
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
        ## Final Accelerator Energy [eV]
        self.E = {"value" : Elast, "units" : "eV", "description" : "Final Accelerator Energy"}

        # Assign Tstep to the global variable
        global Tstep_global
        Tstep_global = self.Tstep['value']

    def __str__(self):
        """Convenient concatenated string output for printout."""

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
        """ Return reference to the SWIG-wrapped C structure. """

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

        ## Pointer to the SWIG-wrapped C structure
        self.C_Pointer = sim

        # Return C Pointer for Simulation
        return sim

    def Get_State_Pointer(self):
        """ Return reference to the SWIG-wrapped State C structure. """
        import accelerator as acc

        # Allocate memory for Noise State
        noise_State_Pointer = self.noise_srcs.Get_State_Pointer()

        sim_state = acc.Simulation_State()
        acc.Sim_State_Allocate(sim_state, self.C_Pointer, noise_State_Pointer)

        ## Pointer to the SWIG-wrapped State C structure
        self.State = sim_state

        return sim_state

class Gun:
    """ Contains parameters specific to an Gun configuration"""

    def __init__(self, confDict):

        gun_entry = confDict["Accelerator"]['gun']

        ## Instance name
        self.name = confDict[gun_entry]['name']
        ## Instance type
        self.type = confDict[gun_entry]['type']

        ## Nominal beam charge [C]
        self.Q = readentry(confDict,confDict[gun_entry]["Q"])
        ## Nominal initial RMS bunch length [m]
        self.sz0 = readentry(confDict,confDict[gun_entry]["sz0"])
        ## Nominal initial incoh. energy spread at nominal gun exit energy [fraction]
        self.sd0 = readentry(confDict,confDict[gun_entry]["sd0"])
        ## Nominal gun exit energy [eV]
        self.E = readentry(confDict,confDict[gun_entry]["E"])

    def __str__(self):
        """Convenient concatenated string output for printout."""

        return ("\n--Gun Object--\n"
        + "name: " + self.name + "\n"
        + "type: " + self.type + "\n"
        + "Q: " + str(self.Q) + "\n"
        + "sz0: " + str(self.sz0) + "\n"
        + "sd0: " + str(self.sd0) + "\n"
        + "E: " + str(self.E) + "\n")

    def Get_C_Pointer(self):
        """ Return reference to the SWIG-wrapped C structure. """

        import accelerator as acc

        gun = acc.Gun()
        acc.Gun_Allocate_In(gun, self.E['value'], self.sz0['value'], self.sd0['value'], self.Q['value']);

        ## Pointer to the SWIG-wrapped C structure
        self.C_Pointer = gun

        return gun

class Noise:
    """ Contains configuration regarding the correlated noise sources in the Accelerator."""

    def __init__(self, confDict):

        noise_entry = confDict['Noise']

        ## Instance name
        self.name = confDict[noise_entry]['name']
        ## Instance type
        self.type = confDict[noise_entry]['type']

        self.dQ_Q = readentry(confDict,confDict[noise_entry]["dQ_Q"])
        self.dtg = readentry(confDict,confDict[noise_entry]["dtg"])
        self.dE_ing = readentry(confDict,confDict[noise_entry]["dE_ing"])
        self.dsig_z = readentry(confDict,confDict[noise_entry]["dsig_z"])
        self.dsig_E = readentry(confDict,confDict[noise_entry]["dsig_E"])
        self.dchirp = readentry(confDict,confDict[noise_entry]["dchirp"])

    def __str__(self):
        """Convenient concatenated string output for printout."""

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
        """ Return reference to the SWIG-wrapped State C structure. """

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
