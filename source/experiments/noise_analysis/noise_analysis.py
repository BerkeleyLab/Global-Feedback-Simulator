"""
A series of exercises to illustrate a noise analysis of the LLRF System.
It includes an open-loop stability analysis, as well as numerical simulations
illustrating the impact of each noise souce in the LLRF system.
"""

import os
import sys
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
if not path in sys.path:
    sys.path.insert(1, path)
del path

import accelerator as acc

import numpy as np
import matplotlib.pylab as plt

def find_zero_crossing(vector, frequency_vector):
    """
    Find the 0-crossing given a vector in frequency domain and the associated frequency vector.
    """
    signs = np.sign(vector)
    signs[signs == 0] = -1
    zero_crossing = np.where(np.diff(signs))[0]
    return vector[zero_crossing], frequency_vector[zero_crossing], zero_crossing

def margins(vector, frequency_vector, configuration, plotit=True):
    """
    Given the open-loop response provided as the vector input,
    find amplitude and phase crossover frequencies and deduce Gain (GM) and Phase (PM) Margins.
    Optionally (enabled by default) plot the results.
    """

    # Calculate amplitude and phase
    vector_amp = 20.0*np.log10(np.abs(vector))  # [dB]
    vector_phase = np.unwrap(np.angle(vector))*180/np.pi     # [degrees]

    # Find amplitude (f_c) and phase (f_180) crossover frequencies
    amp_at_f_c, f_c, f_c_index = find_zero_crossing(vector_amp, frequency_vector)
    phase_at_f_180, f_180, f_180_index = find_zero_crossing(vector_phase+180.0, frequency_vector)

    # Calculate Gain (GM) and Phase (PM) Margins
    GM = -vector_amp[f_180_index]
    PM = 180 + vector_phase[f_c_index]

    if plotit:
        unity = np.ones(len(vector))
        plt.figure(1)

        plt.subplot(211)
        plt.semilogx(frequency_vector, vector_amp, linewidth=4)
        plt.semilogx(frequency_vector, unity-1, linewidth=4, linestyle='dotted')
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        title = 'Open-Loop Transfer Function (%s)' % (configuration)
        plt.title(title, fontsize=40, y=1.02)
        plt.ylabel('Magnitude [dB]', fontsize=30)
        # plt.annotate(s='', xy=(f_180,vector_amp[f_180_index]), xytext=(f_180,unity[f_180_index]-1), arrowprops=dict(arrowstyle='<->') )
        GM_text = 'GM = %.d dB' % int(GM)
        plt.annotate(GM_text, xy=(f_180*1.5, vector_amp[f_180_index]), fontsize=25)
        plt.annotate(s='', xy=(f_180, vector_amp[f_180_index]), xytext=(
            f_180, unity[f_180_index]-1), arrowprops=dict(arrowstyle='<->'))
        f_c_text = r'$f_{\rm c} \approx %d\,{\rm kHz}$' % np.ceil(f_c*1e-3)
        plt.annotate(f_c_text,
                     xy=(f_c, 1), xytext=(f_c*1.2, 25), arrowprops=dict(facecolor='black', shrink=0.05), fontsize=28)
        plt.semilogx(f_c, amp_at_f_c, 'ro')

        plt.subplot(212)
        plt.semilogx(frequency_vector, vector_phase, linewidth=4)
        plt.semilogx(frequency_vector, unity-181, linewidth=4, linestyle='dotted')
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.annotate(s='', xy=(f_c, vector_phase[f_c_index]), xytext=(f_c, unity[f_c_index]-181), arrowprops=dict(arrowstyle='<->'))
        PM_text = 'PM = %d$^{\circ}$' % int(PM)
        plt.annotate(PM_text, xy=(f_c*0.5, -250), fontsize=25)
        f_180_text = r'$f_{\rm 180} \approx %d\,{\rm kHz}$' % int(round(f_180*1e-3))
        plt.annotate(f_180_text,
                     xy=(f_180, -180), xytext=(f_180*1.2, -100), arrowprops=dict(facecolor='black', shrink=0.05), fontsize=28)
        plt.xlabel('Frequency [Hz]', fontsize=30)
        plt.ylabel('Phase [degrees]', fontsize=30)
        plt.semilogx(f_180, phase_at_f_180-180, 'ro')

    return GM, PM, f_c, f_180

def station_analytical(test_file, configuration, margin=False, plot_tfs=False):
    """
    Run the stability analysis for a RF Station configuration.
    """

    # Import JSON parser module
    import numpy as np
    from matplotlib import pyplot as plt
    from get_configuration import Get_SWIG_RF_Station

    # Configuration file for specific test configuration
    # (to be appended to standard test cavity configuration)
    rf_station, Tstep, fund_mode_dict = Get_SWIG_RF_Station(test_file, Verbose=False)

    # Build a linear-scale frequency vector and Laplace operator
    frequency_vector = np.linspace(1, 1e6, 1e6)
    s = 1j*2.0*np.pi*frequency_vector

    # Cavity
    wc = fund_mode_dict['bw']
    cavity = 1.0/(1.0+s/wc)

    # Loop delay
    T = rf_station.loop_delay_size['value']*Tstep
    delay = np.exp(-s*T)

    # Plant (cavity + delay)
    plant = cavity*delay

    # SSA noise-limiting Low-pass filter
    low_pass_pole = 2*np.pi*complex(rf_station.ns_filter_bw['value'])
    lowpass = 1.0/(1.0+s/low_pass_pole)

    # Unity vector
    unity = np.ones(len(s))

    # Controller
    stable_gbw = rf_station.controller.stable_gbw['value']  # [Hz]
    control_zero = rf_station.controller.control_zero['value']  # [Hz]
    # Calculate corresponding values for ki and kp
    Kp = stable_gbw*2*np.pi/wc
    Ki = Kp*(2*np.pi*control_zero)

    # Calculate transfer functions
    PI = (Kp + Ki/s)
    controller = PI*lowpass

    noise_gain = 1.0/(1+plant*controller)
    loop = plant*controller
    closed_loop = cavity/(1+plant*controller)

    # Plot gain and phase margin (open-loop analysis)
    if margin:
        GM, PM, f_c, f_180 = margins(loop, frequency_vector, configuration)
        plt.show()

    # Plot loop transfer functions
    if plot_tfs:
        plt.semilogx(frequency_vector, 20.0*np.log10(np.abs(plant)), linewidth=3, label='Plant')
        plt.semilogx(frequency_vector, 20.0*np.log10(np.abs(controller)), linewidth=3, label='Controller')
        plt.semilogx(frequency_vector, 20.0*np.log10(np.abs(loop)), linewidth=3, label='Loop')
        plt.semilogx(frequency_vector, 20.0*np.log10(np.abs(noise_gain)), linewidth=3, label='Noise Gain')
        plt.semilogx(frequency_vector, 20.0*np.log10(np.abs(unity)), linewidth=3, label='Unity', linestyle='dotted')

        title = 'RF Station Analytical (%s)' % (configuration)
        plt.title(title, fontsize=40, y=1.02)
        plt.ylabel('Magnitude [dB]', fontsize=30)
        plt.xlabel('Frequency [Hz]', fontsize=30)
        plt.legend(loc='upper right')
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.rc('font', **{'size': 20})

        plt.show()

class Signal:
    """ Support class to consolidate signals involved in a numerical simulation of an RF Station.
    """

    def __init__(self, E_fwd, E_reverse, E_probe, cav_v, set_point, beam_current, error, fund_mode_dict, trang):
        """ Constructor:
        Inputs:
            - Signals of interest in an RF Station simulation, including a time vector (trang).
        """

        self.E_fwd = E_fwd
        self.E_reverse = E_reverse
        self.E_probe = E_probe
        self.cav_v = cav_v
        self.set_point = set_point
        self.beam_current = beam_current
        self.error = error
        self.fund_mode_dict = fund_mode_dict
        self.trang = trang

    def plot(self, title, xlim=None, ylim=None, plot_type='amplitude', show_limits=False, show_ssa_limits=False, beam_on=False, annotate=False, beam_start=20, beam_end=22, hide_fwd=False):
        """ Plot signals of interest in RF Station time-series simulation.
        Supports amplitude, phase and cartesian plots with a number of options.
        """
        fund_k_probe = self.fund_mode_dict['k_probe']
        fund_k_drive = self.fund_mode_dict['k_drive']
        fund_k_em = self.fund_mode_dict['k_em']
        fund_k_beam = self.fund_mode_dict['k_beam']

        # Amplitude
        if plot_type == 'amplitude':
            if hide_fwd == False:
                plt.plot(self.trang*1e3, np.abs(self.E_fwd)*fund_k_drive, '-',
                         label=r'Forward $\left(\vec E_{\rm fwd}\right)$', linewidth=2)
                plt.plot(self.trang*1e3, np.abs(self.E_reverse/fund_k_em), '-',
                         label=r'Reverse $\left(\vec E_{\rm reverse}\right)$', linewidth=2)
            plt.plot(self.trang*1e3, np.abs(self.cav_v), '-', label=r'Cavity Field', linewidth=2, color='r')
            plt.plot(self.trang*1e3, np.abs(self.set_point/fund_k_probe), '-',
                     label=r'Set-point $\left(\vec E_{\rm sp}\right)$', linewidth=2, color='c')
            # Y label
            plt.ylabel('Amplitude [V]', fontsize=30)

            if show_limits == True:
                plt.plot(self.trang*1e3, np.abs(self.set_point/fund_k_probe)*(1+1e-4), label='Upper limit', linewidth=2, color='m')
                plt.plot(self.trang*1e3, np.abs(self.set_point/fund_k_probe)*(1-1e-4), label='Lower limit', linewidth=2, color='y')
                plt.axhspan(16e6/1.00005, 16e6*1.00005, color='blue', alpha=0.2)

            if show_ssa_limits == True:
                # If SSA noise were a sine wave, this would be the amplitude limit
                ssa_limit = 4e-2*fund_k_drive*np.sqrt(3.8e3)
                plt.plot(self.trang*1e3, np.abs(self.set_point/fund_k_probe)+ssa_limit /
                         np.sqrt(2)/2, label='SSA Upper (RMS) limit', linewidth=2, color='m')
                plt.plot(self.trang*1e3, np.abs(self.set_point/fund_k_probe)-ssa_limit /
                         np.sqrt(2)/2, label='SSA Lower (RMS) limit', linewidth=2, color='y')

                low_index = np.where(self.trang == 20e-3)[0][0]
                high_index = np.where(self.trang == 49.9e-3)[0][0]

                ssa_std = np.std(self.E_fwd[low_index:high_index]*fund_k_drive)
                ssa_std_percent = 100*ssa_std/(fund_k_drive*np.sqrt(3.8e3))
                ssa_std_text = r'$\sigma_{\rm SSA}\,=\,$' + '%.2f %% RMS' % (np.abs(ssa_std_percent))
                plt.text(35, 13e6, ssa_std_text, verticalalignment='top', fontsize=30)

        # Phase
        if plot_type == 'phase':
            plt.plot(self.trang*1e3, np.angle(self.cav_v, deg=True), '-', label=r'Cavity Field', linewidth=2, color='r')
            plt.plot(self.trang*1e3, np.angle(self.set_point/fund_k_probe, deg=True),
                     label=r'Set-point $\left(\vec E_{\rm sp}\right)$', linewidth=2, color='c')
            if show_limits:
                plt.plot(self.trang*1e3, 1e-2*np.ones(len(self.trang)), label='Upper limit', linewidth=2, color='m')
                plt.plot(self.trang*1e3, -1e-2*np.ones(len(self.trang)), label='Lower limit', linewidth=2, color='y')
                plt.axhspan(-4e-3, 4e-3, color='blue', alpha=0.2)

            # Y label
            plt.ylabel('Phase [degrees]', fontsize=30)

        # Cartesian coordinates
        if plot_type == 'cartesian':
            plt.plot(self.trang*1e3, np.real(self.E_fwd)*fund_k_drive, '-',
                     label=r'Forward $\Re \left(\vec E_{\rm fwd}\right)$', linewidth=2)
            plt.plot(self.trang*1e3, np.imag(self.E_fwd)*fund_k_drive, '-',
                     label=r'Forward $\Im \left(\vec E_{\rm fwd}\right)$', linewidth=2)
            plt.plot(self.trang*1e3, np.real(self.cav_v), '-', label=r'$\Re$ Cavity Field', linewidth=2)
            plt.plot(self.trang*1e3, np.imag(self.cav_v), '-', label=r'$\Im$ Cavity Field', linewidth=2)

        # Add annotation for a very specific plot
        if annotate:
            delta_E_fwd_text = r'$\Delta \left| \vec E_{\rm fwd} \right| \approx \,$' + \
                '%.d MV' % (round(np.abs(fund_k_beam)*100e-6*1e-6))
            plt.annotate(delta_E_fwd_text, xy=(21.05, 18e6), fontsize=25)
            index_of_21 = np.where(self.trang == 21e-3)
            plt.annotate(s='', xy=(21, np.abs(self.E_fwd[index_of_21])*fund_k_drive),
                         xytext=(21, np.abs(self.cav_v[index_of_21])), arrowprops=dict(arrowstyle='<->'))

        # Add clear red shade if beam current is on to highlight the location of the beam train
        if beam_on:
            plt.axvspan(beam_start, beam_end, color='red', alpha=0.2)

        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.title(title, fontsize=40, y=1.02)
        plt.xlabel('Time [ms]', fontsize=30)

        if xlim:
            plt.xlim(xlim)

        if ylim:
            plt.ylim(ylim)

        plt.rc('font', **{'size': 20})
        plt.legend(loc='upper right')

        plt.show()


def run_noise_numerical(Tmax, test_file, beam_configuration, detuning=False):
    """
    Run the time-series simulations to analyze the impact of different noise sources/perturbations in an RF Station.
    """

    # Grab beam configuration parameters
    feedforward = beam_configuration['feedforward']
    beam_charge = beam_configuration['charge']
    n_pulses = beam_configuration['n_pulses']
    beam_start = beam_configuration['beam_start']
    missing_pulse = beam_configuration['missing_pulse']
    charge_noise_percent = beam_configuration['noise']

    # Import JSON parser module
    from get_configuration import Get_SWIG_RF_Station

    # Configuration file for specific test configuration
    # (to be appended to standard test cavity configuration)
    rf_station, Tstep, fund_mode_dict = Get_SWIG_RF_Station(test_file, set_point=16e6, Verbose=False)

    # Create time vector
    trang = np.arange(0, Tmax, Tstep)

    # Number of points
    nt = len(trang)

    # Initialize vectors for test
    cav_v = np.zeros(nt, dtype=np.complex)
    error = np.zeros(nt, dtype=np.complex)
    E_probe = np.zeros(nt, dtype=np.complex)
    E_reverse = np.zeros(nt, dtype=np.complex)
    E_fwd = np.zeros(nt, dtype=np.complex)
    set_point = np.ones(nt, dtype=np.complex)

    # Initialize beam vector with 0s
    beam_current = np.zeros(nt, dtype=np.complex)
    feed_forward = np.zeros(nt, dtype=np.complex)

    # Create beam vector based on configuration parameters
    # Calculate current
    nominal_beam_current = beam_charge/Tstep  # Amps

    # Add charge noise if turned on
    if charge_noise_percent != None:
        beam_current_samples = -np.random.normal(-nominal_beam_current, -nominal_beam_current*charge_noise_percent, n_pulses)
    else:
        beam_current_samples = nominal_beam_current*np.ones(nt, dtype=np.double)

    # Fill in the pulses
    for i in range(n_pulses):
        beam_current[int((beam_start+i*1e-6)/Tstep)] = beam_current_samples[i]

    # Grab some scaling factors
    fund_k_beam = fund_mode_dict['k_beam']
    fund_k_drive = fund_mode_dict['k_drive']
    fund_k_probe = fund_mode_dict['k_probe']

    # Deduce feed-forward signal with some error
    nominal_feed_forward = -0.99*np.sqrt(2)*nominal_beam_current*(Tstep/1e-6)*fund_k_beam/fund_k_drive

    # Generate feed-forward vector if enabled
    if feedforward:
        feed_forward[int((beam_start-0.5e-6)/Tstep):int((beam_start+(n_pulses-0.5)*1e-6)/Tstep)] = nominal_feed_forward

    # Remove a pulse if option is enabled
    if missing_pulse:
        beam_current[int((beam_start+1e-6)/Tstep)] = 0.0
        feed_forward[int((beam_start+0.5e-6)/Tstep):int((beam_start+1.5*1e-6)/Tstep)] = 0.0

    # Introduce a modulation of the detune frequency if detuning is enabled
    if detuning:
        elecMode_state = acc.ElecMode_State_Get(rf_station.State.cav_state, 0)
        delta_omega_vector = 2*np.pi*10*np.sin(2*np.pi*100*trang)
    else:
        delta_omega_vector = np.zeros(nt, dtype=np.double)

    # Record initial set-point
    set_point = set_point*rf_station.C_Pointer.fpga.set_point

    # Run Numerical Simulation
    for i in xrange(0, nt):
        if detuning:
            elecMode_state.delta_omega = delta_omega_vector[i]

        # Adjust the set-point if feedforward is on
        if (i == (beam_start-0.5e-6)/Tstep) & feedforward:
            rf_station.C_Pointer.fpga.set_point = rf_station.C_Pointer.fpga.set_point - 209*fund_k_probe
        if (i == (beam_start+(n_pulses-0.5)*1e-6)/Tstep) & feedforward:
            rf_station.C_Pointer.fpga.set_point = rf_station.C_Pointer.fpga.set_point + 209*fund_k_probe

        # Run RF Station time-series simulation
        cav_v[i] = acc.RF_Station_Step(rf_station.C_Pointer, 0.0, beam_current[i], feed_forward[i], rf_station.State)
        # Record signals of interest
        error[i] = rf_station.State.fpga_state.err
        E_probe[i] = rf_station.State.cav_state.E_probe
        E_reverse[i] = rf_station.State.cav_state.E_reverse
        E_fwd[i] = rf_station.State.cav_state.E_fwd

    # Instantiate Signal object to provide caller with simulation results
    signal = Signal(E_fwd, E_reverse, E_probe, cav_v, set_point, beam_current, error, fund_mode_dict, trang)
    return signal

def run_noise_numerical_tests():
    """
    Set up configuration for different noise analysis tests and call routine to produce the associated time-series data.
    """

    # Simulation time
    Tmax = 25e-3

    # Default beam_configuration
    beam_configuration = {
        'charge': -100e-12,
        'n_pulses': 1,
        'beam_start': 20e-3,
        'feedforward': False,
        'missing_pulse': False,
        'noise': None
    }

    ###### Beam Loading Exercise: Pulse Train @ 100uA ######
    # No Feed-forward
    beam_configuration['feedforward'] = False
    beam_configuration['n_pulses'] = 2000
    test_file = "source/experiments/noise_analysis/nominal.json"
    signal1 = run_noise_numerical(Tmax, test_file, beam_configuration)
    title = 'Beam Train: '+r'$100 \mu A $'+', no feed-forward'
    signal1.plot(title, ylim=[0, 30e6])
    signal1.plot(title, ylim=[0, 30e6], beam_on=True)
    signal1.plot(title, xlim=[19.5, 22.5], ylim=[14e6, 23e6], beam_on=True, annotate=True)
    signal1.plot(title, xlim=[19.98, 20.3], ylim=[16e6-10e3, 16e6+10e3], beam_on=True, show_limits=True)

    # # With Feed-forward
    beam_configuration['feedforward'] = True
    beam_configuration['noise'] = 0.01
    signal1 = run_noise_numerical(Tmax, test_file, beam_configuration)
    title = 'Beam Train: '+r'$100 \mu A $'+', with feed-forward & charge noise'
    signal1.plot(title, xlim=[19.98, 20.3], ylim=[16e6-10e3, 16e6+10e3], beam_on=True, show_limits=True)
    signal1.plot(title, xlim=[19.98, 22], ylim=[16e6-10e3, 16e6+10e3], beam_on=True, show_limits=True)

    ###### Detuning Exercise: With Pulse Train @ 100uA, no charge noise ######
    Tmax = 50e-3
    beam_configuration['n_pulses'] = 15000
    beam_configuration['noise'] = None
    signal1 = run_noise_numerical(Tmax, test_file, beam_configuration, detuning=True)
    title = 'Detuning: '+r'$\rm \Delta f_c = 10\cdot \sin(2\pi 100t)\,Hz $'+', with feed-forward'
    signal1.plot(title, beam_on=True, beam_end=35)
    signal1.plot(title, xlim=[19.9, 35], ylim=[-15e6, 25e6], plot_type='cartesian', beam_on=True, beam_end=35)
    signal1.plot(title, xlim=[19.9, 35], ylim=[16e6-10e3, 16e6+10e3], beam_on=True, show_limits=True, beam_end=35)
    signal1.plot(title, xlim=[19.9, 35], ylim=[-6.25e-2, 6.25e-2], plot_type='phase', beam_on=True, show_limits=True, beam_end=35)

    ###### Measurement Noise Exercise: No Beam ######
    Tmax = 50e-3
    beam_configuration['n_pulses'] = 0
    # Nominal configuration
    # Amplitude
    test_file = ["source/experiments/noise_analysis/nominal.json", "source/experiments/noise_analysis/149_noise.json"]
    signal1 = run_noise_numerical(Tmax, test_file, beam_configuration, detuning=False)
    title = 'Measurement Noise: -149 dBc/Hz'
    signal1.plot(title)

    # Phase
    test_file = ["source/experiments/noise_analysis/nominal.json", "source/experiments/noise_analysis/149_noise.json"]
    signal1 = run_noise_numerical(Tmax, test_file, beam_configuration, detuning=False)
    title = 'Measurement Noise: -149 dBc/Hz'
    signal1.plot(title, xlim=[20, 35], ylim=[16e6-10e3, 16e6+10e3], show_limits=True, hide_fwd=True)
    signal1.plot(title, xlim=[20, 35], ylim=[-6.25e-2, 6.25e-2], plot_type='phase', show_limits=True)

    ###### Measurement Noise Exercise: Impact on SSA Noise, No Beam ######
    # Nominal with -149 dBc/Hz
    title = 'Measurement Noise: Nominal with -149 dBc/Hz'
    signal1.plot(title, xlim=[20, 50], ylim=[12e6, 20e6], show_ssa_limits=True)

    # HoBiCaT with -149 dBc/Hz
    test_file = ["source/experiments/noise_analysis/hobicat.json", "source/experiments/noise_analysis/149_noise.json"]
    signal1 = run_noise_numerical(Tmax, test_file, beam_configuration, detuning=False)
    title = 'Measurement Noise: HoBiCaT with -149 dBc/Hz'
    signal1.plot(title, xlim=[20, 50], ylim=[12e6, 20e6], show_ssa_limits=True)

    # High gain with -149 dBc/Hz
    test_file = ["source/experiments/noise_analysis/high_gain.json", "source/experiments/noise_analysis/149_noise.json"]
    signal1 = run_noise_numerical(Tmax, test_file, beam_configuration, detuning=False)
    title = 'Measurement Noise: High gain with -149 dBc/Hz'
    signal1.plot(title, xlim=[20, 50], ylim=[12e6, 20e6], show_ssa_limits=True)

    # High gain with -138 dBc/Hz
    test_file = ["source/experiments/noise_analysis/high_gain.json", "source/experiments/noise_analysis/138_noise.json"]
    signal1 = run_noise_numerical(Tmax, test_file, beam_configuration, detuning=False)
    title = 'Measurement Noise: High gain with -138 dBc/Hz'
    signal1.plot(title, xlim=[20, 50], ylim=[12e6, 20e6], show_ssa_limits=True)
    signal1.plot(title, xlim=[20, 50], ylim=[16e6-10e3, 16e6+10e3], show_limits=True, hide_fwd=True)


if __name__ == "__main__":
    plt.close('all')

    # # Run Analytical Analysis
    # # Nominal configuration
    test_file = "source/experiments/noise_analysis/nominal.json"
    station_analytical(test_file, 'nominal configuration', margin=True, plot_tfs=True)
    # # High-gain configuration
    test_file = "source/experiments/noise_analysis/high_gain.json"
    station_analytical(test_file, 'high gain', margin=True, plot_tfs=False)

    # Run Numerical simulations
    run_noise_numerical_tests()
