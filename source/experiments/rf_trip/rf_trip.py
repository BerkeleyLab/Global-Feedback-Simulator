"""
Exercise to illustrate system behavior after an RF trip.
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
import scipy.linalg as linalg
from scipy import stats

def run_RF_trip_experiment(Tmax, test_file):

    # Import JSON parser module
    from get_configuration import Get_SWIG_Cryomodule

    cryo_object, Tstep, fund_mode_dicts =  Get_SWIG_Cryomodule(test_file, Verbose=False)

    # Create time vector
    trang = np.arange(0,Tmax,Tstep)

    # Number of points
    nt = len(trang)

    # Initialize vectors for test
    cav_v = np.zeros(nt,dtype=np.complex)   # Overall cavity accelerating voltage
    E_reverse = np.zeros(nt,dtype=np.complex)
    Kg = np.zeros(nt,dtype=np.complex)

    delta_omega = np.zeros(nt,dtype=np.double)

    # Manually set the controller set-point
    # FPGA controller will be operating in open-loop mode,
    # which is equivalent to driving the cavity with the set-point value
    # cryo_object.station_list[0].C_Pointer.fpga.set_point = 30.0
    # cryo_object.station_list[0].State.fpga_state.openloop = 1

    # Get a C pointer to the Electrical Mode State to record its detune frequency
    elecMode_state = acc.ElecMode_State_Get(cryo_object.station_list[0].State.cav_state, 0);

    # Run Numerical Simulation
    for i in xrange(1,nt):
        # Run Cryomodule Step function
        acc.Cryomodule_Step(cryo_object.C_Pointer, cryo_object.State, 0.0, 0.0)

        # Record some waveforms
        # Cavity state
        cav_v[i] = cryo_object.station_list[0].State.cav_state.V
        E_reverse[i] = cryo_object.station_list[0].State.cav_state.E_reverse
        Kg[i] = cryo_object.station_list[0].State.cav_state.Kg

        # pi-mode's detune frequency [rad/s]
        delta_omega[i] = elecMode_state.delta_omega

    # Make some calculations for plots
    # Mode's bandwidth
    omega_f = fund_mode_dicts[0]['bw']
    # Derivative of the cavity field
    cav_v_der = np.diff(cav_v, 1)/Tstep
    factor = 1-1j*delta_omega/omega_f
    # Calculate waveforms with corresponding factors to equate units in order to plot together
    x2 = cav_v_der/omega_f
    x3 = np.multiply(cav_v, factor)

    # Plot cavity signals
    plt.plot(trang*1e6, np.abs(Kg), '-', label=r'Forward ($\vec K_{\rm g}$)', linewidth = 2)
    plt.plot(trang*1e6, np.abs(cav_v)*1e-4, '-', label=r'Cavity ($\vec V_{\pi}\times \, 10^{-4}$)', linewidth = 2)
    plt.plot(trang*1e6, np.abs(E_reverse), '-', label=r'Reverse ($\vec E_{\rm reverse}$)', linewidth = 2)

    plt.title('Cryomodule Test', fontsize=40, y=1.01)
    plt.xlabel(r'Time [$\rm \mu s$]', fontsize=30)
    plt.ylim(0, 60)
    plt.ylabel('Amplitude [V]', fontsize=30)
    plt.legend(loc='upper right')
    plt.rc('font',**{'size':25})

    plt.show()
    plt.figure()

    # Plot detune frequency vs time
    plt.title('Detune Frequency', fontsize=40, y=1.01)
    plt.plot(trang*1e6, 1.9e-7*np.square(np.abs(cav_v))/2.0/np.pi, '-', label=r'$k_{\rm L}|\vec V_{\pi}|^{\rm 2}$', linewidth = 2)
    plt.plot(trang*1e6, -delta_omega/2.0/np.pi, '-', label=r'$-\Delta f_{\pi}$', linewidth = 2)
    plt.ylabel('Frequency [Hz]', fontsize=30)
    plt.xlabel(r'Time [$\rm \mu s$]', fontsize=30)
    plt.legend(loc='upper right')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ylim(0, 8e3)
    plt.show()

    # Plot cavity signals in complex plane
    dpi=100
    f = plt.figure(figsize=(1400/dpi, 1350/dpi), dpi=dpi)
    x = f.gca()
    x.set_aspect("equal")

    plt.title(r'$\pi$-mode in Volts', fontsize=40, y=1.01)

    x.plot(np.real(cav_v), np.imag(cav_v), '-', label=r'Cavity ($\vec V_{\pi}$)', linewidth = 2)
    x.plot(np.real(x2), np.imag(x2), '-', label=r'$\omega_{\rm f}\times \, \frac{d\vec V_{\pi}}{dt}$', linewidth = 2)
    x.plot(np.real(x3), np.imag(x3), '-', label=r'$\left( 1-j\omega_{\rm d}/\omega_{\rm f}\right) \vec V_{\pi}$', linewidth = 2)
    x.set_ylabel(r'$\Im$ [V]', fontsize=30)
    x.set_xlabel(r'$\Re$ [V]', fontsize=30)
    x.set_xlim(-7e5, 7e5)
    x.set_ylim(-7e5, 7e5)
    x.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    x.legend(loc='upper right')

    plt.show()
    plt.figure()

def RF_trip():
    """
    Perform all unit tests for cryomodule.c/h and return PASS boolean (qualitative test).
    """

    # Tmax = 80e-6

    test_file = "source/experiments/rf_trip/rf_trip.json"

    # # This is not a PASS/FAIL test
    # print "\n****\nRunning RF trip experiment..."
    # run_RF_trip_experiment(Tmax, test_file)

    # plt.figure()

    Tmax = 0.05
    # test_file = "source/configfiles/unit_tests/cavity_test_step1.json"

    # This is not a PASS/FAIL test
    print "\n****\nRunning RF Station test..."
    run_RF_Station_test(Tmax, test_file)

    plt.figure()


    return True

def run_RF_Station_test(Tmax, test_file):

    # Import JSON parser module
    from get_configuration import Get_SWIG_Cryomodule

    cryo_object, Tstep, fund_mode_dicts =  Get_SWIG_Cryomodule(test_file, Verbose=False)

    # Create time vector
    trang = np.arange(0,Tmax,Tstep)

    # Number of points
    nt = len(trang)

    # Initialize vectors for test
    cav_v = np.zeros(nt,dtype=np.complex)   # Overall cavity accelerating voltage
    E_reverse = np.zeros(nt,dtype=np.complex)
    Kg = np.zeros(nt,dtype=np.complex)

    delta_omega = np.zeros(nt,dtype=np.double)

    # Initialize vectors for test
    E_probe = np.zeros(nt,dtype=np.complex)
    E_reverse = np.zeros(nt,dtype=np.complex)
    E_fwd = np.zeros(nt,dtype=np.complex)
    set_point = np.zeros(nt,dtype=np.complex)

    # Run Numerical Simulation
    for i in xrange(1,nt):
        # Run Cryomodule Step function
        acc.Cryomodule_Step(cryo_object.C_Pointer, cryo_object.State, 0.0, 0.0)
        set_point[i] = cryo_object.station_list[0].C_Pointer.fpga.set_point
        E_probe[i] = cryo_object.station_list[0].State.cav_state.E_probe
        E_reverse[i] = cryo_object.station_list[0].State.cav_state.E_reverse
        E_fwd[i] = cryo_object.station_list[0].State.cav_state.E_fwd

    fund_mode_dict = fund_mode_dicts[0]
    fund_k_probe = fund_mode_dict['k_probe']
    fund_k_drive = fund_mode_dict['k_drive']
    fund_k_em = fund_mode_dict['k_em']

    plt.plot(trang*1e3,np.abs(E_reverse/fund_k_em), '-', label=r'Reverse $\left(\vec E_{\rm reverse}\right)$', linewidth=2)
    plt.plot(trang*1e3,np.abs(E_fwd)*fund_k_drive, '-',  label=r'Forward $\left(\vec E_{\rm fwd}\right)$', linewidth=2)
    plt.plot(trang*1e3,np.abs(set_point/fund_k_probe), label=r'Set-point $\left(\vec E_{\rm sp}\right)$', linewidth=2)
    plt.plot(trang*1e3,np.abs(E_probe/fund_k_probe), '-', label=r'Probe $\left(\vec E_{\rm probe}\right)$', linewidth=2)

    plt.title('RF Station Test', fontsize=40, y=1.01)
    plt.xlabel('Time [ms]', fontsize=30)
    plt.ylabel('Amplitude [V]', fontsize=30)
    plt.legend(loc='upper right')
    plt.rc('font',**{'size':25})

    # Show plot
    plt.show()


def unit_RF_Station():
    """
    Unit test for rf_station.c/h
    It emulates a cavity fill-up, where forward signal is saturated and
    cavity field reaching steady-state with controller stabilizing the field
    around the set point. Plot is generated for qualitative analysis.
    This is not a PASS/FAIL test but responds to the typical RF Station behavior.
    """
    Tmax = 0.05

    test_file = "source/configfiles/unit_tests/cavity_test_step1.json"

    run_RF_Station_test(Tmax, test_file)

if __name__=="__main__":
    plt.close('all')
    RF_trip()
