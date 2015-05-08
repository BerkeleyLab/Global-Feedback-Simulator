#!/usr/bin/python

#
# A Series of Unit tests for cryomodule.c/h
#

import accelerator as acc

import numpy as np
import matplotlib.pylab as plt
import scipy.linalg as linalg
from scipy import stats

####################################
#
# Unit tests for step functions in Cryomodule
#
####################################

#
# Unit test for Cryomodule
#

def run_Cryomodule_test(Tmax, test_file):

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
    cryo_object.station_list[0].C_Pointer.fpga.set_point = 30.0

    # Get a C pointer to the Electrical Mode State to record its detune frequency
    elecMode_state = acc.ElecMode_State_Get(cryo_object.station_list[0].State.cav_state, 0);

    # Run Numerical Simulation
    for i in xrange(1,nt):
        # Run Cryomodule Step function
        acc.Cryomodule_Step(cryo_object.C_Pointer, cryo_object.State, 0.0, 0.0, 0.0, 0.0, 1)

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
    plt.ylim(0, 55)
    plt.ylabel('Amplitude [V]', fontsize=30)
    plt.legend(loc='upper right')
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

def unit_Cryomodule():

    Tmax = 80e-6
    
    test_file = "source/configfiles/unit_tests/cryomodule_test.json"

    run_Cryomodule_test(Tmax, test_file)

######################################
#
# Now execute the tests...
#
######################################

def perform_tests():

    # This is not a PASS/FAIL test
    print "\n****\nTesting Cryomodule..."
    unit_Cryomodule()
    print ">>> (Visual inspection only)\n" 

    plt.figure()
    
    return True

if __name__=="__main__":
    plt.close('all')
    perform_tests()