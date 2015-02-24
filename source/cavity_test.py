#!/usr/bin/python

#
# A Series of Unit tests for filter.c/h, cavity.c/h
#

import accelerator as acc

import numpy as np
import matplotlib.pylab as plt
import scipy.linalg as linalg

####################################
#
# Unit test for Filter_Step
#
####################################

def unit_filter(dt=0.002):
    tmax = 5.0
    nt = int(tmax/dt)
    st = np.zeros(nt,dtype=np.complex)

    fil = acc.Filter() # Declare one in python and then set it up
    acc.Filter_Allocate_In(fil,3,3)

    # Configure filter for two conjugate poles at (-1+1j) and (-1-1j)
    poles = acc.complexdouble_Array(1)
    poles[0] = -1.0+1.0j
    acc.Filter_Append_Modes(fil,poles,1,dt)
    poles[0] = -1.0-1.0j
    acc.Filter_Append_Modes(fil,poles,1,dt)

    # Allocate Filter State Vector
    filnow = acc.Filter_State()
    acc.Filter_State_Allocate(filnow,fil)

    # Run Numerical Simulation
    for i in xrange(1,nt):
        st[i] = acc.Filter_Step(fil, 1.0+0.0j, filnow)

    # Create time vector
    trang = np.arange(0,tmax,dt)

    # Calculate Analytical step response
    anal = 1.0-np.exp(-trang)*(np.sin(trang)+np.cos(trang))

    # RMS Error between Numerical and Analytical
    error = linalg.norm(st.real-anal)/linalg.norm(anal)

    # Format text to include in plot
    error_text = 'RMS error is {:.2e}'.format(error)

    # Plot
    plt.plot(trang,st.real,'bo', label='Numerical')
    plt.plot(trang,anal,':r', label='Analytical')
    plt.title("Filter step reponse: Numerical vs. Analytical", fontsize=40, y=1.01)
    plt.xlabel('Time [s]', fontsize=30)
    plt.ylabel('Amplitude [Normalized]', fontsize=30)
    plt.legend(loc='upper right')
    plt.text(3.5,0.1, error_text, verticalalignment='top', fontsize=30)
    plt.text(3.5,0.2, r'$TF(s) = \frac{1}{(s+1)^2+1}$', fontsize=40)
    plt.rc('font',**{'size':25})

    # Establish a threshold for RMS error to claim success of numerical simulation
    if error < 10e-3: filter_pass = True
    else: filter_pass = False

    plt.show()

    # Return PASS boolean to escalate diagnostics
    return filter_pass

def cavity_curve_fit(Tstep, drive_in, cav_v, beam_charge):

    from scipy.signal import lfilter

    # First find the derivative of the cavity field,
    # and fit to find bandwidth (pole location)
    cav_v_der = np.diff(cav_v, 1)
    
    # Stack vectors
    A = np.column_stack([drive_in[:-1], beam_charge[:-1], cav_v[:-1]])
    
    # Curve fit and find the three coefficients in the least-squares sense
    # Fit cavity step response according to equation:
    # D[i] = a*Kg[i] + b*Beam[i] - c*V[i]
    # where: Kg is the drive signal, 
    # and D is the discrete derivative of V (cavity accelerating voltage)
    # (all quantities above are complex numbers)
    a, b, c = np.linalg.lstsq(A, cav_v_der)[0]

    # Test fiting
    # Calculate cavity derivative based on fitted coefficients
    cav_v_der =  a*drive_in + b*beam_charge + c*cav_v;
    # Integrate to obtain measured cavity voltage
    cav_v_meas = lfilter(np.array([1.0]),np.array([1.0,-1.0]),cav_v_der, axis=0);
    # Calculate measured cavity open-loop bandwidth based on c coefficient
    bw_meas =  -np.real(c)/(Tstep*np.pi)

    # Return fitted coefficients
    return a, b, c, cav_v_meas, bw_meas, cav_v

def run_cavity_step_test(Tmax, test_file, drive_on = True, beam_on = False):

    # Import JSON parser module
    from get_configuration import Get_SWIG_Cavity

    # Configuration file for specific test configuration
    # (to be appended to standard test cavity configuration)
    cav, cav_state, Tstep, modes_config = Get_SWIG_Cavity(test_file, Verbose=False)
    
    # Create time vector
    trang = np.arange(0,Tmax,Tstep)
    
    # Number of points
    nt = len(trang)

    # Initialize vectors for test
    cav_v_drive = np.zeros(nt,dtype=np.complex)   # Overall cavity accelerating voltage
    cav_v_beam = np.zeros(nt,dtype=np.complex)   # Overall cavity accelerating voltage
    E_probe = np.zeros(nt,dtype=np.complex)   # Cavity field probe
    E_reverse = np.zeros(nt,dtype=np.complex)   # Reverse field port
    
    # Drive signal
    drive_in_d = np.ones(nt,dtype=np.complex)
    drive_in_b = np.zeros(nt,dtype=np.complex)
    
    # Beam charge
    beam_charge_d = np.zeros(nt,dtype=np.complex)
    beam_charge_b = np.ones(nt,dtype=np.complex)

    delta_tz = 0.0  # Timing noise

    # Run Numerical Simulation
    for i in xrange(1,nt):
        cav_v_drive[i] = acc.Cavity_Step(cav, delta_tz, drive_in_d[i], beam_charge_d[i], cav_state);
        E_probe[i] = cav_state.E_probe
        E_reverse[i] = cav_state.E_reverse

    # Fit cavity step response
    drive_step = cavity_curve_fit(Tstep, drive_in_d, cav_v_drive, beam_charge_d)

    # # Clear cavity state for a new run
    acc.Cavity_Clear(cav, cav_state)

    # Run Numerical Simulation
    for i in xrange(1,nt):
        cav_v_beam[i] = acc.Cavity_Step(cav, delta_tz, drive_in_b[i], beam_charge_b[i], cav_state);

    # # Fit cavity step response
    beam_step = cavity_curve_fit(Tstep, drive_in_b, cav_v_beam, beam_charge_b)
    
    # Pass along the 1st mode configuration dictionary (useful for single mode tests)
    mode_dict = modes_config[0]
    # Return results
    # drive_step and beam_step contain: [a, b, c, cav_v_meas, bw_meas, cav_v]
    return trang, drive_step, beam_step, mode_dict, E_probe, E_reverse

def run_cavity_freq_test(Tmax, test_file):

    # Import JSON parser module
    from get_configuration import Get_SWIG_Cavity

    # Configuration file for specific test configuration
    # (to be appended to standard test cavity configuration)
    cav, cav_state, Tstep, modes_config = Get_SWIG_Cavity(test_file, Verbose=False)
    
    # Create time vector
    trang = np.arange(0,Tmax,Tstep)
    
    # Number of points
    nt = len(trang)

    # Initialize vectors for test
    cav_v = np.zeros(nt,dtype=np.complex)   # Overall cavity accelerating voltage

    # Drive signal
    drive_in = np.ones(nt,dtype=np.complex)
    
    # Beam charge
    beam_charge = np.zeros(nt,dtype=np.complex)

    delta_tz = 0.0  # Timing noise

    # Run Numerical Simulation
    for i in xrange(1,nt):
        cav_v[i] = acc.Cavity_Step(cav, delta_tz, drive_in[i], beam_charge[i], cav_state)

    # Pass along the 1st mode configuration dictionary (useful for single mode tests)
    mode_dict = modes_config[0]

    # Fit the step response
    curve_fit = cavity_curve_fit(Tstep, drive_in, cav_v, beam_charge)

    # Measure mode's offset frequency
    foffset_meas = -np.imag(curve_fit[2])/(Tstep*2.0*np.pi)

    # Return results
    # drive_step and beam_step contain: [a, b, c, cav_v_meas, bw_meas, cav_v]
    return trang, cav_v, drive_in, beam_charge, mode_dict, curve_fit, foffset_meas

def show_cavity_step(title):
    plt.title(title, fontsize=40, y=1.02)
    plt.xlabel('Time [s]', fontsize=30)
    plt.ylabel(r'$| \vec V_{\rm acc}|$ [V]', fontsize=30)

    plt.ticklabel_format(style='sci', axis='y', scilimits=(1,0))
    
    plt.rc('font',**{'size':20})
    plt.legend(loc='upper right')
    plt.show()

def cavity_test_step():

    # Configuration file for specific test configuration
    # (to be appended to standard test cavity configuration)
    test_files = ["source/configfiles/unit_tests/cavity_test_step1.json", \
        "source/configfiles/unit_tests/cavity_test_step2.json", \
        "source/configfiles/unit_tests/cavity_test_step3.json"]

    # Total simulation time 
    Tmax = 0.2

    drive_steps = []
    beam_steps = []
    mode_dicts = []
    E_probes = []
    E_reverses = []

    # Pass/FAIL boolean for curve fits
    # True = curve fit RMS error are under threshold
    fit_pass = True
    fit_threshold = 1e-3

    for test_file in test_files:
        trang, drive_step, beam_step, mode_dict, E_probe, E_reverse = run_cavity_step_test(Tmax, test_file)
        drive_steps.append(drive_step)
        beam_steps.append(beam_step)
        mode_dicts.append(mode_dict)
        E_probes.append(E_probe)
        E_reverses.append(E_reverse)

    drive_in = np.ones(len(trang),dtype=np.complex)
    
    print "\n**** Step drive signal at modes' frequencies...\n"
    for idx, mode_dict in enumerate(mode_dicts):
        plt.plot(trang,np.abs(drive_steps[idx][5]),'-', label='Numerical ('+mode_dict['mode_name']+')', linewidth=5)
        plt.plot(trang,np.abs(drive_steps[idx][3]),'-', label='Curve fit ('+mode_dict['mode_name']+')', linewidth=3)
        
        # Print some results for log file
        print "Mode "+mode_dict['mode_name']+":"
        # Measure mode's bandwidth (calculated in run_and_fit_cavity based on curve fit)
        bw_meas = drive_steps[idx][4]
        print "  Bandwidth: Measured = %.2f Hz, Set to = %.2f Hz" %(bw_meas, mode_dict['bw']/np.pi)

        # Measure gain on drive path
        k_drive_meas = max(np.abs(drive_steps[idx][5]))
        print "  Drive coupling: Measured = {:.2e}, Set to = {:.2e}".format(k_drive_meas, np.abs(mode_dict['k_drive']))
        
        # Calculate error RMS on curve fit
        cav_fit_error = linalg.norm(np.abs(drive_steps[idx][5])-np.abs(drive_steps[idx][3]))/linalg.norm(np.abs(drive_steps[idx][5]))
        E_probe_fit_error = linalg.norm(np.abs(E_probes[idx])-mode_dicts[idx]['k_probe']*np.abs(drive_steps[idx][5]))/linalg.norm(np.abs(drive_steps[idx][5]))
        E_reversre_fit_error = linalg.norm(np.abs(E_reverses[idx] - mode_dicts[idx]['k_em']*drive_steps[idx][5] + drive_in))/linalg.norm(np.abs(E_reverses[idx]))
        
        print '  Cavity accelerating field fit RMS error is {:.2e}'.format(cav_fit_error)
        print '  Cavity probe fit RMS error is {:.2e}'.format(E_probe_fit_error)
        print '  Cavity reverse fit RMS error is {:.2e}'.format(E_reversre_fit_error)
        
        # Compare error to threshold and establish PASS/FAIL
        if ((cav_fit_error < fit_threshold) & (E_probe_fit_error < fit_threshold) & (E_reversre_fit_error < 3e-3)): this_fit_pass  = True
        else: this_fit_pass = False
        fit_pass = fit_pass & this_fit_pass

    # Show response to step on drive signal
    show_cavity_step("Cavity Step Response: Drive @ "+r'$\omega_{\mu}$')

    # Plot response to step on beam
    print "\n**** Step beam signal at modes' frequencies...\n"
    for idx, mode_dict in enumerate(mode_dicts):
        plt.plot(trang,np.abs(beam_steps[idx][5]),'-', label='Numerical ('+mode_dict['mode_name']+')', linewidth=5)
        plt.plot(trang,np.abs(beam_steps[idx][3]),'-', label='Curve fit ('+mode_dict['mode_name']+')', linewidth=3)
        
        # Measure gain on beam path
        k_beam_meas = max(np.abs(beam_steps[idx][5]))

        print "Mode "+mode_dict['mode_name']+":"
        print "  Beam coupling: Measured= {:.2e}, Set to = {:.2e}".format(k_beam_meas, np.abs(mode_dict['k_beam'])*1e-12)

        # Calculate error RMS on curve fit
        fit_error = linalg.norm(np.abs(beam_steps[idx][5])-np.abs(beam_steps[idx][3]))/linalg.norm(np.abs(beam_steps[idx][5]))
        print '  Fit RMS error is {:.2e}'.format(fit_error)
        
        # Compare error to threshold and establish PASS/FAIL
        if fit_error < fit_threshold: this_fit_pass = True
        else: this_fit_pass = False
        fit_pass = fit_pass & this_fit_pass

    # Show response to step on drive signal
    show_cavity_step("Cavity Step Response: 1 pC Beam step")
    
    # print PASS/FAIL
    if (fit_pass):
        result = 'PASS' 
    else:
        result = 'FAIL'
    print "\n*** Step test >>> " + result 

    # Return PASS/FAIL boolean
    return fit_pass

def cavity_test_freqs():

    # Configuration file for specific test configuration
    # (to be appended to standard test cavity configuration)
    test_files = ["source/configfiles/unit_tests/cavity_test_freqs1.json", \
        "source/configfiles/unit_tests/cavity_test_freqs2.json", \
        "source/configfiles/unit_tests/cavity_test_freqs3.json"]

    # Total simulation time 
    Tmax = 1

    # Empty lists to append simulation results to
    drive_in_list = []
    cav_v_list = []
    mode_dicts = []

    # Pass/FAIL boolean for curve fits
    # True = curve fit RMS error are under threshold
    fit_pass = True
    fit_threshold = 1e-3

    print "\n**** Test detuning...\n"
    # Iterate over simulation tests
    for idx, test_file in enumerate(test_files):
        # Run numerical simulation
        trang, cav_v, drive_in, beam_charge, mode_dict, curve_fit, foffset_meas = run_cavity_freq_test(Tmax, test_file)
        # Store signals
        drive_in_list.append(drive_in)
        cav_v_list.append(cav_v)
        mode_dicts.append(mode_dict)

        # Calculate error RMS on curve fit
        fit_error = linalg.norm(np.abs(curve_fit[5])-np.abs(curve_fit[3]))/linalg.norm(np.abs(curve_fit[5]))
        print "Mode "+str(idx+1)+":"
        print '  Fit RMS error is {:.2e}'.format(fit_error)
        print '  Offset: Measured = %.2f Hz, Set to = %.2f Hz' %(foffset_meas, mode_dict['foffset'])

        # Compare error to threshold and establish PASS/FAIL
        if fit_error < fit_threshold: this_fit_pass = True
        else: this_fit_pass = False
        fit_pass = fit_pass & this_fit_pass


    # Plot title    
    plt.title("Cavity Step Response: Drive @ "+r'$\omega_0$', fontsize=40, y=1.05)

    # Iterate over modes and plot amplitude
    for idx, mode_dict in enumerate(mode_dicts):
        plt.plot(np.real(cav_v_list[idx]),np.imag(cav_v_list[idx]),'-', label= 'offset = ' + str(mode_dict['foffset']) + ' Hz', linewidth=2)
    
    # Format plot
    plt.ticklabel_format(style='sci', axis='y', scilimits=(1,0))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(1,0))
    plt.ylim([-1e5,9e5])
    plt.xlabel(r'$\Re ( \vec V_{\mu})$ [V]', fontsize=30)
    plt.ylabel(r'$\Im (\vec V_{\mu})$ [V]', fontsize=30)
    plt.rc('font',**{'size':20})
    plt.legend(loc='upper right')

    # Show figure
    plt.show()

    # print PASS/FAIL
    if (fit_pass):
        result = 'PASS' 
    else:
        result = 'FAIL'
    print "\n*** Detuning test >>> " + result 

    # Return PASS/FAIL boolean
    return fit_pass

def unit_cavity():
    print "\n**** Test Cavity step response..."
    test_step_pass = cavity_test_step()
    test_freqs_pass = cavity_test_freqs()

    return test_step_pass & test_freqs_pass

def perform_tests():
    print "\n****\nTesting Filter..."
    filter_pass = unit_filter()
    if (filter_pass):
        result = 'PASS' 
    else:
        result = 'FAIL'
    print ">>> " + result 

    plt.figure()

    print "\n****\nTesting Cavity..."
    cavity_pass = unit_cavity()
    if (cavity_pass):
        result = 'PASS' 
    else:
        result = 'FAIL'
    print "\n>>> " + result 

    plt.figure()

    return filter_pass and cavity_pass

if __name__=="__main__":
    plt.close('all')
    perform_tests()