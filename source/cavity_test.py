"""
Unit tests for cavity.c/h and filter.c/h
"""


import accelerator as acc

import numpy as np
import matplotlib.pylab as plt
import scipy.linalg as linalg

def unit_filter(dt=0.002):
    """
    Unit test for filter.c/h
    Compare numerical simulation with analytical solution, plot and return PASS/FAIL boolean.
    """

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
    """
    Fit cavity field signal to 1st-order exponential response.
    """

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

def run_cavity_step_test(Tmax, test_file):

    # Import JSON parser module
    from get_configuration import Get_SWIG_Cavity

    # Configuration file for specific test configuration
    # (to be appended to standard test cavity configuration)
    cav, Tstep, modes_config = Get_SWIG_Cavity(test_file, Verbose=False)

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
        cav_v_drive[i] = acc.Cavity_Step(cav.C_Pointer, delta_tz, drive_in_d[i], beam_charge_d[i], cav.State);
        E_probe[i] = cav.State.E_probe
        E_reverse[i] = cav.State.E_reverse

    # Fit cavity step response
    drive_step = cavity_curve_fit(Tstep, drive_in_d, cav_v_drive, beam_charge_d)

    # # Clear cavity state for a new run
    acc.Cavity_Clear(cav.C_Pointer, cav.State)

    # Run Numerical Simulation
    for i in xrange(1,nt):
        cav_v_beam[i] = acc.Cavity_Step(cav.C_Pointer, delta_tz, drive_in_b[i], beam_charge_b[i], cav.State);

    # # Fit cavity step response
    beam_step = cavity_curve_fit(Tstep, drive_in_b, cav_v_beam, beam_charge_b)

    # Pass along the 1st mode configuration dictionary (useful for single mode tests)
    mode_dict = modes_config[0]
    # Return results
    # drive_step and beam_step contain: [a, b, c, cav_v_meas, bw_meas, cav_v]
    return trang, drive_step, beam_step, mode_dict, E_probe, E_reverse

def run_cavity_freq_test(Tmax, test_file, delta_omega=0.0):

    # Import JSON parser module
    from get_configuration import Get_SWIG_Cavity

    # Configuration file for specific test configuration
    # (to be appended to standard test cavity configuration)
    cav, Tstep, modes_config = Get_SWIG_Cavity(test_file, Verbose=False)

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

    elecMode_state = acc.ElecMode_State_Get(cav.State,0);
    elecMode_state.delta_omega = delta_omega
    # Run Numerical Simulation
    for i in xrange(1,nt):
        cav_v[i] = acc.Cavity_Step(cav.C_Pointer, delta_tz, drive_in[i], beam_charge[i], cav.State)

    # Pass along the 1st mode configuration dictionary (useful for single mode tests)
    mode_dict = modes_config[0]

    # Fit the step response
    curve_fit = cavity_curve_fit(Tstep, drive_in, cav_v, beam_charge)

    # Measure mode's offset frequency
    foffset_meas = np.imag(curve_fit[2])/(Tstep*2.0*np.pi)

    # Return results
    # drive_step and beam_step contain: [a, b, c, cav_v_meas, bw_meas, cav_v]
    return trang, cav_v, drive_in, beam_charge, mode_dict, curve_fit, foffset_meas

def run_cavity_detune_test(Tmax, test_file, delta_omega_step=0.0):

    # Import JSON parser module
    from get_configuration import Get_SWIG_Cavity

    # Configuration file for specific test configuration
    # (to be appended to standard test cavity configuration)
    cav, Tstep, modes_config = Get_SWIG_Cavity(test_file, Verbose=False)

    # Create time vector
    trang = np.arange(0,Tmax,Tstep)

    # Number of points
    nt = len(trang)

    # Initialize vectors for test
    cav_v = np.zeros(nt,dtype=np.complex)   # Overall cavity accelerating voltage

    # Drive signal
    drive_in = np.ones(nt,dtype=np.complex)

    drive_in[0:int(nt*0.1)] = 0.0
    # Beam charge
    beam_charge = np.zeros(nt,dtype=np.complex)

    # Frequency offset
    w_offset = np.zeros(nt,dtype=np.double)

    delta_tz = 0.0  # Timing noise

    elecMode_state = acc.ElecMode_State_Get(cav.State,0);
    elecMode_state.delta_omega = 0.0
    # Run Numerical Simulation
    for i in xrange(1,nt):
        if(i>int(nt*0.4)): elecMode_state.delta_omega = delta_omega_step
        cav_v[i] = acc.Cavity_Step(cav.C_Pointer, delta_tz, drive_in[i], beam_charge[i], cav.State)
        w_offset[i] = elecMode_state.delta_omega

    # Pass along the 1st mode configuration dictionary (useful for single mode tests)
    mode_dict = modes_config[0]

    # Return results
    return trang, cav_v, drive_in, beam_charge, mode_dict, w_offset


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

        k_probe_meas = max(np.abs(E_probes[idx]))/k_drive_meas

        print '  Cavity accelerating field fit RMS error is {:.2e}'.format(cav_fit_error)
        print '  Cavity probe fit RMS error is {:.2e}'.format(E_probe_fit_error)
        print "  Cavity probe coupling: Measured = {:.2e}, Set to = {:.2e}".format(k_probe_meas, np.abs(mode_dict['k_probe']))
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
        "source/configfiles/unit_tests/cavity_test_freqs3.json", \
        "source/configfiles/unit_tests/cavity_test_freqs4.json"]
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

    print "\n**** Test detuning..."
    print "******* Iterate over basis frequency offset...\n"

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

    # Now exercise the dynamic deturning input
    # List the detune frequencies to test
    delta_fs = [0.0,10.0, 20.0, 100.0]  # [Hz]
    test_file = test_files[0] # Use test file with no basis frequency offset

    # Empty lists to append simulation results to
    drive_in_list2 = []
    cav_v_list2 = []
    mode_dicts2 = []

    print "\n******* Iterate over frequency perturbation (basis offset = 0)...\n"
    # Iterate over detune frequencies
    for idx, delta_f in enumerate(delta_fs):
        # Run numerical simulation
        trang, cav_v, drive_in, beam_charge, mode_dict, curve_fit, foffset_meas = run_cavity_freq_test(Tmax, test_file, 2.0*np.pi*delta_f)
        # Store signals
        drive_in_list2.append(drive_in)
        cav_v_list2.append(cav_v)
        mode_dicts2.append(mode_dict)

        print "Mode "+str(idx+1)+":"
        print '  Fit RMS error is {:.2e}'.format(fit_error)
        print '  Offset: Measured = %.2f Hz, Set to = %.2f Hz' %(foffset_meas, delta_f)

        # Compare error to threshold and establish PASS/FAIL
        if fit_error < fit_threshold: this_fit_pass = True
        else: this_fit_pass = False
        fit_pass = fit_pass & this_fit_pass

    # Plot cavity signals in complex plane
    ## Make some arrangements so axis are equal and size of plot is fixed
    dpi=96
    f = plt.figure(figsize=(1400/dpi, 1350/dpi), dpi=dpi)
    x = f.gca()
    x.set_aspect("equal")

    # Plot title
    plt.title("Cavity Step Response: Drive @ "+r'$\omega_{\rm ref}$', fontsize=40, y=1.05)

    # Iterate over modes and plot amplitude
    # # First, basis frequency offset tests
    for idx, mode_dict in enumerate(mode_dicts):
        x.plot(np.real(cav_v_list[idx]),np.imag(cav_v_list[idx]),'-', label= 'basis offset = ' + str(mode_dict['foffset']) + ' Hz', linewidth=5)

    # # Then perturbation frequency detuning tests
    for idx, mode_dict in enumerate(mode_dicts):
        x.plot(np.real(cav_v_list2[idx]),np.imag(cav_v_list2[idx]),'-', label= 'perturbation offset = ' + str(mode_dict['foffset']) + ' Hz', linewidth=2)

    # Format plot
    x.ticklabel_format(style='sci', axis='y', scilimits=(1,0))
    x.ticklabel_format(style='sci', axis='x', scilimits=(1,0))
    x.set_ylim([-0.5e5,6e5])
    x.set_xlim([-0.5e5,6e5])
    x.set_xlabel(r'$\Re ( \vec V_{\mu})$ [V]', fontsize=30)
    x.set_ylabel(r'$\Im (\vec V_{\mu})$ [V]', fontsize=30)
    plt.rc('font',**{'size':15})
    x.legend(loc='upper right')

    # Show figure
    plt.show()

    # print PASS/FAIL
    if (fit_pass):
        result = 'PASS'
    else:
        result = 'FAIL'
    print "\n*** Detuning tests >>> " + result

    # Return PASS/FAIL boolean
    return fit_pass

def cavity_test_detune():

    # Configuration file for specific test configuration
    # (to be appended to standard test cavity configuration)
    test_file = "source/configfiles/unit_tests/cavity_test_freqs1.json"

    # Total simulation time
    Tmax = 0.5

    # Empty lists to append simulation results to
    drive_in_list = []
    cav_v_list = []
    mode_dicts = []
    w_offset_list = []

    delta_fs = [100.0, 200.0]   # [Hz]

    print "\n**** Test detuning step..."
    print ">>> (Visual inspection only)\n"

    # # First, basis frequency offset tests
    for idx, delta_f in enumerate(delta_fs):
        # Run numerical simulation
        trang, cav_v, drive_in, beam_charge, mode_dict, w_offset = run_cavity_detune_test(Tmax, test_file, 2.0*np.pi*delta_f)
        # Store signals
        drive_in_list.append(drive_in)
        cav_v_list.append(cav_v)
        w_offset_list.append(w_offset)
        mode_dicts.append(mode_dict)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Plot title
    plt.title("Cavity Test: Fill & Step on Detune Frequency", fontsize=40, y=1.05)

    # Format plots
    # # Cavity fields for the two simulation runs
    lns1 = ax.plot(trang,np.abs(cav_v_list[0]),'-', label= 'Cavity Field ('+r'$\Delta f_1$'+')', linewidth=2)
    lns2 = ax.plot(trang,np.abs(cav_v_list[1]),'-', label= 'Cavity Field ('+r'$\Delta f_2$'+')', linewidth=2)
    # # Drive signal is identical for both
    lns3 = ax.plot(trang,np.abs(drive_in_list[0]*mode_dicts[0]['k_drive']),'-', label= 'Drive', linewidth=2)
    ax2 = ax.twinx()
    # # Plot frequency offset as a function of time
    lns4 = ax2.plot(trang, w_offset_list[0]/2.0/np.pi,'-c', label= r'$\Delta f_1$', linewidth=2)
    lns5 = ax2.plot(trang, w_offset_list[1]/2.0/np.pi,'-m', label= r'$\Delta f_2$', linewidth=2)

    # Add all the lines
    lns = lns1 + lns2 + lns3 + lns4 + lns5
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)

    # Format plot
    ax.ticklabel_format(style='sci', axis='y', scilimits=(1,0))
    ax.set_ylim([0,6e5])
    ax.set_xlabel('Time [s]', fontsize=30)
    # Y axis label for voltages
    ax.set_ylabel('Amplitude [V]', fontsize=30)
    # Y axis label for frequency shifts
    ax2.set_ylabel('Frequency [Hz]', fontsize=30)
    ax2.set_ylim(0,delta_f*1.5)
    plt.rc('font',**{'size':20})

    # Show figure
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Plot title
    plt.title("Cavity Test: Fill & Step on Detune Frequency", fontsize=40, y=1.05)

    # # Cavity fields for the two simulation runs

    ln1 = ax.plot(trang, np.unwrap(np.angle(cav_v_list[0])),'-', label= 'Cavity Field ('+r'$\Delta f_1$'+')', linewidth=2)
    ln2 = ax.plot(trang, np.unwrap(np.angle(cav_v_list[1])),'-', label= 'Cavity Field ('+r'$\Delta f_2$'+')', linewidth=2)

    lns = ln1 + ln2
    labs = [l.get_label() for l in lns]
    ax.set_xlabel('Time [s]', fontsize=30)
    ax.set_ylabel(r'$\angle \vec V_{\mu}$ [rad]', fontsize=30)
    plt.rc('font',**{'size':20})
    ax.legend(lns, labs, loc='upper right')

    ax.set_ylim([0,100])

    ax.annotate(r'$d \theta_2/dt=\omega_{\rm d_2}=2\pi \, \Delta f_2\, rad/s$', xy=(0.245, 50), xytext=(0.295, 55),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )

    ax.annotate(r'$d \theta_1/dt=\omega_{\rm d_1}=2\pi \, \Delta f_1\, rad/s$', xy=(0.23, 15), xytext=(0.28, 20),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )


    # Show figure
    plt.show()

    # Return PASS/FAIL boolean (this is not a PASS/FAIL test)
    return True


def unit_cavity():
    """
    Unit test for cavity.c/h
    It performs step responses, measures port couplings. detuning, etc.
    and returns a PASS/FAIL boolean.
    """

    print "\n**** Test Cavity step response..."
    test_step_pass = cavity_test_step()
    test_freqs_pass = cavity_test_freqs()
    test_detune_pass = cavity_test_detune()

    return test_step_pass & test_freqs_pass & test_detune_pass

def perform_tests():
    """
    Perform all unit tests for filter.c/h and cavity.c/h and return a PASS/FAIL boolean.
    """
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
