#!/usr/bin/python

#
# A Series of Unit tests for rf_station.c/h
#

import accelerator as acc

import numpy as np
import matplotlib.pylab as plt
import scipy.linalg as linalg
from scipy import stats

####################################
#
# Unit tests for step functions in RF Station
#
####################################

#
# Unit test for Phase_Shift
#
def unit_phase_shift():

    trang = np.arange(0.0, 1.0, 0.01)

    nt = len(trang)

    signal_in = np.exp(4j*np.pi*trang)
    thetas = np.arange(0,np.pi,np.pi/6.0)

    out = np.zeros(nt,dtype=np.complex)
    out_phase = np.zeros(nt,dtype=np.double)
    in_phase = np.zeros(nt,dtype=np.double)

    # Initialize PASS variable
    phase_shift_pass = True
    
    for theta in thetas:
        for i in xrange(nt):
            out[i] = acc.Phase_Shift(signal_in[i],theta)
            out_phase[i] = np.angle(out[i])
            in_phase[i] =  np.angle(signal_in[i])

        # Calculate phase difference
        # Unwrap phase
        out_phase_unw = np.unwrap(out_phase)
        in_phase_unw = np.unwrap(in_phase)
        # Average phase
        out_phase_avg = np.average(out_phase_unw)
        in_phase_avg = np.average(in_phase_unw)
        # Calculate measured phase difference
        theta_meas = np.average(out_phase_avg - in_phase_avg)
        # Calculate error
        error = theta_meas - theta

        # If error > threshold, turn phase_shift_pass to False
        pass_now = error < 1e-8
        phase_shift_pass = phase_shift_pass & pass_now

        # Plot and record c parameter for label
        label_text = r'$\theta=$%.1f' % theta

        plt.plot(trang,out.real, label=label_text)
    
    plt.title("Phase shift test", fontsize=40, y=1.01)
    plt.xlabel('Time [s]', fontsize=30)
    plt.ylabel('Real Part [Normalized]', fontsize=30)
    plt.legend(loc='upper right')
    plt.show()

    # Return PASS/FAIL
    return phase_shift_pass

#
# Unit test for step_fpga
#

def unit_fpga(Tstep=0.01):
    
    # Get a pointer to the FPGA C-structure
    fpga = acc.FPGA()

    # Initial settings
    kp = 5.0
    ki = 3.0
    set_point_start = 0.0 + 0.0j  # Start with 0
    set_point_step = 1.0 + 0.0j  # Apply step on real component
    cav_in = 0.0 + 0.0j
    open_loop = 0   # Closed-loop setting (apply control to drive signal)

    # Fill in C data structure with settings
    acc.FPGA_Allocate_In(fpga, kp, ki, set_point_start, Tstep)

    # Allocate State structure and initialize
    fpga_state = acc.FPGA_State()
    fpga_state.drive = 0.0 + 0.0j
    fpga_state.state = 0.0 + 0.0j
    
    # Total simulation time for test (seconds)
    Tmax = 2

    # Buld a time axis
    trang = np.arange(0.0,Tmax,Tstep)

    nt = len(trang) # Number of points

    # Initialize complex vectors
    error = np.zeros(nt,dtype=np.complex)
    drive = np.zeros(nt,dtype=np.complex)
    state = np.zeros(nt,dtype=np.complex)

    # Determine when to step set-point
    sp_step = int(nt*0.1)

    # Run time-series simulation
    for i in xrange(nt):
        
        # Apply step on set-point at sp_step simulation time
        if i == sp_step: fpga.set_point = set_point_step

        # Call step_PI_fpga and record signals of interest
        error[i] = acc.step_PI_fpga(fpga, cav_in, fpga_state, open_loop)
        drive[i] = fpga_state.drive
        state[i] = fpga_state.state

    # Measure kp and ki
    kp_measured = (drive[sp_step] - drive[sp_step-1]-set_point_step*ki*Tstep)/set_point_step
    kp_text = r'$k_{\rm p}$'+' set to: %.1f, measured: %.1f'% (np.real(kp_measured), np.real(kp))
    ki_measured,b = np.polyfit(trang[sp_step:-1], drive[sp_step:-1], 1) # Slope
    ki_text = r'$k_{\rm i}$'+' set to: %.1f, measured: %.1f'% (np.real(ki), np.real(ki_measured))

    # Plot
    plt.plot(trang,np.real(drive),'-r', label='Drive')
    plt.plot(trang,np.real(error),'-+', label='Error')
    plt.plot(trang,np.real(state),'-*', label='Integrator state')

    plt.title("FPGA Unit Test", fontsize=40, y=1.01)
    plt.xlabel('Time [s]', fontsize=30)
    plt.ylabel('Amplitude [Unitless]', fontsize=30)
    plt.legend(loc='upper right')
    plt.ylim([0,15])

    # Add text with results on plot
    plt.text(1,6, kp_text, verticalalignment='top', fontsize=30)
    plt.text(1,4.8, ki_text, verticalalignment='top', fontsize=30)
    plt.rc('font',**{'size':15})

    plt.show()

    # Evaluate PASS/FAIL for unit test (error must be 0)
    kp_error = np.abs(kp_measured - kp)
    ki_error = np.abs(ki_measured - ki)

    if (kp_error <1e-12) and (ki_error <1e-12): 
        unit_fpga_pass = True
    else:  
        unit_fpga_pass = False

    # PASS == True
    return unit_fpga_pass

#
# Unit test for saturation
#
def unit_saturate():

    # Iterate over harshness parameter c
    for c in 10.0**np.arange(0.0,2,0.4):
        # Input vector
        inp = np.arange(0.0,10.0,0.1,dtype=np.complex)
        # Output vector
        oup = np.zeros(inp.shape,dtype=np.complex)

        # Sweep input
        for i in xrange(len(inp)):
            oup[i] = acc.Saturate(inp[i],c)

        # Plot and record c parameter for label
        label_text = "c = %.1f" % c
        plt.plot(inp.real,oup.real, label=label_text)

    # Plot
    plt.title("Saturation Test", fontsize=40, y=1.01)
    plt.xlabel(r'$|\vec V_{\rm in}| [\rm V]$', fontsize=30)
    plt.ylabel(r'$|\vec V_{\rm out}| [\rm Normalized]$', fontsize=30)

    # Add equation text on plot
    plt.text(5,0.4, r'$\vec V_{\rm out} = \vec V_{\rm in} \cdot \left(1 + |\vec V_{\rm in}|^c\right)^{-1/c}$', fontsize=30)

    # Adjust legend placement and vertical axis limits
    plt.legend(loc=7)
    plt.ylim([0,1.1])
    
    # Show plot
    plt.show()


######################################
#
# Now execute the tests...
#
######################################

def perform_tests():
    print "\n****\nTesting Phase_Shift..."
    phase_shift_pass = unit_phase_shift()
    if (phase_shift_pass):
        result = 'PASS' 
    else: 
        result = 'FAIL'
    print ">>> " + result 

    plt.figure()
    
    print "\n****\nTesting FPGA PI controller..."
    fpga_pass = unit_fpga()
    if (fpga_pass):
        result = 'PASS' 
    else: 
        result = 'FAIL'
    print ">>> " + result 
    
    plt.figure()

    # This is not a PASS/FAIL test
    print "\n****\nTesting Saturate..."
    unit_saturate()
    print ">>> (Visual inspection only)\n" 

    return fpga_pass & phase_shift_pass

if __name__=="__main__":
    plt.close('all')
    perform_tests()
