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

    cryo, cryo_state, Tstep, fund_mode_dicts, rf_station_pointers, mechMode_pointers =  Get_SWIG_Cryomodule(test_file, Verbose=False)

    # Create time vector
    trang = np.arange(0,Tmax,Tstep)
    
    # Number of points
    nt = len(trang)

    # Initialize vectors for test
    cav_v = np.zeros(nt,dtype=np.complex)   # Overall cavity accelerating voltage
    fpga_drive_out = np.zeros(nt,dtype=np.complex)
    set_point = np.zeros(nt,dtype=np.complex)
    error = np.zeros(nt,dtype=np.complex)
    
    rf_station = acc.Get_RF_Station(cryo, 0)
    rf_state = acc.Get_RF_State(cryo_state, 0)

    # Get Pointers to RF States to access waveforms
    # rf_states = []
    # rf_stations = []
    # for idx in range(cryo.n_rf_stations):
    #     rf_states.append(acc.Get_RF_State(cryo_state, idx))
    #     rf_stations.append(acc.Get_RF_Station(cryo, idx))
    
    # # Get Pointers to MechMode States to access waveforms
    # mechMode_states = []
    # for idx in range(cryo.n_rf_stations):
    #     mechMode_states.append(acc.Get_MechMode_State(cryo_state, idx))

    # Run Numerical Simulation
    for i in xrange(1,nt):
        acc.Cryomodule_Step(cryo, cryo_state, 0.0, 0.0, 0.0, 0.0, 0)
        cav_v[i] = rf_state.cav_state.V
        set_point[i] = rf_station.fpga.set_point
        fpga_drive_out[i] = rf_state.fpga_state.drive
        error[i] = rf_state.fpga_state.err

    fund_k_probe = fund_mode_dicts[0]['k_probe']
    fund_k_drive = fund_mode_dicts[0]['k_drive']

    plt.plot(trang,np.abs(cav_v), '-', label='Cavity voltage')
    plt.plot(trang,np.abs(fpga_drive_out*fund_k_drive), label='Drive')
    plt.plot(trang,np.abs(set_point/fund_k_probe), label='Set-point')
    plt.plot(trang,np.abs(error/fund_k_probe), label='Cavity field error')
    
    plt.title('Cryomodule Test', fontsize=40, y=1.01)
    plt.xlabel('Time [s]', fontsize=30)
    plt.ylabel('Amplitude [V]', fontsize=30)
    plt.legend(loc='upper right')

    # Show plot
    plt.show()

def unit_Cryomodule():

    Tmax = 0.05
    
    # test_file = "source/configfiles/unit_tests/cavity_test_step1.json"
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