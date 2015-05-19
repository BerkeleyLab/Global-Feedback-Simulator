#!/usr/bin/python

#
# A Series of Unit tests for simulation_top.c/h
#

import accelerator as acc

import numpy as np

####################################
#
# Unit tests for step functions in Simulation Top Level
#
####################################

#
# Unit test for Simulation
#

def run_Simulation_test(Tmax, test_files):

    from get_configuration import Get_SWIG_Simulation
    from plotting import plotdata as pd

    # Get Simulation Object, including Pointers to C structures
    sim = Get_SWIG_Simulation(test_files, Verbose=False)

    # Calculate the length of the Linac
    Nlinac = len(sim.linac_list)

    out_filename = "out.dat"

    acc.Simulation_Run(sim.C_Pointer, sim.State, out_filename, 1)

    pd.PlotData('source/plotting/plotdata_test.json', 'source/plotting/columns.json', True)


def unit_Simulation():

    Tmax = 0.5

    test_files = [
        "source/configfiles/unit_tests/doublecompress_test.json",
        "source/configfiles/unit_tests/simulation_test.json"
    ]

    run_Simulation_test(Tmax, test_files)

######################################
#
# Now execute the tests...
#
######################################

def perform_tests():

    # This is not a PASS/FAIL test
    print "\n****\nTesting Simulation Top Level..."
    unit_Simulation()
    print ">>> (Visual inspection only)\n"

    import matplotlib.pylab as plt
    plt.figure()

    return True

if __name__=="__main__":
    plt.close('all')
    perform_tests()
