"""
A Series of Unit tests for simulation_top.c/h.
Unit tests for step functions in Simulation Top Level.
"""

import accelerator as acc

import numpy as np

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
    """
    Unit test for simulation_top.c/h.
    It exercises the full simulation cycle:
        - JSON configuration parsing,
        - Conversion to Python Objects,
        - Run numerical simulation of a complete model,
        - Generation of time-series data,
        - Plotting of time-series data.
    """

    Tmax = 0.5

    test_files = [
        "source/configfiles/unit_tests/doublecompress_test.json",
        "source/configfiles/unit_tests/simulation_test.json"
    ]

    run_Simulation_test(Tmax, test_files)

def perform_tests():
    """
    Perform all unit tests for simulation_top.c/h and return PASS boolean (qualitative test).
    """

    # This is not a PASS/FAIL test
    print "\n****\nTesting Simulation Top Level..."
    unit_Simulation()
    print ">>> (Visual inspection only)\n"

    import matplotlib.pylab as plt
    plt.figure()

    return True

if __name__ == "__main__":
    plt.close('all')
    perform_tests()
