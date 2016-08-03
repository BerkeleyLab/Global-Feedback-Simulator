"""
A Series of Unit tests for ALL of the C back-end components.
It runs all of the unit tests one by one, collects the PASS/FAIL results
and provides a final PASS/FAIL boolean indicating if all tests passed or not.
"""

import matplotlib.pylab as plt

import cavity_test
import rf_station_test
import cryomodule_test
import doublecompress_test
import simulation_test

# import unit_test_bbf

if __name__ == "__main__":

    print "===== Cavity tests: cavity_test.py ====="
    cavity_pass = cavity_test.perform_tests()
    if (cavity_pass):
        result = 'PASS'
    else:
        result = 'FAIL'
    print "===== Cavity tests >>> " + result + " =====\n"

    print "===== RF Station tests: rf_station_test.py ====="
    rf_station_pass = rf_station_test.perform_tests()
    if (rf_station_pass):
        result = 'PASS'
    else:
        result = 'FAIL'
    print "===== RF Station tests >>> " + result + " =====\n"

    print "===== Cryomodule tests: cryomodule_test.py ====="
    cryomodule_pass = cryomodule_test.perform_tests()
    if (cryomodule_pass):
        result = 'PASS'
    else:
        result = 'FAIL'
    print "===== Cryomodule tests >>> " + result + " =====\n"

    print "===== DoubleCompress test: doublecompress_test.py ====="
    doublecompress_pass = doublecompress_test.perform_tests()
    if (doublecompress_pass):
        result = 'PASS'
    else:
        result = 'FAIL'
    print "===== DoubleCompress tests >>> " + result + " =====\n"

    # print "===== Simulation Top Level test: simulation_test.py ====="
    # sim_test_pass = simulation_test.perform_tests()
    # if (sim_test_pass):
    #     result = 'PASS'
    # else:
    #     result = 'FAIL'
    # print "===== Simulation Top Level tests >>> " + result + " =====\n"

    # print "===== LLRF Components: unit_tests_components.py ====="
    # utc = unit_tests_components.perform_tests()
    # plt.figure()
    # print "===== BBF tests: unit_test_bbf.py ====="
    # utb = unit_test_bbf.perform_tests()

    # if ut and utc and utb:
    if rf_station_pass & cavity_pass & cryomodule_pass & doublecompress_pass:
        print "ooooo ALL TESTS PASSED ooooo"
    else:
        print "xxxxx ALL TESTS FAILED xxxxx"
