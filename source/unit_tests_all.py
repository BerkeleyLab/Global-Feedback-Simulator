#!/usr/bin/python

#
# This file collects and performs all of the unit tests in the code.
# Not all of them say passed, some of them only show plots.
#

import matplotlib.pylab as plt

import unit_tests
import unit_tests_components
import unit_test_bbf

#import test_doublecompress.test_doublecompress
import test_doublecompress.test_doublecompress_new
#import test_doublecompress.test_doublecompress_octave

if __name__ == "__main__":
    print "===== Basic tests: unit_tests.py ====="
    ut = unit_tests.perform_tests()
    print "===== LLRF Components: unit_tests_components.py ====="
    plt.figure()
    utc = unit_tests_components.perform_tests()
    print "===== BBF tests: unit_test_bbf.py ====="
    utb = unit_test_bbf.perform_tests()


    if ut and utc and utb:
        print "ooooo ALL TESTS PASSED ooooo"
    else:
        print "xxxxx ALL TESTS FAILED xxxxx"
