#!/usr/bin/python

#
# Unit tests for Beam Based Feedback Routines
#

import numpy as np
from scipy import linalg

import linac
import full_dc_matrix
import bbf_config
from readjson.loadconfig import LoadConfig as LoadConfig


def perform_tests():
    #[ pa, allaccell, linp, gun, bbfconf, noiseconf] = loadall \
    #    ("configfiles/NGLS_Accelerator.cfg", \
    #         "configfiles/bbf_lim.cfg", \
    #         "configfiles/noise_test.cfg", \
    #         default_acc="configfiles/default_accelerator.cfg")
    [ pa, linp_pylist, linp, gun, bbf, nsrc ] = LoadConfig(
        ["configfiles/default_accelerator.cfg", \
             "configfiles/NGLS_Accelerator.cfg", \
             "configfiles/bbf_causal.cfg", \
             "configfiles/noise_test.cfg"]
        )
    #print pa
    bbf,M,U,S,Vh,Mpinv = bbf_config.BBF_Config(pa,gun,linp,len(linp_pylist))

    print Mpinv
    return True

if __name__=="__main__":
    perform_tests()
