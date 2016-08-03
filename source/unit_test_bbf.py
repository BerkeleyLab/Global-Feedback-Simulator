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
    [pa, linp_pylist, linp, gun, bbf, nsrc] = LoadConfig(
        ["source/configfiles/NGLS/default_accelerator.cfg",
         "source/configfiles/NGLS/NGLS_Accelerator.cfg",
         "source/configfiles/NGLS/bbf_causal.cfg",
         "source/configfiles/NGLS/noise_test.cfg"]
    )
    # print pa
    bbf, M, U, S, Vh, Mpinv = bbf_config.BBF_Config(pa, gun, linp, len(linp_pylist))

    print Mpinv
    return True

if __name__ == "__main__":
    perform_tests()
