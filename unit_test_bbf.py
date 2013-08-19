#!/usr/bin/python

#
# Unit tests for Beam Based Feedback Routines
#

import numpy as np
from scipy import linalg

import linac
import full_dc_matrix
import bbf_config
from readjson.readjson import loadall as loadall

[ pa, allaccell, linp, gun, bbfconf, noiseconf] = loadall \
    ("configfiles/NGLS_Accelerator.cfg", \
         "configfiles/bbf_lim.cfg", \
         "configfiles/noise_test.cfg", \
     default_acc="configfiles/default_accelerator.cfg")

bbf,M,U,S,Vh,Mpinv = bbf_config.BBF_Config(pa,gun,linp,len(allaccell),bbfconf)

print Mpinv
