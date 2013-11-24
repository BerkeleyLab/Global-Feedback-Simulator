#!/usr/bin/python


############################
# 
# unit_tests_components.py
#
# Unit tests for bigger components of the LLRF system
#
# Alejandro F Queiruga
# Daniel Scott Driver
# 2013 LBL
#
############################

#
# TODO: Fill in license
#

import linac
from linac_pretty_print import *
from readjson.readjson import *

import numpy as np
import matplotlib.pylab as plt


# Because octave writes sqrt(-1)=i, and 
# numpy expects sqrt(-1)=j... >.< 
# TODO: submit a patch to numpy that can handle "i"s in loading
i2j = lambda x:np.complex(x.replace("i","j"))


[pa,allaccell,linp,gun] = loadaccelerator("readjson/footest.cfg",
                                          defaultfile="readjson/default.cfg")

Nlinac = len(allaccell)
linss_array = linac.allocate_states(linp.cast(), Nlinac, 3)

linac.state_space_top(gun,linp.cast(),Nlinac, linss_array,3, 
                      pa['Simulation']['dt'], 10,0, None)

linac.deallocate_states(linss_array, Nlinac, 3)
