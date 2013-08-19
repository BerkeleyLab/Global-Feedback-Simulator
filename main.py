#!/usr/bin/python


#
# NGLS Accelerator model
#
# main.py
#
# This is the top level script for the NGLS model.
# It loads data from a series of JSON files,
# evaulates the fields in the imported dictionary,
# configures and allocates the associated swig C data structures,
# calls the C code with the swig wrappers, which write to a file
# and deallocates the C datastructures.
#

import sys

import numpy as np

import linac

from readjson.loadconfig import LoadConfig
from readjson.readjson import readentry

def DoSimulation(ConfigFiles,OutputFile):
    #
    # Read in the configuration files
    #  confdict:    is the raw dictionary from the json file
    #  linp_pylist: is a list of pointers to the linarry objects containing 
    #                parameters relavant to an individual linac
    #  linp_array:  is a ptr to an c-object with the pointers in PyListAccel
    #  gun:         is a pointer to a c-object with parameters related to the gun
    #  bbf:         is a ptr to a c-object with information need to do beam 
    #                base feedback. It includes ptrs to Mpinv, input and output
    #                vectors and indices of measurements and inputs used in the
    #                feadback operation.
    #
    confdict, linp_pylist,linp_arr,gun, bbf, nsrc = \
        LoadConfig(ConfigFiles)
    Nlinac = len(linp_pylist)
    
    
    #
    # Allocate the delay buffer
    #
    linss_array = linac.allocate_states(linp_arr.cast(), Nlinac, 3)


    Nstep=int(readentry(confdict,confdict['Simulation']['Nstep'],localdic=confdict['Simulation']))
    dt=readentry(confdict,confdict['Simulation']['dt'],localdic=confdict['Simulation'])
    #
    # Actually run the simulation
    #
    linac.state_space_top(gun,linp_arr.cast(),Nlinac, linss_array,3, 
                          bbf, nsrc, dt, Nstep,0,
                          OutputFile)
    #
    # Deallocate all c-objects needed
    #
    linac.deallocate_states(linss_array, Nlinac, 3)
    for linp in linp_pylist:
        linac.Linac_Deallocate(linp)


#
# If this file is called as main, run it from command line arguments,
# or some defualt settings.
#
if __name__=="__main__":
    outputfile= sys.argv[1] #"outputdata/test.txt"

    ConfigFiles = sys.argv[2:]
#[
#        "configfiles/default_accelerator.cfg",
#        "configfiles/NGLS_Accelerator.cfg",
#        "configfiles/bbf_causal.cfg",
#        "configfiles/noise_test.cfg"
#        ]
    print sys.argv
    DoSimulation(ConfigFiles,outputfile)
