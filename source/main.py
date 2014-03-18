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

import linac, linac_pretty_print

from readjson.loadconfig import LoadConfig
from readjson.readjson import readentry

#
# Call this routine to script!
#
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
    
    #for l in xrange(len(linp_pylist)):
    #    print "*****LINAC ",l,"*****"
    #    print linac_pretty_print.lintostr(linp_pylist[l])

    #
    # Allocate the delay buffer
    #
    # How much history data is needed... 
    # TODO: This should probably be in the JSON parameters
    Nhist = 3
    linss_array = linac.allocate_states(linp_arr.cast(), Nlinac, Nhist)


    #
    # Pull out the number of step, dt and frequency of outputdata from the config file
    #
    Nstep=int(readentry(confdict,confdict['Simulation']['Nstep'],localdic=confdict['Simulation']))
    dt=readentry(confdict,confdict['Simulation']['dt'],localdic=confdict['Simulation'])
    Outputfreq=int(readentry(confdict,confdict['Simulation']['Outputfreq'],localdic=confdict['Simulation']))
    

    #
    # Actually run the simulation
    #
    linac.state_space_top(gun,linp_arr.cast(),Nlinac, linss_array,Nhist, 
                          bbf, nsrc, dt, Nstep,0,
                          OutputFile,Outputfreq)

    #
    # Deallocate all c-objects needed
    #
    linac.deallocate_states(linss_array, Nlinac, Nhist)
    for linp in linp_pylist:
        linac.Linac_Deallocate(linp)
    linac.BBF_Param_Free(bbf)




#
# If this file is called as main, run it from command line arguments,
# or some defualt settings.
#
if __name__=="__main__":
    outputfile= sys.argv[1] #"outputdata/test.dat"

    ConfigFiles = sys.argv[2:]
    
    print sys.argv
    DoSimulation(ConfigFiles,outputfile)
