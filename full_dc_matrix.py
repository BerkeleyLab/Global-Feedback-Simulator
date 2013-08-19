#!/usr/bin/python

import numpy as np
import linac

#
# Functions for building a linearization matrix
#  of the function doublecompress.c
#
# Full_dc_matrix:
#     This function builds a jacobian/sensitivty
#     matrix from all accelerator inputs (Voltage
#     and phase) to all accelrator ouputs (dE_E,sz,dt,sd). 
#     this forms a (4*Nlinac)x(2*Nlinac) matrix
#     The derivatives to calcualte the jacobian
#     are calculated numerically using a midpoint 
#     approximation and evaluations of 
#     doublecompress_new.c
#
#     To make sure everyone is speaking the same language,
#     synonyms for this matrix include:
#       Jacobian, Tangent, Derivative, Sensitivity, Small Signal,
#       Gradient, etc...
#       dc(x) =~ dc(x0) + [M]*(x-x0) + O((x-x0)^2) ...
#     Jacobian matrix is the term used in this file.
#
# Inputs:
#     gun: linac.gun data structure storing parameters
#           of the electron gun. It supplies charge(gun.Q)
#           for use in double comprress. Should be
#           read in using configuation files
#     linp_arr:  Array of points to linac objects
#                which store parameters for the linac
#                This object is passed to doublecompress and
#                should come from reading in configuration
#                files.
#     dv:  This is the deviation by which the derivative 
#          will be calculated for the voltage error/voltage. 
#          double compress is evaluted with dV_Vvr+-dv  [fraction]
#     dphi:This is the deviation by which the derivative 
#          will be calculated for the phase error. 
#          double compress is evaluted with dphivr+-dphi [radians] 
#     Nlinac: Is an integer which is the number of linacs
#             in the accelator
# Output:The output matrix is the sensitivity matrix
#        The rows map to the meausrement by
#          [dE_E[0],sz[0],dt[0],sd[0],dE_E[1],sz[1]...
#            ...dE_E[Nlinac],sz[Nlinac],dt[Nlinac],sd[Nlinac]
#        The columns map as
#          [dv[0],dphi[0],dv[1].....dv[Nlinac],dphi[Nlinac]
#        In other words, by multipying output matrix(M) by a vector
#        arranged as the columns we get an output in the form
#        of a row with output arranged according to the mapping of the rows
#        Units are SI except energy is in eV.
#        doublecompress.c has a more detailed description.
#
#
#     Daniel Driver
#     Summer 2013  
#
# EDIT: afq 7/29/13



#
# This routine packs the outputs from a double_compress state data structure
# into a single numpy array, numpyarray, in the order 
#    [dE_E[0],sz[0],dt[0],sd[0],dE_E[1],sz[1]...
#           ...dE_E[Nlinac],sz[Nlinac],dt[Nlinac],sd[Nlinac]
#    used by full_dc_matrix
#
def dcs_to_numpy(dE_E,sz,dt,sd, numpyarray,Nlinac,Nmeasurements):
    for n in xrange(Nlinac):
        numpyarray[Nmeasurements*n+0]=dE_E[n]
        numpyarray[Nmeasurements*n+1]=sz[n]
        numpyarray[Nmeasurements*n+2]=dt[n]
        numpyarray[Nmeasurements*n+3]=sd[n]

#
# This routine builds the Jacobian matrix to double compress numerically
# around the point in the phase space specified in the parameters.
#        
def full_dc_matrix(gun,linp_arr,dv,dphi,Nlinac):
    Nmeasurements=4
    Ncontrols=2

    #
    # Set up empty stuctures to call double compress
    #
    dynp=linac.Dynamic_Param()
    dcs=linac.Doublecompress_State()
    linac.Doublecompress_State_Alloc(dcs,Nlinac)
    
    #
    # Create Swig C-Array wrappers
    #
    sz= linac.double_Array_frompointer(dcs.sz)
    dE_E = linac.double_Array_frompointer(dcs.dE_E)
    sd= linac.double_Array_frompointer(dcs.sd)
    dt= linac.double_Array_frompointer(dcs.dt)
    dphivr= linac.double_Array(Nlinac)
    dV_Vvr= linac.double_Array(Nlinac)

    #
    # Initialize data structures to zero for good measure.
    #
    for j in xrange(Nlinac):  
        dV_Vvr[j]=0.0
        dphivr[j]=0.0
    dynp.dQ_Q=0.0
    dynp.dtg=0.0
    dynp.dE_ing=0.0
    dynp.dsig_z=0.0
    dynp.dsig_E=0.0
    dynp.dchirp=0.0

    #
    # Initialize M matrix and the output perturbation vectors
    #
    M=np.zeros((Nmeasurements*Nlinac,Ncontrols*Nlinac),
               dtype=np.double)
    dcs_plus_dv=np.zeros(Nlinac*Nmeasurements,dtype=np.double)
    dcs_minus_dv=np.zeros(Nlinac*Nmeasurements,dtype=np.double)
    dcs_plus_dphi=np.zeros(Nlinac*Nmeasurements,dtype=np.double)
    dcs_minus_dphi=np.zeros(Nlinac*Nmeasurements,dtype=np.double)


    #
    # Find the jacobian matrix numerically using midpoint approximation
    # to the first derivative.
    #
    for j in xrange(Nlinac):
        #
        # Perform the calcutions needed for midpoint rule
        #
        # First do the voltage first add dv and then subtract dv
        dV_Vvr[j]=dv
        linac.doublecompress_new(gun,linp_arr.cast(),Nlinac,
                                 dynp,dphivr,dV_Vvr,
                                 dcs)
        dcs_to_numpy(dE_E,sz,dt,sd,dcs_plus_dv,Nlinac,Nmeasurements)

        dV_Vvr[j]=-dv
        linac.doublecompress_new(gun,linp_arr.cast(),Nlinac,
                                 dynp,dphivr,dV_Vvr,
                                 dcs)
        dcs_to_numpy(dE_E,sz,dt,sd,dcs_minus_dv,Nlinac,Nmeasurements)
        dV_Vvr[j]=0.0  #set back to zero
        
        # Do the phase variations first add dphi and do minus dphi
        dphivr[j]=dphi
        linac.doublecompress_new(gun,linp_arr.cast(),Nlinac,
                                 dynp,dphivr,dV_Vvr,
                                 dcs)
        dcs_to_numpy(dE_E,sz,dt,sd,dcs_plus_dphi,Nlinac,Nmeasurements)

        dphivr[j]=-dphi
        linac.doublecompress_new(gun,linp_arr.cast(),Nlinac,
                                 dynp,dphivr,dV_Vvr,
                                 dcs)
        dcs_to_numpy(dE_E,sz,dt,sd,dcs_minus_dphi,Nlinac,Nmeasurements)
        dphivr[j]=0.0  #set bck to zero
        
        # Calculate derivatives and place into the columns of matix M
        M[:,Ncontrols*j+0]=(dcs_plus_dv- dcs_minus_dv)/(2.0*dv) 
        M[:,Ncontrols*j+1]=(dcs_plus_dphi- dcs_minus_dphi)/(2.0*dphi) 

    return M

