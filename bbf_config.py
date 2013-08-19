#!/usr/bin/python

#
# Python routines for configuring a BBF_Param
#


import numpy as np
from scipy import linalg

import linac
import full_dc_matrix

#
# bbf_conf is located in pa['BBF'] in the master dictionary from JSON
# it should have the fields:
#  Measurements: a list of strings of which measurements to 
#    include.
#  Controls: The same thing
#  dv: The step to take in V in the numerical Jacobian
#  dphi: The step to take in Phi
#
# TODO: THIS IS WHERE I WAS WORKING AFQ 7/30/13

def BBF_Config(pa,gun,linp_arr,Nlinac):
    Nmeas = 4
    Ncont = 2
    
    bbf_conf = pa['BBF']

    #
    # Dictionaries that are used by the selectors
    #
    measdict = {
        "dE_E" : 0,
        "sz" : 1,
        "dt" : 2,
        "sd" : 3
        }
    contdict = {
        "dv" : 0,
        "dphi" : 1
        }
    connect = pa['Accelerator']['connect']
    linidx = { connect[i]:i for i in range(len(connect)) }
    

    #
    # Build the list of utilized pairs
    #
    used_measured = []
    used_control = []
    mask = np.zeros([Nmeas*Nlinac, Ncont*Nlinac])
    KIdict = {}

    for (LinName,LinConts) in bbf_conf.iteritems():
        if LinName not in connect:
            continue
        lval = linidx[LinName]
        for (ContName,Affected) in LinConts.iteritems():
            cval = contdict[ContName]
            used_control += [(lval,cval)]
            for (LinMeas,Measured) in Affected.iteritems():
                if LinMeas not in connect:
                    continue
                mval = linidx[LinMeas]
                for mm in Measured:
                    used_measured += [(mval,measdict[mm])]
                    mask[Nmeas*mval+measdict[mm],Ncont*lval+cval] = 1
            # Extract Ki
            KIdict[(lval,cval)] = Affected['ki']
    used_measured = sorted(set(used_measured))
    used_control = sorted(set(used_control))
    print KIdict
    KI = np.array([KIdict[x] for x in used_control])
    Um = len(used_measured)
    Uc = len(used_control)
    
    print KI
    print used_measured
    print used_control
    print mask
    #return



    #
    # Make the Jacobian, take the SVD and construct the Pseudoinv
    #
    M = full_dc_matrix.full_dc_matrix(gun,linp_arr, 0.0005,
                       0.0008,Nlinac)
    Mmasked = (M*mask)
    print Mmasked
    Mmasked = Mmasked[
        [Nmeas*A+i for (A,i) in used_measured]
        ,:][:,[Ncont*A+i for (A,i) in used_control ] ]
    #print Mmasked
    # Chop chop
    #print M
    U, S, Vh = linalg.svd(Mmasked,full_matrices=False)

    Mpinv = Vh.conj().T .dot(np.diag(KI)) .dot(np.diag(1/S)) .dot(U.conj().T)
    
    #print U
    #print S
    #print Vh

    #Mpinv = linalg.pinv2(M) #2 uses the SVD
    #Mpinv = Vh.H * KI * inv(S) * U.H;
    # Cut down the Mpinv matrix
    
    #
    # Allocate the Param data structure
    #
    bbf = linac.BBF_Param()
    linac.BBF_Param_Alloc(bbf, Uc, Um)

    #
    # Pack the data into BBF
    #
    idxC = linac.intArray_frompointer(bbf.idx_control)
    idxM = linac.intArray_frompointer(bbf.idx_measured)
    for k in xrange(Uc):
        idxC[2*k+0] = used_control[k][0]
        idxC[2*k+1] = used_control[k][1]
    for k in xrange(Um):
        idxM[2*k+0] = used_measured[k][0]
        idxM[2*k+1] = used_measured[k][1]

    CM = linac.double_Array_frompointer(bbf.Mpinv)
    for i in xrange(Uc):
        for j in xrange(Um):
            CM[Um*i+j] = Mpinv[i,j]
    return [bbf,M,U,S,Vh,Mpinv]
