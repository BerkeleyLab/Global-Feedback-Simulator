#!/usr/bin/python

####################################################################
#
# unit test for full_dc_matrix
# Uses Oct2Py to call the octave version double_compressxv.m
# and compares the two output for accuracy
# 
# reads inputs from footest.cfg and default.cfg to set up
# data stuctures for the c routine
# double_compress_params_octave.sav has the parameter structure for
# the octave routine
#
# Random values are selected for the input fields testnum times
# and the outputs are compared. If the maximum difference is greater 
# than tol the test is considered a failure
#
# LBL July 2013
# Daniel Driver
#
####################################################################


print "Starting unit test for dc Matrix..."

import sys
sys.path.append("../")
import numpy as np
from numpy.random import rand
from readjson.readjson import loadaccelerator
import linac
import oct2py
from full_dc_matrix import full_dc_matrix

readfile=True
if(readfile):

    a,accelout,linp_arr,gun=loadaccelerator("footest.cfg", 
                                            defaultfile="default.cfg")
    Nlinac=len(accelout)
    
    loadout=oct2py.octave.call('load',
                               "double_compress_params_octave.sav",verbose=False)
    params=loadout.pdc
    Nlinmax=Nlinac

else:#curently not working set readfile to True

    Nlinmax=10
    gun=linac.Gun_Param()
    linp_arr=linac.Linac_Param_Array(Nlinmax)

    for l in range(Nlinmax):
        lin=linac.Linac_Param()
        print 'here'
        linp_arr[l]=lin

    params=oct2py.Struct()
    params.Ev=[None]
    params.lamv=[None]
    params.Lv=[None]
    params.av=[None]
    params.R56v=[None]
    params.T566v=[None]
    params.phiv=[None]
    params.s0v=[None]

    params.Ev[0]=Nlinmax*[None]
    params.lamv[0]=Nlinmax*[None]
    params.Lv[0]=Nlinmax*[None]
    params.av[0]=Nlinmax*[None]
    params.R56v[0]=Nlinmax*[None]
    params.T566v[0]=Nlinmax*[None]
    params.phiv[0]=Nlinmax*[None]
    params.s0v[0]=Nlinmax*[None]

Q=gun.Q
g=1

maxerr=0
testnum=10
for cnt in range(testnum):
    np.random.seed(1301+cnt)

    #find random values
    for l in range(Nlinac):
        params.lamv[0][l]=rand()+.2
        params.Lv[0][l]=100*rand()+10
        params.av[0][l]=(20*rand())+15.0
        params.R56v[0][l]=.002*(rand()-.5)*2.0
        params.T566v[0][l]=.002*(rand()-.5)*2.0
        params.phiv[0][l]=360.0*rand()
        params.s0v[0][l]=3.0*rand()+.5

        #copy and do unit conversion for c version
        linp_arr[l].lam=params.lamv[0][l]
        linp_arr[l].L=params.Lv[0][l]
        linp_arr[l].a=params.av[0][l]/1000.0
        linp_arr[l].R56=params.R56v[0][l]
        linp_arr[l].T566=params.T566v[0][l]
        linp_arr[l].phi=params.phiv[0][l]*np.pi/180.0
        linp_arr[l].s0=params.s0v[0][l]/1000.0



    #find matrix M from Octave
    bbf=oct2py.octave.call('bbf_config',params,Q,g,verbose=False)

    #set up for the python version
    Nlinac=len(accelout)
    ki=np.ones(Nlinac)*g
    dv=.05/100.0
    dphi=.05*np.pi/180.0

    M=full_dc_matrix(gun,linp_arr,dv,dphi,Nlinac)

    #scale the outputs so that derivative units match
    scalenum=[100.0,1000.0,1.0e12,100.0]
    scaleden=[100.0,180.0/np.pi]

    Nbbfval=4;
    Nsetpnt=2;
    for k in range(Nbbfval):
        for j in range(Nsetpnt):
            M[k::Nbbfval,j::Nsetpnt]*=scalenum[k]/scaleden[j]
    #for k in range(Nbbfval*Nlinac):
    #    for j in range(Nsetpnt*Nlinac):
    #         M[k,j]*=scalenum[k%Nbbfval]/scaleden[j%Nsetpnt]

    #compare ouputs to ensure that they are the same
    err=np.max(abs(bbf.M-M))
    maxerr=max(maxerr,err)

tol=1.0e-10   
print "Maximum error: ", maxerr 
if maxerr<tol:
    print "******************************"
    print "*********    PASS    *********"
    print "******************************"
else:
    print "******************************"
    print "*********    FAIL    *********"
    print "******************************"
