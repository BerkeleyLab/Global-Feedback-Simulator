"""
Unit tests for doublecompress.c/h
"""

import numpy as np
from numpy.random import rand

import readjson.parse_simulation as parseSim

import accelerator as acc
import oct2py


def unit_doublecompress(Verbose=False):
    """ Unit test for doublecompress.c
        Uses Oct2Py to call the octave version double_compressxv.m
        and compares the two output for accuracy. It reads inputs from:
            - doublecompress_test.json to set up data structures for the C routine
            - double_compress_params_octave.sav has the parameter structure for the Octave routine.
        Random values are selected for the input fields testnum times and the outputs are compared.
        If the maximum error is greater than a threshold the test FAILS """

    from get_configuration import Get_SWIG_Simulation

    # Get Simulation Object, including Pointers to C structures
    test_config_files = ["source/configfiles/unit_tests/doublecompress_test.json"]
    sim = Get_SWIG_Simulation(test_config_files, Verbose)

    # Calculate the length of the Linac
    Nlinac = len(sim.linac_list)

    # Load parameters from Octave sav file with stored configuration
    oct_config_file = "./source/configfiles/unit_tests/"
    oct2py.octave.addpath(oct_config_file)
    loadout = oct2py.octave.load("double_compress_params_octave.sav")

    # Extract Octave data structure containing Linacs' configuration
    params = loadout.pdc

    # Allocate C Arrays for doublecompress inputs
    dphivr = acc.double_Array(Nlinac)
    dV_Vvr = acc.double_Array(Nlinac)

    # Outputs
    dc_noise_srcs = acc.Noise_Srcs()
    dcs = acc.Doublecompress_State()
    acc.Doublecompress_State_Allocate(dcs,Nlinac)

    # Pointers to outputs
    Ipk = acc.double_Array_frompointer(dcs.Ipk)
    sz = acc.double_Array_frompointer(dcs.sz)
    dE_E = acc.double_Array_frompointer(dcs.dE_E)
    sd = acc.double_Array_frompointer(dcs.sd)
    dt = acc.double_Array_frompointer(dcs.dt)
    sdsgn = acc.double_Array_frompointer(dcs.sdsgn)
    k = acc.double_Array_frompointer(dcs.k)
    Eloss = acc.double_Array_frompointer(dcs.Eloss)
    dE_Ei = acc.double_Array_frompointer(dcs.dE_Ei)
    dE_Ei2 = acc.double_Array_frompointer(dcs.dE_Ei2)
    cor = acc.double_Array_frompointer(dcs.cor)

    # Number of Simulation runs
    testnum = 30

    # Store Maximum error after each run (initialize at 0)
    maxerr = 0.0

    # Run testnum Simulation runs
    for cnt in range(testnum):

        np.random.seed(1301+cnt)

        # Find random values for Octave Linac configuration
        for l in range(Nlinac):
            params.lamv[0][l]=rand()+.2
            params.Lv[0][l]=100*rand()+10
            params.av[0][l]=(20*rand())+15.0
            params.R56v[0][l]=.002*(rand()-.5)*2.0
            params.T566v[0][l]=.002*(rand()-.5)*2.0
            params.phiv[0][l]=360.0*rand()
            params.s0v[0][l]=3.0*rand()+.5

            # Copy and do unit conversion for C version
            sim.linac_list[l].C_Pointer.lam = params.lamv[0][l]
            sim.linac_list[l].C_Pointer.L = params.Lv[0][l]
            sim.linac_list[l].C_Pointer.a = params.av[0][l]/1000.0
            sim.linac_list[l].C_Pointer.R56 = params.R56v[0][l]
            sim.linac_list[l].C_Pointer.T566 = params.T566v[0][l]
            sim.linac_list[l].C_Pointer.phi = params.phiv[0][l]*np.pi/180.0
            sim.linac_list[l].C_Pointer.s0 = params.s0v[0][l]/1000.0

        # Introduce some pseudo-random numbers into noise source inputs
        dN_N=100*(rand()-.5)*2
        dtg=500*(rand()-.5)*2 #[picoseconds]
        dEg=.001*(rand()-.5)*2  #[Gev]
        dsig_z=params.sz0*(rand()-.5)*2 #[mm]
        dsig_E=100*(rand()-.5)*2 #[%]
        chirp = .00001*(rand()-0)*2*0 #[m]

        # Vectors for Octave code
        dphiv_oct,dV_Vv_oct=np.zeros((2,Nlinac),dtype=float)
        for i in range(Nlinac):
            dphiv_oct[i]=180*(rand()-.5)*2 #[deg]
            dV_Vv_oct[i]=100*(rand()-.5)*2 #[%]

        # Copy and do unit conversion for C version
        dc_noise_srcs.dQ_Q=dN_N/100
        dc_noise_srcs.dtg=dtg/1e12
        dc_noise_srcs.dE_ing=dEg*1e9
        dc_noise_srcs.dsig_z=dsig_z/1000
        dc_noise_srcs.dsig_E=dsig_E/100
        dc_noise_srcs.dchirpt=chirp

        # Copy and do unit conversion for C version
        for i in range(Nlinac):
            # Inputs
            dphivr[i]=dphiv_oct[i]*np.pi/180
            dV_Vvr[i]=dV_Vv_oct[i]/100

        # Call the C routine via SWIG
        acc.Doublecompress(sim.gun.C_Pointer, sim.C_Pointer.linac_net, Nlinac,\
            dc_noise_srcs,dphivr,  dV_Vvr, dcs)

        # Call the Octave routine using Oct2py
        Ipk_o,sz_o,dE_E_o,sd_o,dt_o,sdsgn_o,k_o,Eloss_o\
            ,dE_Ei_o,dE_Ei2_o,cor1_o=oct2py.octave.double_compressxv(params,\
                 dN_N, dtg,dEg,dsig_z,\
                 dsig_E,chirp,\
                 dphiv_oct,dV_Vv_oct,sim.gun.C_Pointer.Q,verbose=False)

        # Subroutine for subtracting and scaling outputs
        # for error calcualtion
        def diffoctcarr(inputoct,inputcarr,scale=1.0):
            N=len(inputoct[0])
            out=np.zeros(N,dtype=float)
            for k in range(N):
                out[k]=(inputoct[0][k]-inputcarr[k]*scale)
            return out

        # Get the errors for the different outputs
        errIpk = abs(diffoctcarr(Ipk_o,Ipk))
        errsz = abs(diffoctcarr(sz_o,sz,1.0e3))
        errdE_E = abs(diffoctcarr(dE_E_o,dE_E,100.0))
        errsd = abs(diffoctcarr(sd_o,sd,100.0))
        errdt = abs(diffoctcarr(dt_o,dt,1.0e12))
        errsdsgn = abs(diffoctcarr(sdsgn_o,sdsgn,100.0))
        errk = abs(diffoctcarr(k_o,k))
        errEloss = abs(diffoctcarr(Eloss_o,Eloss,1.0e-9))
        errdE_Ei = abs(diffoctcarr(dE_Ei_o,dE_Ei,100.0))
        errdE_Ei2 = abs(diffoctcarr(dE_Ei2_o,dE_Ei2,100.0))
        errcor = abs(diffoctcarr([cor1_o[0,1:]],cor))

        # Find maximum error for this run
        temp=np.max([errIpk,
                     errsz,
                     errdE_E,
                     errsd,
                     errdt,
                     errsdsgn,
                     errk,
                     errEloss,
                     errdE_Ei,
                     errdE_Ei2,
                     errcor])

        # Compare to maximum of previous runs and keep highest error
        maxerr=max(temp,maxerr)

    print " After {0} runs with random inputs the".format(testnum)
    print " Maximum difference is {0}\n".format(maxerr)

    # Error threshold
    tol = 1e-8

    # Unit passes if maximum error lower than error threshold
    unit_pass = maxerr<tol

    return unit_pass

def perform_tests():
    """
    Perform all unit tests for doublecompress.c/h and return a PASS/FAIL boolean.
    """
    print "\n****\nTesting Doublecompress...\n"
    dc_pass = unit_doublecompress()
    if (dc_pass):
        result = 'PASS'
    else:
        result = 'FAIL'
    print ">>> " + result

    return dc_pass

if __name__=="__main__":
    perform_tests()
