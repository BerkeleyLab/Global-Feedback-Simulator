#def double_compress_unit_test():
####################################################################
#
# unit test for doublecompress.c
# Uses Oct2Py to call the octave version double_compressxv.m
# and compares the two output for accuracy
# 
# reads inputs from footest.cfg and default.cfg to set up
# data stuctures for the c routine
# double_compress_params_octave.sav has the parameter structure for
# the octave routine
#
# Random values are selected for the input fields testnum times
# and the ouputs are compared. If the maximum difference is greater 
# than tol the test is considered a failure
#
#LBL July 2013
#Daniel Driver
####################################################################


print "Starting unit test for double Compress..."

import sys
sys.path.append("../")
import numpy as np
from numpy.random import rand
from readjson.readjson import loadaccelerator
import linac
import oct2py


readfile=True
if(readfile):
    a,accelout,linp_arr,gun=loadaccelerator("footest.cfg", defaultfile="default.cfg")
    Nlinac=len(accelout)

    loadout=oct2py.octave.call('load',"double_compress_params_octave.sav")
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


#allocate the c arrays
dphivr= linac.double_Array(Nlinmax)
dV_Vvr= linac.double_Array(Nlinmax)
    #outputs
Ipk= linac.double_Array(Nlinmax)
sz= linac.double_Array(Nlinmax)
dE_E= linac.double_Array(Nlinmax)
sd= linac.double_Array(Nlinmax)
dt= linac.double_Array(Nlinmax)
sdsgn= linac.double_Array(Nlinmax)
k= linac.double_Array(Nlinmax)
Eloss= linac.double_Array(Nlinmax)
dE_Ei= linac.double_Array(Nlinmax)
dE_Ei2= linac.double_Array(Nlinmax)
cor= linac.double_Array(Nlinmax)

testnum=60
maxerr=0.0
track_worst=False
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

        #copy and do unit conversion for c verion
        linp_arr[l].lam=params.lamv[0][l]
        linp_arr[l].L=params.Lv[0][l]
        linp_arr[l].a=params.av[0][l]/1000.0
        linp_arr[l].R56=params.R56v[0][l]
        linp_arr[l].T566=params.T566v[0][l]
        linp_arr[l].phi=params.phiv[0][l]*np.pi/180.0
        linp_arr[l].s0=params.s0v[0][l]/1000.0


    dN_N=100*(rand()-.5)*2
    dtg=500*(rand()-.5)*2 #[picoseconds]
    dEg=.001*(rand()-.5)*2  #  [Gev]
    dsig_z=params.sz0*(rand()-.5)*2# [mm]
    dsig_E=100*(rand()-.5)*2 # [%]
    chirp = .00001*(rand()-0)*2*0 #[m]

    #vectors to go to the octave code
    dphiv_oct,dV_Vv_oct=np.zeros((2,Nlinac),dtype=float)
    for i in range(Nlinac):
        dphiv_oct[i]=180*(rand()-.5)*2 #[deg]
        dV_Vv_oct[i]=100*(rand()-.5)*2 #[%]

    #N=1.8724e9
    #Q=N/6.241E18 #charge for octave code. C takes in N (number of particles)
    #gun.Q=Q

    #print Q
    #put values in c structures and set outputs to zero
    for i in range(Nlinac):
        #inputs
        dphivr[i]=dphiv_oct[i]*np.pi/180
        dV_Vvr[i]=dV_Vv_oct[i]/100

        Ipk[i]=0.0  
        sz[i]=0.0
        dE_E[i]=0.0
        sd[i]=0.0
        dt[i]=0.0
        sdsgn[i]=0.0 
        k[i]=0.0 
        Eloss[i]=0.0 
        dE_Ei[i]=0.0
        dE_Ei2[i]=0.0  
    for i in range(Nlinac):
        cor[i]=0.0 

    #run the C routine
    linac.doublecompress_octave(gun,linp_arr.cast(),Nlinac,
                           dN_N/100,  dtg/1e12,  dEg*1e9,
                           dsig_z/1000,   dsig_E/100,  chirp,
                           dphivr,  dV_Vvr,
                           Ipk,  sz,  dE_E,
                           sd ,  dt,  sdsgn, 
                           k,  Eloss,  dE_Ei,
                           dE_Ei2,  cor)

    #call the octave routine using Oct2py
    Ipk_o,sz_o,dE_E_o,sd_o,dt_o,sdsgn_o,k_o,Eloss_o\
        ,dE_Ei_o,dE_Ei2_o,cor1_o=oct2py.octave.call\
        ('double_compressxv',params,\
             dN_N, dtg,dEg,dsig_z,\
             dsig_E,chirp,\
             dphiv_oct,dV_Vv_oct,gun.Q,verbose=False)

    #subroutine for subtracting and scaling outputs
    #for error calcualtion
    def diffoctcarr(inputoct,inputcarr,scale=1.0):
        N=len(inputoct[0])
        out=np.zeros(N,dtype=float)
        for k in range(N):
            out[k]=(inputoct[0][k]-inputcarr[k]*scale)
        return out

    #get the errors for the different outputs
    errIpk=abs(diffoctcarr(Ipk_o,Ipk))
    errsz=abs(diffoctcarr(sz_o,sz,1.0e3))
    errdE_E=abs(diffoctcarr(dE_E_o,dE_E,100.0))
    errsd=abs(diffoctcarr(sd_o,sd,100.0))
    errdt=abs(diffoctcarr(dt_o,dt,1.0e12))
    errsdsgn=abs(diffoctcarr(sdsgn_o,sdsgn,100.0))
    errk=abs(diffoctcarr(k_o,k))
    errEloss=abs(diffoctcarr(Eloss_o,Eloss,1.0e-9))
    errdE_Ei=abs(diffoctcarr(dE_Ei_o,dE_Ei,100.0))
    errdE_Ei2=abs(diffoctcarr(dE_Ei2_o,dE_Ei2,100.0))
    errcor=abs(diffoctcarr([cor1_o[0,1:]],cor))

    #find maximum error for this run
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

    #print errIpk
    #print errk
    #print errEloss
    #print errdE_E

    #compare to maximum of previous runs and keep highest error
    maxerr=max(temp,maxerr)

    if maxerr==temp and track_worst:
        print 'new max'
        print cnt
        print temp2


print "After {0} runs with random inputs the".format(testnum)
print "Maximum difference is {0}".format(maxerr)

tol=1e-8
if maxerr<tol:
    print "******************************"
    print "*********    PASS    *********"
    print "******************************"
else:
    print "******************************"
    print "*********    FAIL    *********"
    print "******************************"

   



