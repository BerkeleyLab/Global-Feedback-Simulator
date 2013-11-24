#!/usr/bin/python

#
# A Series of Unit tests for filter.c/h, step_llrf.c/h
#
# Alejandro F Queiruga
# Daniel Scott Driver
# 2013 LBL
#

import linac
import linac_pretty_print


import numpy as np
import matplotlib.pylab as plt
import scipy.linalg as linalg

####################################
#
# Unit test for Filter_Step
#
####################################
def step_response(dt=0.002):
    tmax = 5.0
    nt = int(tmax/dt)
    st = np.zeros(nt,dtype=np.complex)

    fil = linac.Filter() #declare one in python and then set it up
    linac.Filter_Allocate_In(fil,3,3)
    
    poles = linac.complexdouble_Array(1)
    poles[0] = -1.0+1.0j
    linac.Filter_Append_Modes(fil,poles,1,dt)
    poles[0] = -1.0-1.0j
    linac.Filter_Append_Modes(fil,poles,1,dt)

    filnow = linac.Filter_State()
    linac.Filter_State_Allocate(filnow,fil)
    filpast = linac.Filter_State()
    linac.Filter_State_Allocate(filpast,fil)

    for i in xrange(1,nt):
        st[i] = linac.Filter_Step(fil,1.0+0.0j,
                        filnow,filpast)
        tmp = filpast
        filpast=filnow
        filnow=tmp
        
    trang = np.arange(0,tmax,dt)
    anal = 1.0-np.exp(-trang)*(np.sin(trang)+np.cos(trang))
    plt.plot(trang,st.real,'x',trang,anal,'-')
    plt.title("Step reponse compared to analytic solution")
    print "Normalized error is ", linalg.norm(st.real-anal)/linalg.norm(anal)

####################################
#
# Unit tests for step_llrf helpers
#
####################################

#
# Unit test for phase_shift
#
def unit_phase_shift(inp=1.0j):
    thetas = np.arange(0,2.0*np.pi,np.pi/16.0)
    phs = np.zeros(thetas.shape,dtype=np.complex)
    for i in xrange(len(thetas)):
        phs[i] = linac.phase_shift(inp,thetas[i])
    plt.plot(thetas,phs.real,"-x",thetas,phs.imag,"-o")
    plt.title("Phase shift test")
    plt.show()

#
# Unit test for step_fpga
#
def unit_fpga(dt=0.01):
    fpgap = linac.FPGA_Param()
    print "Configuration:"
    linac.FPGA_Config(fpgap,
                      10.0,1.5,1.0)
    print linac_pretty_print.fpgatostr(fpgap)
    
    stnow = linac.FPGA_State()
    stnow.drive = 0.0
    stnow.state = 0.0
    stpast = linac.FPGA_State()
    stpast.drive = 0.0
    stpast.state =0.0
    
    tmax = 10.0
    nt = int(tmax/dt)
    st = np.zeros(nt,dtype=np.complex)
    drv = np.zeros(nt,dtype=np.complex)
    sta = np.zeros(nt,dtype=np.complex)
    trang = np.arange(0.0,tmax,dt)

    plantx = 0.0
    plantv = 0.0
    plantk = 10.0
    plantxa = np.zeros(nt,dtype=np.complex)
    
    for i in xrange(nt):
        plantv += dt*(-plantk*plantx+stnow.drive-1.0*plantv)
        plantx += dt*plantv
        st[i] = linac.step_PI_fpga(fpgap,dt,  plantx,
                          stnow,stpast)
        drv[i] = stnow.drive
        sta[i] = stnow.state
        
        plantxa[i] = plantx
        tmp = stpast
        stpast=stnow
        stnow=tmp

    plt.plot(trang,plantxa.real,'o',trang,drv.real,'-+',trang,sta.real,'-*')
    plt.title("Testing an FPGA")
    plt.show()

#
# Unit test for saturation
#
def unit_saturate():
    for c in 10.0**np.arange(0.0,3.0,0.5):
        inp = np.arange(0.0,10.0,0.1,dtype=np.complex)
        oup = np.zeros(inp.shape,dtype=np.complex)
        for i in xrange(len(inp)):
            oup[i] = linac.saturate(inp[i],c)
        plt.plot(inp,oup.real,'-x',inp,oup.imag,'-+')
    plt.title("Testing saturate()")
    plt.show()



######################################
#
# Now execute the tests...
#
######################################

def perform_tests():
    print "\n****\nPlotting a the step response to the filter (-1-1j),(-1+1j)"
    step_response()
    
    plt.figure()
    print "\n****\nPlotting some phase shifts..."
    unit_phase_shift(1.0j)
    unit_phase_shift(2.0+3.0j)
    
    plt.figure()
    print "\n****\nTesting the FPGA module"
    unit_fpga()
    
    plt.figure()
    print "\n****\nTesting Saturate"
    unit_saturate()
    return True

if __name__=="__main__":
    plt.close('all')
    perform_tests()
