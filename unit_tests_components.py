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
from readjson.loadconfig import *

import numpy as np
import matplotlib.pylab as plt


# Because octave writes sqrt(-1)=i, and 
# numpy expects sqrt(-1)=j... >.< 
# TODO: submit a patch to numpy that can handle "i"s in loading
i2j = lambda x:np.complex(x.replace("i","j"))

##########################
#
# Utility Routine
# inputs a file name pointing to the octave data,
# and a single-argument function wrapping the function
# in question, passing it the first column from the file
#  and comparing the result to the second column
#
##########################
def unit_compare(fname,cfunc,showplots=True,TOL=1.0e-14,title=None):
    octdat=np.loadtxt(fname,
            delimiter=",",dtype=np.complex,converters={0:i2j,1:i2j})
    inps = octdat[:,0] #np.arange(0.0,50.0,0.2,dtype=np.complex)
    outs = np.zeros(inps.shape,dtype=np.complex)
    for i in xrange(len(inps)):
        outs[i] = cfunc(inps[i])
    if showplots:
        plt.plot(inps.real,outs.real,'b-+',label='pyreal')
        plt.plot(inps.real,outs.imag,'g-+',label='pyimag')
        plt.plot(octdat[:,0].real,octdat[:,1].real,'b-x',label='octreal')
        plt.plot(octdat[:,0].real,octdat[:,1].imag,'g-x',label='octrimag')
        if title==None:
            plt.title("Unit test on "+fname)
        else:
            plt.title(title)

        plt.legend()
        plt.show()
    norm = np.linalg.norm(outs-octdat[:,1])/np.linalg.norm(octdat[:,1])
    passed = norm <= TOL

    return [passed,norm]


##########################
#
# Unit Test for the response of the Triode system
# Compares it to data files generated from the Octave model
#
##########################
def unit_triode(showplots=True,TOL=1.0e-14):
    # Load a linac configuration file and select out the first one 
    [pa,allaccell,linp,gun, bbf,nrsc] = LoadConfig("configfiles/all_config.cfg")
    linnow = linac.Linac_State()
    linpast = linac.Linac_State()
    linac.Linac_State_Allocate(linnow,linp[0]);
    linac.Linac_State_Allocate(linpast,linp[0]);
    #print lintostr(linp[0])


    print "  Testing Triode..."
    #
    # First Unit Test: 
    # Keep all of the state variables at zero, 
    # and advance a timestep with a varying input signal
    #
    pass_first,norm_first = unit_compare(
        "unit_test_data/triode_zero_load.csv",
        lambda x: linac.step_triode(linp[0],x,linnow,linpast),
        showplots=showplots, TOL=TOL)
    print "   First test: ||py-oct||_2 = ",norm_first

    #
    # Second Unit Test:
    # Fix the input signal at unity, and let the system "ride out"
    # for some number of timesteps
    #
    if showplots: plt.figure()
    lins = [ linnow,linpast ];
    def ride(x):
        
        out = linac.step_triode(linp[0],1.0,lins[0],lins[1])
        tmp=lins[1]
        lins[1]=lins[0]
        lins[0]=tmp
        return out
    pass_second,norm_second = unit_compare(
        "unit_test_data/triode_unity_ride.csv",
        ride, showplots=showplots, TOL=TOL)
    print "   Second test: ||py-oct||_2 = ",norm_second

    
    if pass_first and pass_second:
        print "o Triode Unit Test Passed."
        return True
    else:
        print "x TRIODE UNIT TEST FAILED!!!!"
        return False


################################
# 
# Unit Tests for the response of the Cavity system, including filter
#
################################
def unit_cavity(showplots=True,TOL=1.0e-8):
    [pa,allaccell,linp,gun, bbf,nrsc] = LoadConfig("configfiles/all_config.cfg")

    linnow = linac.Linac_State()
    linpast = linac.Linac_State()
    linac.Linac_State_Allocate(linnow,linp[0]);
    linac.Linac_State_Allocate(linpast,linp[0]);

    print "  Testing Cavity..."
    
    #
    # First Unit Test: 
    # Keep all of the state variables at zero, 
    # and advance a timestep with a varying drive_in signal
    #
    pass_first,norm_first = unit_compare(
        "unit_test_data/cavity_zero_drive.csv",
        lambda x: linac.step_cavity(linp[0],0.0,x,1.0,
                                    linnow,linpast),
        showplots=showplots, TOL=TOL )
    print "   First test: ||py-oct||_2 = ",norm_first
 

    #
    # Second Unit Test: 
    # Keep all of the state variables at zero, 
    # and advance a timestep with a varying beam_charge signal
    #
    if showplots: plt.figure()
    pass_second,norm_second = unit_compare(
        "unit_test_data/cavity_zero_beam.csv",
        lambda x: linac.step_cavity(linp[0],0.0,1.0,x,
                                    linnow,linpast),
        showplots=showplots, TOL=TOL )
    print "   Second test: ||py-oct||_2 = ",norm_second

    #
    # Third Unit Test: 
    # Keep all of the state variables at zero, 
    # and advance a timestep with a varying delta_tz signal
    #
    if showplots: plt.figure()
    pass_third,norm_third = unit_compare(
        "unit_test_data/cavity_zero_deltatz.csv",
        lambda x: linac.step_cavity(linp[0],1e-12*x.real,1.0,1.0,
                                    linnow,linpast),
        showplots=showplots, TOL=TOL )
    print "   Third test: ||py-oct||_2 = ",norm_third

    #
    # Fourth Unit Test:
    # Maintain unity inputs, and ride out the cavity to see
    # the transient response
    #
    if showplots: plt.figure()
    lins = [linnow,linpast]
    def ride(x):
        out = linac.step_cavity(linp[0],0.0,1.0,0.1,lins[0],lins[1])
        tmp=lins[1]
        lins[1]=lins[0]
        lins[0]=tmp
        return out
    pass_fourth,norm_fourth = unit_compare(
        "unit_test_data/cavity_unity_ride.csv",
        ride, showplots=showplots, TOL=TOL )
    print "   Fourth test: ||py-oct||_2 = ",norm_fourth


    if pass_first and pass_second and pass_third and pass_fourth:
        print "o Cavity Unit Tests Passed."
        return True
    else:
        print "x CAVITY UNIT TESTS FAILED!!"
        return False




################################
#
# Unit Test for the step_llrf routine, the whole system
#
################################

def unit_step_llrf(showplots=True,TOL=1.0e-8):
    # Load a linac configuration file and select out the first one 
    [pa,allaccell,linp,gun, bbf,nrsc] = LoadConfig("configfiles/all_config.cfg")

    dt = pa["Simulation"]["dt"]
    #print lintostr(linp[0])
    # Allocate the history array and each state
    # linss = linac.Linac_State_Array(5);
    # for i in xrange(5):
    #     newlin = linac.Linac_State()

    #     linac.Linac_State_Allocate(newlin,linp[0])
    #     linss[i] = newlin
    #     linss[i].fpga.state = 1.0*i+2.0j
    linssp = linac.make_state_arrays(5)
    linss = linac.Linac_State_Array_frompointer(linssp)
    for i in xrange(5):
        linac.Linac_State_Allocate(linss[i],linp[0])

    #
    # First Unit Test
    # Load on meas_d from a zero state
    #

    pass_first,norm_first = unit_compare(
         "unit_test_data/llrf_step_beam.csv",
         lambda x: linac.step_llrf(linp[0],dt,0.0, x,0.0, 0,
                                   linss.cast()),
         showplots=showplots,TOL=TOL)
    print "   First test: ||py-oct||_2 = ",norm_first

    #
    # Second Unit Test
    # Load on beam_charge
    #
    plt.figure()
    pass_second,norm_second = unit_compare(
         "unit_test_data/llrf_step_deltatz.csv",
         lambda x: linac.step_llrf(linp[0],dt, 1.0e-12*x.real, 0.1,0.1, 0,
                                   linss.cast()),
         showplots=showplots,TOL=TOL)
    print "   Second test: ||py-oct||_2 = ",norm_second

    #
    # Third Unit Test
    # Unity ride 
    #
    plt.figure()
    def ride(x):
        
        out = linac.step_llrf(linp[0],dt,0.1e-12, 0.1, 0.0, 0,
                                   linss.cast())
        
        # tmp = linss[4]
        # for i in xrange(4):
        #     linss[i+1]=linss[i]
        # linss[0]=tmp
        linac.cycle_buffer(linss.cast(),5)
        return out

    pass_third,norm_third = unit_compare(
         "unit_test_data/llrf_unity_ride.csv",
         ride,
         showplots=showplots,TOL=TOL)
    print "   Third test: ||py-oct||_2 = ",norm_third


    #
    # Deallocate data
    #
    for i in xrange(5):
        linac.Linac_State_Deallocate(linss[i])
    for l in allaccell:
        linac.Linac_Deallocate(l)

    
    if pass_first and pass_second and pass_third:
        print "o Step_LLRF Unit Tests Passed."
        return True
    else:
        print "x STEP_LLRF UNIT TESTS FAILED!!"
        return False
    return False
#
# Run all of them
#
def perform_tests():
    tri = unit_triode()
    plt.figure()
    cav = unit_cavity(TOL=1.0e-8)
    plt.figure()
    step = unit_step_llrf()
    if tri and cav and step:
        print "oo PASSED: unit_tests_components.py oo"
        return True
    else:
        print "xx FAILED: unit_tests_components.py xx"
        return False

if __name__ == "__main__":
    plt.close('all')
    perform_tests()
