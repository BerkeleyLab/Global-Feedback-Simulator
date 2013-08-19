#!/usr/bin/python

import linac

#
# Create a Filter and test it out
#
#fil = linac.Filter_Allocate_New(3,3) #this uses malloc, lets avoid this 
fil = linac.Filter() #declare one in python and then set it up
linac.Filter_Allocate_In(fil,3,3)

#push in some random poles

poles = linac.complexdouble_Array(1)
poles[0] = -1.0+1.0j
linac.Filter_Append_Modes(fil,poles,1,0.04)
poles[0] = -1.0-1.0j
linac.Filter_Append_Modes(fil,poles,1,0.04)

# poles = linac.complexdouble_Array(1)
# poles[0] = -1.0-1.5j
# linac.Filter_Append_Modes(fil,poles,1,1.0)
# poles[0] = -1.0-2.5j
# linac.Filter_Append_Modes(fil,poles,1,1.0)

# poles = linac.complexdouble_Array(2)
# poles[0] = -1.0-3.5j
# poles[1] = 2.1
# linac.Filter_Append_Modes(fil,poles,2,1.0)

# Use the swig carrays.i to get handles to the pointer arrays

fil.A_modes = linac.intArray_frompointer(fil.modes)
fil.A_coeff_start = linac.intArray_frompointer(fil.coeff_start)
fil.A_coeffs = linac.complexdouble_Array_frompointer(fil.coeffs)
fil.A_poles = linac.complexdouble_Array_frompointer(fil.poles)



#
# This creates a new linac
#
p_TRF1 = linac.complexdouble_Array(2)
p_TRF1[0] = 1.0-2.0j
p_TRF1[1] = 3.0-4.0;
p_TRF2 = linac.complexdouble_Array(1)
p_TRF2[0] = 5.0+6.0j;
p_RXF = linac.complexdouble_Array(3)
p_RXF[0] = -18+0j
p_RXF[1] = -9.5+15j
p_RXF[2] = -9.5-15j
for i in xrange(3):
    p_RXF[i] *= 1e6
lin = linac.Linac_Param()
linac.Linac_Config(lin,
                   0.5,
                   1.0,2.0,3.0,4.0,
                   5.0,6.0,7.0,8.0,
                   8.5,
                   9.0,
                   p_TRF1,p_TRF2,p_RXF,
                   10, 11.0,
                   12.0,13.0,14.0,
                   15.0,16.0,
                   17.0,18.0,19.0,
                   30.0e6, 50.0e3
                   )

import linac_pretty_print


filter_dotify(fil)


